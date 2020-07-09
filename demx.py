#!/usr/bin/env python3

"""
# Description

Demultiplexing Illumina reads

1. Demultiplexing P7-index, (index sequence saved at the end of name-line of fastq)

eg: 
@ST-E00318:957:H7VYVCCX2:6:1101:4807:1520 1:N:0:GACGACAT+AGATCTCG
NGCTTTATATATCTTGTGGAAAGGACGAAACACCGCTTTCCTGGTGGAGCAGGAGGTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGAGATCGGAAGAGCACACGTCTGA
+
#AAAFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJ7FFJJJJJJ7A-7FAJJJJJFAFJJJAJJJFAF-FJJJJJJJJJJJJJFJJFAFJJJJJJAJJFJJJJJJJJAJ<JJ<A-77AA7F-F7FFFFFAJFJJFA-<-F7
(

in above example: GACGACAT+AGATCTCG are P7-index + P5-index

2. Demultiplexing inline barcode (saved in the beginning of read1 or read2)

# Usage


# Requirements

## Python packages (using conda)
xopen (=0.9.0)
python-levenshtein (=0.12.0)

"""

import os
import sys
import re
import json
import pickle
import shutil
import argparse
import fnmatch
from xopen import xopen
from contextlib import ExitStack
import Levenshtein as lev # distance
from xopen import xopen
import logging
from multiprocessing import Pool
from itertools import combinations

logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    stream=sys.stdout)
log = logging.getLogger(__name__)
log.setLevel('INFO')


def get_args():
    """
    Demultiplexing 
    """
    parser = argparse.ArgumentParser(description='hiseq_demx')
    parser.add_argument('-1', '--fq1', required=True,
        help='read1 in fastq format, gzipped')
    parser.add_argument('-2', '--fq2', 
        help='read2 in fastq format, gzipped, (optional)') 
    parser.add_argument('-o', '--outdir', required=True,
        help='directory to save the reulsts')
    parser.add_argument('-s', '--index-csv', dest='index_csv', required=True,
        help='index list in csv format, [filename,index1,NULL,barcode]')
    parser.add_argument('--demo', action='store_true',
        help='run demo (1M reads) for demostration, default: off')
    parser.add_argument('-m', '--mismatch', type=int, default=0,
        help='mismatches allowed to search index, default: [0]') 
    parser.add_argument('-x', '--barcode-in-read', type=int, dest='barcode_in_read',
        choices=[1, 2], default=2,
        help='barcode in the 5\' end of, 1:read1 or 2:read2, default: [2]')
    parser.add_argument('-l', '--barcode-n-left', type=int, dest='barcode_n_left',
        default=0, help='bases locate on the left of barcode')
    parser.add_argument('-r', '--barcode-n-right', type=int, dest='barcode_n_right',
        default=0, help='bases locate on the right of barcode')
    parser.add_argument('-p', '--threads', type=int, default=1,
        help='number of threads, default: [1]')
    parser.add_argument('-j', '--parallel-jobs', type=int, dest='parallel_jobs',
        default=1, help='number of josb run in parallel, default: [1]')
    parser.add_argument('-w', '--overwrite', action='store_true',
        help='Overwrite exists files, default: off')
    args = parser.parse_args()
    return args
    

class IndexList(object):
    """
    Input the index for demultiplexing
    #file_name, #index1, #index2, #barcode

    criteria:
    1. file_name: unique
    2. index1 + index2 + barcode unique
    3. mismatch required
    """
    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)
        self.init_args()
        self.idx = self.load()
        self.sample = self.idx.get('sample', {})
        self.count = self.idx.get('count', {})
        self.main = self.idx.get('main', {})
        # sequences for each index
        m = [i.split(',') for i in self.sample.keys()]
        self.index1 = list(set([i[0] for i in m if not i == 'NULL'])) #list(self.main.keys()) #
        self.index2 = list(set([i[1] for i in m if not i == 'NULL']))
        self.barcode = list(set([i[2] for i in m if not i == 'NULL']))
        self.status = self.check()


    def init_args(self):
        """
        required arguments: csv/txt/excel
        """
        args_init = {
            'input': None,
            'outdir': os.getcwd(),
            'mismatch': 0
        }
        for k, v in args_init.items():
            if not hasattr(self, k):
                setattr(self, k, v)

        if self.input is None:
            log.error('input, required')


    def load(self):
        """
        Read index info file
        csv file
        index1-index2-barcode
        """
        da = {}
        db = {}
        dx = {}
        with open(self.input) as r:
            for line in r:
                if line.startswith('#'):
                    continue # skip
                s = line.strip().split(',')
                if not len(s) == 4:
                    log.error('input: format unknown: fname,index1,index2,barcode')
                    flag = False
                    break
                fname, ia, ib, bc = s
                ia = ia.upper()
                ib = ib.upper()
                ic = bc.upper()
                n = ','.join([ia, ib, ic])
                da[n] = fname # file name
                db[n] = db.get(n, 0) + 1
                # save all
                dx[ia] = dx.get(ia, {})
                dx[ia][ib] = dx[ia].get(ib, {})
                if ic in dx[ia][ib]:
                    dx[ia][ib][ic] = [dx[ia][ib][ic], fname]
                else:
                    dx[ia][ib][ic] = fname
        # output
        return {
            'sample': da,
            'count': db,
            'main': dx
        }


    def str_distance(self, a, b, partial=True):
        """
        Check index
        """
        if partial:
            a = a[:len(b)]
            b = b[:len(a)]
        return lev.distance(a, b)


    def is_compatible(self, x):
        """
        x list of index
        Check index/barcode, mismatch
        """
        if isinstance(x, str):
            flag = True
        elif isinstance(x, list):        
            # remove null from x
            x = [i.upper() for i in x if not i.upper() == 'NULL']
            flag = True
            for i, k in combinations(x, 2):
                if self.str_distance(i, k) <= self.mismatch:
                    flag = False
            return flag
        elif isinstance(x, dict):
            return self.is_compatible(list(x.keys()))


    def is_valid_index(self, x):
        """
        x, str, list, dict;
        Check the index sequence is valid: ACGTN
        skip: NULL
        """
        if isinstance(x, str):
            x = x.upper()
            return True if x == 'NULL' or re.match('^[ACGTN]*$', x) else False
        elif isinstance(x, list):
            return all([self.is_valid_index(i) for i in x])
        elif isinstance(x, dict):
            return all([self.is_valid_index(i) for i in list(x.keys())])
        else:
            return False


    def check(self):
        """
        Check the index: index1, index2, barcode
        ACGTN only
        barcode, same width
        """
        ## unique all
        m = [i for i in self.count.values() if i > 1] # index unique
        mn = len(self.sample.values()) == len(set(self.sample.values())) # filename unique

        # width
        self.barcode_width = list(set([len(i) for i in self.barcode if not i == 'NULL']))

        # status
        ss = 'ok' if len(m) == 0 else 'failed'
        sn = 'ok' if mn else 'failed'
        s1 = 'ok' if self.is_valid_index(self.index1) else 'failed'
        s2 = 'ok' if self.is_valid_index(self.index2) else 'failed'
        s3 = 'ok' if self.is_valid_index(self.barcode) else 'failed'
        s4 = 'ok' if self.is_compatible(self.index1) else 'failed'
        s5 = 'ok' if self.is_compatible(self.index2) else 'failed'
        s6 = 'ok' if self.is_compatible(self.barcode) else 'failed'
        s7 = 'ok' if len(self.barcode_width) < 2 else 'failed'

        # report
        msg = '\n'.join([
            'Check index status',
            '{0:.<40}: {1:<10}'.format('number of mismatch', self.mismatch),
            '{0:.<40}: {1:<10}'.format('all index unique', ss),
            '{0:.<40}: {1:<10}'.format('filename unique', sn),
            '{0:.<40}: {1:<10}'.format('index-1 [ACGTN]', s1),
            '{0:.<40}: {1:<10}'.format('index-2 [ACGTN]', s2),
            '{0:.<40}: {1:<10}'.format('barcode [ACGTN]', s3),
            '{0:.<40}: {1:<10}'.format('index-1 mismatch', s4),
            '{0:.<40}: {1:<10}'.format('index-2 mismatch', s5),
            '{0:.<40}: {1:<10}'.format('barcode mismatch', s6),
            '{0:.<40}: {1:<10}'.format('barcode width consistent', s7)
        ])

        log.info(msg)
        
        return all([i == 'ok' for i in [ss, sn, s1, s2, s3, s4, s5, s6, s7]])


class Demx(object):
    """
    P7+P5+inline
    SE or PE
    """
    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)
        self.init_args()
        self.save_config()


    def save_config(self):
        """
        Save config to file
        """
        chk0 = args_checker(self.__dict__, self.config_pickle)
        chk1 = args_logger(self.__dict__, self.config_txt)


    def init_args(self):
        """
        Deault arguments
        """
        args_init = {
            'fq1': None,
            'fq2': None,
            'outdir': os.getcwd(),
            'index_csv': None,
            'barcode_in_read': 2,
            'barcode_n_left': 0,
            'barcode_n_right': 0,
            'mismatch': 0,
            'threads': 1,
            'parallel_jobs': 4,
            'overwrite': False,
            'demo': False
        }
        for k, v in args_init.items():
            if not hasattr(self, k):
                setattr(self, k, v)

        # config
        self.config_dir = os.path.join(self.outdir, 'config')
        self.config_pickle = os.path.join(self.config_dir, 'arguments.pickle')
        self.config_txt = os.path.join(self.config_dir, 'arguments.txt')
        check_path(self.config_dir)

        # fastq file
        if self.fq1 is None:
            log.error('fq1, file not exists')
            raise ValueError('fq1, fastq file required')

        # subsample reads
        if self.demo:
            self.sample_fq()

        # sample: csv file
        if self.index_csv is None:
            log.error('index_csv,  file not exists')
            raise ValueError('index_csv, csv file required, filename,index1,index2,barcode')
        else:
            self.idx = IndexList(input=self.index_csv, mismatch=self.mismatch)
            self.sample = self.idx.sample # index1,index2,barcode=sample
            self.index1 = self.idx.index1
            self.index2 = self.idx.index2
            self.barcode = self.idx.barcode
            self.barcode_width = self.idx.barcode_width[0] if len(self.idx.barcode_width) > 0 else 0
            self.idx_main = self.idx.main
            if not self.idx.status:
                raise ValueError('Check out the index list file')

        # mismatch
        if not self.mismatch in range(3):
            log.error('mismatch, [0,1,2,3] expected, {} got'.format(self.mismatch))
            self.mismatch = 0

        # check variables: int
        for i in ['barcode_in_read', 'barcode_n_left', 'barcode_n_right', 'mismatch', 'threads', 'parallel_jobs']:
            ix = getattr(self, i, None)
            if not isinstance(ix, int):
                log.error('{}, int expected, {} got'.format(i, type(i).__name__))

        # resource
        self.init_cpu() 

        ## clean outdir
        # tmp = listdir(self.outdir, include_dir=True)
        # if len(tmp) > 0:
        #     raise Exception('outdir: empty dir required, files detected in: {}'.format(self.outdir))


    def sample_fq(self, smp='demo', sample_size=1000000):
        """
        Create a subsample of input fastq files, default: 1M reads
        Run the whole process for demostration
        
        update: fq1, fq2, outdir
        """
        log.info('Running demo with subsample: {} reads'.format(sample_size))
        # update args
        self.outdir = os.path.join(self.outdir, smp)
        self.data_dir = os.path.join(self.outdir, 'data')
        check_path(self.data_dir)
        
        # subsample
        fq1 = os.path.join(self.data_dir, os.path.basename(self.fq1))
        if os.path.exists(fq1):
            log.info('file eixsts, {}'.format(fq1))
        else:
            with xopen(self.fq1, 'rt') as r1, xopen(fq1, 'wt') as w1:
                i = 0 # counter
                for line in r1:
                    i += 1
                    if i > sample_size*4: # fastq file: 4 lines per read
                        break
                    w1.write(line)

        if isinstance(self.fq2, str):
            fq2 = os.path.join(self.data_dir, os.path.basename(self.fq2))
            if os.path.exists(fq2):
                log.info('file exists, {}'.format(fq2))
            else:
                with xopen(self.fq2, 'rt') as r2, xopen(fq2, 'wt') as w2:
                    i = 0 # counter
                    for line in r2:
                        i += 1
                        if i > sample_size*4: # fastq file: 4 lines per read
                            break
                        w2.write(line)
        else:
            fq2 = None

        # update fq1, fq2
        self.fq1 = fq1
        self.fq2 = fq2


    def init_cpu(self):
        """
        threads, CPUs
        """
        ## check number of threads, parallel_jobs
        ## parallel jobs * threads
        n_cpu = os.cpu_count() # alternative: multiprocessing.cpu_count()

        max_jobs = int(n_cpu / 10.0)
        if max_jobs < 1:
            max_jobs = 1
        ## check parallel_jobs (max: 1/10 of n_cpus)
        if self.parallel_jobs > max_jobs: 
            log.warning('Too large, change parallel_jobs from {} to {}'.format(
                self.parallel_jobs, max_jobs))
            self.parallel_jobs = max_jobs

        ## check threads
        max_threads = int(0.3 * n_cpu / self.parallel_jobs)
        if self.threads * self.parallel_jobs > 0.3 * n_cpu:
            log.warning('Too large, change threads from {} to {}'.format(
                self.threads, max_threads))
            self.threads = max_threads       


    def mission(self):
        """
        Determine the mission of demx
        [pe/se][p7][p5][barcode]
        0101
        0100
        0001
        0000
        1101
        1100
        1001
        1000
        """
        is_pe = 0 if self.fq2 is None else 1 # PE
        is_p7 = 0
        is_p5 = 0
        is_bc = 0
        if not(all([i == 'NULL' for i in self.index1])):
            is_p7 = 1
        if not(all([i == 'NULL' for i in self.index2])):
            is_p5 = 1
        if not(all([i == 'NULL' for i in self.barcode])):
            is_bc = 1

        code = [is_pe, is_p7, is_p5, is_bc]
        log.info('Mission: pe:{} p7:{} p5:{} barcode:{}'.format(is_pe, is_p7, is_p5, is_bc))

        if code == [1, 1, 0, 1]:
            self.run_pe_p7_bc()
        elif code == [1, 1, 0, 0]:
            self.run_pe_p7()
        elif code == [1, 0, 0, 1]:
            self.run_pe_bc()
        elif code == [0, 1, 0, 1]:
            self.run_se_p7_bc()
        elif code == [0, 1, 0, 0]:
            self.run_se_p7()
        elif code == [0, 0, 0, 1]:
            self.run_se_bc()
        else:
            log.error('unknown mission')


    def run_se_p7_bc(self):
        """
        p7 index
        barcode
        """
        # step1. index1
        outdir1 = os.path.join(self.outdir, '_tmp')
        with xopen(self.fq1, 'rt') as r:
            self.index_se(r, outdir1)
        flist1 = self.wrap_dir(outdir1, 'index1') # filenames, not-filenames

        ## parallel-jobs
        n_jobs = min(self.parallel_jobs, len(flist1)) # number of jobs
        with Pool(processes=n_jobs) as pool:
            pool.map(self.run_se_barcode_single, flist1)
        # for fq in flist1:
        #     flist2 = self.run_se_barcode_single(fq)

        # rename files
        self.wrap_read_count()
        self.wrap_file()


    def run_se_p7(self):
        """
        structure: 
        |--_tmp
        |--index1
        |    |--barcode
        |    |    |--file.fq
        """
        # step1. index1
        outdir1 = os.path.join(self.outdir, '_tmp')
        with xopen(self.fq1, 'rt') as r:
            self.index_se(r, outdir1)
        flist1 = self.wrap_dir(outdir1, 'index1') # filenames, not-filenames

        # save file
        self.wrap_read_count()
        self.wrap_file()


    def run_se_bc(self):
        # step1. index1
        outdir1 = os.path.join(self.outdir, '_tmp')
        with xopen(self.fq1, 'rt') as r1:
            self.barcode_se(r1, outdir1)
        # rename files in dirs/
        flist1 = self.wrap_dir(outdir1, mode='barcode')

        # save files
        self.wrap_read_count()
        self.wrap_file()


    def run_pe_p7_bc(self):
        """
        structure: 
        |--_tmp
        |--index1
        |    |--barcode
        |    |    |--file.fq
        """
        # step1. index1
        outdir1 = os.path.join(self.outdir, '_tmp')
        with xopen(self.fq1, 'rt') as r1, xopen(self.fq2, 'rt') as r2:
            self.index_pe(r1, r2, outdir1)
        # rename files in dirs/
        flist1 = self.wrap_dir(outdir1, mode='index1')
        flist1 = sorted(flist1) # sort, read1/read2

        # step2. index2
        # step3. barcode
        # multiple jobs
        n_jobs = min(self.parallel_jobs, len(flist1[0::2])) # number of jobs
        with Pool(processes=n_jobs) as pool:
            pool.map(self.run_pe_barcode_single, flist1[0::2]) # read1
        # for fq1 in flist1[0::2]: # read1
        #     self.run_pe_barcode_single(fq1)

        ## step4. rename files
        self.wrap_read_count()
        self.wrap_file()


    def run_pe_p7(self):
        """
        structure: 
        |--_tmp
        |--index1
        |    |--barcode
        |    |    |--file.fq
        """
        # step1. index1
        outdir1 = os.path.join(self.outdir, '_tmp')
        with xopen(self.fq1, 'rt') as r1, xopen(self.fq2, 'rt') as r2:
            self.index_pe(r1, r2, outdir1)
        # rename files in dirs/
        flist1 = self.wrap_dir(outdir1, mode='index1')

        # save files
        self.wrap_read_count()
        self.wrap_file()
        

    def run_pe_bc(self):
        # step1. index1
        outdir1 = os.path.join(self.outdir, '_tmp')
        with xopen(self.fq1, 'rt') as r1, xopen(self.fq2, 'rt') as r2:
            self.barcode_pe(r1, r2, outdir1)
        # rename files in dirs/
        flist1 = self.wrap_dir(outdir1, mode='barcode')

        # save files
        self.wrap_read_count()
        self.wrap_file()


    def fq_merge(self, fout, qlist):
        """
        Compress, multiple fastq files into single file
        """
        with xopen(fout, 'wb') as w:
            for q in qlist:
                with xopen(q, 'rb') as r:
                    shutil.copyfileobj(r, w)


    def wrap_dir(self, x, mode):
        """
        rename fq files in dir
        """
        dirs = []
        if mode == 'index1':
            for f in listfile(x, "*.fq"):
                fname = fq_name(f, pe_fix=True)
                ix = fname + ',NULL,NULL'
                f_new = os.path.join(x, self.sample.get(ix, fname))
                # check suffix
                f_new += '_1.fq' if f.endswith('_1.fq') else ('_2.fq' if f.endswith('_2.fq') else '.fq')
                # check, index1 contains barcode or not
                if f == f_new:
                    if fname in self.sample.values():
                        # files already renamed 
                        log.info('file exists, {}'.format(f))
                    else:
                        dirs.append(f)
                else:
                    # not in self.sample.keys()
                    if os.path.exists(f_new):
                        log.warning('file exists: {}'.format(f_new))
                    else:
                        os.rename(f, f_new)
        elif mode == 'barcode':
            for f in listfile(x, "*.fq"):
                fname = fq_name(f, pe_fix=True)
                if fq_name(x) == '_tmp': # for bc only
                    ix = 'NULL,NULL,' + fname
                else:
                    ix =  fq_name(x) + ',NULL,' + fname
                # ix =  'NULL,NULL,' + fname
                f_new = os.path.join(x, self.sample.get(ix, fname))
                f_new += '_1.fq' if f.endswith('_1.fq') else ('_2.fq' if f.endswith('_2.fq') else '.fq')
                # check, index1 contains barcode or not
                if f == f_new:
                    # in self.samples.keys()
                    if fname in self.sample.values():
                        # files already renamed
                        log.info('file exists, {}'.format(f))
                    else:
                        dirs.append(f)
                else:
                    # not in self.samples.keys()
                    if os.path.exists(f_new):
                        log.warning('file exists: {}'.format(f_new))
                    else:
                        os.rename(f, f_new)
        else:
            pass

        return dirs


    def wrap_file(self):
        """
        move "_tmp/" files to "outdir", 
        combine undemx file
        """
        f_undemx = []
        f_hits = []
        # level-1: index-1
        dv1 = os.path.join(self.outdir, '_tmp')
        for d1 in listdir(dv1, include_dir=True):
            d1_name = fq_name(d1, pe_fix=True)
            # for exists files
            if d1_name in self.sample.values():
                f_hits.append(d1)
            elif d1_name == 'undemx':
                f_undemx.append(d1)
            elif os.path.isdir(d1): # barcode level
                # level-2: barcode
                for d2 in listfile(d1, "*.fq"):
                    d2_name = fq_name(d2, pe_fix=True)
                    if d2_name in self.sample.values():
                        f_hits.append(d2)
                    elif d2_name == "undemx":
                        f_undemx.append(d2)
                    else:
                        pass
            else:
                pass
        ## save files
        # ## multiple threads
        # n_jobs = 4
        # with Pool(processes=n_jobs) as pool:
        #     pool.map(self.compress_output, f_hits) # read1

        for f in f_hits:
            f_out = os.path.join(self.outdir, os.path.basename(f) + '.gz')
            log.info('Saving file: {}'.format(f_out))
            if os.path.exists(f_out):
                log.warning('file exists, {}'.format(f_out))
            else:
                with xopen(f, 'rb') as r:
                    with xopen(f_out, 'wb') as w:
                        shutil.copyfileobj(r, w)
        
        # save undemx
        log.info('Saving file: undemx')
        f_undemx = sorted(f_undemx)
        if f_undemx[0].endswith('_1.fq'):
            # PE mode
            f_r1_undemx_out = os.path.join(self.outdir, 'undemx_1.fq.gz')
            f_r2_undemx_out = os.path.join(self.outdir, 'undemx_2.fq.gz')
            if os.path.exists(f_r1_undemx_out) and os.path.exists(f_r2_undemx_out):
                log.warning('file exists, merging undemx file skipped...')
            else:
                self.fq_merge(f_r1_undemx_out, f_undemx[0::2]) # read1
                self.fq_merge(f_r2_undemx_out, f_undemx[1::2]) # read2
        else:
            # SE mode
            f_undemx_out = os.path.join(self.outdir, 'undemx.fq.gz')
            if os.path.exists(f_undemx_out):
                log.warning('file exists, merging undemx file skipped...')
            else:            
                self.fq_merge(f_undemx_out, f_undemx)

        # remove temp files "outdir/_tmp"


    def compress_output(self, f_in):
        """
        Compress f_in, save to self.outdir
        """
        f_out = os.path.join(self.outdir, os.path.basename(f_in) + '.gz')
        log.info('Saving file: {}'.format(f_out))
        # pigz faster than gzip
        with xopen(f_in, 'rb') as r:
            with xopen(f_out,'wb') as w:
                shutil.copyfileobj(r, w)


    def wrap_read_count(self):
        """
        parse the "read_number.json" file
        save read counts
        """        
        ##----------------------------------##
        ## for index1
        ## main dir
        f1 = os.path.join(self.outdir, '_tmp', 'read_number.json')
        if os.path.exists(f1):
            da = Json(f1).reader()
        else:
            log.error('file not exists: {}'.format(f1))

        ##----------------------------------##
        ## for barcode (within index)
        db = {} # subdir
        for fv1 in listdir(os.path.join(self.outdir, '_tmp'), include_dir=True):
            fv1_name = fq_name(fv1)
            if not os.path.isdir(fv1):
                continue
            # subdir
            fv2 = os.path.join(fv1, 'read_number.json')
            if os.path.exists(fv2):
                db[fv1_name] = Json(fv2).reader()
            else:
                log.error('file not exists: {}'.format(fv2))

        ##----------------------------------##
        ## sum
        dd = {} # all report


        ##----------------------------------##
        ## index_only / bc_only
        # assign read number for each file
        n_undemx = 0
        dx = {}
        for k, v in da.items():
            ix1 = k + ',NULL,NULL' # for index
            ix2 = 'NULL,NULL,' + k # for bc only
            if ix1 in self.sample:
                dx[self.sample.get(ix1, k)] = v
            if ix2 in self.sample:
                dx[self.sample.get(ix2, k)] = v
            elif k == 'undemx':
                n_undemx += v
            else:
                pass

        ##----------------------------------##
        ## index+bc
        for k1, v1 in db.items():
            dd[k1] = {}
            for k2, v2 in v1.items():
                ix = k1 + ',NULL,' + k2
                if ix in self.sample:
                    smp = self.sample.get(ix, k2)
                    dx[smp] = v2
                    dd[k1][smp] = v2
                elif k2 == 'undemx':
                    n_undemx += v2
                    dd[k1]['undemx'] = v2
                else:
                    pass
        # for undemx
        dx['undemx'] = n_undemx

        # for all
        dd['main'] = {}
        for k, v in sorted(dx.items(), key=lambda x: x[0], reverse=False):
            dd['main'][k] = v

        ##----------------------------------##
        # save to final report
        self.demx_report = dd
        f_out = os.path.join(self.outdir, 'demx_report.json')
        Json(dd).writer(f_out)


    def run_se_barcode_single(self, fq):
        """
        Demultiplex barcode: read again
        """
        fname = fq_name(fq, pe_fix=True)
        if fname == 'undemx':
            return None # skip
        outdir = os.path.join(os.path.dirname(fq), fname)
        with xopen(fq) as r:
            self.barcode_se(r, outdir, index1=fname)
        return self.wrap_dir(outdir, 'barcode')


    def run_pe_barcode_single(self, fq1):
        """
        Demultiplex barcode: read again
        """
        fname = fq_name(fq1, pe_fix=True)
        fq2 = fq1.replace('_1.fq', '_2.fq')
        if fname == 'undemx':
            return None # skip
        outdir = os.path.join(os.path.dirname(fq1), fname)
        with open(fq1, 'rt') as r1, open(fq2, 'rt') as r2:
            self.barcode_pe(r1, r2, outdir, index1=fname)
        return self.wrap_dir(outdir, 'barcode')


    def search(self, x, mode='index1', index1=None):
        """
        Search by index
        support for: wrap dir
        return the sample name
        """
        if mode == 'index1':
            h = [i for i in self.index1 if self.str_distance(i, x) <= self.mismatch]
        elif mode == 'barcode':
            h = [i for i in self.get_barcode(index1) if self.str_distance(i, x) <= self.mismatch]
        else:
            h = []

        return h


    def readfq(self, fh): # this is a generator function
        """
        source: https://github.com/lh3/readfq/blob/master/readfq.py
        processing fastq file
        """
        last = None # this is a buffer keeping the last unprocessed line
        while True: # mimic closure; is it a bad idea?
            if not last: # the first record or a record following a fastq
                for l in fh: # search for the start of the next record
                    if l[0] in '>@': # fasta/q header line
                        last = l[:-1] # save this line
                        break
            if not last: break
            [name, _, comment], seqs, last = last[1:].partition(" "), [], None
            for l in fh: # read the sequence
                if l[0] in '@+>':
                    last = l[:-1]
                    break
                seqs.append(l[:-1])
            if not last or last[0] != '+': # this is a fasta record
                yield name, ''.join(seqs), None, comment # yield a fasta record
                if not last: break
            else: # this is a fastq record
                seq, leng, seqs = ''.join(seqs), 0, []
                for l in fh: # read the quality
                    seqs.append(l[:-1])
                    leng += len(l) - 1
                    if leng >= len(seq): # have read enough quality
                        last = None
                        yield name, seq, ''.join(seqs), comment; # yield a fastq record
                        break
                if last: # reach EOF before reading enough quality
                    yield name, seq, None, comment # yield a fasta record instead
                    break


    def str_distance(self, a, b, partial=True):
        """
        a, string 
        b, string
        partial, bool, whether compare the same length
        The distance between a and b
        """
        if partial:
            a = a[:len(b)]
            b = b[:len(a)]
        return lev.distance(a, b)


    def get_barcode(self, index1=None):
        """
        Return the barcode list
        according to the index1
        """
        if index1 in self.idx_main:
            return list(self.idx_main[index1]['NULL'].keys())
        else:
            return self.barcode


    # p7 index
    def index_se(self, fh, outdir):
        """
        Demultiplex Illumina Hiseq fastq file
        single index or dual index
        1. [ACGTN]{8}
        2. [ACGTN]{8}+[ACGTN]{8} 
        """
        log.info('Running index_se, {}'.format(outdir))
        check_path(outdir)

        # all availabe output files
        fn = {} # fq number
        fnames = self.index1
        fnames.append('undemx')
        fouts = [os.path.join(outdir, i + '.fq') for i in fnames]

        # check status
        fn_json = os.path.join(outdir, 'read_number.json')
        if os.path.exists(fn_json):
            log.info('Skipped index_se, {}, read_number.json found'.format(outdir))
            return fouts        

        # open multiple files at same time
        with ExitStack() as stack:
            fws = [stack.enter_context(open(f, 'wt')) for f in fouts]
            n = 0 # counter
            for name, seq, qual, comment in self.readfq(fh):
                n += 1
                if n%1000000 == 0 :
                    log.info('Processed reads: {}'.format(n))

                fq = '\n'.join(['@' + name + ' ' + comment, seq, '+', qual])
                index = comment.split(':')[-1] # '2:N:0:ATCACGAT+AGATCTCG'
                # index = index.partition('+')[0] # for single index only
                index = index.split('+')[0]
                fhit = self.search(index, 'index1')
                fname = fhit.pop() if len(fhit) == 1 else 'undemx'
                fw = fws[fnames.index(fname)]
                fw.write(fq + '\n')
                fn[fname] = fn.get(fname, 0) + 1
        
        # save dict to json
        Json(fn).writer(fn_json)        

        return fouts


    # p7 index: 
    def index_pe(self, fh1, fh2, outdir):
        """
        Demultiplex Illumina Hiseq fastq file
        single index or dual index
        1. [ACGTN]{8}
        2. [ACGTN]{8}+[ACGTN]{8}
        """
        log.info('Running index_pe, {}'.format(outdir))
        check_path(outdir)

        # all availabe output files
        fn = {} # fq read number
        fnames = self.index1
        fnames.append('undemx')
        r1_fouts = [os.path.join(outdir, i + '_1.fq') for i in fnames]
        r2_fouts = [os.path.join(outdir, i + '_2.fq') for i in fnames]


        # check status
        fn_json = os.path.join(outdir, 'read_number.json')
        if os.path.exists(fn_json):
            log.info('Skipped index_pe, {}, read_number.json found'.format(outdir))
            return (r1_fouts, r2_fouts)

        # open multiple files at same time
        with ExitStack() as stack:
            fws = [stack.enter_context(open(f, 'wt')) for f in r1_fouts + r2_fouts]
            n = 0 # counter
            for r1, r2 in zip(self.readfq(fh1), self.readfq(fh2)):
                n += 1
                if n%1000000 == 0:
                    log.info('Processed reads: {}'.format(n))

                r1_name, r1_seq, r1_qual, r1_comment = r1
                r2_name, r2_seq, r2_qual, r2_comment = r2
                fq1 = '\n'.join(['@' + r1_name + ' ' + r1_comment, r1_seq, '+', r1_qual])
                fq2 = '\n'.join(['@' + r2_name + ' ' + r2_comment, r2_seq, '+', r2_qual])
                # check output
                index = r1_comment.split(':')[-1] # '2:N:0:ATCACGAT+AGATCTCG'
                # index = index.partition('+')[0] # for single index only
                index = index.split('+')[0]
                fhit = self.search(index, 'index1')
                fname = fhit.pop() if len(fhit) == 1 else 'undemx'
                fw1 = fws[fnames.index(fname)]
                fw2 = fws[fnames.index(fname) + len(fnames)]
                fw1.write(fq1 + '\n')
                fw2.write(fq2 + '\n')
                fn[fname] = fn.get(fname, 0) + 1
        
        # save dict to json
        Json(fn).writer(fn_json)


    # barcode
    def barcode_se(self, fh, outdir, index1=None):
        """
        check barcode
        save as files
        """
        log.info('Running barcode_se, {}'.format(outdir))
        check_path(outdir)

        # add undemx
        fn = {} # fq read number
        bc = self.get_barcode(index1)
        bc.append('undemx')
        fouts = [os.path.join(outdir, i + '.fq') for i in bc]

        # check status
        fn_json = os.path.join(outdir, 'read_number.json')
        if os.path.exists(fn_json):
            log.info('Skipped barcode_se, {}, read_number.json found'.format(outdir))
            return fouts

        # open multiple files at same time
        with ExitStack() as stack:
            fws = [stack.enter_context(open(f, 'wt')) for f in fouts]
            n = 0 # counter
            for name, seq, qual, comment in self.readfq(fh):
                n += 1
                if n%1000000 == 0 :
                    log.info('Processed reads: {}'.format(n))
                fq = '\n'.join(['@' + name + ' ' + comment, seq, '+', qual])
                # fq = '@' + name + ' ' + comment + '\n' + seq + '\n+\n' + qual
                s = self.barcode_n_left
                w = self.barcode_width
                index = seq[s:(s+w)]
                fhit = self.search(index, 'barcode', index1) # search names
                fname = fhit.pop() if len(fhit) == 1 else 'undemx'
                fw = fws[bc.index(fname)]
                fw.write(fq + '\n')
                fn[fname] = fn.get(fname, 0) + 1
        
        # save dict to json
        Json(fn).writer(fn_json)
        return fouts


    # barcode
    def barcode_pe(self, fh1, fh2, outdir, index1=None):
        """
        check barcode
        save as files
        """
        log.info('Running barcode_pe, {}'.format(outdir))
        check_path(outdir)

        # add undemx
        fn = {} # fq read number 
        bc = self.get_barcode(index1)
        bc.append('undemx')
        r1_fouts = [os.path.join(outdir, i + '_1.fq') for i in bc]
        r2_fouts = [os.path.join(outdir, i + '_2.fq') for i in bc]

        # check status
        fn_json = os.path.join(outdir, 'read_number.json')
        if os.path.exists(fn_json):
            log.info('Skipped barcode_pe, {}, read_number.json found'.format(outdir))
            return (r1_fouts, r2_fouts)

        # open multiple files at same time
        with ExitStack() as stack:
            fws = [stack.enter_context(open(f, 'wt')) for f in r1_fouts + r2_fouts]
            n = 0 # counter
            for r1, r2 in zip(self.readfq(fh1), self.readfq(fh2)):
                n += 1
                if n%1000000 == 0:
                    log.info('Processed reads: {}'.format(n))

                r1_name, r1_seq, r1_qual, r1_comment = r1
                r2_name, r2_seq, r2_qual, r2_comment = r2
                fq1 = '\n'.join(['@' + r1_name + ' ' + r1_comment, r1_seq, '+', r1_qual])
                fq2 = '\n'.join(['@' + r2_name + ' ' + r2_comment, r2_seq, '+', r2_qual])

                # check barcode
                s = self.barcode_n_left
                w = self.barcode_width
                seq = r1_seq if self.barcode_in_read == 1 else r2_seq # which read1/2
                index = seq[s:(s+w)]
                fhit = self.search(index, 'barcode') # search names
                # print('!AAAA-1', index, fhit)
                fname = fhit.pop() if len(fhit) == 1 else 'undemx'
                fw1 = fws[bc.index(fname)]
                fw2 = fws[bc.index(fname) + len(bc)]
                fw1.write(fq1 + '\n')
                fw2.write(fq2 + '\n')
                fn[fname] = fn.get(fname, 0) + 1
        
        # save dict to json
        Json(fn).writer(fn_json)
        return (r1_fouts, r2_fouts)


    def run(self):
        self.mission()
        total = sum(self.demx_report['main'].values())
        # report
        msg = ['RT {:-^74}'.format('Demx Report: BEGIN')]
        msg.append('RT {0:>3} {1:<50} {2:>10}  {3:>7}'.format('num', 'filename', 'count', 'percent'))
        i = 0
        for k, v in self.demx_report['main'].items():
            i += 1
            msg.append('RT {0:>3} {1:<50} {2:>10} {3:>7.2f}%'.format(i, k, v, v / total * 100))
        # end
        msg.append('RT {0:>3} {1:<50} {2:>10} {3:>7}%'.format('', 'sum', total, '100.00'))
        msg.append('RT {:-^74}'.format('Demx Report: END'))

        msg_log = '\n'.join(msg) #

        # save report to file
        report_txt = os.path.join(self.outdir, 'report.txt')
        with open(report_txt, 'wt') as w:
            w.write(msg_log + '\n')

        print(msg_log)
        log.info('Demultiplexing finished!')


def args_checker(d, x, update=False):
    """Check if dict and x are consitent
    d is dict
    x is pickle file
    """
    assert isinstance(d, dict)
    flag = None
    if os.path.exists(x):
        # read file to dict
        with open(x, 'rb') as fh:
            d_checker = pickle.load(fh)
        if d == d_checker:
            flag = True
        else:
            if update:
                with open(x, 'wb') as fo:
                    pickle.dump(d, fo, protocol=pickle.HIGHEST_PROTOCOL)
    elif isinstance(x, str):
        # save dict to new file
        with open(x, 'wb') as fo:
            pickle.dump(d, fo, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        log.error('illegal x= argument: %s' % x)

    return flag


def args_logger(d, x, overwrite=False):
    """Format dict, save to file
        key: value
    """
    assert isinstance(d, dict)
    n = ['%30s |    %-40s' % (k, d[k]) for k in sorted(d.keys())]
    if os.path.exists(x) and overwrite is False:
        return True
    else:
        with open(x, 'wt') as fo:
            fo.write('\n'.join(n) + '\n')
        return '\n'.join(n)


def check_path(x, show_log=False, create_dirs=True):
    """
    Check if x is path, Create path
    """
    if isinstance(x, str):
        if os.path.isdir(x):
            tag = True
        else:
            if create_dirs is True:
                try:
                    os.makedirs(x)
                    tag = True
                except:
                    tag = False
            else:
                tag = False
        # show log
        flag = 'ok' if tag is True else 'failed'
        if show_log is True:
            log.info('{:<6s} : {}'.format(flag, x))
        return tag
    elif isinstance(x, list):
        return all([check_path(i, show_log, create_dirs) for i in x])
    else:
        log.warning('expect str and list, not {}'.format(type(x)))
        return None


class Json(object):
    """
    Wrapper for *.json file
    1.json to dict
    2.dict to json -> file
    """
    def __init__(self, input, **kwargs):
        self.input = input
        self.mission()


    def mission(self):
        """
        Check what to do, based on the input args
        """
        # json_file to dict
        if input is None:
            log.warning('require, dict or str (file), Nonetype detected')
        elif isinstance(self.input, str) and os.path.exists(self.input):
            self.json_file_in = self.input
            # self.reader() # updated
        elif isinstance(self.input, dict):
            self.d = self.input
            # self.writer()
        else:
            log.warning('expect dict or str, get {}'.format(type(self.input).__name__))


    def _tmp(self):
        """
        Create a tmp file to save json object
        """
        tmp = tempfile.NamedTemporaryFile(prefix='tmp', suffix='.json',
            delete=False)
        return tmp.name


    def to_dict(self):
        # self.json_file_in = self.input
        return self.reader()


    def to_json(self, file=None):
        return self.writer(file)


    def reader(self):
        """
        Read json file as dict
        """
        if isinstance(self.input, dict):
            self.d = self.input
        else:
            try:
                with open(self.json_file_in) as r:
                    if os.path.getsize(self.json_file_in) > 0:
                        self.d = json.load(r)
                    else:
                        self.d = {}
            except:
                log.warning('failed reading json file: {}'.format(self.json_file_in))
                self.d = {}

        return self.d


    def writer(self, json_file_out=None):
        """
        Write d (dict) to file x, in json format
        """
        # save to file
        if json_file_out is None:
            json_file_out = self._tmp()

        try:
            with open(json_file_out, 'wt') as w:
                json.dump(self.d, w, indent=4, sort_keys=True)
        except:
            log.warning('failed saving file: {}'.format(json_file_out))

        return json_file_out


## 1. files and path ##
def listdir(path, full_name=True, recursive=False, include_dir=False):
    """
    List all the files within the path
    """
    out = []
    for root, dirs, files in os.walk(path):
        if full_name:
            dirs = [os.path.join(root, d) for d in dirs]
            files = [os.path.join(root, f) for f in files]
        out += files

        if include_dir:
            out += dirs

        if recursive is False:
            break

    return sorted(out)


def listfile(path='.', pattern='*', full_name=True, recursive=False):
    """
    Search files by the pattern, within directory
    fnmatch.fnmatch()

    pattern:

    *       matches everything
    ?       matches any single character
    [seq]   matches any character in seq
    [!seq]  matches any char not in seq

    An initial period in FILENAME is not special.
    Both FILENAME and PATTERN are first case-normalized
    if the operating system requires it.
    If you don't want this, use fnmatchcase(FILENAME, PATTERN).

    example:
    listfile('./', '*.fq')
    """
    fn_list = listdir(path, full_name, recursive, include_dir=False)
    fn_list = [f for f in fn_list if fnmatch.fnmatch(os.path.basename(f), pattern)]
    return sorted(fn_list)


def fq_name(fq, include_path=False, pe_fix=False):
    """
    parse the name of fastq file:
    .fq.gz
    .fastq.gz
    .fq
    .fastq
    (also for fasta, fa)
    """
    p1 = re.compile('[.](fast|f)[aq](.gz)?$', re.IGNORECASE)
    p2 = re.compile('[._][12]$', re.IGNORECASE)
    if isinstance(fq, str):
        fq = fq if include_path is True else os.path.basename(fq)
        fq_tmp = re.sub(p1, '', fq) # r1_1.fq.gz : r1_1
        if pe_fix is True:
            fq_tmp = re.sub(p2, '', fq_tmp) # r1_1.fq.gz: r1
        return fq_tmp
    elif isinstance(fq, list):
        return [fq_name(x, include_path=include_path, pe_fix=pe_fix) for x in fq]
    else:
        log.warning('unknown type found: {}'.format(type(fq)))
        return fq


def main():
    args = vars(get_args())
    Demx(**args).run()


if __name__ == '__main__':
    main()
