# Hiseq_demx

Demultiplexing Illumina sequencing reads



## Getting Started

### Prerequisites

This script was tested for `python 3.6` only, and the following two modules are required, (recommend `conda` to manage the modules)

+ xopen (=0.9.0)
+ python-levenshtein (=0.12.0)

```
conda install xopen=0.9.0 python-levenshtein=0.12.0
```

### Installing

No need to install the script. Just run the script directly.

1. Clone this repo to your local machine   

```
$ git clone https://github.com/bakerwm/hiseq_demx.git
```

2. It works, if you can see the following message    

```
$ cd hiseq_demx
$ python demx.py -h

usage: demx.py [-h] -1 FQ1 [-2 FQ2] -o OUTDIR -s INDEX_CSV [--demo]
               [-m MISMATCH] [-x {1,2}] [-l BARCODE_N_LEFT]
               [-r BARCODE_N_RIGHT] [-p THREADS] [-j PARALLEL_JOBS] [-w]

hiseq_demx

optional arguments:
  -h, --help            show this help message and exit
  -1 FQ1, --fq1 FQ1     read1 in fastq format, gzipped
  -2 FQ2, --fq2 FQ2     read2 in fastq format, gzipped, (optional)
  -o OUTDIR, --outdir OUTDIR
                        directory to save the reulsts
  -s INDEX_CSV, --index-csv INDEX_CSV
                        index list in csv format,
                        [filename,index1,NULL,barcode]
  --demo                run demo (1M reads) for demostration, default: off
  -m MISMATCH, --mismatch MISMATCH
                        mismatches allowed to search index, default: [0]
  -x {1,2}, --barcode-in-read {1,2}
                        barcode in the 5' end of, 1:read1 or 2:read2, default:
                        [2]
  -l BARCODE_N_LEFT, --barcode-n-left BARCODE_N_LEFT
                        bases locate on the left of barcode
  -r BARCODE_N_RIGHT, --barcode-n-right BARCODE_N_RIGHT
                        bases locate on the right of barcode
  -p THREADS, --threads THREADS
                        number of threads, default: [1]
  -j PARALLEL_JOBS, --parallel-jobs PARALLEL_JOBS
                        number of josb run in parallel, default: [1]
  -w, --overwrite       Overwrite exists files, default: off

```

## Running the tests

It support `fastq` file in the following format:  

+ `P7, (P5, optional)` index saved in name line (**1st line**) 

`@ST-E00318:957:H7VYVCCX2:6:1101:27143:1538 1:N:0:CAGATCAT+CGATCTCG`


+ `inline-barcode` located at the beginning of read1 (or read2), (5' end)

`5'--{GAGGAG}------3'`


### P7 index

For `Paired-end` mode, ONLY `P7-index` from `read1` were checked.

```
$ cd hiseq_demx/test
$ python ../demx.py -1 idx_1.fq.gz -2 idx_2.fq.gz -o results/p7/pe -s info_idx.csv

...

RT ----------------------------Demx Report: BEGIN----------------------------
RT num filename                                                count  percent
RT   1 sample1                                                   536   53.60%
RT   2 sample2                                                   149   14.90%
RT   3 undemx                                                    315   31.50%
RT     sum                                                      1000  100.00%
RT -----------------------------Demx Report: END-----------------------------
```

`-1`   : path to `read1` file  
`-2`   : path to `read2` file   
`-o`   : path to directory, saving the results
`-s`   : path to sample_info file (CSV)
`-x 1` : barcode located in `read1`  
`-l 2` : `2` bp on the left of barcode   
`-r 3` : `3` bp on the right of barcode    
`-m 0` : Number of mismatches allowed, for searching barcode


### Inline Barcode

In this example, (iCLIP reads), barcode were located at the 5' end of `read1`, in the following format: `5'-NN{4bp}NNN---`, so the following arguments are requried:

```
$ python ../demx.py -1 iclip_1.fq.gz -2 iclip_2.fq.gz -o results/bc/pe -s info_iclip.csv -x 1 -l 2 -r 3 -m 0 

...

RT ----------------------------Demx Report: BEGIN----------------------------
RT num filename                                                count  percent
RT   1 sample1                                                   200   40.00%
RT   2 sample2                                                   100   20.00%
RT   3 undemx                                                    200   40.00%
RT     sum                                                       500  100.00%
RT -----------------------------Demx Report: END-----------------------------
```


`-1`   : path to `read1` file  
`-2`   : path to `read2` file   
`-o`   : path to directory, saving the results
`-s`   : path to sample_info file (CSV)
`-x 1` : barcode located in `read1`  
`-l 2` : `2` bp on the left of barcode   
`-r 3` : `3` bp on the right of barcode    
`-m 0` : Number of mismatches allowed, for searching barcode


### Both P7 index and Inline Barcode

In this example, inline-barcode (eCLIP-like reads) were located at the 5' end of `read2`, in the following format: `5'{6bp}---`, so the following arguments are requried:

```
$ python ../demx.py -1 idx_eclip_1.fq.gz -2 idx_eclip_2.fq.gz -o results/p7_bc/pe -s info_idx_eclip.csv -x 2 -l 0 -r 1 -m 0 

...

RT ----------------------------Demx Report: BEGIN----------------------------
RT num filename                                                count  percent
RT   1 sample1                                                   100   20.00%
RT   2 sample2                                                   100   20.00%
RT   3 sample3                                                   100   20.00%
RT   4 sample4                                                   100   20.00%
RT   5 undemx                                                    100   20.00%
RT     sum                                                       500  100.00%
RT -----------------------------Demx Report: END-----------------------------
```

`-1`   : path to `read1` file  
`-2`   : path to `read2` file   
`-o`   : path to directory, saving the results
`-s`   : path to sample_info file (CSV)
`-x 2` : barcode located in `read2`  
`-l 0` : `0` bp on the left of barcode   
`-r 1` : `1` bp on the right of barcode    
`-m 0` : Number of mismatches allowed, for searching barcode



### For your data  

+ Prepare a `sample_info.csv` file, including the following columns   

`sample_name`, `P7_index`, `P5_index`, `barcode`  

```
#filename,index1,index2,barcode
sample1,GCCAAT,NULL,NULL
sample2,CAGATC,NULL,NULL
```

## Authors

+ Ming Wang

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details  


## Acknowledgements  


