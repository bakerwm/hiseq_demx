

## 1. barcode only
## - iCLIP barcode, NNACGTNNN at 5' end of read1
python ../demx.py -1 iclip_1.fq.gz -2 iclip_2.fq.gz -o results/bc/pe -s info_iclip.csv -m 0 -x 1 -l 2 -r 3 
python ../demx.py -1 iclip_1.fq.gz -o results/bc/se -s info_iclip.csv -m 0 -x 1 -l 2 -r 3  

## 2. P7 index only
## - index saved at the tail of name_line:
python ../demx.py -1 idx_1.fq.gz -2 idx_2.fq.gz -o results/p7/pe -s info_idx.csv 
python ../demx.py -1 idx_1.fq.gz -o results/p7/se -s info_idx.csv 

## 3. both P7 index and barcode
## - index saved at the tail of name_line:
## - eCLIP barcode, ACGTAC at 5' end of read2
python ../demx.py -1 idx_eclip_1.fq.gz -2 idx_eclip_2.fq.gz -o results/p7_bc/pe -s info_idx_eclip.csv -m 0 -x 2 -l 0 -r 1
python ../demx.py -1 idx_eclip_2.fq.gz -o results/p7_bc/se -s info_idx_eclip.csv -m 0 -x 1 -l 0 -r 1

## for barcode, SE mode
zcat iclip_1.fq.gz | fastx_trimmer -f 3 | fastx_barcode_splitter.pl --bcfile bc.txt --bol --mismatches 0 --prefix aaaaaa. --suffix .fq
