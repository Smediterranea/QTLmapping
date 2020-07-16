#align multiple fastq files with bwa mem

bwa mem '<zcat R1_001.fastq.gz R1_002.fastq.gz R1_003.fastq.gz' '<zcat R2_001.fastq.gz R2_002.fastq.gz R2_003.fastq.gz' > out.sam

####process radtags, cat fastq first###
process_radtags -1 ./raw/plate1.1.fq.gz   -2 ./raw/plate1.2.fq.gz  -b ./barcodes/index.txt -o ./samples/ -c -q -r -i gzfastq --inline_inline --renz_1 nheI --renz_2 ecoRI 
