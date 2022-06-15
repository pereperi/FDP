#!/bin/bash

kb count --h5ad -i index.idx -g t2g.txt -x 10xv3 -o PNF19 -c1 spliced_t2c.txt -c2 unspliced_t2c.txt --workflow lamanno --filter bustools -t 2 PNF19_S1_L001_R1_001.fastq.gz PNF19_S1_L001_R2_001.fastq.gz PNF19_S1_L002_R1_001.fastq.gz PNF19_S1_L002_R2_001.fastq.gz 

kb count --h5ad -i index.idx -g t2g.txt -x 10xv3 -o PNF20 -c1 spliced_t2c.txt -c2 unspliced_t2c.txt --workflow lamanno --filter bustools -t 2 PNF20_S1_L001_R1_001.fastq.gz PNF20_S1_L001_R2_001.fastq.gz PNF20_S1_L002_R1_001.fastq.gz PNF20_S1_L002_R2_001.fastq.gz 

kb count --h5ad -i index.idx -g t2g.txt -x 10xv3 -o PNF23 -c1 spliced_t2c.txt -c2 unspliced_t2c.txt --workflow lamanno --filter bustools -t 2 PNF23_S1_L001_R1_001.fastq.gz PNF23_S1_L001_R2_001.fastq.gz PNF23_S1_L002_R1_001.fastq.gz PNF23_S1_L002_R2_001.fastq.gz 

