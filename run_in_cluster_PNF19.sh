#!/bin/bash
# request Bourne shell as shell for job
#$ -S /bin/bash
# Name for the script in the queuing system
#$ -N scRNA_PNF19_cellranger
# In order to load the environment variables and your path
# You can either use this or do a : source /etc/profile
#$ -V
# name of the queue you want to use
#$ -q d10imppcv3
# You can redirect the error output to a specific file
#$ -e /imppc/labs/eslab/share/PerePericot/Projects/scRNA_PNF/preprocessed_data/cellranger_cluster_PNF19.log
# You can redirect the output to a specific file
#$ -o /imppc/labs/eslab/share/PerePericot/Projects/scRNA_PNF/preprocessed_data/cellranger_cluster_PNF19.log
# In order to receive an e-mail at the begin of the execution and in the end of it
#$ -m be
# You have to specify an address
#$ -M pere.pericot@alum.esci.upf.edu
#$ -pe smp 27

echo "`date` Launching cell ranger for scRNA_PNF19 data"

cd /imppc/labs/eslab/share/PerePericot/Projects/scRNA_PNF
/imppc/labs/eslab/bgel/Software/CellRanger/cellranger-3.1.0/cellranger count --id=run_count_scRNA_PNF19 --fastqs=/imppc/labs/eslab/share/PerePericot/Projects/scRNA_PNF/data --sample=PNF19 --transcriptome=/imppc/labs/eslab/share/PerePericot/Projects/testing/yard/run_cellranger_count/refdata-cellranger-GRCh38-3.0.0

echo "`date` All done"
 
