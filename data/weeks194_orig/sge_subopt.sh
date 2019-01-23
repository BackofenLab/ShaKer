#!/bin/bash
#$ -cwd
#$ -l h_vmem=8G
#$ -M mautner@cs.uni-freiburg.de
#$ -m as
##$ -pe smp 1
#$ -R y
#$ -o /home/mautner/JOBZ/asd/$JOB_ID.o_$TASK_ID
#$ -e /home/mautner/JOBZ/asd/$JOB_ID.e_$TASK_ID

##mkdir -p /home/mautner/JOBZ/reconstr_o/
##mkdir -p /home/mautner/JOBZ/reconstr_e/
python /scratch/bi01/mautner/ShaKer/ShaKer/data/weeks194_orig/calcsubopt.py $SGE_TASK_ID
#qsub -V -t 1-194 sge_subopt.sh

