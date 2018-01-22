#!/bin/sh
#!/scratch/bi01/mautner/miniconda2/bin/fish


#$ -cwd
##$ -l h_vmem=10G
#$ -pe smp 24
#$ -V    
#$ -R y
#$ -M mautner@cs.uni-freiburg.de
#$ -m a 
#$ -m s 
##$ -o out_stdout
##$ -e out_stderr

# layercomp.py  TASKFILE 
python optimize_forest.py
