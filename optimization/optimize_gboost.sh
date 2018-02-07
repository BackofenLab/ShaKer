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
##$ -o out_grad_std
##$ -e out_grad_err

# layercomp.py  TASKFILE 
python optimize_gradient_boosting.py
