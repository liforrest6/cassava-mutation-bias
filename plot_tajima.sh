#!/bin/bash -l

#SBATCH -D /home/fli21/
#SBATCH -o /home/fli21/slurm-log/output.txt
#SBATCH -e /home/fli21/slurm-log/error.txt
#SBATCH -J slim-dm05
#SBATCH -t 24:00:00
#SBATCH --mem 128G
#SBATCH -n 8
#SBATCH --account=gmonroegrp
#SBATCH --partition=bmm
#SBATCH --mail-type=END
#SBATCH --mail-user=frrli@ucdavis.edu

Rscript /home/fli21/cassava-mutation-bias/tss_tts.tajima.batch.R 