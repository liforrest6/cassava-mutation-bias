#!/bin/bash -l

#SBATCH -D /home/fli21/
#SBATCH -o /home/fli21/slurm-log/output.txt
#SBATCH -e /home/fli21/slurm-log/error.txt
#SBATCH -J slim-dm05
#SBATCH -t 24:00:00
#SBATCH --mem 50G
#SBATCH -n 48
#SBATCH --account=gmonroegrp
#SBATCH --partition=bmm
#SBATCH --mail-type=END
#SBATCH --mail-user=frrli@ucdavis.edu

module load SLiM/3.5
module load vcftools/0.1.13

for dm in 0.5 1
do
	for gds in 0 0.3
	do
		for ids in 0 0.1
		do
			for sf in 0.3 0.98
			do
				for i in {1..50}
				do 
					# slim -d dm=$dm \
					# -d gds=$gds \
					# -d ids=$ids \
					# -d i=$i \
					# -d sf=$sf \
					# -d m=1e-6 \
					# -d n=0.5 \
					# /home/fli21/cassava-mutation-bias/gene_intergene.slim

					vcftools --gzvcf /home/fli21/cassava-mutation-bias/vcf/cassava_dm${dm}_gds${gds}_ids${ids}_sf${sf}_i${i}.vcf.gz \
					--out /home/fli21/cassava-mutation-bias/tajima/cassava_dm${dm}_gds${gds}_ids${ids}_sf${sf}_i${i} --TajimaD 10

					# gzip /home/fli21/cassava-mutation-bias/vcf/cassava_dm${dm}_gds${gds}_ids${ids}_sf${sf}_i${i}.vcf
				 done
			done
		done
	done
done



