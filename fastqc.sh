#!/bin/sh

#PBS -N gc6_fastqc
#PBS -l select=1:ncpus=24,walltime=72:00:00
#PBS -q smp
#PBS -e gc6_fastqc.err
#PBS -o gc6_fastqc.out
#PBS -P CBBI0999

cd $PBS_O_WORKDIR

module load java/17.0.2
module load FastQC/0.11.9

for A in /mnt/lustre/groups/CBBI0999/odile/1_GC6_corrected/*.fq.gz
do
    fastqc -o ~/lustre/msc/analyses/gc6_fastqc/ $A
done
echo 'DONE'

exit
