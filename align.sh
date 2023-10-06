#!/bin/bash

#PBS -N align
#PBS -l select=1:ncpus=24:mpiprocs=4,walltime=48:00:00
#PBS -q smp
#PBS -e align.err
#PBS -o align.out
#PBS -P CBBI0999

cd $PBS_O_WORKDIR

echo $HOSTNAME

path='/mnt/lustre/users/oortell/msc/analyses/alignment/bam'

module load samtools

for f in $path/*.bam
do
    samtools index $f
done