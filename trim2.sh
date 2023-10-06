#!/bin/sh
#PBS -N trim_tandem
#PBS -l select=1:ncpus=24,walltime=3:00:00
#PBS -q smp
#PBS -e trim_tandem.err
#PBS -o trim_tandem.out
#PBS -P CBBI0999

cd $PBS_O_WORKDIR
module load trimmomatic/0.36
module load java/17.0.2

p='/mnt/lustre/users/oortell/msc'
# fq=${p}'/seq_data'
adapt=${p}'/scripts/adapters.fa'
data='/mnt/lustre/groups/HEAL1360/TANDEM_rawdata/rawdata'

THR=24
minlen=15          # seq less than that will be removed
trim5p=3           # nucleotides trimmed on each end
trim3p=3


###############
### PROGRAM
###############

reads=${data}/*.fq.gz

for f in $(<${data}); do
  trimmomatic PE -threads 24 -phred33 $reads \
    ILLUMINACLIP:${adapt}:2:30:10 SLIDINGWINDOW:4:20 LEADING:3 TRAILING:4 MINLEN:5
done

exit