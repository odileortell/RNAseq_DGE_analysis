#!/bin/bash

#PBS -N tstar_map
#PBS -l select=1:ncpus=24
#PBS -l walltime=96:00:00
#PBS -q smp
#PBS -e tstar_map.err
#PBS -o tstar_map.out
#PBS -P CBBI0999


cd $PBS_O_WORKDIR

echo $HOSTNAME

module load chpc/BIOMODULES
module load STAR/2.7.7a

data='/mnt/lustre/users/oortell/msc/analyses/tandem_data/rawdata'
path='/mnt/lustre/users/oortell/msc/analyses/tandem_alignment3'
staridx=$path/tandem_references/star.overlap100.gencode.40

rm -R star
mkdir -p star
cd star

# adjusting memory per thread (in kb)
M=$((120*1000000000))

#STAR \
#--runThreadN 3 \
#    --genomeDir $path/References/star.overlap100.gencode.40 \
#    --outSAMtype BAM SortedByCoordinate \
#    --quantMode GeneCounts \
#    --outFileNamePrefix $path/ \
#    --readFilesCommand zcat \
#    --readFilesIn $data/*.fq.gz


#########################################
### Map each individual to ref genome ###
#########################################

for i in ${data[@]}; do  ### run through each file in folder
    echo $i
    ulimit -v $M
    STAR --genomeDir $staridx --runThreadN 24 --outFileNamePrefix $i \
    --readFilesIn $data/*.fq --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts \
    --limitGenomeGenerateRAM $M
done