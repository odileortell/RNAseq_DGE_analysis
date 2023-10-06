#!/bin/sh
#PBS -N trim_tandem
#PBS -l select=1:ncpus=24,walltime=96:00:00
#PBS -q smp
#PBS -e trim_tandem.err
#PBS -o trim_tandem.out
#PBS -P CBBI0999

cd $PBS_O_WORKDIR
MOD=/mnt/lustre/groups/CBBI0999/mods
module load gnu-parallel/20200322
module load $MOD/perl/5.34.0
module load trimmomatic/0.36
module load java/9.0.1


THR=24
p='/mnt/lustre/groups/HEAL1360/TANDEM_rawdata/rawdata'
P='~/lustre/msc/scripts'

minlen=15
t5p=3
t3p=3

function trim {
 F=$1;
 R=`echo $1 | sed -r 's/_R1/_R2/'`

 Fp=`echo $F | sed -r 's/.fq.gz/.P.fq.gz/'`
 Rp=`echo $R | sed -r 's/.fq.gz/.P.fq.gz/'`

 Fu=`echo $F | sed -r 's/.fq.gz/.U.fq.gz/'`
 Ru=`echo $R | sed -r 's/.fq.gz/.U.fq.gz/'`

 trimmomatic PE -threads $2 -phred33 $F $R $Fp $Fu $Rp $Ru \
                ILLUMINACLIP:adapters.fa:2:30:10 LEADING:$3 TRAILING:$4 SLIDINGWINDOW:4:15 MINLEN:$5
}
export -f trim

cd $p
ls *_R1* > list.lst && chmod +x list.lst
cat list.lst | parallel -j 24 -I% --max-args 1 "trim % 1 ${t5p} ${t3p} ${minlen}"

exit
