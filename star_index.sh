#!/bin/bash
#PBS -N tandem_index
#PBS -l select=1:ncpus=24,walltime=03:00:00
#PBS -q test
#PBS -e tandem_index.err
#PBS -o tandem_index.out
#PBS -P CBBI0999

cd $PBS_O_WORKDIR

echo $HOSTNAME

wrk_dir="/mnt/lustre/users/oortell/msc"
outpath="tandem_references"
mkdir -p $outpath

ref="$wrk_dir/seq_data/references"

cd $ref

# wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh38.primary_assembly.genome.fa.gz
# gunzip GRCh38.primary_assembly.genome.fa.gz
FASTA="GRCh38.primary_assembly.genome.fa"

# wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.basic.annotation.gtf.gz
# gunzip gencode.v40.basic.annotation.gtf.gz
GTF="gencode.v40.basic.annotation.gtf"

cd $wrk_dir/$outpath
mkdir -p star.overlap100.gencode.40
cd star.overlap100.gencode.40

module load chpc/BIOMODULES
module load STAR/2.7.7a

STAR --runThreadN 12 --runMode genomeGenerate --genomeDir . --genomeFastaFiles $ref/$FASTA --sjdbGTFfile $ref/$GTF --sjdbOverhang 100
