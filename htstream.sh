#!/bin/bash

#PBS -N htstream
#PBS -l select=1:ncpus=24:mpiprocs=4
#PBS -P CBBI0999
#PBS -l walltime=48:00:00
#PBS -q smp
#PBS -o htstream2.out
#PBS -e htstream2.err

cd $PBS_O_WORKDIR

module load chpc/BIOMODULES
module load HTStream/1.3.1

wrk_dir="/mnt/lustre/users/oortell/msc"


### bash array with samples
samples=($(ls -d $wrk_dir/analyses/tandem_data/*fq.gz | awk -F 'tandem_data/' '{print $2}'))

inpath="tandem_data"
outpath="hts_preproc"
### create directory if it does not exist
[[ -d $outpath ]] || mkdir -p $wrk_dir/analyses/tandem_alignment2/$outpath

### first element from array
init_idx=0
end_idx=$((${#samples[@]} -1))

for i in $(eval echo "{$init_idx..$end_idx}")
do
target_sample=${samples[i]}
target_sample=$(echo $target_sample | awk -F '.fq.gz' '{print $1}')
### run the pipe for each sample
hts_Stats -L $wrk_dir/analyses/tandem_alignment2/$outpath/$target_sample'.json' -N 'initial stats' -U \
$wrk_dir/analyses/$inpath/$target_sample'.fq.gz' | \
hts_SeqScreener -A $wrk_dir/analyses/tandem_alignment2/$outpath/$target_sample'.json' -N 'screen phix' | \
hts_SeqScreener -A $wrk_dir/analyses/tandem_alignment2/$outpath/$target_sample'.json' -N 'count the number of rRNA reads' -r \
-s $wrk_dir/seq_data/references/human_rrna.fasta | \
hts_AdapterTrimmer -A $wrk_dir/analyses/tandem_alignment2/$outpath/$target_sample'.json' -N 'trim adapters' | \
hts_PolyATTrim -A $wrk_dir/analyses/tandem_alignment2/$outpath/$target_sample'.json' -N 'remove polyAT tails' | \
hts_NTrimmer -A $wrk_dir/analyses/tandem_alignment2/$outpath/$target_sample'.json' -N 'remove any remaining N characters' | \
hts_QWindowTrim -A $wrk_dir/analyses/tandem_alignment2/$outpath/$target_sample'.json' -N 'quality trim the ends of reads' | \
hts_LengthFilter -A $wrk_dir/analyses/tandem_alignment2/$outpath/$target_sample'.json' -N 'remove reads < 45bp' -n -m 45 | \
hts_Stats -A $wrk_dir/analyses/tandem_alignment2/$outpath/$target_sample'.json' -N 'final stats' -f \
$wrk_dir/analyses/tandem_alignment2/$outpath/$target_sample

echo "Finished: " $target_sample "at: " date
date
done


