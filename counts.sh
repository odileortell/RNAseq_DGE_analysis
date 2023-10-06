#!/bin/sh


path='/mnt/lustre/users/oortell/msc/analyses/tandem_data'
group='/mnt/lustre/groups/CBBI0999/odile/alignment'

mkdir -p $path/counts/tmp

cat $path/sample_names_edit.lst | while read sample; do
    echo ${sample}
    cat $group/star/${sample}.ReadsPerGene.out.tab | tail -n +5 | cut -f4 > $path/counts/tmp/${sample}.count
done



#for sample in 'cat samples_edited.lst'; do \
#    echo ${sample}
#    cat star2/${sample}.fq.gz.ReadsPerGene.out.tab | tail -n +5 | cut -f4 > counts/tmp/${sample}.count
#done
