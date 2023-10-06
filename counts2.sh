#!/bin/bash

group='/mnt/lustre/groups/CBBI0999/odile/alignment'
path='/mnt/lustre/users/oortell/msc/analyses/gc6_alignment'

cat $path/sample_names_edit.lst | while read sample; do
    echo ${sample}
    tail -n +5 $group/star/${sample}.ReadsPerGene.out.tab | cut -f1 > $path/counts/tmp/geneids.txt
    head counts/tmp/geneids.txt
done