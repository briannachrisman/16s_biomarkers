#!/bin/bash

cd $MY_HOME/sequence_based_biomarkers


for id in 90 95 97 99 ; do
    for datatype in obesity autism ; do
        echo $datatype $id

        /oak/stanford/groups/dpwall/computeEnvironments/usearch -cluster_smallmem data/${datatype}/seqs.fa -id 0.$id -centroids $SCRATCH/tmp/${datatype}_otus$id.fa --sortedby other

            /oak/stanford/groups/dpwall/computeEnvironments/usearch -otutab data/${datatype}/seqs.fa -otus $SCRATCH/tmp/${datatype}_otus$id.fa -mapout results/generate_biomarkers/otu_asv_map_${id}_${datatype}.txt  -otutabout $SCRATCH/tmp/otutab.txt  -id 0.$id


    done
done