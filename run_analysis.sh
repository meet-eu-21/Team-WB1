#!/bin/bash

python3.6 preprocess_mtx.py -i /mnt/chr12/Data/ZIP/WB1/dataforstudent/HiC -o /mnt/chr12/Data/ZIP/WB1/preprocessed

python3.6 preprocess_epi.py -e ../E116_15_coreMarks_stateno.bed -s ../hg19.chrom.sizes -o Epi_files

python3.6 calc_comps.py -i /mnt/chr12/Data/ZIP/WB1/preprocessed -c ../centromericpos_hg19.txt -e Epi_files -o results_full \
-p results_full/plots --cmp_consensus --cmp_min 2 --cmp_max 15 -a 0.7
