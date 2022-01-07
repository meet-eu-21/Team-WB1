#!/bin/bash

python3.6 preprocess_mtx.py -i HiC -o preprocessed

python3.6 preprocess_epi.py -e E116_15_coreMarks_stateno.bed -s hg19.chrom.sizes -o Epi_files

python3.6 calc_comps.py -i preprocessed -c centromericpos_hg19.txt -e Epi_files -o results_full \
-o results_full -p results_full/plots --cmp_consensus --cmp_min 2 --cmp_max 15 -a 0.7

python3.6 calc_comps.py -i preprocessed -c centromericpos_hg19.txt -e Epi_files -o results_full \
-o results_separate -p results_separate/plots --cmp_min 2 --cmp_max 15 -a 0.7
