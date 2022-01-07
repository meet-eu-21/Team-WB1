# Team-WB1

Meet-EU Team WB1

Topic B : Chromosome compartments
## Introduction
We tried to use HMM for detecting compartments and noticed that when we use division to three compartments, two of the detected compartments have very simmilar epigenetic profiles (which suggest that they shouldn't have beed devided to separate compartments) but lay on different chromosome arms. Then we performed Principal Component Analysis and noticed that often there is one component that corresponds to chromosome arms (it is positive on one arm and negative on the other arm). To remove the arm bias we tried running HMM on 6 first principal components, except the component which corresponds the most to the chromosome arms (so 5 out of 6 components). 

## Basic pipeline
You should first download the data from http://www.lcqb.upmc.fr/meetu/dataforstudent/  
Then preprocess HiC matrices and the file with epigenetic marks using ```preprocess_mtx.py``` and ```preprocess_epi.py```.
Finally you can run the analysis with ```calc_comps.py```.

Check out bash script ```run_analysis.sh``` that we used for obtaining the results shared on github.

## calc_comps.py
This program detects the chromatine compartments. It builds HMM models for different number of compartments and chooses best scoring number. It can either choose a different number of compartments for each chromosome (default mode) or a consensus number that has the highest overall score (```--cmp_consensus``` option).  
For each chromosome the program checks if there is a principal component that corresponds well enough to the chromosome arms (parameter ```min_arm_score```). If there is, then it calculates compartments before and after removing the arm component. Then it chooses the setup that results with the highest score.

```
  -i INPUT_DIR, --input_dir INPUT_DIR
                        Directory that contains the preprocessed HiC data
  -c CENTR_FILE, --centr_file CENTR_FILE
                        file with centromeric positions
  -e EPI_DIR, --epi_dir EPI_DIR
                        Directory that contains preprocessed data with
                        epigenetic marks
  -o OUT_DIR, --out_dir OUT_DIR
                        Directory for saving the results (compartments)
  -p PLOTS_DIR, --plots_dir PLOTS_DIR
                        Directory for saving the plots
  --cmp_consensus       Find best number of compartments for all chromosomes
                        (default: find best number for each chromosome
                        separately)
  --cmp_min CMP_MIN     Minimum number of compartments to consider
  --cmp_max CMP_MAX     Maximum number of compartments to consider
  -a MIN_ARM_SCORE, --min_arm_score MIN_ARM_SCORE
                        If arm_score of a principal component is higher than
                        min_arm_score than it can be removed (between 0 and 1,
                        recommended values: 0.5-0.75
```
