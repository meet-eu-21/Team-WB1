import numpy as np
import hictoolbox
import argparse
import os

coverages = [0.7, 0.5, 0.1, 3.6, 11.6, 0.4, 2.8, 0.2, 2.6, 0.1, 0.1, 0.1, 1.2, 8.3, 67.8]
coverages = np.array([c/100 for c in coverages])


def compute_epi_heights(data, epi, compartments_nr):
    comps = hictoolbox.makecompartimentbyGaussianHMM(data, N=compartments_nr)[1]
    epi_heights = []
    for k in range(compartments_nr):
        epi_k = np.mean(epi[comps==k], axis=0)/coverages   # vectors of normalised markers' values  
        epi_heights.append(epi_k)
    return epi_heights


def score(epi_heights, score_function):
    values = []
    for i in range(len(epi_heights)-1): 
        for j in range(i+1, len(epi_heights)):
            differences = abs(epi_heights[i] - epi_heights[j])
            values.append(score_function(differences))
    return sum(values)/len(values)


def file_exists(path):
    if not os.path.isfile(path):
        raise argparse.ArgumentTypeError(f"Directory {path} does not exist.")
    return path


def dir_exists(path):
    if not os.path.isdir(path):
        raise argparse.ArgumentTypeError(f"Directory {path} does not exist.")
    return path
