from hictoolbox import SCN, observed_expected, filteramat, bin2d
from utils import dir_exists
import os
import numpy as np
from scipy import sparse
import argparse


def process_matrix(mfile, out_mtx, file_ind):
    R = 100000
    A = np.loadtxt(mfile)
    A = np.int_(A)
    A = np.concatenate((A, np.transpose(np.array([A[:, 1], A[:, 0], A[:, 2]]))), axis=0)
    A = sparse.coo_matrix((A[:, 2], (A[:, 0], A[:, 1])))
    A = bin2d(A, R, R)
    A, ind = filteramat(A)  # [0] - mtx , [1] -vector with index
    A = SCN(A.copy())
    A = observed_expected(A)
    A = np.corrcoef(A)

    # save processed matrix
    np.save(out_mtx, A)
    # save indices
    np.save(file_ind, ind)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_dir", type=dir_exists, help="HiC directory (containing cell line directories)",
                        required=True)
    parser.add_argument("-o", "--out_dir", type=str, default="preprocessed",
                        help="Directory for saving the preprocessed data")
    args = parser.parse_args()

    HiC_dir = args.input_dir
    out_dir = args.out_dir

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    # cell line GM12878
    sample = 'GM12878'
    out_path = os.path.join(out_dir, sample)
    if not os.path.isdir(out_path):
        os.makedirs(out_path)

    gm_path = os.path.join(HiC_dir, sample, '100kb_resolution_intrachromosomal')
    for file in os.listdir(gm_path):
        if file.endswith('RAWobserved'):
            file_path = os.path.join(gm_path, file)
            fileout = os.path.join(out_path, file.replace(".RAWobserved", "_processed"))
            fileout_ind = os.path.join(out_path, file.replace(".RAWobserved", "_indices"))
            print('Processing file', file, "\n Results will be saved at:", fileout_ind)
            process_matrix(file_path, fileout, fileout_ind)

    # other cell lines
    for sample in ['HMEC', 'HUVEC', 'IMR90', 'NHEK']:
        inpath = os.path.join(HiC_dir, sample, '100kb_resolution_intrachromosomal')
        out_path = os.path.join(out_dir, sample)
        if not os.path.isdir(out_path):
            os.makedirs(out_path)

        for chromosome in os.listdir(inpath):
            if chromosome.startswith('chr'):
                # now we are iterating through chr dirs
                chromosome = os.path.join(inpath, chromosome, 'MAPQGE30')
                for file in os.listdir(chromosome):
                    if file.endswith('RAWobserved'):
                        file_path = os.path.join(chromosome, file)
                        fileout = os.path.join(out_path, file.replace(".RAWobserved", "_processed"))
                        fileout_ind = os.path.join(out_path, file.replace(".RAWobserved", "_indices"))
                        print('Processing file', file, "\n Results will be saved at:", fileout_ind)
                        process_matrix(file_path, fileout, fileout_ind)
