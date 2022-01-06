import numpy as np
from scipy import sparse
import pandas as pd
import hictoolbox
import os
import argparse
from utils import file_exists


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-e", "--epi_file", type=file_exists, default='E116_15_coreMarks_stateno.bed', 
                        help="bed file with epigenetic marks", required=True)
    parser.add_argument("-s", "--sizes_file", type=file_exists, default='hg19.chrom.sizes', 
                        help="file with chromosome sizes", required=True)
    parser.add_argument("-o", "--out_dir", type=str, default="Epi_files", 
                        help="Directory for saving the preprocessed data")
    args = parser.parse_args()
    
    out_dir = args.out_dir

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)


    R = 100000
    EpiGfilename = '../E116_15_coreMarks_stateno.bed'

    colors = pd.read_csv(args.epi_file,delimiter='\t',header=None,names=[1,2,3,4], dtype={2:int, 3:int, 4:int}, skipfooter=1)
    lengths = pd.read_csv(args.sizes_file, sep='\t', header=None, index_col=0)

    for ch in sorted(list(set(colors[1]))):
        if ch=='chrM':
            continue
        print(ch, 'in progress...')
        print('Result will be saved to:', os.path.join(args.out_dir, 'Epi_'+ch+'.npy'))
        color = colors[colors[1]==ch]
        number = color[4].max() # number of color in the file
        length = int(lengths.loc[ch])
        color_vec = np.zeros((length,number+1)) # build array at pb resolution LENchr * number of color
        i = 0
        while i<np.shape(color)[0]:
            color_vec[color[2].iloc[i]:color[3].iloc[i],color[4].iloc[i]]=1
            i+=1
        color_bins = hictoolbox.bin2d(color_vec,R,1)
        del color_vec

        color_bins = color_bins / np.amax(color_bins)

        color_bins = np.float64(color_bins.todense())
        np.save(os.path.join(out_dir, 'Epi_'+ch+'.npy'), color_bins)
