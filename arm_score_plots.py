import numpy as np
import pandas as pd
import hictoolbox
import matplotlib.pyplot as plt
from sklearn import decomposition
from utils import file_exists, dir_exists
import os
import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_dir", type=dir_exists, default="preprocessed",
                        help="Directory that contains the preprocessed data", required=True)
    parser.add_argument("-c", "--centr_file", type=file_exists, default='centromericpos_hg19.txt', 
                        help="file with centromeric positions", required=True)
    parser.add_argument("-e", "--epi_dir", type=dir_exists, default="Epi_files", 
                        help="Directory that contains preprocessed data with epigenetic marks", required=True)
    parser.add_argument("-p", "--plots_dir", type=str, default="pca_plots", 
                        help="Directory for saving the plots")

    args = parser.parse_args()

    plots_dir = args.plots_dir
    if not os.path.isdir(plots_dir):
        os.makedirs(plots_dir)

    input_dir = args.input_dir
    

    R = 100000
    chromosomes = [str(i) for i in range(1, 23)] + ['X']

    centr = pd.read_csv(args.centr_file, sep='\t', names=['chr', 'start', 'end', 'q', 'acen'])
    centr['start'] = centr['start']/R
    centr['end'] = centr['end']/R
    centr.index = centr['chr'].str[3:]

    for cell_line in ['GM12878', 'HMEC', 'HUVEC', 'IMR90', 'NHEK']:
        print(cell_line)
        
        in_path = os.path.join(input_dir, cell_line)
        plots_path = os.path.join(plots_dir, cell_line)

        if not os.path.isdir(plots_path):
            os.makedirs(plots_path)

        best_scores = np.zeros(len(chromosomes))
        best_pcs = np.zeros(len(chromosomes))

        for c in range(len(chromosomes)):
            ch = chromosomes[c]
            print(ch)
            D = np.load(os.path.join(in_path, 'chr'+ch+'_100kb_processed.npy'))
            binsaved = np.load(os.path.join(in_path, 'chr'+ch+'_100kb_indices.npy'))

            epi = np.load(os.path.join(args.epi_dir, 'Epi_chr'+ch+'.npy'))
            length = len(epi)
            epi = epi[binsaved]

            pca = decomposition.PCA(n_components=6)
            res = pca.fit_transform(D)

            best_arm_score = 0
            best_pc = None

            fig, axes = plt.subplots(2, 3, figsize=(15,8))
            k = 0
            for i in range(2):
                for j in range(3):
                    s = centr.loc[ch, 'start']
                    e = centr.loc[ch, 'end']
                    axes[i, j].axhline(0, c='gray', lw=1)
                    axes[i, j].axvspan(s, e, facecolor='r', alpha=0.2)
                    axes[i, j].scatter(binsaved, res[:, k], s=3)
                    id_centr = len(binsaved[binsaved<s])
                    arm1_frac = id_centr / len(binsaved)
                    if arm1_frac<0.1 or arm1_frac>0.9:  # one of the arms is too short for the analysis
                        arm_score = 0
                    else:
                        arm_scores1 = [np.mean(res[:id_centr, k]>0), np.mean(res[id_centr:, k]<0)]
                        arm_scores2 = [np.mean(res[:id_centr, k]<0), np.mean(res[id_centr:, k]>0)]
                        arm_score = np.abs(max(np.mean(arm_scores1), np.mean(arm_scores2)))
                    axes[i, j].set_title(f'PC{k+1}, arm_score={np.round(arm_score, 3)}')
                    if arm_score > best_arm_score:
                        best_arm_score = arm_score
                        best_pc = k
                    k+=1
            best_scores[c] = best_arm_score
            if best_arm_score > 0:
                best_pcs[c] = best_pc + 1

            plt.suptitle(f'First six principal components for chromosome {ch}', size=14)
            plt.savefig(os.path.join(plots_path, ch+'.png'), dpi=150, facecolor='white', bbox_inches='tight')
            plt.close()

        plt.plot(best_scores)
        plt.xticks(range(len(chromosomes)), chromosomes)
        plt.savefig(os.path.join(plots_path, 'best_scores.png'), dpi=150, facecolor='white', bbox_inches='tight')
        plt.xlabel('Chromosome')
        plt.ylabel('Highest arm score')
        plt.close()

        plt.plot(best_pcs)
        plt.xticks(range(len(chromosomes)), chromosomes)
        plt.yticks(range(7), ['None', 1, 2, 3, 4, 5, 6])
        plt.xlabel('Chromosome')
        plt.ylabel('PC with highest arm score')
        plt.savefig(os.path.join(plots_path, 'best_pcs.png'), dpi=150, facecolor='white', bbox_inches='tight')
        plt.close()
