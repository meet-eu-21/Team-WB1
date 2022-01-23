import numpy as np
import pandas as pd
import hictoolbox
import matplotlib.pyplot as plt
from sklearn import decomposition
from utils import compute_epi_heights, score, file_exists, dir_exists
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
    parser.add_argument("-o", "--out_dir", type=str, default="results", 
                        help="Directory for saving the results (compartments)")
    parser.add_argument("-p", "--plots_dir", type=str, default="results/plots", 
                        help="Directory for saving the plots")
    parser.add_argument("--cmp_consensus", action="store_true", 
                        help="Find best number of compartments for all chromosomes (default: find best number for each chromosome separately)")
    parser.add_argument("--cmp_min", type=int, default=2, 
                        help="Minimum number of compartments to consider")
    parser.add_argument("--cmp_max", type=int, default=15, 
                        help="Maximum number of compartments to consider")
    parser.add_argument("-a", "--min_arm_score", type=float, default=0.7, 
                        help="If arm_score of a principal component is higher than min_arm_score than it can be removed \
                        (between 0 and 1, recommended values: 0.6-0.9")

    args = parser.parse_args()

    out_dir = args.out_dir
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    plots_dir = args.plots_dir
    if not os.path.isdir(plots_dir):
        os.makedirs(plots_dir)

    input_dir = args.input_dir

    cmp_min, cmp_max = args.cmp_min, args.cmp_max
    cmp_max = cmp_max + 1
    assert cmp_min > 1
    assert cmp_max > cmp_min

    R = 100000
    chromosomes = [str(i) for i in range(1, 23)] + ['X']

    centr = pd.read_csv(args.centr_file, sep='\t', names=['chr', 'start', 'end', 'q', 'acen'])
    centr['start'] = centr['start']/R
    centr['end'] = centr['end']/R
    centr.index = centr['chr'].str[3:]

    for cell_line in ['GM12878', 'HMEC', 'HUVEC', 'IMR90', 'NHEK']:
        print(cell_line)
        scores1 = np.zeros(cmp_max-cmp_min)
        scores2 = np.zeros(cmp_max-cmp_min)
        scores_before = np.zeros(cmp_max-cmp_min)
        used_trimmed = np.zeros((len(chromosomes), cmp_max-cmp_min))
        best_pcs = np.zeros(len(chromosomes)) - 1
        cmp_sizes = range(cmp_min, cmp_max)
        
        in_path = os.path.join(input_dir, cell_line)
        out_path = os.path.join(out_dir, cell_line)
        plots_path = os.path.join(plots_dir, cell_line)

        if not os.path.isdir(out_path):
            os.makedirs(out_path)
        if not os.path.isdir(plots_path):
            os.makedirs(plots_path)

        lengths = []
        for c in range(len(chromosomes)):
            ch = chromosomes[c]
            print(ch)
            D = np.load(os.path.join(in_path, 'chr'+ch+'_100kb_processed.npy'))
            binsaved = np.load(os.path.join(in_path, 'chr'+ch+'_100kb_indices.npy'))

            epi = np.load(os.path.join(args.epi_dir, 'Epi_chr'+ch+'.npy'))
            length = len(epi)
            lengths.append(length)
            epi = epi[binsaved]

            pca = decomposition.PCA(n_components=6)
            res = pca.fit_transform(D)

            best_arm_score = 0
            best_pc = None

            k = 0
            for i in range(2):
                for j in range(3):
                    s = centr.loc[ch, 'start']
                    e = centr.loc[ch, 'end']

                    id_centr = len(binsaved[binsaved<s])
                    arm1_frac = id_centr / len(binsaved)
                    if arm1_frac<0.1 or arm1_frac>0.9:  # one of the arms is too short for the analysis
                        arm_score = 0
                    else:
                        arm_scores1 = [np.mean(res[:id_centr, k]>0), np.mean(res[id_centr:, k]<0)]
                        arm_scores2 = [np.mean(res[:id_centr, k]<0), np.mean(res[id_centr:, k]>0)]
                        arm_score = np.abs(max(np.mean(arm_scores1), np.mean(arm_scores2)))

                    if arm_score > best_arm_score:
                        best_arm_score = arm_score
                        best_pc = k
                        best_pcs[c] = best_pc
                    k+=1


            y_pc_epi_heights = [compute_epi_heights(res, epi, i) for i in cmp_sizes]
            
            y_pc = [score(heights, np.mean) for heights in y_pc_epi_heights]
            plt.plot(cmp_sizes, y_pc, label="First 6 components")
            best_max = np.max(y_pc)
            best_argmax = np.argmax(y_pc)
            used = 'y_pc'

            if best_arm_score > args.min_arm_score:
                res_trimmed = np.delete(res, best_pc, 1)  
                y_pc_trimmed_epi_heights = [compute_epi_heights(res_trimmed, epi, i) for i in cmp_sizes]
                y_pc_trimmed = [score(heights, np.mean) for heights in y_pc_trimmed_epi_heights]
                plt.plot(cmp_sizes, y_pc_trimmed, label="First 6 components except the arm component")
                if best_max < max(y_pc_trimmed):
                    best_max = max(y_pc_trimmed)
                    best_argmax = np.argmax(y_pc_trimmed)
                    used = 'y_pc_trimmed'
                print(f'best_arm_score > {args.min_arm_score}')
            else:
                y_pc_trimmed = y_pc

            plt.xticks(range(cmp_min, cmp_max))
            plt.title("Compartmentalization score")
            plt.xlabel("Number of compartments")
            plt.ylabel("Mean score")
            plt.legend(title="Data for HMM")
            plt.savefig(os.path.join(plots_path, ch+'.png'), dpi=150, facecolor='white', bbox_inches='tight')
            plt.close()


            if args.cmp_consensus:
                scores1 += np.maximum(y_pc, y_pc_trimmed)
                scores2 += y_pc_trimmed
                scores_before += y_pc
                used_trimmed[c] = np.array(y_pc) < np.array(y_pc_trimmed)
            else:   
                print('used:', used)
                n_comps = best_argmax+min(cmp_sizes)
                print('best n_comps:', n_comps)

                if used == 'y_pc':
                    data = res

                else:
                    data = res_trimmed
                comps = hictoolbox.makecompartimentbyGaussianHMM(data, N=n_comps)[1]

                result = np.zeros(length) - 1
                result[binsaved] = comps
                np.savetxt(os.path.join(out_path, str(ch)+'.txt'), result, delimiter='\n', fmt='%1.1f')

        if args.cmp_consensus:
            scores1 = scores1 / len(chromosomes)
            scores2 = scores2 / len(chromosomes)
            scores_before = scores_before / len(chromosomes)

            plt.plot(cmp_sizes, scores_before, label="Before")
            plt.plot(cmp_sizes, scores1, label="After")
            plt.title("Compartmentalization score before and after \n removing the arm component (only if it improves the score)")
            plt.xlabel("Number of compartments")
            plt.ylabel("Average score for all chromosomes")
            plt.xticks(range(cmp_min, cmp_max))
            plt.legend()
            plt.savefig(os.path.join(plots_path, 'n_comps_max.png'), dpi=150, facecolor='white', bbox_inches='tight')
            plt.close()

            plt.plot(cmp_sizes, scores_before, label="Before")
            plt.plot(cmp_sizes, scores2, label="After")
            plt.title("Compartmentalization score before and after \n removing the arm component (whenever possible)")
            plt.xlabel("Number of compartments")
            plt.ylabel("Average score for all chromosomes")
            plt.xticks(range(cmp_min, cmp_max))
            plt.legend()
            plt.savefig(os.path.join(plots_path, 'n_comps_trimmed.png'), dpi=150, facecolor='white', bbox_inches='tight')
            plt.close()

            # choosing best number of compartments for all chromosomes
            n_comps = np.argmax(scores1)
            used_trimmed = used_trimmed[:, n_comps].flatten()
            n_comps = n_comps + min(cmp_sizes)
            print('Best number of compartments for all chromosomes: ', n_comps)


            for c in range(len(chromosomes)):   # infer the compartments once more (only for the best number)
                ch = chromosomes[c]
                D = np.load(os.path.join(in_path, 'chr'+ch+'_100kb_processed.npy'))
                binsaved = np.load(os.path.join(in_path, 'chr'+ch+'_100kb_indices.npy'))

                length = lengths[c]

                pca = decomposition.PCA(n_components=6)
                res = pca.fit_transform(D)
                if used_trimmed[c]:
                    res_trimmed = np.delete(res, best_pcs[c], 1)
                    data = res_trimmed
                else:
                    data = res

                comps = hictoolbox.makecompartimentbyGaussianHMM(data, N=n_comps)[1]

                result = np.zeros(length) - 1
                result[binsaved] = comps
                np.savetxt(os.path.join(out_path, str(ch)+'.txt'), result, delimiter='\n', fmt='%1.1f')
