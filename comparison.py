import numpy as np
import os
import matplotlib.pyplot as plt


def epi_score(comps, epi_markers):
    coverages = [0.7, 0.5, 0.1, 3.6, 11.6, 0.4, 2.8, 0.2, 2.6, 0.1, 0.1, 0.1, 1.2, 8.3, 67.8]
    coverages = np.array([c / 100 for c in coverages])

    binsaved = comps > -1
    comps = comps[binsaved]
    epi_markers = epi_markers[binsaved]
    compartments_nr = int(max(comps) + 1)
    epi_heights = []
    for k in range(compartments_nr):
        epi_k = np.mean(epi_markers[comps == k], axis=0) / coverages  # vectors of normalised markers' values  
        epi_heights.append(epi_k)

    values = []
    for i in range(len(epi_heights) - 1):
        for j in range(i + 1, len(epi_heights)):
            differences = abs(epi_heights[i] - epi_heights[j])
            values.append(np.mean(differences))
    return sum(values) / len(values)


scores1 = []
scores2 = []
scores3 = []
ncomps1 = []
ncomps2 = []
ncomps3 = []

chromosomes = [str(i) for i in range(1, 23)] + ['X']

for ch in chromosomes:
    path = 'results_separate/GM12878/' + ch + '.txt'
    with open(path) as f:
        comps1 = f.readlines()
    comps1 = np.array([float(i.strip()) for i in comps1])
    ncomps1.append(max(comps1) + 1)
    epi = np.load('Epi_files/Epi_chr' + ch + '.npy')
    scores1.append(epi_score(comps1, epi))

    path = '../Team-SB2/Results/GM12878/' + ch + '/100/'
    files = os.listdir(path)
    labels_files = [f for f in files if f.startswith('Labels_epi')]
    labels_file = sorted(labels_files)[-1]
    with open(path + labels_file) as f:
        comps2 = f.readlines()
    comps2 = [float(i.strip()) for i in comps2] + [-1]
    ncomps2.append(max(comps2) + 1)
    l = len(epi) - len(comps2)
    comps2 = np.array(comps2 + [-1] * l)
    scores2.append(epi_score(comps2, epi))

    labels_files = [f for f in files if f.startswith('Labels_contacts')]
    labels_file = sorted(labels_files)[-1]
    with open(path + labels_file) as f:
        comps3 = f.readlines()
    comps3 = [float(i.strip()) for i in comps3] + [-1]
    ncomps3.append(max(comps3) + 1)
    comps3 = np.array(comps3 + [-1] * l)
    scores3.append(epi_score(comps3, epi))


plt.figure(figsize=(8, 6))
plt.plot(scores1, label="Our model")
plt.plot(scores2, label="SB2 epi")
plt.plot(scores3, label="SB2 contacts")
plt.title("Comparison of prediction scores")
plt.xlabel("Chromosome")
plt.ylabel("Epi score of best compartmentalization")
plt.xticks(range(len(chromosomes)), chromosomes)
plt.legend()
plt.grid(axis='y')
plt.savefig('comparison.png', dpi=100, facecolor='white', bbox_inches='tight')
plt.close()

plt.figure(figsize=(8, 6))
plt.plot(ncomps1, label="Our model")
plt.plot(ncomps2, label="SB2 epi")
plt.plot(ncomps3, label="SB2 contacts")
plt.title("Comparison of compartment number")
plt.xlabel("Chromosome")
plt.ylabel("Inferred number of compartments")
plt.xticks(range(len(chromosomes)), chromosomes)
plt.legend()
plt.grid(axis='y')
plt.savefig('comparison_ncomps.png', dpi=100, facecolor='white', bbox_inches='tight')
plt.close()
