import numpy as np
from scipy.spatial import ConvexHull
from scipy import sparse
import hictoolbox


def get_color(in_file):
    with open(in_file) as f:
        color = f.readlines()
    color = np.array([int(float(i.strip())) for i in color])
    color = color[color > -1]
    return color


hic_path = 'HiC/GM12878/100kb_resolution_intrachromosomal/chr1_100kb.RAWobserved'
comps = ['results_separate/GM12878/1.txt',
         '../Team-SB2/Results/GM12878/1/100/Labels_epi_2compartments.txt',
         '../Team-SB2/Results/GM12878/1/100/Labels_contacts_7compartments.txt']
names = ['our', 'epi', 'contacts']

alpha = 0.227

R = 100000
A = np.loadtxt(hic_path)
A = np.int_(A)
A = np.concatenate((A, np.transpose(np.array([A[:, 1], A[:, 0], A[:, 2]]))), axis=0)
A = sparse.coo_matrix((A[:, 2], (A[:, 0], A[:, 1])))
A = hictoolbox.bin2d(A, R, R)
A, ind = hictoolbox.filteramat(A)  # [0] - mtx , [1] -vector with index
A = hictoolbox.SCN(A.copy())

colors = []
for comp in comps:
    colors.append(get_color(comp))

# 3D
print('3D')
contact_map = np.asarray(A) ** alpha
dist_matrix = hictoolbox.fastFloyd(1 / contact_map)  # shortest path on the matrix
dist_matrix = dist_matrix - np.diag(np.diag(dist_matrix))  # remove the diagonal
dist_matrix = (dist_matrix + np.transpose(dist_matrix)) / 2  # just to be sure that the matrix is symetric, not really usefull in theory

XYZ, E = hictoolbox.sammon(dist_matrix, 3)

# point rescale
hull = ConvexHull(XYZ)
scale = 100 / hull.area ** (1 / 3)
XYZ = XYZ * scale

for name, color in zip(names, colors):
    hictoolbox.writePDB('3Dcolors_' + name + '.pdb', XYZ, color)
