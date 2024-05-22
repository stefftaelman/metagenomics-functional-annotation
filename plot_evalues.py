import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sb

# set input
#REAL_HITS_DIR = "./virulence_factors/VF-hits/"
REAL_HITS_DIR = "./AMR_genes/AMR-hits/"
#RAND_HITS_DIR = "./virulence_factors/shuffled-hits/"
RAND_HITS_DIR = "./AMR_genes/shuffled-hits/"
#OUT_FIGURE = "./virulence_factors/evalue_histogram_VFs.png"
OUT_FIGURE = "./AMR_genes/evalue_histogram_AMRs.png"


print("Reading in blast results on real and randomized reads...")
list_of_hits = []
for root, dirs, files in os.walk(REAL_HITS_DIR):
	for file in files:
         if 'DS_Store' not in file:
            list_of_hits.append(os.path.join(root,file))
list_of_sims = []
for root, dirs, files in os.walk(RAND_HITS_DIR):
	for file in files:
         if 'DS_Store' not in file:
            list_of_sims.append(os.path.join(root,file))


print("Plotting e-values & calculating median...")
num_files = len(list_of_hits)
for idx, name in enumerate(list_of_hits):
    evals = list(pd.read_csv(name, delimiter="\t", usecols=[10]).to_numpy().transpose()[0])
    alpha = 0.05
    sb.kdeplot(evals, bw_method=0.005, color='#1565A9', alpha=alpha, linewidth=0.5)

simulated_evals = []
for idx, name in enumerate(list_of_sims):
    # in these it can happen that there are no hits, so we need to check for that:
    if os.stat(name).st_size != 0:
        evals = list(pd.read_csv(name, delimiter="\t", usecols=[10]).to_numpy().transpose()[0])
        sb.kdeplot(evals, bw_method=0.1, color='#C00000', alpha=0.3, linewidth=0.5)
        simulated_evals.extend(evals)

med = np.median(simulated_evals)

plt.axvline(med, color='black', linestyle='--')
plt.xlim(0, 0.001)
plt.xlabel("E-values")
plt.ylabel("Count")
plt.savefig(OUT_FIGURE, dpi=300, format="png")


print(f"The median hit e-value of the randomized reads is: {med}.")