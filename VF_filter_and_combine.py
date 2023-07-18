from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import re
from tqdm import tqdm


EVALUE_THRESH_VF = 0.000487

print("Retrieving VFG IDs & names...")
fasta_sequences = SeqIO.parse(open("./virulence_factors/VFDB_setA_pro.fas"),'fasta')
fasta_headers = [f.description for f in fasta_sequences]
vfg_to_vfid = {header.split('(')[0]: re.findall("\(VF\d{4}\)", header)[0][1:-1] for header in fasta_headers}


print("Retrieving files with BLAST results...")
dir = "./virulence_factors/VF-hits/"
list_of_hits = []
for root, dirs, files in os.walk(dir):
	for file in files:
            if file.endswith("tsv"):
                list_of_hits.append(os.path.join(root,file))


print("Retrieving BLAST hits...")
print(f"E-value threshold set at {EVALUE_THRESH_VF}.")
vf_overview = {}
sample_ids = []
for idx, file in enumerate(tqdm(list_of_hits)):
    # read in sample-specific file
    fname = file.split('/')[-1].split('V')[0]
    sample_ids.append(fname)
    vf_hits = pd.read_csv(
        file, 
        delimiter="\t", 
        usecols=[0, 1, 10], 
        header=None, 
        names=["read", "VFG", "eval"]
    )

    # keep only the lowest e-value hit per read
    vf_hits = vf_hits.groupby('read', as_index=False).min()
    vf_hits.drop(columns=["read"], inplace=True)

    # filter out e-values above the threshold
    vir_fax = list(vf_hits.VFG)
    hit_evals = vf_hits['eval'].to_numpy()
    for name, e in zip(vir_fax, hit_evals):
        vf_id = vfg_to_vfid[name.split('(')[0]]
        if e < EVALUE_THRESH_VF:
            if vf_id in vf_overview.keys():
                vf_overview[vf_id][idx] += 1
            else:
                vf_overview[vf_id] = np.zeros(len(list_of_hits))
                vf_overview[vf_id][idx] = 1



print("Exporting...")
export_name = "./virulence_factors/virulence_factors_absolute_overview.csv"
vf_overview_df = pd.DataFrame.from_dict(vf_overview, orient='index', columns=sample_ids)
vf_overview_df.to_csv(export_name)
print(f"Wrote out overview to {export_name}!")
vf_overview_rel_df = vf_overview_df / vf_overview_df.sum(0)
vf_overview_rel_df.to_csv("./virulence_factors/virulence_factors_relative_overview.csv")

print("Retrieving group names and mapping...")
vfs_found_ordered = list(vf_overview_df.index)
vf_metadata = pd.read_excel("./virulence_factors/VFs.xls", index_col=0, header=1)
vf_metadata.drop(index=[i for i in vf_metadata.index if i not in vfs_found_ordered], inplace=True)
vf_metadata.to_csv("./virulence_factors/virulence_classes.csv")

print("concatenating categories...")
vf_metadata.drop(columns=[i for i in vf_metadata.columns if i != "VFcategory"], inplace=True)
vf_overview_df  = pd.merge(vf_overview_df, vf_metadata, left_index=True, right_index=True)

vf_classes_df = vf_overview_df.groupby('VFcategory').sum()
vf_classes_df.to_csv("./virulence_factors/virulence_factors_classes_absolute.csv")
vf_classes_rel_df = vf_classes_df / vf_classes_df.sum(0)
vf_classes_rel_df.to_csv("./virulence_factors/virulence_factors_classes_relative.csv")