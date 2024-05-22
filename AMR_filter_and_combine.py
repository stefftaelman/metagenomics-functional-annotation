from Bio import SeqIO
from matplotlib import use
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from tqdm import tqdm


EVALUE_THRESH_AMR = 0.000432

print("Retrieving AMR genbank IDs & names...")
acc_to_gene = pd.read_csv('./AMR_genes/aro_categories_index.tsv', sep="\t", index_col=0, usecols=[0, 2])["AMR Gene Family"].to_dict()
missing_amr_names = [k for (k,v) in acc_to_gene.items() if pd.isnull(v)]
found_missing = []

print("Retrieving files with BLAST results...")
dir = "./AMR_genes/AMR-hits/"
list_of_hits = []
for root, dirs, files in os.walk(dir):
    for file in files:
        if file.endswith("tsv"):
            list_of_hits.append(os.path.join(root,file))


print("Retrieving BLAST hits...")
print(f"E-value threshold set at {EVALUE_THRESH_AMR}.")
amr_overview = {}
amr_hit_coverage = {}
amr_hit_identity = {}
sample_ids = []
for idx, file in enumerate(tqdm(list_of_hits)):
    # read in sample-specific file
    fname = file.split('/')[-1].split('m')[0][:-1]
    sample_ids.append(fname)
    amr_hits = pd.read_csv(
        file, 
        delimiter="\t", 
        usecols=[0, 1, 2, 4, 6], 
        header=None, 
        names=["read", "AMR", "identity", "eval", "coverage"]
    )

    # keep only the lowest e-value hit per read
    amr_hits = amr_hits.groupby('read', as_index=False).min()
    amr_hits.drop(columns=["read"], inplace=True)
    amr_hits.dropna(how='all', inplace=True)

    # filter out e-values above the threshold
    amr_names = list(amr_hits.AMR)
    hit_evals = amr_hits['eval'].to_numpy()
    hit_pid = amr_hits['identity'].to_numpy()
    hit_cov = amr_hits['coverage'].to_numpy()
    for name, e, pid, cov in zip(amr_names, hit_evals, hit_pid, hit_cov):
        if name.split('|')[1] in missing_amr_names:
            found_missing.append(name)
        amr_gene = acc_to_gene[name.split('|')[1]]
        if e < EVALUE_THRESH_AMR:
            if amr_gene in amr_overview.keys():
                amr_overview[amr_gene][idx] = 1
            else:
                amr_overview[amr_gene] = np.zeros(len(list_of_hits))
                amr_overview[amr_gene][idx] = 1
            amr_hit_identity[sample_ids[idx]] = [amr_gene, pid]
            amr_hit_coverage[sample_ids[idx]] = [amr_gene, cov]



print("Exporting...")
export_name = "./AMR_genes/AMR_elements_absolute_overview.csv"
amr_overview_df = pd.DataFrame.from_dict(amr_overview, orient='index', columns=sample_ids)
amr_overview_df.to_csv(export_name)
print(f"Wrote out overview to {export_name}!")
amr_overview_rel_df = amr_overview_df / amr_overview_df.sum(0)
amr_overview_rel_df.to_csv("./AMR_genes/AMR_elements_relative_overview.csv")


amr_pid_df = pd.DataFrame.from_dict(amr_hit_identity, orient='index', columns=["gene", "identity"])
amr_pid_df.identity.describe()
amr_cov_df = pd.DataFrame.from_dict(amr_hit_coverage, orient='index', columns=["gene", "coverage"])
amr_cov_df.coverage.describe()

print("Retrieving group names and mapping...")
amrs_found_ordered = list(amr_overview_df.index)
amr_metadata = pd.read_csv("./AMR_genes/aro_categories_index.tsv", sep="\t", header=0, usecols=[2, 3, 4])
amr_metadata.drop_duplicates(subset=["AMR Gene Family"], inplace=True)
amr_metadata.set_index("AMR Gene Family", inplace=True)
amr_metadata.drop(index=[i for i in amr_metadata.index if i not in amrs_found_ordered], inplace=True)
amr_metadata.to_csv("./AMR_genes/AMR_classes.csv")

print("concatenating categories...")
amr_mech_metadata = amr_metadata.drop(columns=[i for i in amr_metadata.columns if i != "Resistance Mechanism"])
amr_mech_overview_df  = pd.merge(amr_overview_df, amr_mech_metadata, left_index=True, right_index=True)
amr_mechanisms_df = amr_mech_overview_df.groupby("Resistance Mechanism").sum()
amr_mechanisms_df.to_csv("./AMR_genes/AMR_mechanism_classes_absolute.csv")
amr_mechanisms_rel_df = amr_mechanisms_df / amr_mechanisms_df.sum(0)
amr_mechanisms_rel_df.to_csv("./AMR_genes/AMR_mechanism_classes_relative.csv")

amr_drugs_metadata = amr_metadata.drop(columns=[i for i in amr_metadata.columns if i != "Drug Class"])
amr_drugs_overview_df  = pd.merge(amr_overview_df, amr_drugs_metadata, left_index=True, right_index=True)
amr_drugs_df = amr_drugs_overview_df.groupby("Drug Class").sum()
amr_drugs_df.to_csv("./AMR_genes/AMR_drug_classes_absolute.csv")
amr_drugs_rel_df = amr_drugs_df / amr_drugs_df.sum(0)
amr_drugs_rel_df.to_csv("./AMR_genes/AMR_drug_classes_relative.csv")