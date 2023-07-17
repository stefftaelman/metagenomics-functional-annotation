# metagenomics-functional-annotation
This repo contains the code to perform the annotation of virulence factors (VFs) and antimicrobial resistance (AMR) genes to whole genome metagenomic sequencing data. These scripts were used in preparation of the following publications:
- TBA
- TBA

# Methodology
## Reference data

In this repo, two public databases are used to respectively annotate virulence factors ([VFDB](http://www.mgc.ac.cn/VFs/main.htm)) and antimicrobial resistance elements ([CARD]https://card.mcmaster.ca/(https://card.mcmaster.ca/)). While reference data from these databases is publically available for download, CARD only allows its use for academic purposes, which is why they are not included in this repository. To be able to run the code in the scripts above, please download the following files from the latest releases of CARD and VFDB:
- "aro_categories_index.tsv" (from CARD; accessed 28/03/2023)
- "protein_fasta_protein_homolog_model.fasta" (from CARD; accessed 28/03/2023)
- "VFDB_setA_pro.fas" (from VFDB; accessed 28/03/2023)
- "VFs.xls" (from VFDB; accessed 28/03/2023)

## Metagenomic read annotation

To annotate the reads with VF and AMR functions, the following steps were taken:
1. [Diamond blastx]() was run on each individual read against the protein reference dataset
2. All reads for a random subset of samples were shuffled as to retrieve a mock-set of reads with the same distributions of nucleotide fractions.
3. Diamond blastx was also run on these shuffled reads to retrieve an empirical E-value cut-off (to distinguish true VF/AMR gene matches from random hits).
4. Matches with higher E-values than the empirical cut-off were removed and for each read, the most likely true hit was annotated and logged into a count table.
5. Extra count tables were also generated aggregating the counts on prespecified VF- and AMR--groups, as specified in the source databases.

# References
- CARD (to add)
- VFDB (to add)
- DIAMOND (to add)
