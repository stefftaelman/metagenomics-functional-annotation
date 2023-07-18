# metagenomics-functional-annotation
This repo contains the code to perform the annotation of virulence factors (VFs) and antimicrobial resistance (AMR) genes to whole genome metagenomic sequencing data. These scripts were used in preparation of the following publications:
- TBA
- TBA

# Methodology

These scripts use BLAST to annotate known VFs/AMRs to metagenomic sequencing reads. When BLASTing sequences, there will also be a chance of subsequences randomly matching each other. In BLAST, this is quantified in the e-value, and the e-value cutoff below which you should assume every match to be a true hit will depend on library sizes of both reference database and your samples.To fix this, we’ll empirically determine the cutoff e-value by randomizing a subset of our fasta files (thus preserving library size and individual nucleotide proportions), and running the annotation against those as well. We can then plot the e-values for all our different reads, including the randomized ones, and threshold the true hits on the median e-value of the hits on the randomized reads.

# Reference data

In this repo, two public databases are used to respectively annotate virulence factors ([VFDB](http://www.mgc.ac.cn/VFs/main.htm)) and antimicrobial resistance elements ([CARD](https://card.mcmaster.ca/)). While reference data from these databases is publically available for download, CARD only allows its use for academic purposes, which is why they are not included in this repository. To be able to run the code in the scripts above, please download the following files from the latest releases of CARD and VFDB:
- "aro_categories_index.tsv" (from CARD; accessed 28/03/2023)
- "protein_fasta_protein_homolog_model.fasta" (from CARD; accessed 28/03/2023)
- "VFDB_setA_pro.fas" (from VFDB; accessed 28/03/2023)
- "VFs.xls" (from VFDB; accessed 28/03/2023)

## Metagenomic read annotation

To annotate the reads with VF and AMR functions, the following steps were taken:
1. Shuffling the reads of a subset of fasta files (see `shuffle_reads.py`).
2. Run [Diamond BLASTx](https://github.com/bbuchfink/diamond) on each fasta file against the reference dataset.
	a. Install diamond:
		```
		wget http://github.com/bbuchfink/diamond/releases/download/v2.1.6/diamond-linux64.tar.gz
		tar xzf diamond-linux64.tar.gz
		```
	b. Make a diamond compatible reference database out of the reference data fasta file
		```
		./diamond makedb --in VFDB_setA_pro.fas -d vf_reference
		```
	c. Run diamond (see `get_blastx_hits.sh`)
3. Extract an empirical threshold e-value (see `plot_evalues.py`) and discard hits with higher e-values (see `*_filter_and_combine.py`).

# References
- Alcock, Brian P et al. “CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database.” Nucleic acids research vol. 51,D1 (2023). D690-D699. doi:10.1093/nar/gkac920
- Chen, Lihong et al. “VFDB: a reference database for bacterial virulence factors.” Nucleic acids research vol. 33,Database issue (2005). D325-8. doi:10.1093/nar/gki008
- Buchfink B, Reuter K, Drost HG, "Sensitive protein alignments at tree-of-life scale using DIAMOND", Nature Methods 18, 366–368 (2021). doi:10.1038/s41592-021-01101-x
