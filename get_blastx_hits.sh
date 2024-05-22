export COUNTER=0
export NO_SAMPLES=$(cat samples | wc -l)
for s in $(cat samples); 
	do \
	let COUNTER++
	echo "On sample: ${s} (${COUNTER}/${NO_SAMPLES})"
	./diamond blastx -d vf_reference -q ${s}.fasta -o ${s}_matches.tsv --outfmt 6 qseqid sseqid pident length evalue bitscore qcovhsp; 
	done
