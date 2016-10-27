
clean:
	find . -type f -name '*~' -delete
	find . -type f -name '.*~' -delete

# parameters
cores=20
rsemDir=/u/w/e/welch/Desktop/Docs/Code/lib/RSEM-1.3.0
idxDir=../refs
hg19idx=$(idxDir)/hg19/hg19-rsem
akataidx=$(idxDir)/AKATA-GFP/AKATA_GFP_RSEM

# RSEM evaluation
data/RSEM/hg19/%.genes.results:data/FASTQ/%.fastq
	$(rsemDir)/rsem-calculate-expression -p $(cores) --estimate-rspd --append-names $^ $(hg19idx) $(@D)/$(basename $(basename $(@F)))

data/RSEM/AKATA-GFP/%.genes.results:data/FASTQ/%.fastq
	$(rsemDir)/rsem-calculate-expression -p $(cores) --estimate-rspd --append-names $^ $(akataidx) $(@D)/$(basename $(basename $(@F)))


# parse RSEM-bowtie into table
manuscript/%.txt:manuscript/logs/%/*
	rscripts/create_aligned_reads_report.R --verbose TRUE --directory $(<D) --outputfile $@ --aligner bowtie
