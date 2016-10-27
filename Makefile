
clean:
	find . -type f -name '*~' -delete
	find . -type f -name '.*~' -delete

# parameters
cores=20
rsemDir=/u/w/e/welch/Desktop/Docs/Code/lib/RSEM-1.3.0
idxDir=../refs
akataEBVidx=$(idxDir)/AKATA-EBV/AKATA_EBV-RSEM
hg19idx=$(idxDir)/hg19/hg19-rsem

# data/SAM/%.sam:data/FASTQ/%.fastq
# 	bowtie -q -S --best -p $(cores) -v 2 -m 1 -a --strata $(index) $^ $@ 

# data/BAM/%.bam:data/SAM/%.sam
# 	samtools view -bS $^ > $@

# data/RSEM/%.genes.results:data/BAM/%.bam
# 	$(rsemDir)/rsem-calculate-expression -p $(cores) --bam --estimate-rspd --append-names $^ $(index) $(@D)/$(basename $(basename $(@F)))

data/RSEM/%.genes.results:data/FASTQ/%.fastq
	$(rsemDir)/rsem-calculate-expression -p $(cores) --estimate-rspd --append-names $^ $(hg19idx) $(@D)/$(basename $(basename $(@F)))

