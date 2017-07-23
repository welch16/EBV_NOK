

# Directories

BASE=/store01/Collab_EJohannsen_RWelch

RSEMDIR:=/u/w/e/welch/Desktop/Docs/Code/lib/RSEM-1.3.0
INDEXDIR:=$(BASE)/refs

INDEX:=$(INDEXDIR)/hg19/hg19-rsem
BASE0=$(BASE)/EBV_NOK
DATADR=$(BASE0)/data
RDR=$(BASE0)/rscripts
FIGSDR=$(BASE0)/figs/diff_expression/DESeq2_marginal/dRdZ
CODEDR=/p/stat/genomics/bin


# Parameters

CORES:=12
IDXSIZES=/p/keles/SOFTWARE/hg19.chrom.sizes

file1=$(DATADR)/RSEM/hg19/RNAseq-akata-noks-no_treatment-clone1.genes.results
file2=$(DATADR)/RSEM/hg19/RNAseq-akata-noks-methyl_cell-clone1.genes.results
file3=$(DATADR)/RSEM/hg19/RNAseq-akata-noks-no_treatment-clone2.genes.results
file4=$(DATADR)/RSEM/hg19/RNAseq-akata-noks-methyl_cell-clone2.genes.results

all:$(file1) $(file2) $(file3) $(file4)

# Diff expression analysis with DESeq2
dRdZ_analysis:$(file1) $(file2) $(file3) $(file4)
	$(RDR)/perform_differential_expression_marginal_contrasts.R \
            --mono_files "$(DATADR)/RSEM/hg19/RNAseq-akata-noks-no_treatment-clone?.genes.results" \
            --treat_files "$(DATADR)/RSEM/hg19/RNAseq-akata-noks-methyl_cell-clone?.genes.results" \
            --treatment 'none,MC' \
            --outfile $(DATADR)/Diff.Genes/hg19/DESeq2_marginal/dRdZ_MC_diff_genes_EBV.tsv \
            --var 'clone' --plot_title 'dRdZ : MC vs No treatment' --figs $(FIGSDR)/dRdZ_EBV

$(DATADR)/metadata/dRdZ_TPM_matrix.tsv:$(file1) $(file2) $(file3) $(file4)
	$(RDR)/create_TPM_matrix.R --all_files "$(DATADR)/RSEM/hg19/RNAseq-akata-noks-*clone*.genes.results" --outfile $@

#
bwig_tracks:$(DATADR)/BW/hg19/RNAseq-akata-noks-no_treatment-clone1.bw \
            $(DATADR)/BW/hg19/RNAseq-akata-noks-methyl_cell-clone1.bw \
            $(DATADR)/BW/hg19/RNAseq-akata-noks-no_treatment-clone2.bw \
            $(DATADR)/BW/hg19/RNAseq-akata-noks-methyl_cell-clone2.bw

# RSEM evaluation
$(DATADR)/RSEM/hg19/%.genes.results:$(DATADR)/FASTQ/%.fastq
	$(RSEMDIR)/rsem-calculate-expression -p $(CORES) --estimate-rspd --append-names $^ $(INDEX) $(@D)/$(basename $(basename $(@F)))

# Sort bam files from transcripts
$(DATADR)/BAM/hg19/%.genome.bam:$(DATADR)/RSEM/hg19/%.transcript.bam
	$(RSEMDIR)/rsem-tbam2gbam $(INDEX) $^ $@ 

$(DATADR)/BAM/hg19/%.genome.sort.bam:$(DATADR)/BAM/hg19/%.genome.bam
	$(CODEDR)/bamtools sort -in $^ -out $@

$(DATADR)/WIG/hg19/%.wig:$(DATADR)/BAM/hg19/%.genome.sort.bam
	$(RSEMDIR)/rsem-bam2wig $^ $@ $(@F)

$(DATADR)/BW/hg19/%.bw:$(DATADR)/WIG/hg19/%.wig
	wigToBigWig $^ $(IDXSIZES) $@


# Basic QC with FASTQC
# fastqc_analysis:
# 	$(CODEDR)/fastqc -o $(DATADR)/quality_control/FASTQC $(DATADR)/FASTQ/RNAseq*clone*.fastq


# # parse RSEM-bowtie into table
# manuscript/%.txt:manuscript/logs/%/*
# 	rscripts/create_aligned_reads_report.R --verbose TRUE --directory $(<D) --outputfile $@ --aligner bowtie

