# Directories

BASE=/store01/Collab_EJohannsen_RWelch

RSEMDIR:=/u/w/e/welch/Desktop/Docs/Code/lib/RSEM-1.3.0
INDEXDIR:=$(BASE)/refs

INDEX:=$(INDEXDIR)/hg19/hg19-rsem
SIZES:=/p/keles/SOFTWARE/hg19.chrom.sizes

BASE0=$(BASE)/EBV_NOK
DATADR=$(BASE0)/data
RDR=$(BASE0)/rscripts
FIGSDR=$(BASE0)/figs/diff_expression/DESeq2_marginal/dRdZ
CODEDR=/p/stat/genomics/bin


# Parameters

CORES:=16
IDXSIZES=/p/keles/SOFTWARE/hg19.chrom.sizes

# RSEM evaluation
$(DATADR)/RSEM/hg19/Sept17/%.genes.results:$(DATADR)/FASTQ/Sept17/%.fastq
	$(RSEMDIR)/rsem-calculate-expression  \
		-p $(CORES) \
		--estimate-rspd \
		--append-names $^ \
		$(INDEX) $(@D)/$(basename $(basename $(@F)))

$(DATADR)/RSEM/AKATA-GFP/Sept17/%.genes.results:$(DATADR)/FASTQ/Sept17/%.fastq
	$(RSEMDIR)/rsem-calculate-expression \
		-p $(CORES) \
		--estimate-rspd \
		--append-names $^ \
		$(INDEXDIR)/AKATA-GFP/AKATA_GFP_RSEM \
		$(@D)/$(basename $(basename $(@F)))

# Build TPM
$(DATADR)/TPM_matrices/Sept17/Genes_TPM_matrix_newBatch.tsv:$(DATADR)/RSEM/hg19/Sept17/*.genes.results
	Rscript rscripts/generate_TPM_matrices.R \
		--all_files "data/RSEM/hg19/Sept17/*.genes.results" \
		--out_file $@

$(DATADR)/TPM_matrices/Sept17/Isoforms_TPM_matrix_newBatch.tsv:$(DATADR)/RSEM/hg19/Sept17/*.isoforms.results
	Rscript rscripts/generate_TPM_matrices.R \
		--all_files "data/RSEM/hg19/Sept17/*.isoforms.results" \
		--out_file $@

# Tracks from RSEM
$(DATADR)/BAM/hg19/Sept17/%.bam:$(DATADR)/RSEM/hg19/Sept17/%.transcript.bam
	$(RSEMDIR)/rsem-tbam2gbam $(INDEX) $^ $@

$(DATADR)/BAM/hg19/Sept17/transcriptome/%.bam:$(DATADR)/BAM/hg19/Sept17/%.bam
	samtools sort $^ $(@:.bam=) ; samtools index $@

$(DATADR)/WIG/hg19/Sept17/transcriptome/%.wig:$(DATADR)/BAM/hg19/Sept17/%.bam
	$(RSEMDIR)/rsem-bam2wig $^ $@ $(@F)

$(DATADR)/BW/hg19/Sept17/transcriptome/%.bw:$(DATADR)/WIG/hg19/Sept17/transcriptome/%.wig
	wigToBigWig $^ $(SIZES) $@

# Marginal differential expression

MARGINAL:$(DATADR)/Diff.Genes/hg19/Sept17/EBV_dRdZ_diff_genes.tsv \
	$(DATADR)/Diff.Genes/hg19/Sept17/EBV_dRdZ_diff_isoforms.tsv \
	$(DATADR)/Diff.Genes/hg19/Sept17/EBV_diff_genes.tsv \
	$(DATADR)/Diff.Genes/hg19/Sept17/EBV_diff_isoforms.tsv \
	$(DATADR)/Diff.Genes/hg19/Sept17/NOKS_diff_genes.tsv \
	$(DATADR)/Diff.Genes/hg19/Sept17/NOKS_diff_isoforms.tsv


$(DATADR)/Diff.Genes/hg19/Sept17/EBV_dRdZ_diff_genes.tsv:$(DATADR)/TPM_matrices/Sept17/Genes_TPM_matrix_newBatch.tsv
	Rscript rscripts/differential_analysis_marginal.R \
		--no_treat_files "data/RSEM/hg19/Sept17/RNAseq-akata-noks-no_treatment-clone?.genes.results" \
		--treat_files "data/RSEM/hg19/Sept17/RNAseq-akata-noks-methyl_cell-clone?.genes.results" \
		--treats "no_treatment,methyl_cell" \
		--iso "gene" \
		--tpm_file $^ --out_file $@

$(DATADR)/Diff.Genes/hg19/Sept17/EBV_dRdZ_diff_isoforms.tsv:$(DATADR)/TPM_matrices/Sept17/Isoforms_TPM_matrix_newBatch.tsv
	Rscript rscripts/differential_analysis_marginal.R \
		--no_treat_files "data/RSEM/hg19/Sept17/RNAseq-akata-noks-no_treatment-clone?.isoforms.results" \
		--treat_files "data/RSEM/hg19/Sept17/RNAseq-akata-noks-methyl_cell-clone?.isoforms.results" \
		--treats "no_treatment,methyl_cell" \
		--iso "isoform" \
		--tpm_file $^ --out_file $@

$(DATADR)/Diff.Genes/hg19/Sept17/EBV_diff_genes.tsv:$(DATADR)/TPM_matrices/Sept17/Genes_TPM_matrix_newBatch.tsv
	Rscript rscripts/differential_analysis_marginal.R \
		--no_treat_files "data/RSEM/hg19/Sept17/RNAseq-akata-noks-no_treatment-rep?.genes.results" \
		--treat_files "data/RSEM/hg19/Sept17/RNAseq-akata-noks-methyl_cell-rep?.genes.results" \
		--treats "no_treatment,methyl_cell" \
		--iso "gene" \
		--tpm_file $^ --out_file $@

$(DATADR)/Diff.Genes/hg19/Sept17/EBV_diff_isoforms.tsv:$(DATADR)/TPM_matrices/Sept17/Isoforms_TPM_matrix_newBatch.tsv
	Rscript rscripts/differential_analysis_marginal.R \
		--no_treat_files "data/RSEM/hg19/Sept17/RNAseq-akata-noks-no_treatment-rep?.isoforms.results" \
		--treat_files "data/RSEM/hg19/Sept17/RNAseq-akata-noks-methyl_cell-rep?.isoforms.results" \
		--treats "no_treatment,methyl_cell" \
		--iso "isoform" \
		--tpm_file $^ --out_file $@

$(DATADR)/Diff.Genes/hg19/Sept17/NOKS_diff_genes.tsv:$(DATADR)/TPM_matrices/Sept17/Genes_TPM_matrix_newBatch.tsv
	Rscript rscripts/differential_analysis_marginal.R \
		--no_treat_files "data/RSEM/hg19/Sept17/RNAseq-noks-no_treatment-rep?.genes.results" \
		--treat_files "data/RSEM/hg19/Sept17/RNAseq-noks-methyl_cell-rep?.genes.results" \
		--treats "no_treatment,methyl_cell" \
		--iso "gene" \
		--tpm_file $^ --out_file $@

$(DATADR)/Diff.Genes/hg19/Sept17/NOKS_diff_isoforms.tsv:$(DATADR)/TPM_matrices/Sept17/Isoforms_TPM_matrix_newBatch.tsv
	Rscript rscripts/differential_analysis_marginal.R \
		--no_treat_files "data/RSEM/hg19/Sept17/RNAseq-noks-no_treatment-rep?.isoforms.results" \
		--treat_files "data/RSEM/hg19/Sept17/RNAseq-noks-methyl_cell-rep?.isoforms.results" \
		--treats "no_treatment,methyl_cell" \
		--iso "isoform" \
		--tpm_file $^ --out_file $@

Full_gene_tests:
	Rscript rscripts/differential_expression_analysis.R \
		--sample_dir $(DATADR)/RSEM/hg19/Sept17 \
		--samples_file $(DATADR)/Diff.Genes/hg19/Sept17/full_model/Sept17_Genes_samples_full.tsv \
		--contrast_file $(DATADR)/Diff.Genes/hg19/Sept17/full_model/Sept17_contrasts_full.tsv \
		--tpm_file $(DATADR)/TPM_matrices/Sept17/Genes_TPM_matrix_newBatch.tsv \
		--iso "gene" \
		--figs_dir figs/diff_expression/Sept17/full_model \
		--cores 20 \
		--out_dir $(DATADR)/Diff.Genes/hg19/Sept17/full_model

Full_isoform_tests:
	Rscript rscripts/differential_expression_analysis.R \
		--sample_dir $(DATADR)/RSEM/hg19/Sept17 \
		--samples_file $(DATADR)/Diff.Genes/hg19/Sept17/full_model/Sept17_Isoforms_samples_full.tsv \
		--contrast_file $(DATADR)/Diff.Genes/hg19/Sept17/full_model/Sept17_contrasts_full.tsv \
		--tpm_file $(DATADR)/TPM_matrices/Sept17/Isoforms_TPM_matrix_newBatch.tsv \
		--iso "isoform" \
		--figs_dir figs/diff_expression/Sept17/full_model \
		--cores 20 \
		--out_dir $(DATADR)/Diff.Genes/hg19/Sept17/full_model

# Alignment to genome
# $(DATADR)/BAM/hg19/bowtie_genome/Sept17/%.sam:$(DATADR)/FASTQ/Sept17/%.fastq
# 	/p/stat/genomics/bin/bowtie -q -n 2 \
# 		-m 200 -e 99999999 -p $(CORES) --best -S $(INDEX) $^ $@

# $(DATADR)/BAM/hg19/bowtie_genome/Sept17/%.bam:$(DATADR)/BAM/hg19/bowtie_genome/Sept17/%.sam
# 	samtools view -bS $^ | samtools sort - $(@:.bam=) ; samtools index $@

# data/BAM/hg19/%.genome.sort.bam:data/BAM/hg19/%.genome.bam
# 	bamtools sort -in $^ -out $@




# $(DATADR)/WIG/hg19/Sept17/%.wig:$(DATADR)/BAM/hg19/bowtie_genome/Sept17/%.bam
# 	Rscript rscripts/create_bw_file.R 

# $(DATADR)/BW/hg19/Sept17/%.bw:$(DATADR)/WIG/hg19/Sept17/%.wig
# 	wigToBigWig $^ $(IDXSIZES) $@

# $(DATADR)/WIG/hg19/Sept17/%.wig:$(DATADR)/BAM/hg19/%.genome.sort.bam
# 	$(rsemDir)/rsem-bam2wig $^ $@ $(@F)


# hg19.sizes=/p/keles/SOFTWARE/hg19.chrom.sizes
# data/BW/hg19/%.bw:data/WIG/hg19/%.wig
# 	wigToBigWig $^ $(hg19.sizes) $@

