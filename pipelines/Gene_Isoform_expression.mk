
## Directories

genedr=data/Diff.Genes/hg19/Gene_Isoform

data/TPM_matrices/Genes_TPM_matrix.tsv:
	Rscript rscripts/generate_TPM_matrices.R \
		--all_files "data/RSEM/hg19/*genes.results"  \
		--out_file $@

data/TPM_matrices/Isoforms_TPM_matrix.tsv:
	Rscript rscripts/generate_TPM_matrices.R \
		--all_files "data/RSEM/hg19/*isoforms.results" \
		--out_file $@

Scott_EBV_vs_NOKS_Genes:data/TPM_matrices/Genes_TPM_matrix.tsv
	Rscript rscripts/differential_analysis_by_cell.R \
		--A_noTr "data/RSEM/hg19/RNAseq-Noks-mono-rep?.genes.results" \
		--B_noTr "data/RSEM/hg19/RNAseq-Noks_EBV-mono-rep2.genes.results,data/RSEM/hg19/RNAseq-Noks_EBV-mono-rep3.genes.results,data/RSEM/hg19/RNAseq-Noks_EBV-mono-rep4.genes.results" \
		--A_Tr "data/RSEM/hg19/RNAseq-Noks-MC-rep?.genes.results" \
		--B_Tr "data/RSEM/hg19/RNAseq-Noks_EBV-MC-rep?.genes.results" \
		--cells "NOKS,EBV" \
		--treats "mono,MC" \
		--iso "gene" \
		--tpm_file $^ \
		--out_file_suff "$(genedr)/Scott_Genes_"
