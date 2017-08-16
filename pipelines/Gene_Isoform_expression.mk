
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

Scott_EBV_vs_NOKS_Isoforms:data/TPM_matrices/Isoforms_TPM_matrix.tsv
	Rscript rscripts/differential_analysis_by_cell.R \
		--A_noTr "data/RSEM/hg19/RNAseq-Noks-mono-rep?.isoforms.results" \
		--B_noTr "data/RSEM/hg19/RNAseq-Noks_EBV-mono-rep2.isoforms.results,data/RSEM/hg19/RNAseq-Noks_EBV-mono-rep3.isoforms.results,data/RSEM/hg19/RNAseq-Noks_EBV-mono-rep4.isoforms.results" \
		--A_Tr "data/RSEM/hg19/RNAseq-Noks-MC-rep?.isoforms.results" \
		--B_Tr "data/RSEM/hg19/RNAseq-Noks_EBV-MC-rep?.isoforms.results" \
		--cells "NOKS,EBV" \
		--treats "mono,MC" \
		--iso "isoform" \
		--tpm_file $^ \
		--out_file_suff "$(genedr)/Scott_Isoforms_"

dRdZ_vs_NOKS_Genes:data/TPM_matrices/Genes_TPM_matrix.tsv
	Rscript rscripts/differential_analysis_by_cell.R \
		--A_noTr "data/RSEM/hg19/RNAseq-Noks-mono-rep?.genes.results" \
		--B_noTr "data/RSEM/hg19/RNAseq-akata-noks-no_treatment-clone?.genes.results" \
		--A_Tr "data/RSEM/hg19/RNAseq-Noks-MC-rep?.genes.results" \
		--B_Tr "data/RSEM/hg19/RNAseq-akata-noks-methyl_cell-clone?.genes.results" \
		--cells "NOKS,EBV.dRdZ" \
		--treats "mono,MC" \
		--iso "gene" \
		--tpm_file $^ \
		--out_file_suff "$(genedr)/dRdZ_Genes_"

dRdZ_vs_NOKS_Isoforms:data/TPM_matrices/Isoforms_TPM_matrix.tsv
	Rscript rscripts/differential_analysis_by_cell.R \
		--A_noTr "data/RSEM/hg19/RNAseq-Noks-mono-rep?.isoforms.results" \
		--B_noTr "data/RSEM/hg19/RNAseq-akata-noks-no_treatment-clone?.isoforms.results" \
		--A_Tr "data/RSEM/hg19/RNAseq-Noks-MC-rep?.isoforms.results" \
		--B_Tr "data/RSEM/hg19/RNAseq-akata-noks-methyl_cell-clone?.isoforms.results" \
		--cells "NOKS,EBV.dRdZ" \
		--treats "mono,MC" \
		--iso "isoform" \
		--tpm_file $^ \
		--out_file_suff "$(genedr)/dRdZ_Isoforms_"

dRdZ_vs_EBV_Genes:data/TPM_matrices/Genes_TPM_matrix.tsv
	Rscript rscripts/differential_analysis_by_cell.R \
		--A_noTr "data/RSEM/hg19/RNAseq-Noks_EBV-mono-rep2.genes.results,data/RSEM/hg19/RNAseq-Noks_EBV-mono-rep3.genes.results,data/RSEM/hg19/RNAseq-Noks_EBV-mono-rep4.genes.results" \
		--B_noTr "data/RSEM/hg19/RNAseq-akata-noks-no_treatment-clone?.genes.results" \
		--A_Tr "data/RSEM/hg19/RNAseq-Noks-MC-rep?.genes.results" \
		--B_Tr "data/RSEM/hg19/RNAseq-Noks_EBV-MC-rep?.genes.results" \
		--cells "EBV,EBV.dRdZ" \
		--treats "mono,MC" \
		--iso "gene" \
		--tpm_file $^ \
		--out_file_suff "$(genedr)/dRdZ_Genes_"

dRdZ_vs_EBV_Isoforms:data/TPM_matrices/Isoforms_TPM_matrix.tsv
	Rscript rscripts/differential_analysis_by_cell.R \
		--A_noTr "data/RSEM/hg19/RNAseq-Noks_EBV-mono-rep2.isoforms.results,data/RSEM/hg19/RNAseq-Noks_EBV-mono-rep3.isoforms.results,data/RSEM/hg19/RNAseq-Noks_EBV-mono-rep4.isoforms.results" \
		--B_noTr "data/RSEM/hg19/RNAseq-akata-noks-no_treatment-clone?.isoforms.results" \
		--A_Tr "data/RSEM/hg19/RNAseq-Noks-MC-rep?.isoforms.results" \
		--B_Tr "data/RSEM/hg19/RNAseq-Noks_EBV-MC-rep?.isoforms.results" \
		--cells "EBV,EBV.dRdZ" \
		--treats "mono,MC" \
		--iso "isoform" \
		--tpm_file $^ \
		--out_file_suff "$(genedr)/dRdZ_Isoforms_"

Scott_reports:
	{ \
	Rscript reports/compile_Diff_Gene_Isoforms.R \
		--report_name "Scott_EBV_NOKS"\
		--gene_file "../$(genedr)/Scott_Genes_EBV_vs_NOKS.tsv" \
		--iso_file "../$(genedr)/Scott_Isoforms_EBV_vs_NOKS.tsv" \
		--test_name "EBV_vs_NOKS" ;\
	Rscript reports/compile_Diff_Gene_Isoforms.R \
		--report_name "Scott_EBV_NOKS_mono"\
		--gene_file "../$(genedr)/Scott_Genes_mono:EBV_vs_NOKS.tsv" \
		--iso_file "../$(genedr)/Scott_Isoforms_mono:EBV_vs_NOKS.tsv" \
		--test_name "mono:EBV_vs_NOKS" ;\
	Rscript reports/compile_Diff_Gene_Isoforms.R \
		--report_name "Scott_EBV_NOKS_MC"\
		--gene_file "../$(genedr)/Scott_Genes_MC:EBV_vs_NOKS.tsv" \
		--iso_file "../$(genedr)/Scott_Isoforms_MC:EBV_vs_NOKS.tsv" \
		--test_name "MC:EBV_vs_NOKS" ;\
	}

dRdZ_vs_NOKS_reports:
	{ \
	Rscript reports/compile_Diff_Gene_Isoforms.R \
		--report_name "dRdZ_vs_NOKS"\
		--gene_file "../$(genedr)/dRdZ_Genes_EBV.dRdZ_vs_NOKS.tsv" \
		--iso_file "../$(genedr)/dRdZ_Isoforms_EBV.dRdZ_vs_NOKS.tsv" \
		--test_name "EBV.dRdZ_vs_NOKS" ;\
	Rscript reports/compile_Diff_Gene_Isoforms.R \
		--report_name "dRdZ_vs_NOKS_mono"\
		--gene_file "../$(genedr)/dRdZ_Genes_mono:EBV.dRdZ_vs_NOKS.tsv" \
		--iso_file "../$(genedr)/dRdZ_Isoforms_mono:EBV.dRdZ_vs_NOKS.tsv" \
		--test_name "mono:EBV.dRdZ_vs_NOKS" ;\
	Rscript reports/compile_Diff_Gene_Isoforms.R \
		--report_name "dRdZ_NOKS_MC"\
		--gene_file "../$(genedr)/dRdZ_Genes_MC:EBV.dRdZ_vs_NOKS.tsv" \
		--iso_file "../$(genedr)/dRdZ_Isoforms_MC:EBV.dRdZ_vs_NOKS.tsv" \
		--test_name "MC:EBV.dRdZ_vs_NOKS" ;\
	}

dRdZ_vs_EBV_reports:
	{ \
	Rscript reports/compile_Diff_Gene_Isoforms.R \
		--report_name "dRdZ_vs_EBV"\
		--gene_file "../$(genedr)/dRdZ_Genes_EBV.dRdZ_vs_EBV.tsv" \
		--iso_file "../$(genedr)/dRdZ_Isoforms_EBV.dRdZ_vs_EBV.tsv" \
		--test_name "EBV.dRdZ_vs_EBV" ;\
	Rscript reports/compile_Diff_Gene_Isoforms.R \
		--report_name "dRdZ_vs_EBV_mono"\
		--gene_file "../$(genedr)/dRdZ_Genes_mono:EBV.dRdZ_vs_EBV.tsv" \
		--iso_file "../$(genedr)/dRdZ_Isoforms_mono:EBV.dRdZ_vs_EBV.tsv" \
		--test_name "mono:EBV.dRdZ_vs_EBV" ;\
	Rscript reports/compile_Diff_Gene_Isoforms.R \
		--report_name "dRdZ_EBV_MC"\
		--gene_file "../$(genedr)/dRdZ_Genes_MC:EBV.dRdZ_vs_EBV.tsv" \
		--iso_file "../$(genedr)/dRdZ_Isoforms_MC:EBV.dRdZ_vs_EBV.tsv" \
		--test_name "MC:EBV.dRdZ_vs_EBV" ;\
	}

