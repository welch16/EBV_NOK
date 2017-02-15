
clean:
	find . -type f -name '*~' -delete
	find . -type f -name '.*~' -delete

# parameters
cores=20
rsemDir=/u/w/e/welch/Desktop/Docs/Code/lib/RSEM-1.3.0
idxDir=../refs
hg19idx=$(idxDir)/hg19/hg19-rsem
# akataidx=$(idxDir)/AKATA-GFP/AKATA_GFP_RSEM

# # RSEM evaluation
# data/RSEM/hg19/%.genes.results:data/FASTQ/%.fastq
# 	$(rsemDir)/rsem-calculate-expression -p $(cores) --estimate-rspd --append-names $^ $(hg19idx) $(@D)/$(basename $(basename $(@F)))

# data/RSEM/AKATA-GFP/%.genes.results:data/FASTQ/%.fastq
# 	$(rsemDir)/rsem-calculate-expression -p $(cores) --estimate-rspd --append-names $^ $(akataidx) $(@D)/$(basename $(basename $(@F)))

# Sort bam files from trasncripts
data/BAM/hg19/%.genome.bam:data/RSEM/hg19/%.transcript.bam
	$(rsemDir)/rsem-tbam2gbam $(idxDir)/hg19/hg19-rsem $^ $@ 

data/BAM/hg19/%.genome.sort.bam:data/BAM/hg19/%.genome.bam
	bamtools sort -in $^ -out $@

data/WIG/hg19/%.wig:data/BAM/hg19/%.genome.sort.bam
	$(rsemDir)/rsem-bam2wig $^ $@ $(@F)


hg19.sizes=/p/keles/SOFTWARE/hg19.chrom.sizes
data/BW/hg19/%.bw:data/WIG/hg19/%.wig
	wigToBigWig $^ $(hg19.sizes) $@

# # parse RSEM-bowtie into table
# manuscript/%.txt:manuscript/logs/%/*
# 	rscripts/create_aligned_reads_report.R --verbose TRUE --directory $(<D) --outputfile $@ --aligner bowtie

#############################################################################################################################################################################################

# genome alignment
hg19idx=$(idxDir)/hg19_bowtie/hg19

data/BAM/hg19/bowtie_genome/%_usual_param.sam:data/FASTQ/%.fastq
	/p/stat/genomics/bin/bowtie -q -v 2 -m 1 -p $(cores) --best -S $(hg19idx) $^ $@

data/BAM/hg19/bowtie_genome/%_usual_param1.bam:data/BAM/hg19/bowtie_genome/%_usual_param.sam
	samtools view -bS -o $@ $^ 

data/BAM/hg19/bowtie_genome/%_usual_param.bam:data/BAM/hg19/bowtie_genome/%_usual_param1.bam
	/p/stat/genomics/bin/bamtools sort -in $^ -out $@

data/BAM/hg19/bowtie_genome/%_rsem_default.sam:data/FASTQ/%.fastq
	/p/stat/genomics/bin/bowtie -q -n 2 -m 200 -e 99999999 -p $(cores) --best -S $(hg19idx) $^ $@

data/BAM/hg19/bowtie_genome/%_rsem_default1.bam:data/BAM/hg19/bowtie_genome/%_rsem_default.sam
	samtools view -bS -o $@ $^ 

data/BAM/hg19/bowtie_genome/%_rsem_default.bam:data/BAM/hg19/bowtie_genome/%_rsem_default1.bam
	/p/stat/genomics/bin/bamtools sort -in $^ -out $@




# dip data
data/BAM/hg19/bowtie_dip/%.sam:data/FASTQ/DIP/%.fastq
	/p/stat/genomics/bin/bowtie -q -v 2 -m 1 --best -p $(cores) -S $(hg19idx) $^ $@ 

data/BAM/hg19/bowtie_dip/%.unsorted.bam:data/BAM/hg19/bowtie_dip/%.sam
	samtools view -bS -o $@ $^

data/BAM/hg19/bowtie_dip/%.sort.bam:data/BAM/hg19/bowtie_dip/%.unsorted.bam
	/p/stat/genomics/bin/bamtools sort -in $^ -out $@ 

#############################################################################################################################################################################################

# DIP data analysis

## quality control of fastq sequences with fastqc

codedr=/p/stat/genomics/bin
indir=data/FASTQ/DIP
outdr=data/quality_control/FASTQC

FASTQC:
	$(codedr)/fastqc -o $(outdr) $(indir)/*.fastq.gz

#############################################################################################################################################################################################

exploratory_plots:
	make figs/exploratory/hg19/ben_CaFBS_NOK_vs_EBV.counts_hexbin_plot.pdf &
	make figs/exploratory/hg19/ben_MC_NOK_vs_EBV.counts_hexbin_plot.pdf &
	make figs/exploratory/hg19/scott_MC_NOK_vs_EBV.counts_hexbin_plot.pdf & 
	make figs/exploratory/hg19/ben_CaFBS_NOK_vs_MC_NOV.counts_hexbin_plot.pdf &

# exploratory plots with count and abundancy data
dataDr=data/RSEM/hg19
figs/exploratory/hg19/ben_CaFBS_NOK_vs_EBV.counts_hexbin_plot.pdf:data/RSEM/hg19/*.genes.results
	rscripts/create_data_summary_plots.R --A '$(dataDr)/RNAseq-noks-no_treatment-rep?.genes.results' --A_diff '$(dataDr)/RNAseq-noks-CaFBS-rep?.genes.results' --B '$(dataDr)/RNAseq-akata-noks-no_treatment-rep?.genes.results' --B_diff '$(dataDr)/RNAseq-akata-noks-CaFBS-rep?.genes.results' --type rsem --outputfile $(@D)/$(basename $(basename $(@F))) --xlab 'log2FC(CaFBS)' --ylab 'EBV:log2FC(CaFBS)'

figs/exploratory/hg19/ben_MC_NOK_vs_EBV.counts_hexbin_plot.pdf:data/RSEM/hg19/*.genes.results
	rscripts/create_data_summary_plots.R --A '$(dataDr)/RNAseq-noks-no_treatment-rep?.genes.results' --A_diff '$(dataDr)/RNAseq-noks-methyl_cell-rep?.genes.results' --B '$(dataDr)/RNAseq-akata-noks-no_treatment-rep?.genes.results' --B_diff '$(dataDr)/RNAseq-akata-noks-methyl_cell-rep?.genes.results' --type rsem --outputfile $(@D)/$(basename $(basename $(@F))) --xlab 'log2FC(MC)' --ylab 'EBV:log2FC(MC)'

figs/exploratory/hg19/scott_MC_NOK_vs_EBV.counts_hexbin_plot.pdf:data/RSEM/hg19/*.genes.results
	rscripts/create_data_summary_plots.R --A '$(dataDr)/RNAseq-Noks-mono-rep?.genes.results' --A_diff '$(dataDr)/RNAseq-Noks-MC-rep?.genes.results' --B '$(dataDr)/RNAseq-Noks_EBV-mono-rep?.genes.results' --B_diff '$(dataDr)/RNAseq-Noks_EBV-MC-rep?.genes.results' --type rsem --outputfile $(@D)/$(basename $(basename $(@F))) --xlab 'log2FC(MC)' --ylab 'EBV:log2FC(MC)'

figs/exploratory/hg19/ben_CaFBS_NOK_vs_MC_NOV.counts_hexbin_plot.pdf:data/RSEM/hg19/*.genes.results
	rscripts/create_data_summary_plots.R --A '$(dataDr)/RNAseq-noks-no_treatment-rep?.genes.results' --A_diff --A_diff '$(dataDr)/RNAseq-noks-CaFBS-rep?.genes.results' --B '$(dataDr)/RNAseq-noks-no_treatment-rep?.genes.results' --B_diff '$(dataDr)/RNAseq-noks-methyl_cell-rep?.genes.results' --type rsem --outputfile $(@D)/$(basename $(basename $(@F))) --xlab 'log2FC(CaFBS)' --ylab 'log2FC(MC)' 

exploratory_EBV_plots:
	make figs/exploratory/hg19/ben_EBV_notr_vs_CaFBS.counts_hexbin_plots.pdf & 
	make figs/exploratory/hg19/ben_EBV_notr_vs_MC.counts_hexbin_plots.pdf  &
	make figs/exploratory/hg19/scott_EBV_notr_vs_MC.counts_hexbin_plots.pdf &
	make figs/exploratory/hg19/ben_EBV_CaFBS_vs_MC.counts_hexbin_plots.pdf & 

figs/exploratory/hg19/ben_EBV_notr_vs_CaFBS.counts_hexbin_plots.pdf:data/RSEM/hg19/*.genes.results
	rscripts/create_data_summary_plots.R --A '$(dataDr)/RNAseq-noks-no_treatment-rep?.genes.results' --A_diff '$(dataDr)/RNAseq-akata-noks-no_treatment-rep?.genes.results' --B '$(dataDr)/RNAseq-noks-CaFBS-rep?.genes.results' --B_diff '$(dataDr)/RNAseq-akata-noks-CaFBS-rep?.genes.results' --type rsem --outputfile $(@D)/$(basename $(basename $(@F))) --xlab 'log2FC(notreatment)' --ylab 'log2FC(CaFBS)'

figs/exploratory/hg19/ben_EBV_notr_vs_MC.counts_hexbin_plots.pdf:data/RSEM/hg19/*.genes.results
	rscripts/create_data_summary_plots.R --A '$(dataDr)/RNAseq-noks-no_treatment-rep?.genes.results' --A_diff '$(dataDr)/RNAseq-akata-noks-no_treatment-rep?.genes.results' --B '$(dataDr)/RNAseq-noks-methyl_cell-rep?.genes.results' --B_diff '$(dataDr)/RNAseq-akata-noks-methyl_cell-rep?.genes.results' --type rsem --outputfile $(@D)/$(basename $(basename $(@F))) --xlab 'log2FC(notreatment)' --ylab 'log2FC(MC)'

figs/exploratory/hg19/scott_EBV_notr_vs_MC.counts_hexbin_plots.pdf:data/RSEM/hg19/*.genes.results
	rscripts/create_data_summary_plots.R --A '$(dataDr)/RNAseq-Noks-mono-rep?.genes.results' --A_diff '$(dataDr)/RNAseq-Noks_EBV-mono-rep?.genes.results' --B '$(dataDr)/RNAseq-Noks-MC-rep?.genes.results' --B_diff '$(dataDr)/RNAseq-Noks_EBV-MC-rep?.genes.results' --type rsem --outputfile $(@D)/$(basename $(basename $(@F))) --xlab 'log2FC(notreatment)' --ylab 'log2FC(MC)'

figs/exploratory/hg19/ben_EBV_CaFBS_vs_MC.counts_hexbin_plots.pdf:data/RSEM/hg19/*.genes.results
	rscripts/create_data_summary_plots.R --A '$(dataDr)/RNAseq-noks-CaFBS-rep?.genes.results' --A_diff '$(dataDr)/RNAseq-akata-noks-CaFBS-rep?.genes.results' --B '$(dataDr)/RNAseq-noks-methyl_cell-rep?.genes.results' --B_diff '$(dataDr)/RNAseq-akata-noks-methyl_cell-rep?.genes.results' --type rsem --outputfile $(@D)/$(basename $(basename $(@F))) --xlab 'log2FC(CaFBS)' --ylab 'log2FC(MC)'

exploratory_plots_EBVbase:
	make figs/exploratory/hg19/ben_EBVbase_notr_vs_MC.counts_hexbin_plots.pdf &
	make figs/exploratory/hg19/ben_EBVbase_notr_vs_CaFBS.counts_hexbin_plots.pdf &
	make figs/exploratory/hg19/scott_EBVbase_notr_vs_MC.counts_hexbin_plots.pdf & 

figs/exploratory/hg19/ben_EBVbase_notr_vs_MC.counts_hexbin_plots.pdf:data/RSEM/hg19/*.genes.results
	rscripts/create_data_summary_plots.R --A '$(dataDr)/RNAseq-noks-no_treatment-rep?.genes.results' --A_diff '$(dataDr)/RNAseq-akata-noks-no_treatment-rep?.genes.results' --B '$(dataDr)/RNAseq-noks-no_treatment-rep?.genes.results' --B_diff '$(dataDr)/RNAseq-noks-methyl_cell-rep?.genes.results' --type rsem --outputfile $(@D)/$(basename $(basename $(@F))) --xlab 'log2FC(notreatment)' --ylab 'log2FC(MC)'


figs/exploratory/hg19/ben_EBVbase_notr_vs_CaFBS.counts_hexbin_plots.pdf:data/RSEM/hg19/*.genes.results
	rscripts/create_data_summary_plots.R --A '$(dataDr)/RNAseq-noks-no_treatment-rep?.genes.results' --A_diff '$(dataDr)/RNAseq-akata-noks-no_treatment-rep?.genes.results' --B '$(dataDr)/RNAseq-noks-no_treatment-rep?.genes.results' --B_diff '$(dataDr)/RNAseq-noks-CaFBS-rep?.genes.results' --type rsem --outputfile $(@D)/$(basename $(basename $(@F))) --xlab 'log2FC(notreatment)' --ylab 'log2FC(CaFBS)'

figs/exploratory/hg19/scott_EBVbase_notr_vs_MC.counts_hexbin_plots.pdf:data/RSEM/hg19/*.genes.results
	rscripts/create_data_summary_plots.R --A '$(dataDr)/RNAseq-Noks-mono-rep?.genes.results' --A_diff '$(dataDr)/RNAseq-Noks_EBV-mono-rep?.genes.results' --B '$(dataDr)/RNAseq-Noks-mono-rep?.genes.results' --B_diff '$(dataDr)/RNAseq-Noks-MC-rep?.genes.results' --type rsem --outputfile $(@D)/$(basename $(basename $(@F))) --xlab 'log2FC(notreatment)' --ylab 'log2FC(MC)'




#############################################################################################################################################################################################

dataDr=data/RSEM/hg19
outDr=data/Diff.Genes/hg19/DESeq
figsDr=figs/diff_expression/DESeq
fdr=0.05

DESeq_diff_expression:
	make $(outDr)/ben_CaFBS_NOK_diff_exp_genes.tsv   &
	make $(outDr)/ben_CaFBS_EBV_NOK_diff_exp_genes.tsv & 
	make $(outDr)/ben_MC_NOK_diff_exp_genes.tsv  &
	make $(outDr)/ben_MC_EBV_NOK_diff_exp_genes.tsv & 
	make $(outDr)/scott_MC_NOK_diff_exp_genes.tsv &
	make $(outDr)/scott_MC_EBV_NOK_diff_exp_genes.tsv &

gene_overlap:
	make $(outDr)/comparison_MC_old_vs_scott.tsv & 

$(outDr)/ben_CaFBS_NOK_diff_exp_genes.tsv:$(dataDr)/*.genes.results
	rscripts/perform_differential_expression_analysis.R --base '$(dataDr)/RNAseq-noks-no_treatment-rep?.genes.results' --cond '$(dataDr)/RNAseq-noks-CaFBS-rep?.genes.results' --outfile $@ --type rsem --package deseq --figs $(figsDr)/ben_CaFBS_NOK --fdr $(fdr)

$(outDr)/ben_CaFBS_EBV_NOK_diff_exp_genes.tsv:$(dataDr)/*.genes.results
	rscripts/perform_differential_expression_analysis.R --base '$(dataDr)/RNAseq-akata-noks-no_treatment-rep?.genes.results' --cond '$(dataDr)/RNAseq-akata-noks-CaFBS-rep?.genes.results' --outfile $@ --type rsem --package deseq --figs $(figsDr)/ben_CaFBS_EBV_NOK --fdr $(fdr)

$(outDr)/ben_MC_NOK_diff_exp_genes.tsv:$(dataDr)/*.genes.results
	rscripts/perform_differential_expression_analysis.R --base '$(dataDr)/RNAseq-noks-no_treatment-rep?.genes.results' --cond '$(dataDr)/RNAseq-noks-methyl_cell-rep?.genes.results' --outfile $@ --type rsem --package deseq --figs $(figsDr)/ben_MC_NOK --fdr $(fdr)

$(outDr)/ben_MC_EBV_NOK_diff_exp_genes.tsv:$(dataDr)/*.genes.results
	rscripts/perform_differential_expression_analysis.R --base '$(dataDr)/RNAseq-akata-noks-no_treatment-rep?.genes.results' --cond '$(dataDr)/RNAseq-Noks_EBV-MC-rep?.genes.results' --outfile $@ --type rsem --package deseq --figs $(figsDr)/ben_MC_EBV_NOK --fdr $(fdr) 

$(outDr)/scott_MC_NOK_diff_exp_genes.tsv:$(dataDr)/*.genes.results
	rscripts/perform_differential_expression_analysis.R --base '$(dataDr)/RNAseq-Noks-mono-rep?.genes.results' --cond '$(dataDr)/RNAseq-Noks-MC-rep?.genes.results' --outfile $@ --type rsem --package deseq --figs $(figsDr)/scott_MC_NOK --fdr $(fdr)

$(outDr)/scott_MC_EBV_NOK_diff_exp_genes.tsv:$(dataDr)/*.genes.results
	rscripts/perform_differential_expression_analysis.R --base '$(dataDr)/RNAseq-Noks_EBV-mono-rep?.genes.results' --cond '$(dataDr)/RNAseq-Noks_EBV-MC-rep?.genes.results' --outfile $@ --type rsem --package deseq --figs $(figsDr)/scott_MC_EBV_NOK --fdr $(fdr)

# ******

$(outDr)/comparison_MC_old_vs_scott.tsv:$(outDr)/*.tsv
	rscripts/perform_gene_overlap_analysis.R --data_1A $(outDr)/ben_MC_NOK_diff_exp_genes.tsv --data_2A $(outDr)/scott_MC_NOK_diff_exp_genes.tsv --data_1B $(outDr)/ben_MC_EBV_NOK_diff_exp_genes.tsv --data_2B $(outDr)/scott_MC_EBV_NOK_diff_exp_genes.tsv --outfile $@ --cells 'NOK,EBV_NOK' --datasets 'Johannsen,Scott' --figs $(figsDr)/gene_overlap_MC_NOK_vs_EBV


#############################################################################################################################################################################################

dataDr=data/RSEM/hg19
outDr=data/Diff.Genes/hg19/DESeq2_contrasts
figsDr=figs/diff_expression/DESeq2

DESeq2_contrast_diff_expression:
	make $(outDr)/ben_CaFBS_diff_genes.tsv &
	make $(outDr)/ben_MC_diff_genes.tsv &
	make $(outDr)/scott_MC_diff_genes.tsv &
	make $(outDr)/scott_MC_diff_genes_wo_noTr_EBV-NOK_rep1.tsv &

$(outDr)/ben_CaFBS_diff_genes.tsv:$(dataDr)/*.genes.results
	rscripts/perform_differential_expression_contrast_analysis.R --A_noTr '$(dataDr)/RNAseq-noks-no_treatment-rep?.genes.results' --B_noTr '$(dataDr)/RNAseq-akata-noks-no_treatment-rep?.genes.results' --A_Tr '$(dataDr)/RNAseq-noks-CaFBS-rep?.genes.results' --B_Tr '$(dataDr)/RNAseq-akata-noks-CaFBS-rep?.genes.results' --cells 'NOK,EBV' --treatments 'none,CaFBS' --outfile $@ --type rsem --figs $(figsDr)/ben_CaFBS

$(outDr)/ben_MC_diff_genes.tsv:$(dataDr)/*.genes.results
	rscripts/perform_differential_expression_contrast_analysis.R --A_noTr '$(dataDr)/RNAseq-noks-no_treatment-rep?.genes.results' --B_noTr '$(dataDr)/RNAseq-akata-noks-no_treatment-rep?.genes.results' --A_Tr '$(dataDr)/RNAseq-noks-methyl_cell-rep?.genes.results' --B_Tr '$(dataDr)/RNAseq-akata-noks-methyl_cell-rep?.genes.results' --cells 'NOK,EBV' --treatments 'none,MC' --outfile $@ --type rsem --figs $(figsDr)/ben_MC

$(outDr)/scott_MC_diff_genes.tsv:$(dataDr)/*.genes.results
	rscripts/perform_differential_expression_contrast_analysis.R --A_noTr '$(dataDr)/RNAseq-Noks-mono-rep?.genes.results' --B_noTr '$(dataDr)/RNAseq-Noks_EBV-mono-rep?.genes.results' --A_Tr '$(dataDr)/RNAseq-Noks-MC-rep?.genes.results' --B_Tr '$(dataDr)/RNAseq-Noks_EBV-MC-rep?.genes.results' --cells 'NOK,EBV' --treatments 'none,MC' --outfile $@ --type rsem --figs $(figsDr)/scott_MC

$(outDr)/scott_MC_diff_genes_wo_noTr_EBV-NOK_rep1.tsv:$(dataDr)/*.genes.results
	rscripts/perform_differential_expression_contrast_analysis.R --A_noTr '$(dataDr)/RNAseq-Noks-mono-rep?.genes.results' --B_noTr '$(dataDr)/RNAseq-Noks_EBV-mono-rep2.genes.results,$(dataDr)/RNAseq-Noks_EBV-mono-rep3.genes.results,$(dataDr)/RNAseq-Noks_EBV-mono-rep4.genes.results' --A_Tr '$(dataDr)/RNAseq-Noks-MC-rep?.genes.results' --B_Tr '$(dataDr)/RNAseq-Noks_EBV-MC-rep?.genes.results' --cells 'NOK,EBV' --treatments 'none,MC' --outfile $@ --type rsem --figs $(figsDr)/scott_MC_wo_noTr_EBV-NOK_rep1

#############################################################################################################################################################################################

# Meta analysis and related stuff as summary tables

metadir=data/metadata

$(metadir)/PCA_definition.tsv:
	R CMD BATCH --vanilla rscripts/create_metadata_for_PCA.R

$(metadir)/gene_count_matrix.tsv:$(metadir)/PCA_definition.tsv
	rscripts/create_counts_matrix.R --metadatafile $^ --outfile $@

figsDr=figs/PCA_MDS

$(figsDr)/PCA_all_samples.pdf:$(metadir)/PCA_definition.tsv
	rscripts/create_PCA_plot_with_lab_effect.R --metadatafile $^ --ntop 500 --figsfile $@

$(figsDr)/MDS_all_samples.pdf:$(metadir)/PCA_definition.tsv
	rscripts/create_MDS_plot_with_lab_effect.R --metadatafile $^ --ntop 500 --figsfile $@

#############################################################################################################################################################################################

# DIP-analysis

### First step, compare for each treatment block: All Inputs, All mc and all hmc

dataDr=data/BAM/hg19/bowtie_dip
figsDr=figs/methylation/correlation
fragLen=200
binSize=200
log=TRUE
sizeFile=/p/keles/SOFTWARE/hg19.chrom.sizes

DIP_step1_mono:
	make NOKS_mono_Input & 
	make NOKS_mono_mc   &
	make NOKS_mono_hmc &
	make NOKS_akata_mono_Input & 
	make NOKS_akata_mono_mc  &
	make NOKS_akata_mono_hmc &

DIP_step1_CaFBS:
	make NOKS_CaFBS_Input & 
	make NOKS_CaFBS_mc  &
	make NOKS_CaFBS_hmc &
	make NOKS_akata_CaFBS_Input & 
	make NOKS_akata_CaFBS_mc  &
	make NOKS_akata_CaFBS_hmc  &

data/bins/MeDIPseq/BinMatrix_binsize$(binSize)_fragLen$(fragLen).tsv:data/metadata/MeDIPseq_definition.tsv
	rscripts/create_bin_matrix.R --design_file $^ --bin_size $(binSize) --frag_len $(fragLen) --size_file /p/keles/SOFTWARE/hg19.chrom.sizes --outfile $@ --cores 12

NOKS_mono_Input:$(dataDr)/*sort*
	./rscripts/compare_MeDIP_bins.R --samples '$(dataDr)/MeDIPseq-NOKS-mono-Input-rep?.sort.bam' --bin_size $(binSize) --frag_len $(fragLen) --size_file $(sizeFile) --figs $(figsDr)/MeDIPseq-NOKS-mono-Input --use_log $(log) --cores 2

NOKS_mono_mc:$(dataDr)/*sort*
	./rscripts/compare_MeDIP_bins.R --samples '$(dataDr)/MeDIPseq-NOKS-mono-mC-rep?.sort.bam' --bin_size $(binSize) --frag_len $(fragLen) --size_file $(sizeFile) --figs $(figsDr)/MeDIPseq-NOKS-mono-mc --use_log $(log) --cores 2

NOKS_mono_hmc:$(dataDr)/*sort*
	./rscripts/compare_MeDIP_bins.R --samples '$(dataDr)/MeDIPseq-NOKS-mono-hmC-rep?.sort.bam' --bin_size $(binSize) --frag_len $(fragLen) --size_file $(sizeFile) --figs $(figsDr)/MeDIPseq-NOKS-mono-hmc --use_log $(log) --cores 2

NOKS_akata_mono_Input:$(dataDr)/*sort*
	./rscripts/compare_MeDIP_bins.R --samples '$(dataDr)/MeDIPseq-NOKS-akata-mono-Input-rep?.sort.bam' --bin_size $(binSize) --frag_len $(fragLen) --size_file $(sizeFile) --figs $(figsDr)/MeDIPseq-NOKS-akata-mono-Input --use_log $(log) --cores 2

NOKS_akata_mono_mc:$(dataDr)/*sort*
	./rscripts/compare_MeDIP_bins.R --samples '$(dataDr)/MeDIPseq-NOKS-akata-mono-mC-rep?.sort.bam' --bin_size $(binSize) --frag_len $(fragLen) --size_file $(sizeFile) --figs $(figsDr)/MeDIPseq-akata-NOKS-mono-mc --use_log $(log) --cores 2

NOKS_akata_mono_hmc:$(dataDr)/*sort*
	./rscripts/compare_MeDIP_bins.R --samples '$(dataDr)/MeDIPseq-NOKS-akata-mono-hmC-rep?.sort.bam' --bin_size $(binSize) --frag_len $(fragLen) --size_file $(sizeFile) --figs $(figsDr)/MeDIPseq-akata-NOKS-mono-hmc --use_log $(log) --cores 2

NOKS_CaFBS_Input:$(dataDr)/*sort*
	./rscripts/compare_MeDIP_bins.R --samples '$(dataDr)/MeDIPseq-NOKS-CaFBS-Input-rep?.sort.bam' --bin_size $(binSize) --frag_len $(fragLen) --size_file $(sizeFile) --figs $(figsDr)/MeDIPseq-NOKS-CaFBS-Input --use_log $(log) --cores 2

NOKS_CaFBS_mc:$(dataDr)/*sort*
	./rscripts/compare_MeDIP_bins.R --samples '$(dataDr)/MeDIPseq-NOKS-CaFBS-mC-rep?.sort.bam' --bin_size $(binSize) --frag_len $(fragLen) --size_file $(sizeFile) --figs $(figsDr)/MeDIPseq-NOKS-CaFBS-mc --use_log $(log) --cores 2

NOKS_CaFBS_hmc:$(dataDr)/*sort*
	./rscripts/compare_MeDIP_bins.R --samples '$(dataDr)/MeDIPseq-NOKS-CaFBS-hmC-rep?.sort.bam' --bin_size $(binSize) --frag_len $(fragLen) --size_file $(sizeFile) --figs $(figsDr)/MeDIPseq-NOKS-CaFBS-hmc --use_log $(log) --cores 2

NOKS_akata_CaFBS_Input:$(dataDr)/*sort*
	./rscripts/compare_MeDIP_bins.R --samples '$(dataDr)/MeDIPseq-NOKS-akata-CaFBS-Input-rep?.sort.bam' --bin_size $(binSize) --frag_len $(fragLen) --size_file $(sizeFile) --figs $(figsDr)/MeDIPseq-NOKS-akata-CaFBS-Input --use_log $(log) --cores 2

NOKS_akata_CaFBS_mc:$(dataDr)/*sort*
	./rscripts/compare_MeDIP_bins.R --samples '$(dataDr)/MeDIPseq-NOKS-akata-CaFBS-mC-rep?.sort.bam' --bin_size $(binSize) --frag_len $(fragLen) --size_file $(sizeFile) --figs $(figsDr)/MeDIPseq-akata-NOKS-CaFBS-mc --use_log $(log) --cores 2

NOKS_akata_CaFBS_hmc:$(dataDr)/*sort*
	./rscripts/compare_MeDIP_bins.R --samples '$(dataDr)/MeDIPseq-NOKS-akata-CaFBS-hmC-rep?.sort.bam' --bin_size $(binSize) --frag_len $(fragLen) --size_file $(sizeFile) --figs $(figsDr)/MeDIPseq-akata-NOKS-CaFBS-hmc --use_log $(log) --cores 2






#############################################################################################################################################################################################




