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
	rscripts/create_data_summary_plots.R \
		--A '$(dataDr)/RNAseq-Noks-mono-rep?.genes.results' \
		--A_diff '$(dataDr)/RNAseq-Noks-MC-rep?.genes.results' \
		--B '$(dataDr)/RNAseq-Noks_EBV-mono-rep?.genes.results' \
		--B_diff '$(dataDr)/RNAseq-Noks_EBV-MC-rep?.genes.results' --type rsem --outputfile $(@D)/$(basename $(basename $(@F))) --xlab 'log2FC(MC)' --ylab 'EBV:log2FC(MC)'

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

dataDr=data/RSEM/hg19
outDr=data/Diff.Genes/hg19/DESeq2_marginal
figsDr=figs/diff_expression/DESeq2_marginal

DESeq2_marginal_diff_expression:
	make $(outDr)/ben_CaFBS_diff_genes_NOKS.tsv & 
	make $(outDr)/ben_CaFBS_diff_genes_EBV.tsv & 
	make $(outDr)/ben_MC_diff_genes_NOKS.tsv &
	make $(outDr)/ben_MC_diff_genes_EBV.tsv &
	make $(outDr)/scott_MC_diff_genes_NOKS.tsv &
	make $(outDr)/scott_MC_diff_genes_EBV.tsv &

$(outDr)/ben_CaFBS_diff_genes_NOKS.tsv:$(dataDr)/*.genes.results
	rscripts/perform_differential_expression_marginal_contrasts.R --mono_files 'data/RSEM/hg19/RNAseq-noks-no_treatment-rep?.genes.results' --treat_files 'data/RSEM/hg19/RNAseq-noks-CaFBS-rep?.genes.results' --treatment 'none,CaFBS' --outfile $@ --type rsem --figs figs/diff_expression/DESeq2_marginal/ben_CaFBS_NOKS --plot_title 'Ben NOKS: CaFBS vs. mono'

$(outDr)/ben_CaFBS_diff_genes_EBV.tsv:$(dataDr)/*.genes.results
	rscripts/perform_differential_expression_marginal_contrasts.R --mono_files 'data/RSEM/hg19/RNAseq-akata-noks-no_treatment-rep?.genes.results' --treat_files 'data/RSEM/hg19/RNAseq-akata-noks-CaFBS-rep?.genes.results' --treatment 'none,MC' --outfile $@ --type rsem --figs figs/diff_expression/DESeq2_marginal/ben_CaFBS_EBV --plot_title 'Ben EBV: CaFBS vs. mono'

$(outDr)/ben_MC_diff_genes_NOKS.tsv:$(dataDr)/*.genes.results
	rscripts/perform_differential_expression_marginal_contrasts.R --mono_files 'data/RSEM/hg19/RNAseq-noks-no_treatment-rep?.genes.results' --treat_files 'data/RSEM/hg19/RNAseq-noks-methyl_cell-rep?.genes.results' --treatment 'none,MC' --outfile $@ --type rsem --figs figs/diff_expression/DESeq2_marginal/ben_MC_NOKS --plot_title 'Ben NOKS: MC vs. mono'

$(outDr)/ben_MC_diff_genes_EBV.tsv:$(dataDr)/*.genes.results
	rscripts/perform_differential_expression_marginal_contrasts.R --mono_files 'data/RSEM/hg19/RNAseq-akata-noks-no_treatment-rep?.genes.results' --treat_file 'data/RSEM/hg19/RNAseq-akata-noks-methyl_cell-rep?.genes.results' --treatment 'none,MC' --outfile $@ --type rsem --figs figs/diff_expression/DESeq2_marginal/ben_MC_EBV --plot_title 'Ben EBV: MC vs. mono'

$(outDr)/scott_MC_diff_genes_NOKS.tsv:$(dataDr)/*.genes.results
	rscripts/perform_differential_expression_marginal_contrasts.R --mono_files 'data/RSEM/hg19/RNAseq-Noks-mono-rep?.genes.results' --treat_files 'data/RSEM/hg19/RNAseq-Noks-MC-rep?.genes.results' --treatment 'none,MC' --outfile $@ --type rsem --figs figs/diff_expression/DESeq2_marginal/scott_MC_NOKS --plot_title 'Scott NOKS: MC vs. mono'

$(outDr)/scott_MC_diff_genes_EBV.tsv:$(dataDr)/*.genes.results
	rscripts/perform_differential_expression_marginal_contrasts.R --mono_files 'data/RSEM/hg19/RNAseq-Noks_EBV-mono-rep2.genes.results,data/RSEM/hg19/RNAseq-Noks_EBV-mono-rep3.genes.results,data/RSEM/hg19/RNAseq-Noks_EBV-mono-rep4.genes.results' --treat_files 'data/RSEM/hg19/RNAseq-Noks_EBV-MC-rep?.genes.results' --treatment 'none,MC' --outfile $@ --type rsem --figs figs/diff_expression/DESeq2_marginal/scott_MC_EBV --plot_title 'Scott EBV: MC vs. mono'


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

data/quality_control/MEDIPS/%_medips.tsv:$(dataDr)/%.sort.bam
	rscripts/MEDIPS_QC_metrics.R --readsfile $^ --outfile $@ --shift 0 --extend 200 --window_size 200 --figs ./figs/methylation/medips/$(@F:_medips.tsv=) --cores 12

data/bins/MeDIPseq/BinMatrix_binsize$(binSize)_fragLen$(fragLen).tsv:data/metadata/MeDIPseq_definition.tsv
	rscripts/create_bin_matrix.R --design_file $^ --bin_size $(binSize) --frag_len $(fragLen) --size_file /p/keles/SOFTWARE/hg19.chrom.sizes --outfile $@ --cores 12

data/bins/MeDIPseq/BinMatrix_binsize$(binSize)_fragLen$(fragLen)_InputPooled.tsv:data/metadata/MeDIPseq_definition_InputPooled.tsv
	rscripts/create_bin_matrix.R --design_file $^ --bin_size $(binSize) --frag_len $(fragLen) --size_file /p/keles/SOFTWARE/hg19.chrom.sizes --outfile $@ --cores 12

#############################################################################################################################################################################################

figsDr=./figs/methylation/peaks
binMat=./data/bins/MeDIPseq/BinMatrix_binsize200_fragLen200_InputPooled.tsv
peakDir=./data/peaks/hg19/MeDIPseq

call_peaks_EBV_CaFBS:$(binMat)
	make binSize=200 fragLen=200 peakSample='EBV_CaFBS_hmC_Rep1' fdr=5 ./data/peaks/hg19/MeDIPseq/Peaks_EBV_CaFBS_hmC_Rep1_binsize200_fragLen200_FDR5.tsv
	make binSize=200 fragLen=200 peakSample='EBV_CaFBS_hmC_Rep2' fdr=5 ./data/peaks/hg19/MeDIPseq/Peaks_EBV_CaFBS_hmC_Rep2_binsize200_fragLen200_FDR5.tsv
	make binSize=200 fragLen=200 peakSample='EBV_CaFBS_hmC_Rep3' fdr=5 ./data/peaks/hg19/MeDIPseq/Peaks_EBV_CaFBS_hmC_Rep3_binsize200_fragLen200_FDR5.tsv
	make binSize=200 fragLen=200 peakSample='EBV_CaFBS_mC_Rep1' fdr=5 ./data/peaks/hg19/MeDIPseq/Peaks_EBV_CaFBS_mC_Rep1_binsize200_fragLen200_FDR5.tsv
	make binSize=200 fragLen=200 peakSample='EBV_CaFBS_mC_Rep2' fdr=5 ./data/peaks/hg19/MeDIPseq/Peaks_EBV_CaFBS_mC_Rep2_binsize200_fragLen200_FDR5.tsv
	make binSize=200 fragLen=200 peakSample='EBV_CaFBS_mC_Rep3' fdr=5 ./data/peaks/hg19/MeDIPseq/Peaks_EBV_CaFBS_mC_Rep3_binsize200_fragLen200_FDR5.tsv

call_peaks_EBV_mono:$(binMat)
	make binSize=200 fragLen=200 peakSample='EBV_mono_hmC_Rep1' fdr=5 ./data/peaks/hg19/MeDIPseq/Peaks_EBV_mono_hmC_Rep1_binsize200_fragLen200_FDR5.tsv
	make binSize=200 fragLen=200 peakSample='EBV_mono_hmC_Rep2' fdr=5 ./data/peaks/hg19/MeDIPseq/Peaks_EBV_mono_hmC_Rep2_binsize200_fragLen200_FDR5.tsv
	make binSize=200 fragLen=200 peakSample='EBV_mono_hmC_Rep3' fdr=5 ./data/peaks/hg19/MeDIPseq/Peaks_EBV_mono_hmC_Rep3_binsize200_fragLen200_FDR5.tsv
	make binSize=200 fragLen=200 peakSample='EBV_mono_mC_Rep1' fdr=5 ./data/peaks/hg19/MeDIPseq/Peaks_EBV_mono_mC_Rep1_binsize200_fragLen200_FDR5.tsv
	make binSize=200 fragLen=200 peakSample='EBV_mono_mC_Rep2' fdr=5 ./data/peaks/hg19/MeDIPseq/Peaks_EBV_mono_mC_Rep2_binsize200_fragLen200_FDR5.tsv
	make binSize=200 fragLen=200 peakSample='EBV_mono_mC_Rep3' fdr=5 ./data/peaks/hg19/MeDIPseq/Peaks_EBV_mono_mC_Rep3_binsize200_fragLen200_FDR5.tsv

call_peaks_NOK_CaFBS:$(binMat)
	make binSize=200 fragLen=200 peakSample='NOK_CaFBS_hmC_Rep1' fdr=5 ./data/peaks/hg19/MeDIPseq/Peaks_NOK_CaFBS_hmC_Rep1_binsize200_fragLen200_FDR5.tsv
	make binSize=200 fragLen=200 peakSample='NOK_CaFBS_hmC_Rep2' fdr=5 ./data/peaks/hg19/MeDIPseq/Peaks_NOK_CaFBS_hmC_Rep2_binsize200_fragLen200_FDR5.tsv
	make binSize=200 fragLen=200 peakSample='NOK_CaFBS_hmC_Rep3' fdr=5 ./data/peaks/hg19/MeDIPseq/Peaks_NOK_CaFBS_hmC_Rep3_binsize200_fragLen200_FDR5.tsv
	make binSize=200 fragLen=200 peakSample='NOK_CaFBS_mC_Rep1' fdr=5 ./data/peaks/hg19/MeDIPseq/Peaks_NOK_CaFBS_mC_Rep1_binsize200_fragLen200_FDR5.tsv
	make binSize=200 fragLen=200 peakSample='NOK_CaFBS_mC_Rep2' fdr=5 ./data/peaks/hg19/MeDIPseq/Peaks_NOK_CaFBS_mC_Rep2_binsize200_fragLen200_FDR5.tsv
	make binSize=200 fragLen=200 peakSample='NOK_CaFBS_mC_Rep3' fdr=5 ./data/peaks/hg19/MeDIPseq/Peaks_NOK_CaFBS_mC_Rep3_binsize200_fragLen200_FDR5.tsv

call_peaks_NOK_mono:$(binMat)
	make binSize=200 fragLen=200 peakSample='NOK_mono_hmC_Rep1' fdr=5 ./data/peaks/hg19/MeDIPseq/Peaks_NOK_mono_hmC_Rep1_binsize200_fragLen200_FDR5.tsv
	make binSize=200 fragLen=200 peakSample='NOK_mono_hmC_Rep2' fdr=5 ./data/peaks/hg19/MeDIPseq/Peaks_NOK_mono_hmC_Rep2_binsize200_fragLen200_FDR5.tsv
	make binSize=200 fragLen=200 peakSample='NOK_mono_hmC_Rep3' fdr=5 ./data/peaks/hg19/MeDIPseq/Peaks_NOK_mono_hmC_Rep3_binsize200_fragLen200_FDR5.tsv
	make binSize=200 fragLen=200 peakSample='NOK_mono_mC_Rep1' fdr=5 ./data/peaks/hg19/MeDIPseq/Peaks_NOK_mono_mC_Rep1_binsize200_fragLen200_FDR5.tsv

$(peakDir)/Peaks_$(peakSample)_binsize$(binSize)_fragLen$(fragLen)_FDR$(fdr).tsv:$(binMat)
	rscripts/MOSAiCS_call_MeDIPseq_peaks.R --bindata_file $^ --peakfile $@  --peak_sample $(peakSample) --read_dir $(dataDr) --fdr $(fdr) --figs $(figsDr)/$(peakSample) --cores 20

#############################################################################################################################################################################################


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




