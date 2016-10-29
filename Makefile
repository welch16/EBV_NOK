
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

##############################################################################################################################################################
exploratory_plots:
	make figs/exploratory/hg19/ben_CaFBS_NOK_vs_EBV.counts_hexbin_plot.pdf
	make figs/exploratory/hg19/ben_MC_NOK_vs_EBV.counts_hexbin_plot.pdf
	make figs/exploratory/hg19/scott_MC_NOK_vs_EBV.counts_hexbin_plot.pdf

# exploratory plots with count and abundancy data
dataDr=data/RSEM/hg19
figs/exploratory/hg19/ben_CaFBS_NOK_vs_EBV.counts_hexbin_plot.pdf:data/RSEM/hg19/*.genes.results
	rscripts/create_data_summary_plots.R --A '$(dataDr)/RNAseq-noks-no_treatment-rep?.genes.results' --A_diff '$(dataDr)/RNAseq-noks-CaFBS-rep?.genes.results' --B '$(dataDr)/RNAseq-akata-noks-no_treatment-rep?.genes.results' --B_diff '$(dataDr)/RNAseq-akata-noks-CaFBS-rep?.genes.results' --type rsem --outputfile $(@D)/$(basename $(basename $(@F))) --xlab 'log2FC(CaFBS)' --ylab 'EBV:log2FC(CaFBS)'

figs/exploratory/hg19/ben_MC_NOK_vs_EBV.counts_hexbin_plot.pdf:data/RSEM/hg19/*.genes.results
	rscripts/create_data_summary_plots.R --A '$(dataDr)/RNAseq-noks-no_treatment-rep?.genes.results' --A_diff '$(dataDr)/RNAseq-noks-methyl_cell-rep?.genes.results' --B '$(dataDr)/RNAseq-akata-noks-no_treatment-rep?.genes.results' --B_diff '$(dataDr)/RNAseq-akata-noks-methyl_cell-rep?.genes.results' --type rsem --outputfile $(@D)/$(basename $(basename $(@F))) --xlab 'log2FC(MC)' --ylab 'EBV:log2FC(MC)'

figs/exploratory/hg19/scott_MC_NOK_vs_EBV.counts_hexbin_plot.pdf:data/RSEM/hg19/*.genes.results
	rscripts/create_data_summary_plots.R --A '$(dataDr)/RNAseq-Noks-mono-rep?.genes.results' --A_diff '$(dataDr)/RNAseq-Noks-MC-rep?.genes.results' --B '$(dataDr)/RNAseq-Noks_EBV-mono-rep?.genes.results' --B_diff '$(dataDr)/RNAseq-Noks_EBV-MC-rep?.genes.results' --type rsem --outputfile $(@D)/$(basename $(basename $(@F))) --xlab 'log2FC(MC)' --ylab 'EBV:log2FC(MC)'
