#!/bin/sh


base=/store01/Collab_EJohannsen_RWelch
indir=$base/mark_RNA
outdir=$base/EBV_NOK/data/FASTQ

echo Files taken from the directory: $indir
echo To the directory: $outdir

treat1="methyl_cell"
treat2="no_treatment"

cell1="noks"
cell2="akata-noks"

ln -s $indir/MO21_AR005_ACAGTG_L005_R1_001.fastq.gz $outdir/RNAseq-$cell2-$treat1-clone3.fastq.gz
ln -s $indir/MO22_AR006_GCCAAT_L005_R1_001.fastq.gz $outdir/RNAseq-$cell2-$treat2-clone3.fastq.gz
ln -s $indir/MO23_AR012_CTTGTA_L005_R1_001.fastq.gz $outdir/RNAseq-$cell2-$treat1-clone4.fastq.gz
ln -s $indir/MO24_AR019_GTGAAA_L005_R1_001.fastq.gz $outdir/RNAseq-$cell2-$treat2-clone4.fastq.gz


ln -s $indir/MO37_AR005_ACAGTG_L007_R1_001.fastq.gz $outdir/RNAseq-$cell1-$treat1-rep1.fastq.gz
ln -s $indir/MO38_AR006_GCCAAT_L007_R1_001.fastq.gz $outdir/RNAseq-$cell1-$treat2-rep1.fastq.gz
ln -s $indir/MO39_AR012_CTTGTA_L007_R1_001.fastq.gz $outdir/RNAseq-$cell1-$treat1-rep2.fastq.gz
ln -s $indir/MO40_AR019_GTGAAA_L007_R1_001.fastq.gz $outdir/RNAseq-$cell1-$treat2-rep2.fastq.gz
ln -s $indir/MO41_AR002_CGATGT_L007_R1_001.fastq.gz $outdir/RNAseq-$cell2-$treat1-rep1.fastq.gz
ln -s $indir/MO42_AR004_TGACCA_L007_R1_001.fastq.gz $outdir/RNAseq-$cell2-$treat2-rep1.fastq.gz
ln -s $indir/MO43_AR007_CAGATC_L007_R1_001.fastq.gz $outdir/RNAseq-$cell2-$treat1-rep2.fastq.gz
ln -s $indir/MO44_AR014_AGTTCC_L007_R1_001.fastq.gz $outdir/RNAseq-$cell2-$treat2-rep2.fastq.gz

