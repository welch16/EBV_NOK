#!/bin/sh

# The function of this script is to rename all the dip data files to
# a more readable format.

# The rule that we followed is taken from Mark's email:

# The first character is the replicate number.
# The second character is either U or I.  U is uninduced (no Ca/FBS), I is induced (with Ca/FBS).
# If the third character is A, then the sample has Akata (EBV); if not, then it is uninfected.
# Finally, the name ends in either I, mC, or hmC.  I is input (the unenriched control),
#    mC is enriched for methyl cytosine, and hmC is enriched for hydroxymethyl cytosine. 
#
# So:
#
# 1UAI is untreated, EBV-infected, input control
# 3IhmC is Ca/FBS treated, uninfected, hmC enriched
#
# 

base=/store01/Collab_EJohannsen_RWelch
indir=$base/mark_dip
outdir=$base/EBV_NOK/data/FASTQ/DIP

echo Files taken from the directory: $indir
echo To the directory: $outdir

cell=NOKS
cond1=CaFBS
cond2=mono

suff1=Input
suff2=mC
suff3=hmC

ln -s $indir/1IAhmC_TGACCA_L002_R1_001.fastq.gz $outdir/MeDIPseq-$cell-akata-$cond1-$suff3-rep1.fastq.gz
ln -s $indir/1UAhmC_ATCACG_L002_R1_001.fastq.gz $outdir/MeDIPseq-$cell-akata-$cond2-$suff3-rep1.fastq.gz
ln -s $indir/2IAhmC_CGATGT_L001_R1_001.fastq.gz $outdir/MeDIPseq-$cell-akata-$cond1-$suff3-rep2.fastq.gz
ln -s $indir/2UAhmC_CGATGT_L002_R1_001.fastq.gz $outdir/MeDIPseq-$cell-akata-$cond2-$suff3-rep2.fastq.gz
ln -s $indir/3IAhmC_GCCAAT_L001_R1_001.fastq.gz $outdir/MeDIPseq-$cell-akata-$cond1-$suff3-rep3.fastq.gz
ln -s $indir/3UAhmC_TTAGGC_L002_R1_001.fastq.gz $outdir/MeDIPseq-$cell-akata-$cond2-$suff3-rep3.fastq.gz
ln -s $indir/1IAI_TGACCA_L002_R1_001.fastq.gz $outdir/MeDIPseq-$cell-akata-$cond1-$suff1-rep1.fastq.gz
ln -s $indir/1UAI_ATCACG_L002_R1_001.fastq.gz $outdir/MeDIPseq-$cell-akata-$cond2-$suff1-rep1.fastq.gz
ln -s $indir/2IAI_ACAGTG_L001_R1_001.fastq.gz $outdir/MeDIPseq-$cell-akata-$cond1-$suff1-rep2.fastq.gz
ln -s $indir/2UAI_CGATGT_L002_R1_001.fastq.gz $outdir/MeDIPseq-$cell-akata-$cond2-$suff1-rep2.fastq.gz
ln -s $indir/3IAI_GCCAAT_L002_R1_001.fastq.gz $outdir/MeDIPseq-$cell-akata-$cond1-$suff1-rep3.fastq.gz
ln -s $indir/3UAI_TTAGGC_L001_R1_001.fastq.gz $outdir/MeDIPseq-$cell-akata-$cond2-$suff1-rep3.fastq.gz
ln -s $indir/1IAmC_TGACCA_L001_R1_001.fastq.gz $outdir/MeDIPseq-$cell-akata-$cond1-$suff2-rep1.fastq.gz
ln -s $indir/1UAmC_ATCACG_L001_R1_001.fastq.gz $outdir/MeDIPseq-$cell-akata-$cond2-$suff2-rep1.fastq.gz
ln -s $indir/2IAmC_ACAGTG_L002_R1_001.fastq.gz $outdir/MeDIPseq-$cell-akata-$cond1-$suff2-rep2.fastq.gz
ln -s $indir/2UAmC_ACAGTG_L002_R1_001.fastq.gz $outdir/MeDIPseq-$cell-akata-$cond2-$suff2-rep2.fastq.gz
ln -s $indir/3IAmC_GCCAAT_L002_R1_001.fastq.gz $outdir/MeDIPseq-$cell-akata-$cond1-$suff2-rep3.fastq.gz
ln -s $indir/3UAmC_TTAGGC_L002_R1_001.fastq.gz $outdir/MeDIPseq-$cell-akata-$cond2-$suff2-rep3.fastq.gz
ln -s $indir/1IhmC_TGACCA_L001_R1_001.fastq.gz $outdir/MeDIPseq-$cell-$cond1-$suff3-rep1.fastq.gz
ln -s $indir/1UhmC_ATCACG_L001_R1_001.fastq.gz $outdir/MeDIPseq-$cell-$cond2-$suff3-rep1.fastq.gz
ln -s $indir/2IhmC_ACAGTG_L001_R1_001.fastq.gz $outdir/MeDIPseq-$cell-$cond1-$suff3-rep2.fastq.gz
ln -s $indir/2UhmC_ACAGTG_L002_R1_001.fastq.gz $outdir/MeDIPseq-$cell-$cond2-$suff3-rep2.fastq.gz
ln -s $indir/3IhmC_GCCAAT_L001_R1_001.fastq.gz $outdir/MeDIPseq-$cell-$cond1-$suff3-rep3.fastq.gz
ln -s $indir/3UhmC_TTAGGC_L002_R1_001.fastq.gz $outdir/MeDIPseq-$cell-$cond2-$suff3-rep3.fastq.gz
ln -s $indir/1II_TGACCA_L001_R1_001.fastq.gz $outdir/MeDIPseq-$cell-$cond1-$suff1-rep1.fastq.gz
ln -s $indir/1UI_ATCACG_L001_R1_001.fastq.gz $outdir/MeDIPseq-$cell-$cond2-$suff1-rep1.fastq.gz
ln -s $indir/2II_ACAGTG_L001_R1_001.fastq.gz $outdir/MeDIPseq-$cell-$cond1-$suff1-rep2.fastq.gz
ln -s $indir/2UI_CGATGT_L002_R1_001.fastq.gz $outdir/MeDIPseq-$cell-$cond2-$suff1-rep2.fastq.gz
ln -s $indir/3II_GCCAAT_L002_R1_001.fastq.gz $outdir/MeDIPseq-$cell-$cond1-$suff1-rep3.fastq.gz
ln -s $indir/3UI_TTAGGC_L001_R1_001.fastq.gz $outdir/MeDIPseq-$cell-$cond2-$suff1-rep3.fastq.gz
ln -s $indir/1ImC_TGACCA_L002_R1_001.fastq.gz $outdir/MeDIPseq-$cell-$cond1-$suff2-rep1.fastq.gz
ln -s $indir/1UmC_ATCACG_L002_R1_001.fastq.gz $outdir/MeDIPseq-$cell-$cond2-$suff2-rep1.fastq.gz
ln -s $indir/2ImC_CGATGT_L001_R1_001.fastq.gz $outdir/MeDIPseq-$cell-$cond1-$suff2-rep2.fastq.gz
ln -s $indir/2UmC_CGATGT_L001_R1_001.fastq.gz $outdir/MeDIPseq-$cell-$cond2-$suff2-rep2.fastq.gz
ln -s $indir/3ImC_GCCAAT_L001_R1_001.fastq.gz $outdir/MeDIPseq-$cell-$cond1-$suff2-rep3.fastq.gz
ln -s $indir/3UmC_TTAGGC_L001_R1_001.fastq.gz $outdir/MeDIPseq-$cell-$cond2-$suff2-rep3.fastq.gz
