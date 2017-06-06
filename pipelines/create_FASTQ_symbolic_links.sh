

INDIR=/store01/Collab_EJohannsen_RWelch/mark_RNAseq_MC_1st
OUTDIR=/store01/Collab_EJohannsen_RWelch/EBV_NOK/data/FASTQ

genome="akata-noks"
treat1="methyl_cell"
treat2="no_treatment"

ln -s $INDIR/MO9_CGATGT_L003_R1_001.fastq $OUTDIR/RNAseq-"$genome"-"$treat1"-clone1.fastq
ln -s $INDIR/MO10_TGACCA_L003_R1_001.fastq $OUTDIR/RNAseq-"$genome"-"$treat2"-clone1.fastq
ln -s $INDIR/MO11_CTTGTA_L003_R1_001.fastq $OUTDIR/RNAseq-"$genome"-"$treat1"-clone2.fastq
ln -s $INDIR/MO12_GTGAAA_L003_R1_001.fastq $OUTDIR/RNAseq-"$genome"-"$treat2"-clone2.fastq
