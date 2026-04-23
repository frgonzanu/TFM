THREADS=6
OUTDIR=~/fgonzaleznun/TFM/data/trimmed/

for R1 in SRR*_1.fastq.gz
do
    SRR=$(basename $R1 | cut -d'_' -f1)
    R2=${SRR}_2.fastq.gz

    echo "Procesando SmartSeq2 - $SRR"

    cutadapt \
      -b AGATCGGAAGAGC -B AGATCGGAAGAGC \
      -g AAGCAGTGGTATCAACGCAGAGTACATGGGAAGCAGTGGTATCAACGCAGAGTACATGGG \
      -g AAGCAGTGGTATCAACGCAGAGTACATGGG \
      -g AAGCAGTGGTATCAACGCAGAGCACACGTCTGAACTCCAGTCAC \
      -G AAGCAGTGGTATCAACGCAGAGCACACGTCTGAACTCCAGTCAC \
      -G AAGCAGTGGTATCAACGCAGAGT \
      -G T{20} \
      -a CCCATGTACTCTGCGTTGATACCACTGCTT \
      -A CCCATGTACTCTGCGTTGATACCACTGCTT \
      -a TCAGACGTGTGCTCTTCCGATCT \
      -A TCAGACGTGTGCTCTTCCGATCT \
      -q 20,20 \
      -e 0.2 \
      --overlap 5 \
      --minimum-length 25 \
      -j $THREADS \
      -o $OUTDIR/${SRR}_1.trim.fastq.gz \
      -p $OUTDIR/${SRR}_2.trim.fastq.gz \
      $R1 $R2
done
