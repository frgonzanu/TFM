STAR \
  --runThreadN 16 \
  --runMode genomeGenerate \
  --genomeDir star_index_hg38 \
  --genomeFastaFiles genome.fa \
  --sjdbGTFfile genes_ann.gtf \
  --sjdbOverhang 149 \
  --limitGenomeGenerateRAM 50000000000
