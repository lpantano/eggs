mkdir -p star
STAR   --runMode genomeGenerate   --genomeDir star   --genomeFastaFiles genome.fa --sjdbGTFfile transcript.gtf    --genomeSAindexNbases 9 --sjdbOverhang 99

STAR --genomeDir star --readFilesIn ../trimmed_aTc_R1.fastq ../trimmed_aTc_R2.fastq --runThreadN 1 --outFileNamePrefix output_ --outFilterMultimapNmax 10 

mkdir -p star_ann140
STAR   --runMode genomeGenerate   --genomeDir star_ann140   --genomeFastaFiles genome.fa --sjdbGTFfile transcript.gtf    --genomeSAindexNbases 9 --sjdbOverhang 140

STAR --genomeDir star_ann140 --readFilesIn ../trimmed_aTc_R1.fastq ../trimmed_aTc_R2.fastq --runThreadN 1 --outFileNamePrefix output_ --outFilterMultimapNmax 10 

mkdir -p star_alone
STAR   --runMode genomeGenerate   --genomeDir star_alone   --genomeFastaFiles genome.fa --genomeSAindexNbases 9

STAR --genomeDir star_alone --readFilesIn ../trimmed_aTc_R1.fastq ../trimmed_aTc_R2.fastq --runThreadN 1 --outFileNamePrefix output_ --outFilterMultimapNmax 10 
