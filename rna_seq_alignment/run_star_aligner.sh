
STAR=$1
INDEX_DIR=$2
FQ1=$3
FQ2=$4
OUTDIR=$5

NUM_THREAD=4

# --quantMode TranscriptomeSAM # This mode will generate Transcriptome mapped BAM file 
#							     will only work if STAR index is generated with Annotation

$STAR --runThreadN $NUM_THREAD \
      --genomeDir $INDEX_DIR \
      --readFilesIn $FQ1 $FQ2 \
      --outSAMtype BAM SortedByCoordinate \
      --quantMode TranscriptomeSAM \
      --outFileNamePrefix $OUTDIR \
      --readFilesCommand gunzip -c 


# Sort and Index Transcriptome BAM
transcriptome_bam="${OUTDIR}/Aligned.toTranscriptome.out.bam"
sorted_transcriptome_bam="${OUTDIR}/Aligned.toTranscriptome.out.sorted.bam"

samtools sort "${transcriptome_bam} > ${sorted_transcriptome_bam}"
samtools index "${sorted_transcriptome_bam}"

