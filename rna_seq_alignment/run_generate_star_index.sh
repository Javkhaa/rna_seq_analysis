
STAR=$1 # STAR executable
REF_FASTA=$2 # Reference FASTA file
ANNOTATIONS=$3 # Annotated GTF file
INDEX_DIR=$4 # STAR index directory
READLEN=$5 
NUM_THREAD=4

$STAR --runThreadN $NUM_THREAD \
	  --runMode genomeGenerate \
	  --genomeDir $INDEX_DIR \
	  --genomeFastaFiles $REF_FASTA \
	  --genomeSAindexNbases 10 \
	  --sjdbGTFfile $ANNOTATIONS \
	  --sjdbOverhang $READLEN

