# This code is for read alignment using STAR bu running array jobs.

#PBS -J 1-16

### Define variables ###
working_dir=<PATH/TO/WORKING_DIRECTORY/>

# Only used for naming
gtf_type="Default"
data=<PATH/TO/DATA/>
ref=<PATH/TO/REF_OUTPUT_DIRECTORY/>
STAR_index=<PATH/TO/REF_OUTPUT_DIRECTORY/INDEX_NAME>
genome=<PATH/TO/GENOME_FILE>
annotation=<PATH/TO/ANNOTATOIN_FILE>
STAR_aligned=<PATH/TO/OUTPUT_DIRECTORY/>

### Create/confirm directories ###
mkdir $STAR_aligned

### Star mapping ###

echo -e "\n\n"

FILES=($(ls $data/*_1.fq.gz))
f=${FILES[$PBS_ARRAY_INDEX-1]} 

bname=$(basename "$f" _1.fq.gz)
dname=$(dirname "$f")

read1=$dname/$bname"_1.fq.gz"
read2=$dname/$bname"_2.fq.gz"

echo "Input path for read 1 to STAR: " $read1
echo -e "\n"
echo "Input path for read 2 to STAR: " $read2
echo -e "\n\n"

### Make directories for each read
reads_out=$STAR_aligned/STARmapped_$bname
mkdir $reads_out

### Star commands for aliging the reads to a reference. zcat - to uncompress .gz files.

STAR_commands=" --runMode alignReads \
--genomeDir $STAR_index/ \
--runThreadN 30 \
--outTmpKeep None \
--readFilesCommand zcat \
--readFilesIn $read1 $read2 \
--sjdbGTFfile $annotation \
--outReadsUnmapped Fastx \
--outFileNamePrefix $reads_out/$bname"_"$gtf_type"_" \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts"

echo -e "Running STAR with commands:\n"
echo $STAR_commands
echo -e "\n\n"

STAR $STAR_commands

