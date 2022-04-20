# This code is for RSEM quantification by running array jobs.

#PBS -J 1-16


### Define variables ###
working_dir=<PATH/TO/WORKING_DIRECTORY/>

# Only used for naming
gtf_type="Default"
genome=<PATH/TO/GENOME_FILE>
annotation=<PATH/TO/ANNOTATOIN_FILE>
STAR_aligned=<PATH/TO/ALIGNED_READ_DIRECTORY/>
RSEM_path=<PATH/TO/RSEM_OUTPUT_DIRECTORY/>
RSEM_out=$RSEM_path/$gtf_type
RSEM_ref=$RSEM_path/ref


echo "This is the raw path for RSEM "$RSEM_path
FILES=($(ls $STAR_aligned/*/*Aligned.toTranscriptome.out.bam))
f=${FILES[$PBS_ARRAY_INDEX-1]}

# Check so that the basename extension makes sense
bname=$(basename "$f" _Aligned.toTranscriptome.out.bam)

RSEM_commands="--bam --no-bam-output \
-p 30 --paired-end --forward-prob 0 --calc-pme \
$f \
$RSEM_ref/$gtf_type \
$RSEM_out/Quant_$bname "



echo -e "RSEM Counting for file: " $f
echo -e "RSEM input commands: "$RSEM_commands

rsem-calculate-expression $RSEM_commands >& "$RSEM_out/$bname".log

