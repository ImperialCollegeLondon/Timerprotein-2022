# This is for creating an index in STAR
### Define variables ###
working_dir=<PATH/TO/WORKING_DIRECTORY/>

# Only used for naming
gtf_type="Default"
data=<PATH/TO/DATA/>
ref=<PATH/TO/REF_OUTPUT_DIRECTORY/>
STAR_index=<PATH/TO/REF_OUTPUT_DIRECTORY/INDEX_NAME>
genome=<PATH/TO/GENOME_FILE>

annotation=<PATH/TO/ANNOTATOIN_FILE>

## Make sure directories exist ##
mkdir $ref $STAR_index $STAR_index/$gtf_type

### STAR index ###


STARindex_COMMANDS=" --runMode genomeGenerate \
--runThreadN 30 \
--outFileNamePrefix $STAR_index/$gtf_type/$gtf_type \
--genomeDir $STAR_index/$gtf_type/ \
--genomeFastaFiles $genome \
--sjdbGTFfile $annotation "


echo -e "Creating STAR index with commands:\n"
echo $STARindex_COMMANDS
echo -e "\n\n"

STAR $STARindex_COMMANDS

