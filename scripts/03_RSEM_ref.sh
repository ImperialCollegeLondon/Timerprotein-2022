# This code is for RSEM reference creation
### Define variables ###
working_dir=<PATH/TO/WORKING_DIRECTORY/>

# Only used for naming
gtf_type="Default"
genome=<PATH/TO/GENOME_FILE>
annotation=<PATH/TO/ANNOTATOIN_FILE>
RSEM_path=<PATH/TO/RSEM_OUTPUT_DIRECTORY/>

RSEM_out=$RSEM_path/$gtf_type
RSEM_ref=$RSEM_path/ref

### Create/confirm directories ###
mkdir $RSEM_path $RSEM_out $RSEM_ref


echo "Creating RSEM index to $RSEM_ref"
RSEM_index_commands=" --gtf $annotation \
$genome \
$RSEM_ref/$gtf_type"
echo "Using RSEM prepare reference commands \
$RSEM_index_commands"
rsem-prepare-reference $RSEM_index_commands

