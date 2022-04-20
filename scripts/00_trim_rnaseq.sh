# Trimming reads using PBS array indexing

##### Defining the "raw" path to the folders that have data files ######## 
RAWPATH=<PATH/TO/RAWDATA/>



##### Defining the files inside these folders #######

FILES=($(ls $RAWPATH/*/*_1.fastq.gz))
f=${FILES[$PBS_ARRAY_INDEX-1]} 

bname=$(basename "$f" _1.fastq.gz)
dname=$(dirname "$f")

read1=$dname/$bname"_1.fastq.gz"
read2=$dname/$bname"_2.fastq.gz"

cp $read1 <PATH/TO/COPIED_DATA/>
cp $read2 <PATH/TO/COPIED_DATA/>

newdname=<PATH/TO/COPIED_DATA/>
copied_read1=$newdname/$bname"_1.fastq.gz"
copied_read2=$newdname/$bname"_2.fastq.gz"



trim_galore --paired --length 30 --fastqc --fastqc_args "-d . -o <PATH/TO/FASTQC_OUTPUT/>" --stringency 3 -o <PATH/TO/TRIMMED_DATA_OUTPUT/> $copied_read1 $copied_read2


