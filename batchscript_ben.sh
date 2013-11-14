#!/bin/bash

# File names and path locations
root="~/data/RAD-seq/gazelles-zebras/"
raw_data="raw-data/new-zebra/*"
barcodes="barcodes/lane3_newzebra.txt"
cdhit_path="~/bin/cd-hit-v4.6.1/"
database="new_zebra.db"

fasta_filepath="processed-data/new_zebra.fasta"
filter_expression=""

cluster_filepath="clusters/new_zebra"

# Prefix for created cluster and members tables in database
table_prefix=""

# Main Script
cd $root

# Filter and clean Reads, then pipe to database 
filter_and_clean_reads.py -i $raw-data -b $barcodes -f 20 -n 0.1 -c TGCAGG -r 2 -g 2  \
| load_reads.py -b $barcodes -d $database

# Extract reads and save to fasta file
extract_reads.py -i $database -o $fasta_filepath -e $filter_expression

# Run CDHIT Clustering
# -n kmer length to use
# -s similarity at which to threshold
# -g compare each read with all other cluster seeds and cluster to nearest (instead of first one above threshold)
# -m Limit the memory used, 0 is unlimited.
# -t Number of threads to use.    
nice cluster_CDHIT.py -i $fasta_filepath -o $cluster_filepath -p $cdhit_path -n 8 -s 0.95 -g -m 0 -t 20 

# Import cluster file into database
# min and max are the minimum and maximum cluster sizes to load into database
load_CDHIT_cluster.py -i $cluster_filepath -o $database -p $table_prefix -n --min 100 --max 1000 