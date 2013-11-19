#!/bin/bash
#
# Example of how to run the various stages of the Rapier pipe line on the accompanying test datasets in testdata/
#

# File names and path locations 
root="."                                         # Define the root directory where project data is stored  
raw_data="testdata/testset_100k.fastq"           # Path after 'root' for raw data files. Can accept a glob to pass multiple files 
barcodes="testdata/barcodes/lane3_newzebra.txt"  # Path after 'root' for barcodes associated with raw input files
cdhit_path="lib/bin/"                            # Location of binaries for cd-hit-est clustering program. cd-hit-v4.6.1 is included in repository
database="test100k.db"                           # Output name for central database.
cluster_filepath="test"            # 

fasta_filepath="testdata/testset100k_precluster.fasta"  # Location to write extracted fasta data for input to cd-hit clustering program. 
filter_expression=""                                    # SQL expression to select reads to extract. 

cluster_filepath="clusters/new_zebra"

# Prefix for created cluster and members tables in database
table_prefix=""    # default is no prefix. only used if multiple sequence and cluster tables are used pre database

###############
# Main Script #
###############

cd $root

# Filter and clean Reads, then pipe to database 
# -i 
#
#
#
#
filter_and_clean_reads.py -i $raw-data -b $barcodes -f 20 -n 0.1 -c TGCAGG -r 2 -g 2  \
| load_reads.py -b $barcodes -d $database

# Extract reads and save to fasta file
# -i input database
# -o output fasta file 
# -e SQL expression to filter the query which selects the sequences in the database. Default is to export all sequences in database.
#    Basic query is:
#           SELECT seqid, seq, phred FROM seqs INNER JOIN samples ON seqs.sampleId=samples.sampleId 
#                   WHERE <filter_expression> 
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
