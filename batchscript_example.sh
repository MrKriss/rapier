#!/bin/bash
#
# Example of how to run the various stages of the Rapier pipe line on the accompanying test datasets in testdata/
#

# File names and path locations 
outpath="output"                                 # Root directory to write output to.  
raw_data="testdata/testset_100k.fastq"           # Path after 'root' for raw data files. Can accept a glob to pass multiple files 
barcodes="testdata/barcodes/lane3_newzebra.txt"  # Path after 'root' for barcodes associated with raw input files
cdhit_path="lib/bin/"                            # Location of binaries for cd-hit-est clustering program. cd-hit-v4.6.1 is included in repository
database="$outpath/test100k.db"                     # Output name for central database.
fasta_filepath="$outpath/testset100k_precluster.fasta"  # Location to write extracted fasta data for input to cd-hit clustering program. 
cluster_filepath="$outpath/testset100k_cluster"    # Location and file name to write CDHIT clustering output to. 

###############
# Main Script #
###############

mkdir -p $outpath

# Filter and clean Reads, then pipe to database 
# -i Input files
# -b barcodes for those input files 
# -f threshold for the mean pred filter 
# -n threshold for the N proportion filter
# -c restriction enzymne cutsite
# -r restriction enzyme overhang in number of bases
# -g max number of errors to correct in  MID tag and cutsite, compared to true values.
# -d file path to central SQL database

python cl_scripts/filter_and_clean_reads.py -i $raw_data -b $barcodes -f 20 -n 0.1 -c TGCAGG -r 2 -g 2 -p $outpath \
| python cl_scripts/load_reads.py -b $barcodes -d $database

# Extract reads and save to fasta file
# -i input database
# -o output fasta file 
python cl_scripts/extract_reads.py -i $database -o $fasta_filepath 

# Run CDHIT Clustering
# -n kmer length to use for the filter
# -s fraction of similarity at which to threshold cluster membership
# -g compare each read with all other cluster seeds and cluster to nearest (instead of first one above threshold)
# -m Limit the memory used, 0 is unlimited.
# -t Number of threads to use.    
nice python cl_scripts/cluster_CDHIT.py -i $fasta_filepath -o $cluster_filepath -p $cdhit_path -n 8 -s 0.95 -g -m 0 -t 2 

# Import cluster file into database
# min and max are the minimum and maximum cluster sizes to load into database. Typically this would be roughly --min 100 --max 1000
# so as to save time not loading a large number of small clusters. Though as this is a demo script, all clusters are loaded
python cl_scripts/load_CDHIT_cluster.py -i $cluster_filepath -o $database -n --min 0 --max 1000 
