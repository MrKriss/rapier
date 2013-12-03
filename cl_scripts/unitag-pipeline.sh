#!/bin/bash
#
# Example of how to run the various stages of the Salmon RAD-seq pipeline. 

# Summary specs of bens machine: 
#     198 GB RAM
#     32 core Intel(R) Xeon(R) CPU E5-2470 0 @ 2.30GHz


# File names and path locations 
RDFS_PATH="/home/musselle/salmon"     # home/musselle is my home on bens machine and the rdf storage is mounted with the name 'salmon'
PROJECT_ROOT="$RDFS_PATH/chris"       # Root directory for project on RDFS Storage. 


outpath="output"                                 # Root directory to write output to.  
raw_data="testdata/testset_100k.fastq"           # Path after 'root' for raw data files. Can accept a glob to pass multiple files 
barcodes="testdata/barcodes/lane3_newzebra.txt"  # Path after 'root' for barcodes associated with raw input files
cdhit_path="lib/bin/"                            # Location of binaries for cd-hit-est clustering program. cd-hit-v4.6.1 is included in repository
database="$outpath/test100k.db"                     # Output name for central database.
fasta_filepath="$outpath/testset100k_precluster.fasta"  # Location to write extracted fasta data for input to cd-hit clustering program. 
cluster_filepath="$outpath/testset100k_cluster"    # Location and file name to write CDHIT clustering output to. 


######################
# Data Preprocessing #
######################
# Assuming that stacks is installed and available on the system path. 
# The raw-data files in $PROJECT_ROOT/raw-data have already had the dummy barcodes added, and are located under $PROJECT_ROOT/barcodes

echo "\nAbout to Run Preprocessing Steps"

process_radtags -p $PROJECT_ROOT/raw-data \                      
				-b $PROJECT_ROOT/barcodes/barcodes_only.txt  \
				-o $PROJECT_ROOT/processed-data/ \
				-e sbfI \    # Specify enzyme for cut site
				-q -c \      # Filter based on phred quality (using sliding window) and discard reads with uncalled bases
				-r \         # Correct RAD tags 
				--index_null \   # Dummy barcode is in sequence header not in sequence itself. 
				
echo "\nPreprocessing Steps Complete"				 
echo "\nAbout to Construct Unitag Reference Sequence"

# Parameters for Unitag Reference Construction

# Experiement with Multiple min max values using parallel execution 
MIN_DEPTHS=2 5 10 15 20
MAX_DEPTHS=500 

parallel "~/bin/make_unitag.py -i $PROJECT_ROOT/processed-data/sample_GCAGGC.fq \
                               -o $PROJECT_ROOT/processed-data/unitagref-m{1}-M{2}.fq \
                               -m {1} -M {2}" ::: $MIN_DEPTHS ::: $MAX_DEPTHS
# Statistics for all runs are logged in unitag-logfile.log

echo "\nBuilding Unitags Complete"		

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