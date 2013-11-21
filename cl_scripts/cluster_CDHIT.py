#! /usr/bin/env python 

'''
Created on 17 Jul 2013

@author: musselle
'''

import os
import sys
import time
import argparse
import shlex
import subprocess
from collections import Counter, defaultdict

import numpy as np

import _addpaths

from fileIO import inputfile_check, outputfile_check


class CDHIT_ClusteringClass(object):
    ''' Class to act as a holder of all wrappers for all CDHIT clustering methods 
    '''
    def __init__(self, args):

        self.args = args

    def run(self, infile_handle):
        ''' Run CD-HIT in parallel on list of fasta files. Each file is clustered seperately.
        
        Other flags used:
        -d 0   --> No limit on description written to cluster file (goes to first space in seq ID). 
        -r 0   --> DO Only +/+ and not -/+ alignment comparisons as reads are done in both directions but on different strands. 
        -s 0.8 --> If shorter sequence is less than 80% of the representative sequence, dont cluster. 
        
        infile_handle -- Takes file object or string of the file path/filename
        ouput -- Defines output location and file name.
        
        Writes stdout to console.
        Counter dictionary and summary logfile are generated after each run. 
        
        '''

        # input checks    
        infile_handle = inputfile_check(infile_handle)
        
        logfile_path = os.path.join(os.path.split(self.args.output)[0], 'clusterfile.log')
        infile_path = os.path.abspath(infile_handle.name)
            
        logfile = outputfile_check(logfile_path)

        # setup internal vars        
        start_time = time.time()
    
        #=======================================================================
        # Run CDHIT
        #=======================================================================
    
        cmd = ('cd-hit-est -i {0} -o {1} -c {2} -n {3} -d 0 -r 0 -s 0.8 -M {4} '
            '-T {5}').format(infile_path, self.args.output, 
                             self.args.similarity, 
                             self.args.n_gram, 
                             self.args.maxmemory, 
                             self.args.threads)   
    
        if self.args.maskN:
            cmd = cmd + ' -mask N'
        if self.args.allvall:
            cmd = cmd + ' -g 1'
        
        cdhitpath = os.path.expanduser(self.args.cdhitpath)
        
        # Spawn Process to run CD-HIT
        subprocess.check_call(shlex.split(os.path.join(cdhitpath, cmd)))
        
        finish_time = time.time()
        
        
        #=======================================================================
        # Generate a summary log file 
        #=======================================================================
        
        # Get cluster size summary counter 
        total_counter, by_seqlen_counter = self.cluster_summary_counter(infile_path=self.args.output,
                                                                        mode='both', report=True)    
        st_idx = cmd.find('-c ')
        CDHIT_parameters = cmd[st_idx:]
        
        # Write summary logfile 
        with logfile as f:
            program_name = os.path.join(self.args.cdhitpath, cmd).split(' -i ')[0]
            f.write('=========================================================\n')
            f.write('Program     : {0}\n'.format(program_name))
            f.write('Input File  : {0}\n'.format(infile_path))
            f.write('Output File : {0}\n'.format(self.args.output))
            f.write('Commands    : {0}\n'.format(CDHIT_parameters))
            f.write('\n')
            f.write('Started     : {0}\n'.format(time.strftime('%a, %d %b %Y, %H:%M:%S', 
                                                    time.gmtime(start_time))))
            f.write('=========================================================\n')
            f.write('\n')
            f.write('                       Report Log\n')
            f.write('---------------------------------------------------------\n')
            
            reads_per_cluster = {key: int(key)*value for key, value in total_counter.iteritems()}
            total_reads = sum(reads_per_cluster.values())
            total_clusters = sum(total_counter.values())
            f.write('Total number of reads     : {0}\n'.format(total_reads))
            f.write('Total number of clusters  : {0}\n'.format(total_clusters))
            read_lengths = [int(key) for key in by_seqlen_counter.keys()]
            f.write('Read length Min and Max    : {0} and {1}\n'.format(min(read_lengths), max(read_lengths)))
            f.write('Time taken                 : {0}\n'.format(time.strftime('%H:%M:%S', 
                                                    time.gmtime(finish_time - start_time))))
            f.write('\n')
            f.write('Top 20 Percentage Reads per cluster \n')
            f.write('---------------------------------------------------------\n')
            f.write('Cluster Size    No. Clusters    Total Reads         %    \n')
            f.write('---------------------------------------------------------\n')
            top_reads_per_cluster = sorted(reads_per_cluster.iteritems(), 
                                           key=lambda tup: int(tup[1]), reverse=True)[:20]
            for tup in top_reads_per_cluster:
                if total_reads == 0:
                    perc = 0.0
                else:
                    perc = float(tup[1]) / total_reads
                
                f.write("{clust_size: <16}{num_clust: <16}{total_reads: <18d}{percentage:.2%}\n".format(
                      clust_size=tup[0], num_clust=total_counter[tup[0]], total_reads=tup[1], 
                      percentage=perc))

        cluster_file_handle = open(self.args.output, 'rb')
        
        return cluster_file_handle, total_counter


    def cluster_summary_counter(self, infile_path, mode='total', report=True):
        ''' Takes cluster file output by CD-Hit and produces two Counter for the 
        
        modes:
        counter_per_sequence_length = { 'sequence_length' : Counter(cluster_size) }
        total = Counter(cluster_sizes_for _all_sequences)
        
        '''
    
        path, infile = os.path.split(infile_path)
    
        if not infile.endswith('.clstr'):
            infile = infile + '.clstr'
        
        # Data structure to store cluster size info is a DefaultDictionary of Counter dictionaries.
        # ds = { 'seq_len' : Counter(cluster_size)  }
        # empty keys of ds are initialised with a Counter dictionary. 
        
        ds = defaultdict(Counter)
    
        seq_in_cluster = 0
        rep_length = 0
    
        print 'Generating cluster summary for  %s ...' % (infile)
    
        try:
            with open(os.path.join(path,infile), 'rb')  as cluster_file:   
                
                for line in cluster_file:              
                    line = line.strip()
                    
                    if line.startswith('>'):
                        # This is start of new cluster
                        if seq_in_cluster and rep_length:
                            # This number is the size of last cluster
                            # Store result
                            ds[str(rep_length)][str(seq_in_cluster)] += 1
                            seq_in_cluster = 0
                            rep_length = 0
                            
                    elif line.endswith('*'): 
                        # This is the representative sequence for the cluster
                        rep_length = int(line.split()[1].strip('nt,'))
                        seq_in_cluster += 1
                    else:
                        seq_in_cluster += 1
                
                # Got to end of file but still one more cluster to add
                ds[str(rep_length)][str(seq_in_cluster)] += 1
        
        except IOError:
            print "Error: can\'t find file or read data"
        else:
            print "Finished Scanning cluster file."
        
        # Construct total cluster size counter 
        total_cluster_size_counter = Counter()
        for v in ds.itervalues():
            total_cluster_size_counter.update(v)
            
        # Construct representative sequence length counter 
        seq_len_counter = Counter()
        for k, v in ds.iteritems():
            seq_len_counter[k] += sum(v.values()) 
        
        if report:
            print 'Top 5 Cluster Sizes: ', total_cluster_size_counter.most_common()[:5]
            print 'Top 5 Sequence Lengths: ', seq_len_counter.most_common()[:5]
        
        # Decide what to output    
        if mode == 'total':
            return total_cluster_size_counter
        elif mode == 'by_seqlen':
            return ds
        elif mode =='both':
            return total_cluster_size_counter, ds
        
    
    def hist_counters(self, counters, labels=None, **kwargs):
        ''' Construct a series of histograms from a list of Counter Dictionarys '''
        
        import matplotlib.pyplot as plt
        
        if type(counters) is not list and type(counters) is not tuple:
            counters = [counters]
        
        if labels is not None:
            assert len(labels) == len(counters), "Number of labels must match number of counters."
        
        
        for i in range(len(counters)):
            data = np.array(list(counters[i].elements()), dtype = np.int)
    
            if labels:
                plt.hist(data, histtype='step', label=labels[i], **kwargs)
            else:
                
                plt.hist(data, histtype='step', label='Counter-'+str(i), **kwargs)
            plt.title("Cluster Size Distribution")
            plt.xlabel("Value")
            plt.ylabel("Frequency")

        plt.legend()
        plt.show()


if __name__ == '__main__':
    
    # Parse arguments
    toc = time.time()

    parser = argparse.ArgumentParser(description='Run Simple clustering on reads in the database.')
    
    # IO parameters
    parser.add_argument('-i',  dest='input', 
                        required=True,
                        help='Fasta file where reads are stored (/path/filename)')
    parser.add_argument('-o',  dest='output', default='clusterfile',
                        help='Filename for output clusters (/path/filename).')
    
    # Clustering parameters
    parser.add_argument('-s',  dest='similarity', required=True, type=float,
                        help='Threshold for percentage similarity between clusters.')
    parser.add_argument('-n',  dest='n_gram', required=True, type=int,
                        help='N gram pattern to use for similarity calculation between reads.')
    parser.add_argument('-m',  dest='maxmemory', default=0, type=int,
                        help='Maximum memory to use for the hash table. Default = 0 (unlimited)')
    parser.add_argument('-t',  dest='threads', default=1, type=int,
                        help='Number of threads to use for the clustering. Default = 1')
    parser.add_argument('-g',  dest='allvall', action = 'store_true', 
                        help='Whether to compare a new read to all previous seeds when growing clusters. Default = False')
    
    parser.add_argument('-p',  dest='cdhitpath', default='~/bin/cdhit',
                        help='Path to where cd-hit is installed. Default is ~/bin/cdhit/')
    parser.add_argument('--maskN',  dest='maskN', action = 'store_true', 
                        help='Whether to mask any Ns in the DNA sequence. Default = False')
    
    print sys.argv
    args = parser.parse_args()
    
    # check if cdhit is in default path 
    if not os.path.exists(args.cdhitpath):
        raise Exception('CDHIT program not found on the default path of {0}'.format(args.cdhitpath))
    
    # Setup and Run clustering 
    clustering = CDHIT_ClusteringClass(args)
    outfile_handle, total_cluster_counter = clustering.run(args.input)
    
    total_records = 0
    for clustsize, number in total_cluster_counter.iteritems():
        total_records += int(clustsize) * number
        
    total_clusters =  sum(total_cluster_counter.values())
    
    total_t = time.time() - toc    
    print >> sys.stderr, 'Clustered {0} records into {1} clusters in {2}'.format(
              total_records, total_clusters,
              time.strftime('%H:%M:%S', time.gmtime(total_t)))
