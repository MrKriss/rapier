'''
Created on 24 Oct 2012

@author: musselle
'''

import os 

import numpy as np 
import matplotlib.pyplot as plt 

from collections import Counter, defaultdict

def cluster_summary_plot(infile, cluster_length_bins=None, mincutoff=10, bins=5000, report=True, plot_hist=True):
    ''' Takes cluster file output by CD-Hit and produces a histergram plot
    
    cluster_length_bins = a list of tuples dictating the left and right most ranges
    of the sequences lengths for which the Histograms are plotted. e.g.
    
    cluster_length_bins = [ (0,30), (50,), (51,np.inf) ] 
    
    would plot three histograms. One for sequences of length 0 to 30; one for sequences 
    of length 50; and one for those 51 upwards.    
    
    '''

    if not infile.endswith('.clstr'):
        infile = infile + '.clstr'
        
    if cluster_length_bins is not None:
        # Create counters for each bin in cluster_length_bins
        for tup in cluster_length_bins:   
            name = 'Length-' + str(tup)
            vars()[name] = Counter()

    # Data structure to store cluster size info
    # Default Dictionary of Counter dictionaries!
    # ds = { 'cluster_size' : Counter(length of representative sequences)  }
    # empty keys of ds are initialised with a Counter dictionary. 
        
    ds = defaultdict(Counter)

    # Helper vars
    cluster_size_counter = Counter()
    seq_length_counter = Counter()

    seq_in_cl_counter = 0
    rep_length = 0

    print 'Scanning %s ...' % (infile)

    try:
        with open(infile, 'rb')  as cluster_file:   
            
            for line in cluster_file:              
                line = line.strip()
                
                if line.startswith('>'):
                    # This is start of new cluster
                    if seq_in_cl_counter and rep_length:
                        # This number is the size of last cluster
                        # Store result
                        ds[str(seq_in_cl_counter)][str(rep_length)] += 1
                        
                        cluster_size_counter[str(seq_in_cl_counter)] +=1
                        seq_length_counter[str(rep_length)] += 1
                        
                        seq_in_cl_counter = 0
                        rep_length = 0
                elif line.endswith('*'): 
                    # This is the representative sequence for the cluster
                    rep_length = int(line.split()[1].strip('nt,'))
                    seq_in_cl_counter += 1
                else:
                    seq_in_cl_counter += 1
            
            # Got to end of file but still one more cluster to add
            ds[str(seq_in_cl_counter)][str(rep_length)] += 1
            cluster_size_counter[str(seq_in_cl_counter)] += 1
            seq_length_counter[str(rep_length)] += 1
    
    except IOError:
        print "Error: can\'t find file or read data"
    else:
        print "Finished Scanning cluster file."
    
    if report:
        print 'Top 5 Cluster Sizes: ', cluster_size_counter.most_common()[:5]
        print 'Top 5 Sequence Lengths: ', seq_length_counter.most_common()[:5]
    
    # Process ds to extract needed data for histergrams
    if cluster_length_bins is not None:
        # Update counters for each bin in cluster_length_bins
        for cluster_size, seq_len_c in ds.iteritems():

            for tup in cluster_length_bins:   
                name = 'Length-' + str(tup)
                 
                if len(tup) == 1:
                    # single number bin
                    target_key = str(tup[0])
                    if seq_len_c.has_key(target_key):
                        vars()[name][cluster_size] += seq_len_c[target_key]
        
                elif len(tup) == 2:
                    # range bin
                    for seq_length in seq_len_c.iterkeys():
                        if int(seq_length) >= tup[0] and int(seq_length) <= tup[1]:
                            vars()[name][cluster_size] += seq_len_c[seq_length]
    
    #===============================================================================
    # Plot results
    #===============================================================================
    
    # TODO Add a level of input checking so if Counter is empty, an empty histogram is plotted
    
    if plot_hist:
    
        plt.figure()
        if cluster_length_bins is not None:   
            for tup in cluster_length_bins:   
                name = 'Length-' + str(tup)
                
                hist_counter(vars()[name], bins=bins, label=name, range=(mincutoff, 10000))
         
        hist_counter(cluster_size_counter, bins=bins, label='All Seq Lengths', range=(mincutoff, 10000))   
        plt.title("Cluster Size Distribution")
        plt.xlabel("Value")
        plt.ylabel("Frequency")
        plt.legend()
        plt.show()

    return ds, cluster_size_counter, seq_length_counter
    
def hist_counter(counter, **kwargs):
    ''' Construct a histogram from a Counter Dictionary '''
    
    data = np.array(list(counter.elements()), dtype = np.int)

    plt.hist(data, histtype='step', **kwargs)
    plt.title("Cluster Size Distribution")
    plt.xlabel("Value")
    plt.ylabel("Frequency")
#    plt.legend()
    plt.show()
    

def makeHist(files , data_inpath, bins = 50):
    ''' Function to plot histograms of stored data files'''

    os.chdir(data_inpath)
 
    data = [0] * len(files)
    
    for i, f in enumerate(files):
        
        name = 'f' + str(i)
        vars()[name] = np.load(f)
        data[i] = vars()[name] 

    plt.figure()
    plt.hist(data, bins = bins)
    plt.ylabel('Frequency')

if __name__ == '__main__' : 
    
    # input vars
#    infile = '/Users/chris/Dropbox/work/code/popGen/utils/test_seqs_clustered'
    infile = '/Users/chris/Downloads/sb_all_clustered'
#    bins = [(0,99),(100,),(101, np.inf)]
    bins = [(25,),(42,)]
    cluster_summary_plot(infile, cluster_length_bins=bins, mincutoff=10, bins=5000, report=True)


#    data_inpath = '/space/musselle/datasets/gazellesAndZebras'
#    fs1 = ['lane6_meanPhred.npy', 'lane8_meanPhred.npy']
#    fs2 = ['lane6_propN.npy', 'lane8_propN.npy']
#    
#    makeHist(fs1, data_inpath)
#    plt.title('Mean Phred Score')
#    makeHist(fs2, data_inpath)
#    plt.title('Proportion of Ns')
