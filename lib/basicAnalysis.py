'''
Created on Oct 9, 2012

@author: chris
'''

from Bio import SeqIO
import numpy as np
import time
import os
import sys
import matplotlib.pyplot as plt
from utils import pklsave

from utils import Cycler

def calc_propN_meanphred(infiles = None, filepattern = '', data_inpath = '', 
                         out_filename='stats', plothist=False, png_filename=''):
    ''' Calculate the proportion of Ns and the mean Phred scores per read for
    the files given. 
    the files given.     
    
    OPTIONS
    OPTIONS    
    - out_filename --> name of .npy file for output
    - plothist -->  if ture plots the histogram to screen.
    - png_filename  -->  saves figure of histogram to file.
    - 
    No longer returns stats obj, as too much mem usage. Now just writes to .npy file 
    
    '''

    RecCycler = Cycler(infiles = infiles, filepattern = filepattern, data_inpath = data_inpath)
    
    print '\nCalculating Proportion of Ns and Mean Phred Score per read.\n'
    
    # Define vars and outputs
    numfiles = RecCycler.numfiles
    output_list = [0] * numfiles
    rec_counts = [0] * numfiles

    # Find the number of records in each file.
    toc = time.time()
    cum_t = 0
    for seqrecgen in RecCycler.seqfilegen:
    
        filename = RecCycler.curfilename
        filenum = RecCycler.curfilenum
        
        # Check filename is a .idx file extention or that one exists for it.
        # If not, create one.
        if os.path.isfile(filename.split('.')[0] + '.idx'):
            idx_filename = filename.split('.')[0] + '.idx'
            S = SeqIO.index_db(idx_filename, format='fastq')
                 
        elif not os.path.isfile(filename.split('.')[0] + '.idx'):
            
            print 'Index file for {0} does not exist. Creating.....'.format(filename)
            tak = time.time()
            idx_filename = filename.split('.')[0] + '.idx'
            S = SeqIO.index_db(idx_filename, filename , 'fastq')
            print '{0} written successfully'.format(idx_filename)
            idx_t = time.time() - tak
            print 'Finished Indexing to {0}\n after {1}\n'.format(idx_filename, 
                                time.strftime('%H:%M:%S', time.gmtime(idx_t)))
            
        # count number of Records 
        rec_counts[filenum] = len(S)
        del S
        
        output_list[filenum] = {'propN'  :  np.zeros(rec_counts[filenum]),
                        'meanPhred' : np.zeros(rec_counts[filenum])}

        # Cycle through all records 
        for num, seqRec in enumerate(seqrecgen):
            # Proportion of Ns
            output_list[filenum]['propN'][num] = float(seqRec.seq.count('N')) / len(seqRec)
            # Mean Phred Score
            output_list[filenum]['meanPhred'][num] = np.array(
                            seqRec.letter_annotations['phred_quality']).mean()
        
        loop_t = time.time() - toc - cum_t
        cum_t += loop_t
        print 'Finished {0} after {1}'.format(filename, time.strftime('%H:%M:%S', time.gmtime(loop_t)))

    total_records = np.sum(rec_counts)
    total_stats = {'propN':  np.zeros(total_records)  ,
                  'meanPhred': np.zeros(total_records) }
    start_idx = 0
    end_idx = 0
    for file_idx, stat_obj in enumerate(output_list):
        end_idx += rec_counts[file_idx]
        total_stats['propN'][start_idx:end_idx] = stat_obj['propN']
        total_stats['meanPhred'][start_idx:end_idx] = stat_obj['meanPhred']
        start_idx += rec_counts[file_idx]

    # Saving data
    if out_filename:
        print 'Saving Data......'
        np.save(out_filename + '_propN', total_stats['propN'])
        np.save(out_filename + '_meanPhred', total_stats['meanPhred'])
        print 'DONE! Saved .npy files to {0}'.format(os.getcwd())

    if plothist or png_filename:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.hist(total_stats['propN'])
        ax.set_title('''Frequency Count for the Proportion of 'N's''')
        ax.set_ylable('Frequency')
        ax.set_ylable('''Fraction of 'N's in read''')
        if png_filename:
            fig.savefig(png_filename + '_propN.png')
            fig.close()
        else:
            fig.show()
            
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.hist(total_stats['meanPhred'])
        ax.set_title('''Frequency Count for the Mead Phred Score''')
        ax.set_ylable('Frequency')
        ax.set_ylable('''Mead Phred Score in read''')
        if png_filename:
            fig.savefig(png_filename + '_meanPhred.png')
            fig.close()
        else:
            fig.show()
            
    total_t = time.time() - toc
    print '\nProcessed {0} files in {1}'.format(numfiles, time.strftime('%H:%M:%S', time.gmtime(total_t)))

def calc_readlengths(infiles = None, filepattern = '', data_inpath = '', 
                     out_filename='stats', plothist=False, png_filename=''):
    ''' Calculate the readlengths per read for
    the files given. 
    
    OPTIONS
    - out_filename --> name of .npy file for output
    - plothist -->  if ture plots the histogram to screen.
    - png_filename  -->  saves figure of histogram to file.
    - 
    No longer returns stats obj, as too much mem usage. Now just writes to .npy file 
    '''     
    
    #Generator for Sequence Record files
    RecCycler = Cycler(infiles = infiles, filepattern = filepattern, data_inpath = data_inpath)
    
    print '\nCalculating length per read ...\n'
    
    # Define vars and outputs
    numfiles = RecCycler.numfiles
    output_list = [0] * numfiles
    rec_counts = [0] * numfiles

    # Find the number of records in each file.
    toc = time.time()
    cum_t = 0
    for seqrecgen in RecCycler.seqfilegen:
    
        filename = RecCycler.curfilename
        filenum = RecCycler.curfilenum
        
        # Check filename is a .idx file extention or that one exists for it.
        # If not, create one.
        if os.path.isfile(filename.split('.')[0] + '.idx'):
            idx_filename = filename.split('.')[0] + '.idx'
            S = SeqIO.index_db(idx_filename, format='fastq')
                 
        elif not os.path.isfile(filename.split('.')[0] + '.idx'):
            
            print 'Index file for {0} does not exist. Creating.....'.format(filename)
            tak = time.time()
            idx_filename = filename.split('.')[0] + '.idx'
            S = SeqIO.index_db(idx_filename, filename , 'fastq')
            print '{0} written successfully'.format(idx_filename)
            idx_t = time.time() - tak
            print 'Finished Indexing to {0}\n after {1}\n'.format(idx_filename, 
                                time.strftime('%H:%M:%S', time.gmtime(idx_t)))
            
        # count number of Records 
        rec_counts[filenum] = len(S)
        del S
        
        output_list[filenum] = {'length'  :  np.zeros(rec_counts[filenum]) }
        
        for num, seqR in enumerate(seqrecgen):
            output_list[filenum]['length'][num] = len(seqR.seq)
        
        loop_t = time.time() - toc - cum_t
        cum_t += loop_t
        print 'Finished {0} after {1}'.format(filename, time.strftime('%H:%M:%S', time.gmtime(loop_t)))

    total_records = np.sum(rec_counts)
    total_stats = {'length':  np.zeros(total_records)}
    start_idx = 0
    end_idx = 0
    for file_idx, statObj in enumerate(output_list):
        end_idx += rec_counts[file_idx]
        total_stats['length'][start_idx:end_idx] = statObj['length']
        start_idx += rec_counts[file_idx]

    # Saving data
    if out_filename:
        print 'Saving Data......'
        np.save(out_filename + '_propN', total_stats['propN'])
        np.save(out_filename + '_meanPhred', total_stats['meanPhred'])
        print 'DONE! Saved .npy files to {0}'.format(os.getcwd())

    if plothist or png_filename:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.hist(total_stats['length'])
        ax.set_title('''Frequency Count for the Read Lengths''')
        ax.set_ylable('Frequency')
        ax.set_ylable('''Read Length''')
        if png_filename:
            fig.savefig(png_filename + '_readlength.png')
            fig.close()
        else:
            fig.show()

    total_t = time.time() - toc
    print 'Processed {0} files in {1}'.format(RecCycler.numfiles, 
                                              time.strftime('%H:%M:%S', time.gmtime(total_t)))

def calc_phredperbase_boxplot(infiles = None, filepattern = '', data_inpath = '', 
                              saveprefix = '', png_filename=''):
    ''' Find the median, upper and lower quartile for the Phred score per base 
    
    Returns the stats and the counter dictionary. 
    
    Counter dictionary may become standard way to store mass Phred/seq bases data.  '''
    from collections import Counter

    RecCycler = Cycler(infiles = infiles, filepattern = filepattern, data_inpath = data_inpath)
    
    print '\nCalculating Box plot stats of phred scores per base position.\n'
    
    # Define vars and outputs
    numfiles = RecCycler.numfiles

    toc = time.time()
    cum_t = 0
    
    counter_list = [0] * 101
    for i in range(len(counter_list)):
        counter_list[i] = Counter()
    
    for seqrecgen in RecCycler.seqfilegen:
        
        filename = RecCycler.curfilename
        filenum = RecCycler.curfilenum
        
        for rec in seqrecgen:
            for basenum, phred in enumerate(rec.letter_annotations['phred_quality']):
                counter_list[basenum][phred] += 1
                
        loop_t = time.time() - toc - cum_t
        cum_t += loop_t
        print 'Finished {0} \nfile {1} of {2} after {3}'.format(filename, filenum, numfiles, time.strftime('%H:%M:%S', time.gmtime(loop_t))) 

    # Calculate min, max Q1, Q2, Median and Average
    stats = getStatsFT(counter_list)

    total_t = time.time() - toc
    print 'Processed all files in {0}'.format(time.strftime('%H:%M:%S', time.gmtime(total_t)))
    
    pklfilename = data_inpath.split('/')[-1]
            
    pklsave(counter_list, '_'.join([pklfilename, saveprefix , 'phredCount']))
    np.save( '_'.join([pklfilename, saveprefix , 'phredStats.npy']) , stats)
    
    plotFTstats(stats, png_filename)
    
    return stats, counter_list

def getStatsFT(counter_list):
    ''' Calculate min, max, Q1, Q2, Median and Average for a list of frequency table of data 
    
    Returns a structured array the length of the input list '''
    
    # Input check    
    if type(counter_list) != list:
        counter_list = [counter_list] 

    # Define data type    
    dt = np.dtype([ ('min', 'u2'), 
                ('Q1', 'f2'), 
                ('median','f2'), 
                ('ave', 'f2'), 
                ('Q2', 'f2'), 
                ('max', 'u2') ])
    
    numCounters = len(counter_list)
    stats = np.zeros(numCounters, dtype = dt)
    
    # Calculate min, max Q1, Q2, Median and Average
    for i, cnt in enumerate(counter_list):
        keys = cnt.keys()
        vals = percentileFT(cnt, [0.25, 0.5, 0.75])
        stats[i]['Q1'] = vals[0]
        stats[i]['median'] = vals[1]
        stats[i]['Q2'] = vals[2]
        stats[i]['min'] = min(keys)
        stats[i]['max'] = max(keys)
        ave = float(0)
        for k in keys:
            ave = ave + (k * cnt[k])
        ave = ave / sum(cnt.values()) 
        stats[i]['ave'] = ave

    return stats

def percentileFT(hashCount, percents):
    ''' Calculate the required percentile from the hash count or frequency table given. 
    
    Can take a list of values for percentile, and will return an output for each'''
        
    import math
    
    if type(percents) != list:
        percents = [percents]
    
    N = sum(hashCount.values())
    sorted_keys = hashCount.keys()
    sorted_keys.sort()
    
    outValues = np.zeros(len(percents), dtype = np.uint8)
    
    for i, per in enumerate(percents):
         
        k = (N-1) * per + 1 # f and c must count from 1 not zero
        f = math.floor(k)  
        c = math.ceil(k) 
        
        cum_count = 0
        if f == c: # No remainder so no need to interpolate
            for key in sorted_keys:    
                cum_count += hashCount[key]
                if cum_count >= k:
                    # Found the right one 
                    break
            outValues[i] = key
    
        else: # Must now check to see if we need to interpolate
            floor_key = 0
            ceil_key = 0
            found_floor_val = 0
            for key in sorted_keys:    
                cum_count += hashCount[key]
                if not found_floor_val:
                    if cum_count >= f:  
                        # Found the right one 
                        floor_key = key
                        found_floor_val = 1
                if cum_count >= c: 
                    ceil_key = key
                    break
            
            if floor_key == ceil_key:
                outValues[i] = key
            else:
                # Interpolate
                val1 = floor_key * (c-k)
                val2 = ceil_key * (k-f)
                outValues[i] = val1 + val2     
    
    if len(outValues) == 1:
        return outValues[0]
    else:                        
        return outValues
               
def plotFTstats(stats, png_filename=''):
    ''' Plot the min, max, Q1, Q2, Median and Average from a stats record structure '''
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(stats['max'], 'k-', label='max')
    ax.plot(stats['Q2'], 'b.', label='Q2')
    ax.plot(stats['median'], 'r.', label='median')
    ax.plot(stats['ave'], 'g-', label='average')
    ax.plot(stats['Q1'], 'b.', label='Q1')
    ax.plot(stats['min'], 'k-', label='min')
    ax.set_ylabel('Phred Score') 
    ax.set_xlabel('Base Position')
    ax.legend(loc=3)
    ax.set_title('Statistics for Phred Scores per Base Position')
    if png_filename:
        fig.savefig(png_filename + '_readlength.png')
    else:
        fig.show()
     
def getMeanPhredPerBase(infiles = None, filepattern = '', data_inpath = ''):
    ''' Find the mean Phred score per base '''
        
    #Generator for Sequence Records
    RecCycler = Cycler(infiles = infiles, filepattern = filepattern, data_inpath = data_inpath)
    
    print '\nCalculating mean phred score per base position.\n'

    # Make Data Structure
    running_aves = np.zeros(101)
    
    toc = time.time()
    cum_t = 0
    count = 1
    for num, seqRecs in enumerate(RecCycler):
        for rec in seqRecs.itervalues():
            xi = np.array(rec.letter_annotations['phred_quality'])
            running_aves = running_aves + ((xi - running_aves) / (count)) 
            count += 1
        loop_t = time.time() - toc - cum_t
        cum_t += loop_t
        print 'Finished file {0} after {1}'.format(num, time.strftime('%H:%M:%S', time.gmtime(loop_t))) 
    
    total_t = time.time() - toc
    print 'Processed all files in {0}'.format(time.strftime('%H:%M:%S', time.gmtime(total_t)))
    
    return running_aves

if __name__ == '__main__':
    pass


