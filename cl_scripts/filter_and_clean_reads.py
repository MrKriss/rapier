#! /usr/bin/env python 

'''
Created on 3 Jul 2013

@author: musselle
'''
import os
import sys
import time
import glob
import gzip
import argparse
import StringIO
import pprint
import editdist
from collections import Counter

import numpy as np
from Bio import SeqIO, bgzf

from release.lib.fileIO import SeqRecCycler


class RecordPreprocessor(object):
    
    def __init__(self, args):
        
        # Set counting variables 
        self.filterfail_counter = Counter()
        self.total_read_count = 0
        
        self.skipped_count = 0
        self.MIDtag_corrected_count = 0
        self.cutsite_corrected_count = 0
        self.read_corrected_count = 0
        
        self.total_read_passes = 0
            
        # Load barcodes 
        if 'barcodes' in args:
            # input checks
            
            if type(args.barcodes) is str:
            
                # Check for path head included in data_files
                if not os.path.split(args.barcodes)[0]:
                    # Append the abs path to current directory 
                    current_dir = os.getcwd()
                    args.barcodes = os.path.join(current_dir, args.barcodes)
                
                # Glob the file pattern 
                barcode_files = glob.glob(args.barcodes)
                assert barcode_files, "No barcode files found at destination"
                barcode_files.sort()
            
            elif type(args.barcodes) is list or type(args.barcodes) is tuple:
            
                if not os.path.split(args.barcodes[0])[0]:
                    current_dir = os.getcwd()
                    args.barcodes = [os.path.join(current_dir, x) for x in args.barcodes]
            else:
                raise Exception('Invalid entry for data_files.')

            # Store barcode dictionary            
            MIDs = []
            individuals = []
            for barcode_file in args.barcodes:
                with open(barcode_file, 'rb') as f: 
                    for line in f:
                        parts = line.strip().split()
                        MIDs.append(parts[0])
                        individuals.append(parts[1])
                    diff = len(MIDs) - len(set(MIDs))
                    if diff > 0:
                        raise Exception('MID tags in barcode files are not unique.\n{0} duplicates found'.format(diff))
                    else:
                        self.barcode_dict = dict(zip(MIDs, individuals))
                        
        # Check length of MIDs
        lengths = [len(mid) for mid in self.barcode_dict.keys()]
        assert all([x == lengths[0] for x in lengths]), 'Error: Lengths of MIDs are not equal'
        self.MIDlength = lengths[0]

        # Setup filter functions 
        self.filter_functions = []
        
        if args.illumina_filter:
            # Ilumina Machine filtering 
            self.filter_functions.append(self.set_illumina_filter())
        if 'n_thresh' in args:
            # Filter by proportion of Ns
            self.filter_functions.append(self.set_propN_filter(args.n_thresh))
        if 'phred_thresh' in args:
            # Filter by mean phred score
            self.filter_functions.append(self.set_phred_filter(args.phred_thresh))
        
        if args.cutsites:
            # Filter by target Cutsites
            self.cutsites = args.cutsites
            
            # Check all lengths are equal
            cutsite_lengths = [len(cutsite) for cutsite in self.cutsites]
            assert all([x == cutsite_lengths[0] for x in cutsite_lengths]), 'Error: Lengths of cutsites are not equal'
            self.cutsite_length = cutsite_lengths[0]

            self.filter_functions.append(self.set_cutsite_filter(target_cutsites=self.cutsites, 
                                                                 max_edit_dist=args.cutsite_editdist, 
                                                                 midtag_length=self.MIDlength))
        if 'overhang_idx' in args:
            # Filter by target Cutsites
            assert 'cutsite_editdist' in args, 'Edit distance must be specified to filter by overhang'
            self.filter_functions.append(self.set_overhang_filter(target_cutsites=self.cutsites,
                                                                    overhang=args.overhang_idx, 
                                                                    max_edit_dist=args.cutsite_editdist, 
                                                                    midtag_length=self.MIDlength)) 
        self.args = args
    
    def write_summary_output(self, path, overwrite=False):
        ''' Write the results of the filtering and cleaning to a summary output file'''
    
        filepath = os.path.join(path, "filtering_summary.log")
        
        if overwrite:
            file_handle = open(filepath, 'wb')
        else:
            # Take care of files already present
            if os.path.exists(filepath):
                count = 1
                while os.path.exists(filepath):
                    name = "filtering_summary{0:d}.log".format(count)
                    filepath = os.path.join(path, name)
                    count += 1
                print >> sys.stderr, 'Filter summary file already present. Saving as ', filepath
            
            file_handle = open(filepath, 'wb')
            
        # Write the filtering part of summary to a file 
        with file_handle as f:
            f.write("Filtering and Cleaning Parameters:\n") 
            f.write("------------------\n")
            
            args_dict = vars(self.args)
            
            temp_str = StringIO.StringIO()
            pprint.pprint(args_dict, stream=temp_str)
                
            for x in temp_str.getvalue().split(','):
                f.write(x)
                
            f.write("\n")
            f.write('Filter stats:\n')
            f.write("-------------------------------------------\n")
            f.write('Filter\t\t\tHits\n')    
            f.write("-------------------------------------------\n")
            for i, x in enumerate(self.filter_functions):
                percent = self.filterfail_counter[i] / float(sum(self.filterfail_counter.values())) 
                f.write('%s\t\t%s\t(%.2f %%)\n' % (x.__name__, self.filterfail_counter[i], percent * 100))
            f.write('\nTotal No. Reads Processed:  \t\t%s\n' % self.total_read_count)
            f.write('Total reads filtered out:  \t\t%s (%.2f %%)\n' % (sum(self.filterfail_counter.values()), 
                                                        100 * (sum(self.filterfail_counter.values())/ 
                                                               float(self.total_read_count))))
            f.write('\n')
            f.write('\nTotal reads skipped in cleaning: \t{0} ({1:.2%})\n'.format(
                    self.skipped_count, float(self.skipped_count) / self.total_read_count))
            
            f.write('\nTotal reads corrected: \t\t\t{0} ({1:.2%})\n'.format(
                    self.read_corrected_count, float(self.read_corrected_count) / self.total_read_count))
                    
            f.write('\nTotal reads with MIDtags corrected: \t{0} ({1:.2%})\n'.format(
                    self.MIDtag_corrected_count, float(self.MIDtag_corrected_count) / self.total_read_count))

            f.write('\nTotal reads with cutsite corrected: \t{0} ({1:.2%})\n'.format(
                    self.cutsite_corrected_count, float(self.cutsite_corrected_count) / self.total_read_count))
            
            f.write('\n')  
            f.write('\nTotal No. reads passed : \t\t{0} ({1:.2%})\n'.format(
                    self.total_read_passes, float(self.total_read_passes) / self.total_read_count))
                              
    def make_processing_gen(self, recgen):
        ''' A generator to only yield records that pass all filter functions
         specified
         
         Plus if error_corrected_dist argument is set, the MIDtag and cutsite will 
         be error corrected up this number of edits. 
         '''
        
        filterfuncs = self.filter_functions
        
        MIDtagslist = self.barcode_dict.keys()
        
        goto_next_rec = False
        
        for rec in recgen:
        
            for n, func in enumerate(filterfuncs):
                #===============================================================
                # Run filter functions on record       
                #===============================================================
                if func(rec) == False:
                    self.filterfail_counter.update([n])
                    self.total_read_count += 1 
                    goto_next_rec = True
                    break # move on to next record and don't do else:
            
            if goto_next_rec:
                goto_next_rec = False
                continue
                
            else: # if runs through all filterfuncs and pased, run this too. 
                self.total_read_count += 1 
                    
                # RECORD HAS PASSED ALL FILTERS
            
                if args.err_correct:
                    
                    MID_corrected = False
            
                    #===============================================================
                    # Analise MID tag section
                    #===============================================================
                    recMID = str(rec.seq[:self.MIDlength])
                    
                    if recMID not in self.barcode_dict:
                        # Sequencing error in the tag. Work out nearest candidate.
                        distvec = [editdist.distance(recMID, tag) for tag in MIDtagslist]
                        
                        distvec_min = min(distvec)
                        count_distvec_min = distvec.count(distvec_min)
                        
                        if distvec_min > self.args.err_correct:
                            self.skipped_count += 1
                            continue
                            
                        count_distvec_min = distvec.count(distvec_min)
                        if count_distvec_min > 1:
                            # Muliple candidates. True MID is Ambiguous 
                            self.skipped_count += 1
                            continue
                    
                        elif count_distvec_min == 1:
                            
                            correct_MIDtag = MIDtagslist[distvec.index(distvec_min)]
                            
                            # Correct the erroneous tag with the candidate.
                            # Letter annotations must be removed before editing rec.seq
                            temp_var = rec.letter_annotations
                            rec.letter_annotations = {}
                            # Change seq to mutableseq
                            rec.seq = rec.seq.tomutable()
                            rec.seq[:self.MIDlength] = correct_MIDtag
                            rec.seq = rec.seq.toseq()
                            rec.letter_annotations.update(temp_var)
                            self.MIDtag_corrected_count += 1
                            self.read_corrected_count += 1
                            MID_corrected = True
                            
                    #===============================================================
                    # Analise Cut site section
                    #===============================================================
                    # Must allow for multiple possible cutsites
                    rec_cutsite = str(rec.seq[self.MIDlength: self.MIDlength + self.cutsite_length])
                    if rec_cutsite not in self.cutsites:
                        # Sequencing error in the cutsite. Correct if less than max_edit_dist
                        
                        cutsite_dists = [editdist.distance(rec_cutsite, cutsite) for cutsite in self.cutsites]
    
                        min_cutsite_dist = min(cutsite_dists)
                        
                        if min_cutsite_dist > self.args.err_correct:
                            # Amount of error in cut site is too large to correct 
                            # May also be from contaminants. So read is skipped altogether 
                            self.skipped_count += 1
                            if MID_corrected:
                                # Roll back counters
                                self.MIDtag_corrected_count -= 1
                                self.read_corrected_count -= 1
                            continue
                        
                        min_cutsite_dist_count = cutsite_dists.count(min_cutsite_dist)
                        
                        if cutsite_dists.count(min_cutsite_dist_count) > 1:
                            # Muliple candidates. True cutsite is Ambiguous 
                            self.skipped_count += 1
                            if MID_corrected:
                                # Roll back counters
                                self.MIDtag_corrected_count -= 1
                                self.read_corrected_count -= 1
                            continue
                            
                            
                        elif cutsite_dists.count(min_cutsite_dist_count) == 1:
                            
                            corrected_cutsite = self.cutsites[cutsite_dists.index(min_cutsite_dist)]
                            
                            # Correct the erroneous cutsite with the actual cutsite.
                            # Letter annotations must be removed before editing rec.seq
                            temp_var = rec.letter_annotations
                            rec.letter_annotations = {}
                            # Change seq to mutableseq
                            rec.seq = rec.seq.tomutable()
                            rec.seq[self.MIDlength: self.MIDlength + len(corrected_cutsite)] = corrected_cutsite
                            rec.seq = rec.seq.toseq()
                            rec.letter_annotations.update(temp_var)
                            self.cutsite_corrected_count += 1
                            if not MID_corrected:
                                self.read_corrected_count += 1
                        else:
                            # Amount of error in cut site is too large to correct 
                            # May also be from contaminants. So read is skipped altogether 
                            self.skipped_count += 1
                            if MID_corrected:
                                # Roll back counters
                                self.MIDtag_corrected_count -= 1
                                self.read_corrected_count -= 1
                            continue
        
                # Note: rec only yielded if the read passes all filters, and if specified,
                # only if the MIDtag and cutsites are cleaned sucessfully.   
                yield rec
                
    def set_illumina_filter(self):
        ''' Returns filtering function based on illumina machine filter         
        N = was not pickup up by machine filter i.e. passed
        Y = was flagged by machine filter i.e. fail 
        '''
        # Define filter function
        def illumina_filter(rec):
            ''' filter function for phred '''
            return rec.description.split()[1].split(':')[1] == 'N'
        return illumina_filter  
    
    def set_propN_filter(self, value):
        ''' Returns a filter function based on the proportion of Ns in the read.
    
        If a read has a higher proportion of Ns than the given value, the filter 
        returns false.
        '''
        # Define filter function
        def propN_filter(rec):
            ''' Filter function for propN'''
            return float(rec.seq.count('N')) / len(rec.seq) < value
        return propN_filter
    
    def set_phred_filter(self, value):
        ''' Returns a filtering function based on the mean phred of each read.
        
        If a read has a mean phred of less than the given value, the filter returns
         false.   
        '''
        # Define filter function
        def phred_filter(rec):
            ''' filter function for phred '''
            return np.array(rec.letter_annotations['phred_quality']).mean() > value
        return phred_filter  
    
    def set_cutsite_filter(self, target_cutsites=None, max_edit_dist=None, midtag_length=None):
        ''' Returns a filter function based on the match of the read cutsite to the 
        target_cutsite given.
        
        Reads that differ in edit distance by more than mindist, cause filter to 
        return false 
        
        '''
        cutsite_length = len(target_cutsites[0])
    
        # Define filterfunc
        def cutsite_filter(rec):
            ''' Filter function for cutsite '''
            
            cutsite = rec.seq[midtag_length: midtag_length + cutsite_length].tostring()
            
            for target_site in target_cutsites:
                cutsite_dist = editdist.distance(target_site, cutsite)
                if cutsite_dist <= max_edit_dist:
                    return True
                
            return False
        
        return cutsite_filter
        
    def set_overhang_filter(self, target_cutsites=None, overhang=None, midtag_length=None, max_edit_dist=0):
        ''' Returns a filter function based on the overhang part of the cutsite. 
        
        The cut site should end with the specified overhang. Those that dont are likely 
        to be genetic contaminants which have been inadvertantly sequenced, and 
        therefore should be discarded. 
           
        Reads that mismatch in the overhang region by more than mindist, cause the 
        filter to return false. 
        
        '''   
        overhang_patterns = [i[-overhang:] for i in target_cutsites]
        cutsite_length = len(target_cutsites[0])
        
        # Define filterfunc
        def overhang_filter(rec):
            ''' Filter function for cutsite'''
            
            cutsite = rec.seq[midtag_length: midtag_length + cutsite_length].tostring()

            for i, pat in enumerate(overhang_patterns):
                
                dist = editdist.distance(target_cutsites[i], cutsite)
                if dist <= max_edit_dist:
                    if cutsite.endswith(pat):
                        return True
            
            return False
            
        return overhang_filter


#===============================================================================

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Filter and clean up FastQ files.')

    parser.add_argument('input', nargs='+',
                        help='Input file(s) to process (/path/filename). Will accept a glob')
    parser.add_argument('-b', '--barcodes', required=True, nargs='+',
                        help='Barcodes accociated with input file(s) (/path/filename). Will accept a glob')
    
    # Output parameters
    parser.add_argument('-o', dest='out_postfix',
                        help='Output file postfix to use when writing to file.')

    parser.add_argument('-p',  dest='output_path', default='.',
                        help='Output path to write reads to. Missing dirs will be created.')
    
    # Filter Parameters
    parser.add_argument('-n', dest='n_thresh', type=float,
                        help='Threshold maximum for filtering by proportion of Ns in the read. Default = 0.1. Set to 0 to skip.')

    parser.add_argument('-f',  dest='phred_thresh', type=int,
                        help='Threshold minimum for filtering by mean phred of the read. Default = 20. Set to 0 to skip.')

    parser.add_argument('-l', action='store_false', dest='illumina_filter', default=True,
                        help='Do not use the Illumina machine filter. On by Default.')

    parser.add_argument('-c', dest='cutsites', nargs='*',  default=[],
                        help='Filter reads that are not within -e edits of the cutsites specified.')

    parser.add_argument('-e', dest='cutsite_editdist', default=1,
                        help='Max edit distance allowed between target cutsite and read.')
    
    # Cleanup parameters
    parser.add_argument('-r', dest='overhang_idx', type=int,
                        help=('Number of bases in cutsite that make up the overhang. Reads are filtered out which' 
                        'have errors in the overhang of the cut site.'))
    
    parser.add_argument('-g', '--err_correct', type=int,
                        help='Max edit distance that is corrected between target MIDtag/cutsite and actual read.'
                        'If matched to more than one candidate barcode, the read is discarded due to abiguity of identity.')
    
    # For displaying filter results to stderr 
    parser.add_argument('-v', '--verbose', action='store_true', default=False, 
                        help='Whether to print out the filtering summary to stdout.')
    
    #===============================================================================
    # Main Fucntion loop
    #===============================================================================
    
    # Parse args and set defaults 
    args = parser.parse_args()
    
    starting_dir = os.getcwd()
            
    # Generator to cycle through files
    reads_generator = SeqRecCycler(data_files=args.input)
    
    # Initalise Class to Filter and cleanup reads
    preprocessor = RecordPreprocessor(args)
    
    # Output Path and file variables
    args.output_path = os.path.abspath(args.output_path)
    outputnames = []
    if not os.path.exists(args.output_path):
        os.makedirs(args.output_path) 
    
    cum_t = 0
    toc = time.time()
    
    total_read_passes = 0
    
    for recgen in reads_generator.seqfilegen:
    
        print >> sys.stderr, '\nProcessing {0}'.format(reads_generator.curfilename)
        
        # Initialise generator
        passes = preprocessor.make_processing_gen(recgen)
        
        # DEfine output location
        # Output Path and file variables
        if args.out_postfix:
    
            # Construct file names
            tail, head =  os.path.split(reads_generator.curfilename)
            name = head.split('.')  
            pass_filename = '.'.join([name[0] + args.out_postfix] + name[1:]) 
            pass_file = os.path.join(args.output_path, pass_filename)
            name = '.'.join(name)
            outputnames.append(pass_filename) 
            
            # Writes to the same filetype as the input
            if name.endswith('.bgzf'):
                pass_filehdl = bgzf.BgzfWriter(pass_file)
            elif name.endswith('.fastq'):
                pass_filehdl = open(pass_file, 'wb')
            elif name.endswith('.gz'):
                pass_filehdl = gzip.open(pass_file, 'wb')
            else:
                print >> sys.stderr, 'Input file format not supported: %s' % name
                sys.exit()
            
        else:
            # Output is written to std out
            pass_filehdl = sys.stdout
        
        if args.out_postfix:          
            print >> sys.stderr, 'Writing passes to \n{0} ....'.format(pass_filename)
        else:
            print >> sys.stderr, 'Writing to stdout ....'
        numwritten = SeqIO.write(passes , pass_filehdl , 'fastq')
        preprocessor.total_read_passes += numwritten
        
        if pass_filehdl == sys.stdout:
            pass_filehdl.flush()
        else:
            pass_filehdl.close()
        print >> sys.stderr, '{0} records Preprocessed'.format(numwritten)
        
        loop_t = time.time() - toc - cum_t
        cum_t += loop_t
        print >> sys.stderr, 'Finished file {0} after {1}'.format(reads_generator.curfilenum, 
                                    time.strftime('%H:%M:%S', time.gmtime(loop_t))) 
    
    if args.verbose:
        print >> sys.stderr, '\nFilter stats'
        print >> sys.stderr, '\nFilter No.\tHits'    
        for n in range(len(preprocessor.filter_functions)):
            print >> sys.stderr, '%s\t\t%s' % (n,preprocessor.filterfail_counter[n])
    
        print >> sys.stderr, '\nTotal No. Reads Processed:  %s' % preprocessor.total_read_count
        
        print >> sys.stderr, '\nTotal No. filtered out:  %s (%.2f %%)' % (sum(preprocessor.filterfail_counter.values()), 
                                                        100 * (sum(preprocessor.filterfail_counter.values())/ 
                                                                float(preprocessor.total_read_count)))
        
        print >> sys.stderr, '\nTotal reads skipped in cleaning: \t{0} ({1:.2%})'.format(
                        preprocessor.skipped_count, float(preprocessor.skipped_count) / preprocessor.total_read_count)
        
        print >> sys.stderr, '\nTotal reads corrected: \t\t\t{0} ({1:.2%})'.format(
                        preprocessor.read_corrected_count, float(preprocessor.read_corrected_count) / preprocessor.total_read_count)
        
        print >> sys.stderr, '\nTotal No. reads passed: \t\t{0} ({1:.2%})'.format(
                        preprocessor.total_read_passes, float(preprocessor.total_read_passes) / preprocessor.total_read_count)
        
    # Write summary file
    preprocessor.write_summary_output(args.output_path)
    
    total_t = time.time() - toc    
    print >> sys.stderr, 'Processed all files in {0}'.format(time.strftime('%H:%M:%S', 
                                                                time.gmtime(total_t)))
    os.chdir(starting_dir)
    
    
        
        