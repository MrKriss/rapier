'''
Created on 1 Jul 2013

@author: musselle
'''

import os
import sys
import cPickle as pkl
from subprocess import PIPE, Popen
import glob
import time

import gzip

from Bio import SeqIO, bgzf

class SeqRecCycler(object):
    ''' Object to hold generators that yield Sequence record objects or
    generators for all sequences records from the given file list

    Takes care of file input checks. Each next call returns the next file
    containing the next set of Sequence record objects.
    
    INPUTS
    data_files - Single file as string or list of files to process as strings with
              full extensions. If data_files is a string, it is treated as a glob to 
              the specified path. e.g. glob.glob('*.fastq') for all files ending in .fastq
              data files contain the full path to the file + filename
              
    METHODS
    __init_files_gen - initiate files generator
    __init_rec_gen - initiate record generator

    ATTRIBUTES
    .seqfilegen - a generator that returns a generator per file to iterate over
                its sequence records.
    .recgen - a generator that iterates over all records serially for all files
                in list.
    .curfilename - Current file name in list of given files (infiles) being
                iterated over.
    .curfilenum - Current file number in list of given files (infiles).
    .maxNoSeq - Number of sequence records to return from each file. Default is
                all. Useful for testing.

    '''

    def __init__(self, data_files=None, maxnumseq=None):
        ''' Constructor '''

        # Handle multiple types of input for data_files
        if type(data_files) is str:
            
            # Check for path head included in data_files
            if not os.path.split(data_files)[0]:
                # Append the abs path to current directory 
                current_dir = os.getcwd()
                data_files = os.path.join(current_dir, data_files)
            
            # Glob the file pattern 
            files = glob.glob(data_files)
            assert files, "No files returned from glob"
            files.sort()
            self.input_files = files
            
        elif type(data_files) is list or type(data_files) is tuple:
            
            if not os.path.split(data_files[0])[0]:
                current_dir = os.getcwd()
                data_files = [os.path.join(current_dir, x) for x in data_files]
            self.input_files = data_files
        else:
            raise Exception('Invalid entry for data_files.')

        print '\nGenerator initiated to process the following files...'
        for f in self.input_files:
            print f

        self.numfiles = len(self.input_files)
        
        # Max number of record to run the generator for (default = all)
        self.maxnumseq = maxnumseq

        ''' May need to process files individually or all at once, so two types
        of generators needed '''
        # Generator that yields a sequence record iterator (generator)
        # for the next file in the list
        self.seqfilegen = self.__init_files_gen()

        # Generator that yields the next individual record by running though
        # all files given in the list.
        self.recgen = self.__init_rec_gen()        

    def __init_files_gen(self):
        ''' Constructs the next file generator object '''
        # Generic file handling code
        for filenum, filename in enumerate(self.input_files):

            self.curfilenum = filenum
            self.curfilename = filename

            # Check file extensions
            if (filename.endswith('.gz') or filename.endswith('.bgzf') or
                                            filename.endswith('.fastq') or filename.endswith('.fq')):
                try:
                    yield SeqIO.parse(smartopen(filename), format='fastq')
                except IOError as e:
                    print e
                    raise Exception(('Invalid file name {0},'
                     ' or path {1} may not exist.').format(filename, os.path.basename(filename)))
            else:
                print 'File extension for {0} currently unsupported.'.format(filename) 
                print 'Accepted formats for cycling through all records = .gz .bgzf and .fastq'
                raise Exception

    def __init_rec_gen(self):
        ''' Return next Sequence Record '''

        if self.maxnumseq:
            for rec_file in self.seqfilegen:
                count = 0
                for record in rec_file:
                    count += 1
                    if count <= self.maxnumseq: 
                        yield record
                    else:
                        break
        else:
            for rec_file in self.seqfilegen:
                for record in rec_file:
                    yield record

def pklsave(obj, filename):
    ''' Pickle the given object '''
    if not filename.endswith('pkl'):
        filename = filename + '.pkl'
    with open(filename, 'wb') as f:
        pkl.dump(obj, f)
        
def smartopen(filename, *args, **kwargs):
    '''opens with open unless file ends in .gz or .bgzf, then use gzip.open

    in theory should transparently allow reading of files regardless of compression
    
    gzip can read in both .gz and .bgzf files.
    '''
    if filename.endswith('.gz') or filename.endswith('.bgzf'):
        return gzip.open(filename, *args, **kwargs)
    else:
        return open(filename, *args, **kwargs)
    
def inputfile_check(handle):
    """ Check if input is a valid file object or a string for a file name.
    
    Always returns a handle for an open file.
    
    Usage:
    handle = input_check(handle)
    
    """
    
    if type(handle) is file:
        # Check if handle is stdin
        if handle.name == '<stdin>' :
            return handle
        if handle.closed:
            handle = open(handle.name, 'rb')
        if handle.mode not in ['r', 'rb']:
            print >> sys.stderr, "File not opened for reading. Mode = ".format(handle.mode)
    
    elif type(handle) is str:
        # Check file exists and open it
        try:
            handle = open(handle, 'rb')
        except IOError as e: 
            print  >> sys.stderr, "I/O error({0}): {1}".format(e.errno, e.strerror)
            raise
    else:
        raise Exception('Invalid file handle.')

    return handle 

def outputfile_check(handle, overwrite=False, mode='wb'):
    ''' Checks for valid file object or if a string, checks if a file already exists
     at that location and generates a new file name or overwrites if so.
    
    overwrite -- If set to true will overwrite any previouse files. Default is to 
                 append an integer at the end of main filename. 
    '''

    if type(handle) is file:
        # Check if handle is stdin
        if handle.name in ['<stdout>', '<stderr>'] :
            return handle
        if handle.closed:
            handle = open(handle.name, mode)
        if handle.mode not in ['a', 'a+', 'w', 'wb']:
            print >> sys.stderr, "File not opened for Writing. Mode = ".format(handle.mode)

    elif type(handle) is str:
        # Check if a file already exists with the name
        if os.path.exists(handle):
            if overwrite:
                print >> sys.stderr, 'File {0} already present.\nOverwriting... '.format(handle)
            else:
                count = 1
                path, name = os.path.split(handle)
                name_parts = name.split('.')
                while os.path.exists(handle):
                    name_parts[0] = name_parts[0].rstrip('0123456789')
                    new_name = name_parts[0] + "{0:d}".format(count)
                    name_parts[0] = new_name
                    new_name = '.'.join(name_parts)
                    handle = os.path.join(path, new_name)
                    count += 1
                print >> sys.stderr, 'File already present. Saving as ', handle

        handle = open(handle, mode)
    else:
        raise Exception('Invalid file handle.')

    return handle 
    

def find_numrec(filename):
    ''' Find number of records in a fastq file. 
    
    Makes assumption that all records come from same machine and have a fastq format
    file that shares the filename given. 
    
    Uses subprocess to call 'grep -c 'filename'
    '''
    
    if filename.endswith('.idx') or filename.endswith('.gz') \
                            or filename.endswith('.bgzf') or filename.endswith('.fastq'):
        # Check fastq file exists.
        if filename.endswith('.fastq'):
            fastq_filename = filename
        else:
            fastq_filename = filename.split('.')[0] + '.fastq' 
        try:
            with open(fastq_filename, 'rb') as f:
                # Load in machine name
                sequencerID = f.readline().strip().split(':')[0]
        except IOError as e:
                print e
                print 'Invalid file name, or {0} may not exist.'.format(fastq_filename)
        filename = fastq_filename
    else:
        print 'File extension not recognised.'
        sys.exit()
    
#    cmd = 'cat {0} | parallel --pipe --block 2M grep -c "{1}"'.format(filename, sequencerID)
    cmd = 'grep -c "{1}" {0}'.format(filename, sequencerID)
    
    process = Popen(cmd, stdout=PIPE, shell=True)

    return int(process.communicate()[0].strip())


#===============================================================================
# File conversions
#===============================================================================

def file2bgzf(infiles=None, filepattern=False, data_inpath='', SQLindex=True):
    ''' Convert given list of files from .gz or .fastq to .bgzf,
    And also produce an SQL index if needed. 
    
    infiles accepts str of file name of list of str for filenames. 
    If not specified will look at file type and glob the result to infiles. 
  
    '''
    if data_inpath:
        os.chdir(data_inpath)
  
    # Handle multiple types of input for infiles
    assert infiles is not None, 'No files listed or file pattern specified.'         
    if filepattern:
        # Fetch files by file types using glob
        import glob 
        infiles = glob.glob(infiles)
    elif type(infiles) == str:
        # Convert to list
        infiles = [infiles]
  
    start_time = time.time() 
    
    for filename in infiles: 
        toc = time.time()
        
        f = smartopen(filename)
        
        # Checks for type of input
        if filename.endswith('.gz'):
            # Drop .gz and append .bgzf
            bgzfFileName = '.'.join(filename.split('.')[:-1]) + '.bgzf'  
        elif filename.endswith('.fastq'):
            # Append .bgzf
            bgzfFileName = filename + '.bgzf'

        print "Producing BGZF output from {0}...".format(filename)
        w = bgzf.BgzfWriter(bgzfFileName, 'wb')
        while True:
            data = f.read(65536)
            w.write(data)
            if not data:
                break
        w.close()
        print '{0} written successfully'.format(bgzfFileName)

        conv_t = time.time() - toc 
        print 'Finished converting {0}\n after {1}\n'.format(filename, time.strftime('%H:%M:%S', time.gmtime(conv_t)))
        
        if SQLindex == True:
            makeSQLindex(bgzfFileName)
      
    total_t = time.time() - start_time
    print 'Finished all processing {0} files in {1}'.format(len(infiles), time.strftime('%H:%M:%S', time.gmtime(total_t)))
 
def makeSQLindex(infiles=None, data_inpath='', mode='grouped', outname=None):
    ''' Creates an SQL index out of either an uncompressed file or a compressed .bgzf file 
    
    if infiles is a string it is interpreted as a glob
    
    if infiles is list, goes through all file names in list.
    
     - mode  - grouped: all files are indexed to a single index file, specified by outname 
    
    '''
    
    starting_dir = os.getcwd()
    
    if data_inpath:
        os.chdir(data_inpath)
        
    if outname is None:
        outname = 'reads.idx'
  
    if type(infiles) is str:
        # Fetch files by file types using glob
        import glob 
        infiles = glob.glob(infiles)
    elif type(infiles) is not list and type(infiles) is not tuple:
        raise Exception("Invalid input files specified.")

    assert infiles, 'No files found, or no files passed.'

    # Handle multiple types of input for infiles
    if mode == 'grouped':
        idx_filename = outname
        tak = time.time()
        print 'Writing {0} files to SQL index ...'.format(len(infiles))
        SeqIO.index_db(idx_filename, infiles , 'fastq')
        idx_t = time.time() - tak
        print 'Finished Indexing to {0}\n after {1}\n'.format(idx_filename, time.strftime('%H:%M:%S', time.gmtime(idx_t)))

    elif mode == 'separate':
    
        for filename in infiles: 
            tak = time.time()
            print 'Writing SQL index file for {0} ...'.format(filename)
            idx_filename = filename.split('.')[0] + '.idx'
            SeqIO.index_db(idx_filename, filename , 'fastq')
            print '{0} written successfully'.format(idx_filename)
            idx_t = time.time() - tak
            print 'Finished Indexing after {1}\n'.format(time.strftime('%H:%M:%S', time.gmtime(idx_t)))

    if os.getcwd() != starting_dir:
        os.chdir(starting_dir) 


def file2fasta(filename):
    ''' Convert fastq file to fasta file '''
    
    handle = smartopen(filename)
    out_filename = filename.split('.')[0] + '.fasta'
    
    count = SeqIO.convert(handle, 'fastq', out_filename, 'fasta')
    
    print 'Converted {0} records to file\n{1}'.format(count, out_filename)

    
def bgzf2fastq(infiles=None, data_inpath=''):
    """ Unzip a list of files in .bgzf 
    
    if infiles is a string it is interpreted as a glob
    
    if infiles is list, goes through all file names in list.
    
    """
    
    starting_dir = os.getcwd()
    
    if data_inpath:
        os.chdir(data_inpath)
        
    if type(infiles) is str:
        # Fetch files by file types using glob
        import glob 
        infiles = glob.glob(infiles)
    elif type(infiles) is not list and type(infiles) is not tuple:
        raise Exception("Invalid input files specified.")

    assert infiles, 'No files found, or no files passed.'

    start_time = time.time() 
    for filename in infiles: 
        toc = time.time()
        
        infile_handle = smartopen(filename)
        
        # Checks for type of input
        if filename.endswith('.bgzf'):
            # Drop .bgzf and append .fastq if necessary
            fastqFileName = '.'.join(filename.split('.')[:-1])
            if not fastqFileName.endswith('.fastq'):
                fastqFileName +=  '.fastq'
              
        print "Producing FastQ output from {0}...".format(filename)
        
        outfile_handle = open(fastqFileName, 'wb')
        
        for rec in SeqIO.parse(infile_handle, 'fastq'):
            SeqIO.write(rec, outfile_handle, 'fastq')
        
        print '{0} written successfully'.format(fastqFileName)

        conv_t = time.time() - toc 
        print 'Finished converting {0}\n after {1}\n'.format(filename, time.strftime('%H:%M:%S', time.gmtime(conv_t)))
        
    total_t = time.time() - start_time
    print 'Finished all processing {0} files in {1}'.format(len(infiles), time.strftime('%H:%M:%S', time.gmtime(total_t)))
    

if __name__ == '__main__':
    
    os.chdir('/space/musselle/datasets/gazellesAndZebras/lane6')
    out = find_numrec('lane6_NoIndex_L006_R1_001.fastq')