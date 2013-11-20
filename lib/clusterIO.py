"""
Created on 24 Apr 2013

@author: musselle

Defines a cluster class to hold all processing methods, as well as IO methods to extract information about sequence
and sampleid from the database.

Also includes utility funtions to provide IO support for CDHIT clustering files in python, such as:
    * Iterate through Clusters in a file
    * Return counters for cluster sizes in a file
    * Analise Distribution of sequences within cluster
    (counts dictionary for how far sequences are away from the representative one)

Sample Cluster File:

>Cluster 0
0    88nt, >2096... *
1    88nt, >25159... at +/100.00%
2    88nt, >52830... at +/98.86%
3    88nt, >76269... at +/100.00%
4    88nt, >12309... at +/100.00%
5    88nt, >31554... at +/100.00%
6    88nt, >11162... at +/100.00%
7    88nt, >35059... at +/98.86%
>Cluster 1
0    88nt, >96... *
1    88nt, >55944... at +/100.00%
2    88nt, >22466... at +/98.86%
3    88nt, >45463... at +/97.73%

"""
import os
import shlex
import time
from collections import defaultdict, Counter
from editdist import distance
import subprocess
from StringIO import StringIO

import pandas as pd

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import AlignIO, SeqIO
import numpy as np


class ClusterObj(object):
    """ Holds all cluster based information. """
    
    def __init__(self):
        
        # Cluster vars
        self.rep_seq_id = ""      
        self.rep_seq = ""
        self.rep_phred = None
        self.rep_sample_id = None
        self.members_id = []    
        self.members_seq = []
        self.members_phred = []
        self.members_sample_id = []
        self.size = 0  
        self.id = 0
        self.edit_dists = []
        
        # Start and end Locations on the cluster file in bytes from the beginning  
        self.start_loc = 0
        self.end_loc = 0
        
    def getfromdb(self, items, target, db):
        """ Lookup to main db to retrieve and store the specified items.
        
        items - string of fields to fetch for the cluster from the database. 
                Option of 'seq', 'phred' and 'sampleId' 
        
        target - Whether to fetch the data for the rep seq ('rep'), or all cluster members ('all'). 
        
        db    - Reference to a Reads_db database object.
        
        Note: This method can be slow due to the number of lookups required. 
        """
        
        # Setup query
        sql_query = """ SELECT {0} FROM seqs WHERE seqid = ? """.format(','.join(items))
        multi_sql_query = """ SELECT {0} FROM seqs WHERE seqid IN """.format(','.join(items))
        
        get_seq = 0
        get_phred = 0
        get_sampleid = 0
        if 'seq' in items:
            get_seq =  True
        if 'phred' in items:
            get_phred = True
        if 'sampleId' in items:
            get_sampleid = True
            
        if target == 'rep':
            # Get for rep seq 
            record_curs = db.con.execute(sql_query, (self.rep_seq_id,))
            record = record_curs.fetchone()
            
            if get_seq: 
                self.rep_seq = record['seq']
            if get_phred:
                phred_ascii = record['phred']
                phred_list = [ord(c) - 33 for c in phred_ascii]
                self.rep_phred = np.array(phred_list)
            if get_sampleid:
                self.rep_sample_id = int(record['sampleId'])
            
        elif target == 'all':
            
            # Get for rep seq 
            record_curs = db.con.execute(sql_query, (self.rep_seq_id,))
            record = record_curs.fetchone()
            
            if get_seq: 
                self.rep_seq = record['seq']
            if get_phred:
                phred_ascii = record['phred']
                phred_list = [ord(c) - 33 for c in phred_ascii]
                self.rep_phred = np.array(phred_list)
            if get_sampleid:
                self.rep_sample_id = int(record['sampleId'])
            
            # get members seq and phred
            
            # Optimise sql query             
            # records are returned in batch sorted by seqid
            self.members_id.sort()
            
            s = [str(i) for i in self.members_id]
            multi_sql_query += '({0})'.format(','.join(s))
            record_curs = db.con.execute(multi_sql_query)
            
            for record in record_curs:
                
                if get_seq: 
                    self.members_seq.append(record['seq'])
                if get_phred:
                    phred_ascii = record['phred']
                    phred_list = [ord(c) - 33 for c in phred_ascii]
                    self.members_phred.append(np.array(phred_list))
                if get_sampleid:
                    self.members_sample_id.append(record['sampleId'])

    def align(self, db=None, start_idx=6, muscle_exec_path='~/bin/muscle'):
        """ Create an Allignment of the sequences in the cluster so as to accomodate indels. """

        # Assert seq data is present to use
        if not self.members_seq:
            assert db, 'Sequence data not present and no lookup database specified.'
            print "Sequence data not present in cluster. Retrieving from data base..."
            self.getfromdb(items=['seq', 'seqid'], target='all', db=db)
        if not self.rep_seq:
            assert db, 'Sequence data not present and no lookup database specified.'
            print "Sequence data not present in cluster. Retrieving from data base..."
            self.getfromdb(items=['seq', 'seqid'], target='rep', db=db)

        # Get list of sequence Records in fasta format.
        allSeqRecs = [
            SeqRecord(Seq(self.rep_seq[start_idx:]), id=str(self.rep_seq_id), description=str(self.rep_sample_id))
        ]

        for i in range(len(self.members_seq)):
            rec = SeqRecord(Seq(self.members_seq[i][start_idx:]), id=str(self.members_id[i]),
                            description=str(self.members_sample_id[i]))
            allSeqRecs.append(rec)

        # Write to a tempfile
        temp_file = open('temp_file_in.fasta', 'w')
        SeqIO.write(allSeqRecs, temp_file, format='fasta')
        temp_file.flush()
        temp_file.close()

        temp_filepath = os.path.join(os.path.abspath(os.getcwd()), temp_file.name)

        # Align with MUSCLE
        cmds = os.path.expanduser(muscle_exec_path) + ' -in ' + temp_filepath + ' -quiet '

        cmds = shlex.split(cmds)

        child = subprocess.Popen(cmds, stdout=subprocess.PIPE)

        stdout = child.communicate()[0]

        align = AlignIO.read(StringIO(stdout), 'fasta')

        return align

    def align2(self, db=None, start_idx=6, muscle_exec_path='~/bin/muscle'):
        """ Create an Allignment of just the unique sequences in the cluster so as to accomodate indels.

        The id in the sequence is the index for the sorted unique sequences table
        """

        # Optimised version
        # Get all unique seqs, then find alignment just for them.

        # Assert seq data is present to use
        if not hasattr(self, 'uniqueseqs_table'):
            assert db, 'Sequence data not present and no lookup database specified.'
            # Calculate uniq seq for each individual
            self.get_unique_seq_by_individual(db=db)

        allSeqRecs = []

        # Get list of sequence Records in fasta format.
        for i, seq in enumerate(self.uniqueseqs_table.columns):
            rec = SeqRecord(Seq(seq[start_idx:]), id=str(i))
            allSeqRecs.append(rec)

        # Write to a tempfile
        temp_file = open('temp_file_in.fasta', 'w')
        SeqIO.write(allSeqRecs, temp_file, format='fasta')
        temp_file.flush()
        temp_file.close()

        temp_filepath = os.path.join(os.path.abspath(os.getcwd()), temp_file.name)

        # Align with MUSCLE
        cmds = os.path.expanduser(muscle_exec_path) + ' -in ' + temp_filepath + ' -quiet '

        cmds = shlex.split(cmds)

        child = subprocess.Popen(cmds, stdout=subprocess.PIPE)

        stdout = child.communicate()[0]

        alignment = AlignIO.read(StringIO(stdout), 'fasta')

        # Update columns in uniqueseqs_table to be aligned sequences
        new_idx = [0] * len(self.uniqueseqs_table.columns)

        for seqrec in alignment:
            new_idx[int(seqrec.id)] = seqrec.seq.tostring()

        self.uniqueseqs_table.columns = pd.Index(new_idx)


    def align2seq(self, alignedseq):
        """ Utility function to convert aligned sequences back to original seq by stripping '-' characters """

        import Bio
        from Bio import SeqRecord, Seq

        if type(alignedseq) is Bio.SeqRecord.SeqRecord:
            alignedseq = alignedseq.seq.tostring().replace('-', '')

        elif type(alignedseq) is Bio.Seq.Seq:
            alignedseq = alignedseq.tostring().replace('-', '')

        return alignedseq


    def genotype(self, db=None):
        """ Call genotypes on the cluster for each individual """

        # Get all necessary data: seq, phred,
        items = []
        if not self.members_seq:
            items.append('seq')
        if not self.members_phred:
            items.append('phred')
        if not self.members_sample_id:
            items.append('seqid')
        if items:
            assert db, 'Sequence data not present and no lookup database specified.'
            print "Sequence data not present in cluster. Retrieving from data base..."
            self.getfromdb(items=items, target='all', db=db)

        # Assert seq data is present to use
        if not hasattr(self, 'uniqueseqs_table'):
            assert db, 'Sequence data not present and no lookup database specified.'
            # Calculate uniq seq for each individual
            self.get_unique_seq_by_individual(db=db)

        # Set columns in uniqueseqs_table to multiple sequence alignment.
        self.align2(start_idx=0)

        # Use most frequent seq in cluster as reference genotype
        useqs_total = self.uniqueseqs_table.sum()
        refseq = useqs_total.index[0]
        n = len(refseq)
        m = len(self.uniqueseqs_table.index)

        # List the Commonest Nucleotide not in ref seq per base position
        # Done for each individual separately

        # D = { 'individual' :  defaultd{ baseposition : Counter }

        ds = defaultdict(lambda: defaultdict(Counter))

        # for each inidvidual
        for indiv in self.uniqueseqs_table.index:
            for bp in range(n):
                indiv_total_useqs = self.uniqueseqs_table.ix[indiv]
                for seq in indiv_total_useqs.index:
                    # Count up base paires for each individual and position
                    if indiv_total_useqs[seq] > 0:
                        ds[indiv][bp][seq[bp]] += indiv_total_useqs[seq]

                # Delete those that are the same as the reference genome.
                del ds[indiv][bp][refseq[bp]]

        next_common_nuc = defaultdict(lambda: [''] * n)

        # for each inidvidual
        for indiv in self.uniqueseqs_table.index:
            for bp in range(n):
                if ds[indiv][bp]:
                    next_common_nucleotide_tuple = ds[indiv][bp].most_common()[0]
                    if next_common_nucleotide_tuple[1] != 0:
                        next_common_nuc[indiv][bp] = next_common_nucleotide_tuple




        # # List the Commonest Nucleotide not in ref seq per base position
        # # Done over all
        #
        # ds = defaultdict(Counter)
        # for bp in range(n):
        #     for seq in useqs_total.index:
        #         ds[bp][seq[bp]] += useqs_total[seq]
        #
        #     del ds[bp][refseq[bp]]
        #
        # next_common_nuc = ['-'] * len(refseq)
        # for bp in range(n):
        #     if ds[bp]:
        #         next_common_nuc[bp] = ds[bp].most_common()[0]
        #
        # return ds, refseq, next_common_nuc, useqs_total

        # -----------------------
        # Genotyping Calculations
        # -----------------------

        # predefined error rates
        epsilon_insert = 0.001
        epsilon_deletion = 0.001

        m = len(refseq)

        #map from seq to aligned seq
        seq2aligned = {}
        for aligned_seq in self.uniqueseqs_table.columns:
            seq2aligned[aligned_seq.replace('-', '')] = aligned_seq

        # For each individual
        for k, indiv in enumerate(self.uniqueseqs_table.index):

            # Fetch all reads for the individual
            idx = np.where(np.array(self.members_sample_id + [self.rep_sample_id]) == self.desc2id[indiv])
            allreads = np.array(self.members_seq + [self.rep_seq])[idx]
            allphreds = np.array(self.members_phred + [self.rep_phred])[idx]

            m = len(refseq)
            n = len(allreads)

            # Reset model probabilities
            # True Nuc in ref Nuc
            pm1_nu = np.zeros((m, n))
            # True Nuc is alternative Nuc
            pm2_nu = np.zeros((m, n))
            # Individual is heterozygos
            pm3_nu = np.zeros((m, n))

            # for all reads for that individual
            for j, read in enumerate(allreads):

                # Get aligned seq
                alignedread = seq2aligned[read]

                # Get phred score for read,
                phred = allphreds[idx]

                # find inserted dash positions
                insert_idxs = []
                start_idx = 0
                value_idx = 0
                while value_idx != -1:
                    value_idx = alignedread.find('-', start_idx)
                    if value_idx != -1:
                        insert_idxs.append(value_idx)
                        start_idx = value_idx + 1

                # insert -1 where there is a dash in phred score
                phred = phred.tolist()
                for position in insert_idxs:
                    phred.insert(-1, position)

                # For each base position
                for bp in range(len(refseq)):

                    # 8 possibilities
                    #
                    # b = base
                    # - = no base
                    #                            seq = ref      seq = Nmcn    seq = ref or Nmcn
                    #     Ref nmun seq  epsilon   pm1.1    1.2   pm2.1   2.2   pm3.1   3.2    3.3
                    # 1    b   b    b    phred
                    # 2    b   b    -    edel      edel             edel         2 * edel / 3
                    # 3    b   -    b    phred     1 - e             e           (1 - e) / 2
                    # 4    b   -    -
                    # 5    -   b    b
                    # 6    -   b    -
                    # 7    -   -    b
                    # 8    -   -    -



                    # Deal with insertions and deletions
                    if refseq[bp] == '-' and next_common_nuc[indiv][bp] == '-' :


                        if read[bp] == '-':
                            # Condition 8: Treat '-' as a 5th base

                            pm1_nu[bp, j] =  1 - epsilon
                            pm2_nu[bp, j] =  1 - epsilon
                            pm3_nu[bp, j] = (2 * epsilon) / 3.




                        else:
                            # Condition 7:
                            pass








                    if read[bp] == '-' and refseq[bp] == '-':

                        # check if next most common is not a gap
                        if next_common_nuc[indiv][bp] == '-':






                            continue
                        else:

                            # TODO: Find the phred for the alternative sequences

                            epsilon = 10 ** () #

                            # bp != either ref or next nuc
                            pm1_nu[bp, j] = epsilon
                            pm2_nu[bp, j] = epsilon
                            pm3_nu[bp, j] = (2 * epsilon) / 3.

                    elif read[bp] == '-' and refseq[bp] != '-':
                        # There has been a deletion in read relative to refseq
                        # no phred will be present, so assume epsilon of 0.01
                        epsilon = epsilon_deletion

                        if read[bp] == next_common_nuc[indiv][bp]:
                            # bp = next most common nuc for that individual
                            pm1_nu[bp, j] = epsilon
                            pm2_nu[bp, j] = 1 - epsilon
                            pm3_nu[bp, j] = ((1 - epsilon) / 2.) + (epsilon / 6.)
                        else:
                            # bp != either ref or next nuc
                            pm1_nu[bp, j] = epsilon
                            pm2_nu[bp, j] = epsilon
                            pm3_nu[bp, j] = (2 * epsilon) / 3.


                    elif refseq[bp] == '-' and read[bp] != '-':
                        # There has been an insertion in read relative to refseq

                        epsilon = epsilon_insert

                        if read[bp] == next_common_nuc[indiv][bp]:
                            # bp = next most common nuc for that individual
                            pm1_nu[bp, j] = epsilon
                            pm2_nu[bp, j] = 1 - epsilon
                            pm3_nu[bp, j] = ((1 - epsilon) / 2.) + (epsilon / 6.)











                    else:
                        # Check for SNPs

                        epsilon = 10 ** (-phred[idx][bp] / 10.)

                        # Record probabilities for models
                        if read[bp] == refseq[bp]:
                            # bp = refseq
                            pm1_nu[bp, j] = 1 - epsilon
                            pm2_nu[bp, j] = epsilon
                            pm3_nu[bp, j] = ((1 - epsilon) / 2.) + (epsilon / 6.)

                        elif read[bp] == next_common_nuc[indiv][bp]:
                            # bp = next most common nuc for that individual
                            pm1_nu[bp, j] = epsilon
                            pm2_nu[bp, j] = 1 - epsilon
                            pm3_nu[bp, j] = ((1 - epsilon) / 2.) + (epsilon / 6.)

                        else:
                            # bp != either ref or next nuc
                            pm1_nu[bp, j] = epsilon
                            pm2_nu[bp, j] = epsilon
                            pm3_nu[bp, j] = (2 * epsilon) / 3.







                    if read[bp] == refseq[bp]:
                        pass




        return ds, refseq, next_common_nuc, useqs_total














    # def correct_repseq(self, db=None):
    #     """" Examine whether representative sequences is truely the most common for the cluster
    #     and correct it in the database if necessary.
    #
    #     For faster computation, the editdist_counter dict should be used. This must be read from the
    #     CDHIT output file as it is scanned to get the cluster info. To get the dictionary set the
    #     'similarity_count' flag to True in the 'parse' function.
    #
    #
    #     """
    #
    #     # Fetch all unique seq data
    #     self.get_unique_seq(seq_start_idx=6, db=db)
    #     most_common_seq = self.unique_seqs.most_common()[0]
    #
    #     if most_common_seq != self.rep_seq:
    #
    #         # Find next matching sequenceid to the genuine rep seq
    #         idx = self.members_seq.index(most_common_seq)
    #
    #         # Store old values
    #         old_rep_seq = self.rep_seq
    #         old_rep_seq_id = self.rep_seq_id
    #         old_rep_sample_id = self.rep_sample_id
    #         old_rep_phred = self.rep_phred
    #
    #         # Update rep seq
    #         self.rep_seq = self.members_seq[idx]
    #         self.rep_seq_id = self.members_id[idx]
    #         self.rep_sample_id = self.members_sample_id[idx]
    #
    #         # Remove id from members
    #         del self.members_seq[idx]
    #         del self.members_id[idx]
    #         self.members_sample_id[idx]
    #
    #         # Add rep seq id and Find the index to insert rep data to
    #         self.members_id.append(old_rep_seq_id)
    #         x = np.array(self.members_id)
    #         sortidx = int([i for i, elem in enumerate(x.argsort()) if elem == len(x) -1])
    #
    #         # Insert old rep data into members
    #         self.members_seq.insert(sortidx, old_rep_seq)
    #         self.members_sample_id.insert(sortidx, old_rep_sample_id)
    #         self.members_phred.insert(sortidx, old_rep_phred)
    #         self.members_id.sort()
    #
    #         return True
    #     else:
    #         return False
        
    def get_unique_seq(self, seq_start_idx=6, db=None):
        """ Work out the counts for unique reads within a cluster """
        
        if not self.members_seq:
            assert db, 'Sequence data not present and no lookup database specified.'
            print "Sequence data not present in cluster. Retrieving from data base..."
            self.getfromdb(items=['seq'], target='all', db=db)
        if not self.rep_seq:
            assert db, 'Sequence data not present and no lookup database specified.'
            print "Sequence data not present in cluster. Retrieving from data base..."
            self.getfromdb(items=['seq'], target='rep', db=db)
        
        unique_seq_counter = Counter()
        unique_seq_counter[self.rep_seq[seq_start_idx:]] += 1
        
        seqs = [s[seq_start_idx:] for s in self.members_seq]        
        unique_seq_counter.update(Counter(seqs))

        return unique_seq_counter
        
    def get_unique_seq_by_individual(self, db):
        """ Show the breakdown on sequence counts per individual in the cluster """
            
        # ds = { sampleId : { unique_seq : count } }
        ds = defaultdict(Counter)
        
        if not self.members_seq or not self.members_sample_id:
            print "Sequence data not present in cluster. Retrieving from data base..."
            self.getfromdb(['seq', 'sampleId'], target='all', db=db)

        # Add representative seq
        ds[self.rep_sample_id][self.rep_seq] += 1

        max_i = 0
        # Run over all memebrs        
        for i in range(len(self.members_seq)):
            ds[self.members_sample_id[i]][self.members_seq[i]] += 1
            if i+1 > max_i: max_i = i+1
        
        # Rank Unique sequences
        all_seq_counter = Counter()
        for c in ds.itervalues():
            all_seq_counter.update(c)
            
        ranked_seqs = all_seq_counter.most_common()

        # print out results
        n = len(ds)
        m = len(ranked_seqs)
        freq_matrix = np.zeros([n, m], dtype=int)

        # Store seqs
        seqs = []
        first = 1

        # Fetch all sampleid present in cluster inascending order
        sampleids = sorted(ds.keys())

        id2desc = {}

        for i, sid in enumerate(sampleids):

            for j, seq_count_tup in enumerate(ranked_seqs):
                if first:
                    seqs.append(seq_count_tup[0])
                freq_matrix[i, j] = ds[sid][seq_count_tup[0]]
            first = 0

        # Get actual description of individuals
        sampledescriptions = []
        for x in sampleids:
            c = db.con.execute('select description from samples where sampleid = ?', (x,))
            desc = c.fetchone()['description']
            sampledescriptions.append(desc)
            id2desc[x] = desc

        desc2id = {v : k for k, v in id2desc.items()}

        import pandas as pd

        df = pd.DataFrame(data=freq_matrix, index=sampledescriptions, columns=seqs, dtype=int)

        # save as attributes as well as return values
        self.uniqueseqs_table = df
        self.ds = ds
        self.id2desc = id2desc
        self.desc2id = desc2id

        return df, ds, id2desc, desc2id


    
    def get_basefraction(self, db=None):
        """ Calculate the fraction of nucleotide bases per base location """
        
        # Make sure seq data is available
        if not self.members_seq:
            print "Sequence data not present in cluster. Retrieving from data base..."
            self.getfromdb(['seq'], db=db)
        if not self.rep_seq:
            print "Sequence data not present in cluster. Retrieving from data base..."
            self.getfromdb(['rep'], db=db)
        
        self.basefrac = {} 
        self.basefrac['A'] = np.zeros(len(self.rep_seq)) 
        self.basefrac['T'] = np.zeros(len(self.rep_seq)) 
        self.basefrac['G'] = np.zeros(len(self.rep_seq)) 
        self.basefrac['C'] = np.zeros(len(self.rep_seq)) 
        self.basefrac['N'] = np.zeros(len(self.rep_seq)) 
        
        # Cheack representative sequence
        for i in range(len(self.rep_seq)):
            self.basefrac[self.rep_seq[i]][i] += 1
            
        # Check Cluster members
        for seq in self.members_seq:
            for i in range(len(self.rep_seq)):  
                self.basefrac[seq[i]][i] += 1
        
        for i in range(len(self.rep_seq)):
            self.basefrac['A'][i] /= float(self.size)
            self.basefrac['T'][i] /= float(self.size)
            self.basefrac['G'][i] /= float(self.size)
            self.basefrac['C'][i] /= float(self.size)
            self.basefrac['N'][i] /= float(self.size)
            
    def plot_basefrac(self, seq_start_idx=6, xlab="", ylab="", title="", **kwargs):
        """ Visualise the fraction of bases """
        
        import matplotlib.pyplot as plt
    
        if not hasattr(self, 'basefrac'):
            self.get_basefraction()
    
        # Default plot options
        if 'ms' not in kwargs:
            kwargs['ms'] = 10.0
        if 'marker' not in kwargs:
            kwargs['marker'] = '.'
        if 'mew' not in kwargs:
            kwargs['mew'] = 0.0
    
        plt.figure()
        for k,v in self.basefrac.iteritems():
        
            data_xs = range(1,len(v[seq_start_idx:])+1) 
            data_ys = v[seq_start_idx:]
            label =  k 
            
            vars()['line' + k], = plt.plot(data_xs, data_ys, label=label, ls='', **kwargs)
    
        # Shink current axis's height by 10% on the bottom
        ax = plt.gca()
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])
    
        # Put a legend below current axis
        
        ax.legend([vars()['lineA'], vars()['lineT'], vars()['lineG'], vars()['lineC'], vars()['lineN']],
                  ['A', 'T', 'G', 'C', 'N'], loc='upper center', bbox_to_anchor=(0.5, -0.05),
                  fancybox=True, shadow=True, ncol=5, numpoints=1, markerscale=3)
        
        plt.ylim(-0.1, 1.1)
        plt.xlim(1, 101)
        plt.title(title)
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.show()
        
    def plot_ATGC(self, db=None):
        ''' Visualise the presence of nucliotides at each base position.

        TODO: Broken, needs to be fixed to work with new get_unique_seq
        '''
        
        import matplotlib.pyplot as plt
        

        # get the unique sequences 
        if not hasattr(self, 'unique_seqs'):
            self.get_unique_seq(seq_start_idx=6, db=db)


        
        # Get dimentions
        n = len(self.unique_seqs)
        m = len(self.unique_seqs[self.unique_seqs.keys()[0]])
        
        y_lab = range(n)
        x_lab = range(m)
        
        fig, ax = plt.subplots()

        data = np.zeros([n, m], dtype = int)

        for i, (seq, count) in enumerate(self.unique_seqs.most_common()):
            
            for j, base in enumerate(seq):
                
                if base == 'A':
                    data[i,j] = 1 


#         ax.imshow(image, cmap=plt.cm.gray, interpolation='nearest')
#         ax.set_title('dropped spines')
# 
#         
#         # Move left and bottom spines outward by 10 points
#         ax.spines['left'].set_position(('outward', 10))
#         ax.spines['bottom'].set_position(('outward', 10))
#         # Hide the right and top spines
#         ax.spines['right'].set_visible(False)
#         ax.spines['top'].set_visible(False)
#         # Only show ticks on the left and bottom spines
#         ax.yaxis.set_ticks_position('left')
#         ax.xaxis.set_ticks_position('bottom')
#         
        plt.show()
        
              
    def write2file(self, handle, db=None, format='fasta', start_idx=6, only_unique=False):
        """ Write cluster to one of the following formats:
              1.  .clstr file format like CD-HIT.
              2.  fasta or fastq

        CHHIT cluster file example:

        >Cluster 0
        0    88nt, >2096... *
        1    88nt, >25159... at +/100.00%
        2    88nt, >52830... at +/98.86%
        3    88nt, >76269... at +/100.00%
        4    88nt, >12309... at +/100.00%
        
        """

        # Input Checks
        if type(handle) is file:
            if handle.closed:
                    handle = open(handle.name, 'a')
        elif type(handle) is str:
            if format == 'cdhit':
                if not handle.endswith('.clstr'):
                    filename = handle + '.clstr'
            elif format == 'fasta':
                if not handle.endswith('.fasta'):
                    filename = handle + '.fasta'
            elif format == 'fastq':
                if not handle.endswith('.fatstq'):
                    filename = handle + '.fastq'
            else:
                raise Exception('Invalide file format specified.')

            handle = open(handle, 'a')

        # Make sure data is available to write
        if format == 'cdhit' or format == 'fasta':
            if not self.members_seq or not self.members_sample_id:
                print "Sequence data not present in cluster. Retrieving from data base..."
                self.getfromdb(target='all', items=['seq', 'seqid'], db=db)
        elif format == 'fastq':
            if not self.members_seq or not self.members_sample_id or not self.members_phred:
                print "Sequence data not present in cluster. Retrieving from data base..."
                self.getfromdb(target='all', items=['seq', 'seqid', 'phred'], db=db)

        # Derived vars
        seq_len = str(len(self.rep_seq[start_idx:]))

        if format == 'cdhit':
            with handle as clust_file:
                # Write Header
                clust_file.write('>Cluster {0}\n'.format(self.id))

                # Write representative sequence
                clust_file.write('0\t{length}nt, >{rep_seq_id}... *\n'.format(length=seq_len ,
                                                                              rep_seq_desc=self.rep_seq_id))
                count = 0
                # Write rest of cluster members
                for idx, desc in enumerate(self.members_id):

                    count += 1

                    # Calculate percentage similarity
                    mismatches = distance(self.members_seq[idx][12:], self.rep_seq[12:])
                    percentage = 100.00 - mismatch2percentage(mismatches, seq_len)
                    clust_file.write('{count}\t{length}nt, >{seq_desc}... at +/{percentage:.2f}%\n'.format(
                        count=str(count), length=seq_len, seq_desc=desc, percentage=percentage))

        elif format == 'fastq' or format == 'fasta':

            # Get list of sequence Records in fasta format.
            allSeqRecs = [
                SeqRecord(Seq(self.rep_seq[start_idx:]), id=str(self.rep_seq_id), description=str(self.rep_sample_id))
            ]

            for i in range(len(self.members_seq)):
                rec = SeqRecord(Seq(self.members_seq[i][start_idx:]), id=str(self.members_id[i]),
                                description=str(self.members_sample_id[i]))
                allSeqRecs.append(rec)

            SeqIO.write(allSeqRecs, handle, format=format)

                
class ClusterFilter(object):
    """ Holds all methods for parsing and filtering a cluster file from CD-HIT
    
    Impliments a two-pass generator to parse reads:
    1. Get the next cluster size, and store start and end places. 

    2. If size passes, go back and rescan file to get detailed info for that cluster.
    """
    
    def __init__(self, handle, idx_file_path, db, filter_params,
                 reads_per_cluster_size_counter, output_dirpath=None):

        # IO
        handle = input_check(handle)
        self.handle4parser = handle
        self.handle4get_desc = open(handle.name, 'rb')
        
        if output_dirpath is None:
            self.output_dirpath = os.getcwd()
        else:
            if not os.path.exists(output_dirpath):
                os.makedirs(output_dirpath)
            self.output_dirpath = output_dirpath
            
        # Store internal params        
        self.idx_file_path = idx_file_path
        self.idx_file_dir = os.path.split(idx_file_path)[0]
        self.db = db
        self.filter_params = filter_params
        self.reads_per_cluster_size_counter = reads_per_cluster_size_counter
        # Setup generator
        self.cluster_size_gen = self.__setup_cluster_size_parser()
        
    def __setup_cluster_size_parser(self):
        
        first = True
        # Setup Data structure
        cluster = ClusterObj()        
        # Cycle through file 
        seq_count = 0    
        try:
            with self.handle4parser as cluster_file:  
                
                while True:
                    line = cluster_file.readline()  # Reading one line from file including end-of-line char '\n'
                    if not line: 
                        break
                
                    if line.startswith('>'):
                        # This is start of new cluster
                        # yield last cluster if not first pass
                        if not first:
                            # Record Cluster info 
                            
                            start_of_line_position = cluster_file.tell() - len(line)
                            cluster.size = seq_count
                            cluster.end_loc = start_of_line_position
                            
                            yield cluster
                            
                            # Reset Cluster and counter
                            cluster = ClusterObj()
                            seq_count = 0
                            cluster.start_loc = start_of_line_position
                        else: 
                            first = False
                    else:
                        seq_count += 1
                            
                # Got to end of file but still one more cluster to add
                cluster.size = seq_count
                yield cluster
        
        except IOError:
            print "Error: can\'t find file or read data"
        else:
            print "Finished Scanning cluster file."    

    def get_desc(self, cluster):
        """ Scan back over cluster and get sequence descriptions """
        
        cluster_file = self.handle4get_desc

        # Update info 
        cluster.idx_file_path = self.idx_file_path

        # Load cluster 
        cluster_file.seek(cluster.start_loc)
        
        first = True
        for line in cluster_file:              
            line = line.strip() # Remove white space at start and end
            
            if line.startswith('>'):
                # This is start of new cluster
                
                # yield this cluster
                if first:
                    cluster.id = line.split()[1]
                    first = False
                else: 
                    break
                
            elif line.endswith('*'): 
                # This is the representative sequence for the cluster
                cluster.rep_seq_id = line.split()[2].strip('>.')
            else: 
                line_parts = line.split()
                next_desc = line_parts[2].strip('>.')
                cluster.members_id.append(next_desc)
                
        if os.getcwd() != self.output_dirpath:
            os.chdir(self.output_dirpath)
                
        return cluster    

    def run_filter(self):
        """ Run through file and write desired clusters to fastq file """
 
        size_counter = Counter()
 
        count = 0
 
        for cluster in self.cluster_size_gen:
            
            count += 1 
            
            if cluster.size >= self.filter_params['min_size'] and cluster.size <= self.filter_params['max_size']:
                
                # Get seq descriptions in cluster 
                cluster = self.get_desc(cluster)
                
                # Write to output file
                size_counter[str(cluster.size)] += 1
                
                # Get the sequence records for the cluster 
                seqs = []
                if os.getcwd() != self.idx_file_dir:
                    os.chdir(self.idx_file_dir)
                seqs.append(self.db[cluster.rep_seq_id])
                for member in cluster.members_id:
                    seqs.append(self.db[member])
                
                if os.getcwd() != self.output_dirpath:
                    os.chdir(self.output_dirpath)
                    
                # Write cluster to a file 
                fname = "clustersize{0}-No{1}.fastq".format(str(cluster.size), str(size_counter[str(cluster.size)]))
                output_handle = open(fname, "wb")
                SeqIO.write(seqs, output_handle, "fastq")
                
            elif self.reads_per_cluster_size_counter[str(cluster.size)] < self.filter_params['min_reads']:
                break
        
        print count


#==============================================================================
# Utility functions
#==============================================================================

def input_check(handle):
    """ Check if input is valid and or a string for a file name.
    Always returns a ahndle for an an open file.
    
    handle = input_check(handle)
    """
    if type(handle) is str:
        handle = open(handle, 'rb')
    else:
        assert type(handle) is file, 'Invalid file handle.'
        if handle.closed:
            handle = open(handle.name, 'rb')

    return handle 

def getfromdb(cluster_list, items, target, db=None):
    """ Adds the specified items to each cluster in list, looked up from the main db.
    
    items - list of things to fetch for the cluster. Option of 'seq', 'phred' or 'seqid'

    target - 'rep' will fetch just the items for the representative sequence,
             'all' will fetch the items for all cluster members includeing representative sequence.
            
    """
    for cluster in cluster_list:
        cluster.getfromdb(items=items, target=target, db=db)
        
    return cluster_list

def mismatch2percentage(mismatches, seq_length):
    ''' Convert a number of mismatches to a percentage rounded to the nearest 2 dp '''
    mismatches = int(mismatches)
    seq_length = int(seq_length)
    return  round((float(mismatches) / seq_length) * 10000) / 100 

def percentage2mismatch(percentage, seq_length):
    ''' Convert a percentage similarity to a number of mismatches (rounded to integer)'''
    percentage = float(percentage)
    seq_length = int(seq_length)
    return  int(round( (percentage / 100.) * seq_length))
                   
def parse(handle, db=None, edit_dist=False, editdist_count=False):
    """ Reads in a CDHIT cluster file and returns a generator for the clusters. 
    
     - handle   - handle to the file, or the filename as a string
     - db       - a Reads_db database object. If present 
     - edit_dist - Only calculated if set to true.
     - similarity_count - Will return a counter for similarity to rep seq 

    Currently iterates in the order the file is read. 
    Typical usage, opening a file to read in, and looping over all record(s):
    """
    
    # input checks    
    handle = input_check(handle)
        
    first = True
    # Setup Data structure
    cluster = ClusterObj()
    
    if editdist_count:
        cluster.editdist_counter = Counter()
    
    # Cycle through file     
    try:
        with handle as cluster_file:   
            
            for line in cluster_file:              
                line = line.strip() # Remove white space at start and end
                
                if line.startswith('>'):
                    # This is start of new cluster
                    
                    # Claculate stats for last cluster
                    cluster.size = len(cluster.members_id) + 1
                    
                    # Store number for next yield 
                    next_cluster_num = line.split()[1]
                    
                    # yeild this cluster
                    if first:
                        cluster.id = line.split()[1]
                        first = False
                    else: 
                        yield cluster
                        # Reset Cluster
                        cluster = ClusterObj()
                        cluster.id = next_cluster_num
                        if editdist_count:
                            cluster.editdist_counter = Counter()
                    
                elif line.endswith('*'): 
                    # This is the representative sequence for the cluster
                    cluster.rep_seq_id = int(line.split()[2].strip('>.'))
                else:
                    
                    line_parts = line.split()
                    
                    cluster.members_id.append(int(line_parts[2].strip('>.')))
                    if edit_dist or editdist_count:
                        similarity = line_parts[4].strip('+/%')
                        seq_len = line_parts[1].strip('nt,')
                        edit_dist = percentage2mismatch( 100 - float(similarity), seq_len)
                        cluster.edit_dists.append(edit_dist)
                        if editdist_count:
                            cluster.editdist_counter[edit_dist] += 1
                    
            # Got to end of file but still one more cluster to add
            cluster.size = len(cluster.members_id) + 1
            yield cluster
    
    except IOError:
        print "Error: can\'t find file or read data"
    else:
        print "Finished Scanning cluster file."    

    
def sortby(handle, reverse=True, mode='cluster_size', outfile_postfix=None, 
           cutoff=None, clustersize_min=None, clustersize_max=None):
    """ Reads in a CDHIT cluster file and writes a sorted file for the clusters.
    by their size.  
    
     - handle   - handle to the file, or the filename as a string.
     - reverse  - Python sorts in ascending order by default, if reverse is true
                  the sort will instead be in descending order.  
     - mode     - "cluster_size" = Sort by the cluster sizes
                  "reads_per_cluster" = Sort by total reads in clusters of each size

     - cutoff   - Minimum threshold for sorted attribute. Any clusters with smaller values 
                  are not written to file.

    Works by scanning file once through to build an index, sorting the index, then 
    rewriting file accordingly. 
    """

    if mode == 'reads_per_cluster':
        reads_per_cluster_counter = summary_counter(handle, mode='reads_per_cluster', report=True)
        
    # Sorting key functions
    key_ops = {
               'cluster_size': lambda x : x[0],
               'reads_per_cluster': lambda x : reads_per_cluster_counter[str(x[0])], 
               }
    # Input checks    
    handle = input_check(handle)

    # Setup Data structure and utility vars
    # cluster_idxs = [ (size, start_bytes, end_bytes ), ... ]
    cluster_idxs = []
    start_time = time.time()
    first = True
    # First cycle through file and store all cluster start and end positions 
    try:
        with handle as cluster_file:   
            
            size = 0
            clust_end = 0
            clust_start = cluster_file.tell()
            
            while True:
                line = cluster_file.readline()  # Reading one line from file including end-of-line char '\n'
                if not line: 
                    break
                
                if first:
                    first = False
                    continue
                
                if line.startswith('>'):
                    # This is start of a new cluster
                    
                    # Record results for last cluster 
                    cluster_idxs.append((size, clust_start, clust_end))
                    # Reset cluster start and size counter
                    clust_start = clust_end
                    size = 0
                else:
                    size += 1
                    clust_end = cluster_file.tell()
                    
            # Got to end of file but still one more cluster to add
            cluster_idxs.append((size, clust_start, clust_end)) 
            
            # Sort index in the desired fashion 
            sorted_cluster_idxs = sorted(cluster_idxs, key=key_ops[mode], reverse=reverse)    
                
            # Rewrite file
            if outfile_postfix is None: 
                outfile_postfix = '-sortedby_' + mode
                
            old_filename_parts = cluster_file.name.split('.')
            old_filename_parts[0]  += outfile_postfix
            new_filename = '.'.join(old_filename_parts)
                
            with open(new_filename, 'wb') as sorted_cluster_file :
                
                if cutoff and clustersize_min and clustersize_max:
                    # Filter the clusters that are written 
                    for tup in sorted_cluster_idxs:
                        if key_ops[mode](tup) >= cutoff:
                            if tup[0] >= clustersize_min and tup[0] <= clustersize_max:
                                cluster_file.seek(tup[1])
                                cluster_text = cluster_file.read(tup[2] - tup[1])
                                sorted_cluster_file.write(cluster_text)
                        else:
                            break
                elif cutoff:
                    for tup in sorted_cluster_idxs:
                        # Stop writing after minimum cutoff is reached
                        if key_ops[mode](tup) >= cutoff:
                            cluster_file.seek(tup[1])
                            cluster_text = cluster_file.read(tup[2] - tup[1])
                            sorted_cluster_file.write(cluster_text)
                        else:
                            break
                else:
                    for tup in sorted_cluster_idxs:
                        # Write all clusters to file
                        cluster_file.seek(tup[1])
                        cluster_text = cluster_file.read(tup[2] - tup[1])
                        sorted_cluster_file.write(cluster_text)
                        
    except IOError:
        print "Error: can\'t find file or read data"
    else:
        t = time.time() - start_time
        print "Finished sorting cluster file after {0}\n".format(time.strftime('%H:%M:%S', time.gmtime(t)))    

    return sorted_cluster_file


def summary_counter(handle, mode='cluster_size', report=True):
    ''' Takes cluster file output by CD-Hit and produces two Counter dictionaries 
    output depends on mode specified:
    
    - handle   - handle to the file, or the filename as a string.
    - mode:
        'by_seqlen'         = { 'sequence_length' : Counter(cluster_size) }
        'cluster_size'      = Counter(cluster_sizes_for_all_sequence_lengths)
        'reads_per_cluster' = { 'cluster size'  : cluster_size * Num clusters }
    '''

    # input checks    
    if type(handle) is str:
        if not handle.endswith('.clstr'):
            handle = handle + '.clstr'
        handle = open(handle, 'rb')
    else:
        assert type(handle) is file, 'Invalid file handle.'
        if handle.closed:
            handle = open(handle.name, handle.mode)
            
    assert mode in ['by_seqlen', 'cluster_size', 'reads_per_cluster'], 'Mode not recognised.'
            
    # Data structure to store cluster size info is a DefaultDictionary of Counter dictionaries.
    # ds = { 'seq_len' : Counter(cluster_size)  }
    # empty keys of ds are initialised with a Counter dictionary. 
    
    ds = defaultdict(Counter)

    seq_in_cluster = 0
    rep_length = 0

    print 'Generating cluster summary for  %s ...' % (handle.name)

    try:
        with handle  as cluster_file:   
            
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
    if mode == 'cluster_size':
        return total_cluster_size_counter
    elif mode == 'by_seqlen':
        return ds
    elif mode =='both':
        return total_cluster_size_counter, ds
    elif mode == 'reads_per_cluster':
        reads_per_cluster = {}
        for k,v in total_cluster_size_counter.iteritems():
            reads_per_cluster[k] = int(k) * v
        return reads_per_cluster
    

def plot_counters_scatter(counters, labels=None, log='xy', xlab="", ylab="", title="", figsize=(6,6), **kwargs):
    ''' Construct a series of scatter plots from a list of Counter Dictionaries '''
    
    import matplotlib.pyplot as plt
    
    # Default plot options
    if 'ms' not in kwargs:
        kwargs['ms'] = 4.0
    if 'marker' not in kwargs:
        kwargs['marker'] = '.'
    if 'mew' not in kwargs:
        kwargs['mew'] = 0.0
    
    if type(counters) is not list and type(counters) is not tuple:
        counters = [counters]
    
    if labels is not None:
        assert len(labels) == len(counters), "Number of labels must match number of counters."

    plt.figure(figsize=figsize)

    for i in range(len(counters)):
        
        data_xs = [int(k) for k in counters[i].keys()]
        data_ys = counters[i].values()

        if labels:
            plt.plot(data_xs, data_ys, label=labels[i], ls='', **kwargs)
        else:
            plt.plot(data_xs, data_ys, label='Counter-'+str(i), ls='', **kwargs)
        
        ax = plt.gca()
        if 'x' in log:
            ax.set_xscale('log')
        if 'y' in log:
            ax.set_yscale('log')
        
        plt.title(title)
        plt.xlabel(xlab)
        plt.ylabel(ylab)

    plt.legend(numpoints=1, markerscale=3)
    plt.show()
    

def plot_counters_hist(counters, bin_width=100, labels=None, log='xy', xlab="", ylab="", title="", calc_rpc=False, **kwargs):
    ''' Construct a series of histogram plots from a list of Counter Dictionaries.  
    
    counter dictionaries should be number of total reads per cluster size'''
    
    import matplotlib.pyplot as plt
    
    # Default plot options
#     if 'ms' not in kwargs:
#         kwargs['ms'] = 4.0
#     if 'marker' not in kwargs:
#         kwargs['marker'] = '.'
#     if 'mew' not in kwargs:
#         kwargs['mew'] = 0.0
    
    if type(counters) is not list and type(counters) is not tuple:
        counters = [counters]
    
    if labels is not None:
        assert len(labels) == len(counters), "Number of labels must match number of counters."
    
    for i, counter in enumerate(counters):   
        
        # Compact to hist counter  
        hist_counter = Counter()
        for size, freq in counter.iteritems():
            bin = np.floor(float(size)/bin_width)
            hist_counter[bin] += freq
        
        # Extract data 
        bags = sorted(hist_counter.items())
        x_data, y_data = zip(*bags)
        x_data = np.array(x_data) * bin_width  
        y_data = np.array(y_data)
        
        plt.step(x_data, y_data)
        
    # Plot formating 
    ax = plt.gca()
    if 'x' in log:
        ax.set_xscale('log')
    if 'y' in log:
        ax.set_yscale('log')
    
    plt.title(title)
    plt.xlabel(xlab)
    plt.ylabel(ylab)

    plt.legend(numpoints=1, markerscale=8)
    plt.show()

    return x_data, y_data


if __name__ == '__main__':
    gaz_clustf = '/space/musselle/data/RAD-seq/gazelles-zebras/clusters/gz_allg_allz_95g1_clustered_reads_c95_g1/gz_allg_allz-gazelle-clustered_c95_g1.clstr'
    zeb_clustf = '/space/musselle/data/RAD-seq/gazelles-zebras/clusters/gz_allg_allz_95g1_clustered_reads_c95_g1/gz_allg_allz-zebra-clustered_c95_g1.clstr'
    idx_file = '/space/musselle/data/RAD-seq/gazelles-zebras/processed-data/all_reads_fastq.idx'

    output_dir = '/space/musselle/data/RAD-seq/gazelles-zebras/clusters/gz_allg_allz_95g1_clustered_reads_c95_g1/gazellles_filtered/'

