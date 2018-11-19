'''
A class to handle collection of isolates

This class has methods to output the nexus file
'''

from Bio import AlignIO
import collections
import pandas as pd
import numpy as np
import time

from beastify.nexus.Isolate import Isolate
from beastify import __VERSION__ as version_string

class Collection:
    def __init__(self):
        self.isolates = collections.OrderedDict({})
        self.aln_file = None
        self.seq_path = None
    def __getitem__(self, key):
        return self.isolates[key]
    def __str__(self):
        n_isolats = len( list(self.isolates.keys()) )
        if n_isolates  == 0:
            print( 'An empty collection of isolates' )
        else:
            print(( 'A collection of {} isolates.'.format( n_isolates )))
    def load_isolates(self, path, seq_file, ignore = None):
        '''
        Takes a path, and searches for snippy output files with the
        following name: snps.consensus.subs.fa
        Assumes that each isolate has its own directory named with
        its ID
        '''
        #default_file = "snps.consensus.subs.fa"
        default_file = seq_file
        list_reads = os.walk(path)
        for root, dirs, files in list_reads:
            for d in dirs:
                if d in ignore:
                    print(("Skipping {}".format(d)))
                    continue
                print(("Trying to load {} at {}".format(d, os.path.join(path, d, default_file))))
                tmp_fn = os.path.join(path, d, default_file)
                self.isolates[d] = Isolate()
                self.isolates[d].load_fasta(tmp_fn, d)
                print("Successfully loaded: ")
                print((self.isolates[d]))
            break
        self.seq_path = path
        print("Finished loading isolates!")
        return
    def load_alignment( self, aln_file, aln_format ):
        '''
        Loads an alignment file. Tries to figure out whether the sizes are correct,
        and how many polymorphic sites there.

        Supported formats, any of those supported by BioPython:AlignIO

        Tested formats, at the moment only multifasta
        '''
        # try to open alignment file
        try:
            aln_open = open( aln_file, 'rt' )
        except IOError:
            print(( "Could not open file {}.".format( aln_file )))
            raise
        except:
            print(( "Something went wrong while trying to open file {}.".format( aln_file )))
            raise
        # try to parse alignment file
        try:
            aln = AlignIO.read( aln_open, aln_format )
        except ValueError:
            print(( "Opened the file {}, but could not load the alignment. Is the format {} correct?".format( aln_file, aln_format )))
            raise
        except:
            print(( "Something happend while trying to parse {}.".format( aln_file )))
            raise
        self.aln_file = aln_file
        aln_list = aln
        for s in aln_list:
            print(( "Trying to load isolate {} into collection...".format( s.id )), end=' ')
            self.isolates[s.id] = Isolate()
            self.isolates[s.id].load_seqRec( s, s.id )
            print(( "\033[92m" + "OK" + "\033[0m"))
        print(( "Successfully loaded {} sequences from file {}.".format( len( list(self.isolates.keys()) ), aln_file ) ))
        aln_open.close()
        return
    def load_reference(self, reference):
        '''
        Loading the reference to the collection to make sure it is included in
        the final file
        '''
        self.isolates["Reference"] = Isolate()
        self.isolates["Reference"].load_seqRec(reference, "Reference")
        return
    def add_info(self, info_file):
        '''
        A function to add dates and other ifnormation to the sample ids
        '''
        fn = open(info_file, 'r')
        info = [l.strip().split("\t") for l in fn]
        fn.close()
        self.info = {}
        for i in info:
            isolate = i[0]
            if len(i) != 2:
                iso_info = ":".join(i[1:])
            else:
                iso_info = i[1]
            self.info[isolate] = iso_info

    def _aln2df(self):
        '''
        Transform the data from different isolates in to a pandas.DataFrame
        '''
        self.aln_pd = pd.DataFrame(np.array([list(self.isolates[rec].seq) for rec in self.isolates], order="F"),
                     index=[self.isolates[rec].id for rec in self.isolates])


    def make_nexus( self, outfile, gene_obj, subsample = None, subsample_seed = 42, partitions = [1,2,3,4,5], mask = None ):
        '''
        This function will take the indexed positions, and make
        an alignment pandas DataFrame.
        Ideally, indexed positions DataFrame would also include a locus_tag column, and an include column, which would have True/False for each position, depicting wether that position should be included or not. Positions not to be included might be recombinant sites, mobile element sites, or any other sites that should be masked.

        Pseudocode:
        1. Create alignment pandas DataFrame
        2. At this point, indexed positions DataFrame would only contain those positions to be included, no further processing required here. --- need to add methods to filter positions
        3. Pick from the alignment all positions grouped by codon_pos --- here we sort positions by codon_pos, and just use this as an index on the columns of 1.
        4. Generate the individual sequences for printing to nexus file
        5. Figure out min and max index for each codon_pos type after ordering in order to print out the charset block
        '''
        # position index
        #partitions = [1, 2, 3, 4, 5]
        pos_index = gene_obj.indexed_positions
        # create a pandas DataFrame of the alignment
        self._aln2df()
        # sort the columns so codon position categories are contiguous
        #import pdb; pdb.set_trace()
        if mask != None:
            print('Masking sites')
            positions = pd.read_table(mask, comment = '#',header=None)
            exclude = []
            for i in positions.iterrows():
                exclude.extend(np.arange(i[1][1],i[1][2]+1))
            print(('Found {} sites to mask.'.format( len(exclude) ) ))
            pos_index = pos_index.loc[~pos_index.positions.isin(exclude)]
            print(('New alignment has {} sites.'.format(pos_index.shape[0])))
        pos_index = pos_index.sort_values( ['codon_pos', 'positions'])
        new_column_order = pos_index.index
        self.aln_pd = self.aln_pd.iloc[:, new_column_order ]

        # figure out the ranges for each partition
        pos_index = pos_index.reset_index() # this generates a new index for the sorted DataFrame
        if partitions != [1,2,3,4,5]:
            #import pdb; pdb.set_trace()
            ix = pos_index['codon_pos'].isin(partitions)
            pos_index = pos_index.loc[ix,:]
            self.aln_pd = self.aln_pd.loc[:,pos_index['index']]
            #pos_index = pos_index.reset_index()
        if subsample != None:
            # randomly sample columns within partitions
            # there is currently no way to check whether the requested number
            # of columns exceeds the number of columns in a particular partition
            np.random.seed(subsample_seed)
            fn = lambda obj: obj.loc[np.random.choice(obj.index, subsample, False),:]
            pos_index = pos_index.groupby( 'codon_pos', as_index = False).apply(fn)
            pos_index = pos_index.reset_index()
            self.aln_pd = self.aln_pd.loc[:,pos_index['index']]
            pos_index = pos_index.groupby( 'codon_pos' )
        else:
            pos_index = pos_index.groupby( 'codon_pos' ) # now, grouping by codon_pos will generate list
                                                     # index positions
        partition_ranges = {}
        for partition in partitions:
            partition_ranges[ partition ] = {}
            partition_ranges[ partition ]['min'] = min( pos_index.groups[ partition ] ) + 1
            partition_ranges[ partition ]['max'] = max( pos_index.groups[ partition ] ) + 1
        # prep nexus file
        sp = "    "
        out = "#NEXUS\n"
        out += "[Data from:\n"
        out += "beastify.py version {}\n".format( version_string )
        out += "Date: {}\n".format( time.strftime("%d/%m/%Y") )
        out += "Reference: {}\n".format( gene_obj.reference.id )
        if (subsample != None):
            out += "Random seed for subsampling columns: {}\n".format(subsample_seed)
        if (self.aln_file != None):
            out += "Isolate data was taken from alignment file: {}\n".format(self.aln_file)
        else:
            out += 'Isolate data was parsed from folder: {}\n'.format(self.seq_path)
        out += "]\n\n"
        out+= "begin taxa;\n"
        out += sp + "dimensions ntax={};\n".format(len(list(self.isolates.keys())))
        out += sp + "taxlabels\n"
        try:
            for i in list(self.isolates.keys()):
                out += i + ":" + self.info[i] + "\n"
        except:
            for i in list(self.isolates.keys()):
                out += i + "\n"
        out += sp + ";\n"
        out += "end;\n\n"
        out += "begin characters;\n"
        out += sp + "dimensions nchar={};\n".format( self.aln_pd.shape[ 1 ])
        out += sp + "format missing=? gap=- datatype=dna;\n"
        out += sp + "gapmode=missing;\n"
        out += sp + "matrix\n"
        #same as above when outputting the tax labels
        #import pdb; pdb.set_trace()
        for i in self.aln_pd.index:
            try:
                ident = i + ":" + self.info[i]
            except:
                ident = i
            out += "{:20} {}\n".format(ident, ''.join( self.aln_pd.loc[ i, : ] ) )
        out += sp + ";\n"
        out += "end;\n\n"
        out += "begin assumptions;\n"
        # added the sort to make sure the genes are in the same order in which
        # they were concatenated
        for partition in sorted(partition_ranges.keys()):
            out += sp + "charset partition_{} = {}-{};\n".format( partition, partition_ranges[ partition ][ 'min' ], partition_ranges[ partition ][ 'max' ])
        out += "end;\n"
        fn = open(outfile, 'w')
        fn.write(out)
        fn.close()
        print("All done! Happy BEASTing!")
    def gen_align(self, outfile, gene_obj):
        '''
        iterates over the keys in the self.__dict__ to return
        individual gen alignments, and format the output so it is
        suitable for BEAST.
        '''
        #import pdb; pdb.set_trace()
        concat_seqs = {}
        gene_lens = {}
        for i in list(self.isolates.keys()):
            genome = self.isolates[i]
            genome.get_genes(gene_obj)
            concat_seqs[i] = ""
            for g in sorted(genome.genes.keys()):
                concat_seqs[i] += genome.genes[g]
                try:
                    gene_lens[g]
                    if (gen_lens[g] != len(genome.genes[g])):
                        print(("Gene {} has different length in isolate {} than in previously recorded: {}.".format(g, len(genome.genes[g], gen_lens[g]))))
                except:
                    gene_lens[g] = len(genome.genes[g])
                    pass
        seq_len = len(concat_seqs[list(concat_seqs.keys())[0]])
        var_sites = 0
        for n in range(seq_len):
            tmp = []
            for i in concat_seqs:
                tmp.append(concat_seqs[i][n])
            tmp = [t for t in tmp if t not in ['-']]
            tmp = set(tmp)
            if len(tmp) > 1:
                var_sites += 1
                # print("Site {} is variable.".format(n))
                # print(tmp)
        print(("Total variable sites found: {}.".format(var_sites)))
        # outputting nexus file
        sp = "    "
        out = "#NEXUS\n"
        out += "[Data from:\n"
        out += "beastify.py\n"
        out += "]\n\n"
        out+= "begin taxa;\n"
        out += sp + "dimensions ntax={};\n".format(len(list(self.isolates.keys())))
        out += sp + "taxlabels\n"
        # if self.dates exists, then add it to the id name
        # otherwise, ignore
        # for i in self.dates:
        #     print(i, self.dates[i])
        try:
            for i in list(self.isolates.keys()):
                out += i + ":" + self.info[i] + "\n"
        except:
            for i in list(self.isolates.keys()):
                out += i + "\n"
        out += sp + ";\n"
        out += "end;\n\n"
        out += "begin characters;\n"
        out += sp + "dimensions nchar={};\n".format(len(concat_seqs[list(concat_seqs.keys())[0]]))
        out += sp + "format missing=? gap=- datatype=dna;\n"
        out += sp + "gapmode=missing;\n"
        out += sp + "matrix\n"
        # same as above when outputting the tax labels
        try:
            for i in concat_seqs:
                ident = i + ":" + self.info[i]
                out += "{:20} {}\n".format(ident, concat_seqs[i])
        except:
            for i in concat_seqs:
                out += "{:20} {}\n".format(i, concat_seqs[i])
        out += sp + ";\n"
        out += "end;\n\n"
        # out += "begin assumptions;\n"
        # start = 1
        # end = 0
        # # added the sort to make sure the genes are in the same order in which
        # # they were concatenated
        # for g in sorted(gene_lens.keys()):
        #     end += gene_lens[g]
        #     out += sp + "charset {} = {}-{};\n".format(g, start, end)
        #     start = end + 1
        # out += "end;\n"
        fn = open(outfile, 'w')
        fn.write(out)
        fn.close()
        print("All done. Happy BEASTing!")
        # creating an nexus file
        return
