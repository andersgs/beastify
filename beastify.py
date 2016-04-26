'''
beastify.py

Anders Goncalves da Silva
October 2015

Updated on 7 January 2016
'''

import click
import collections
from Bio import SeqIO
from Bio import AlignIO
import os
import random
import re
import sys
import unittest
import warnings
import pandas as pd
from pandas.util.testing import assert_frame_equal as pd_assert_df_equal
import numpy as np
import pdb
import time

VERSION=0.1

class Genes:
    '''
    This class holds information about features.

    It provides functions to index positions in the genome according to four different partitions:
        1. First codon position
        2. Second codon position
        3. Third codon position
        4. Shared codon position --- for positions that are shared with multiple genes
        5. Other positions --- positions not annotated, or with annotations other than CDS
    '''
    def __init__(self):
        self.features = {}
        self.snippy_list = []
    def load_genome(self, path ):
        '''
        This functions loads the genome data from a Genbank file

        Currently, it will only process a single locus. It should allow for multiple loci in the future.
        '''
        try:
            genome = SeqIO.read(path, 'genbank')
        except IOError:
            print("Could not open file {}".format(path))
            raise
        except ValueError:
            # introduced this to mainly deal with the case when there a GENBANK
            # file has a Chromosome and one or more plamids
            # possibly not the most general approach but it solves the issue
            # at hand.
            print("##### WARNING ######")
            print("Found more than record in the Genbank file.")
            print("Picking the first one, and discarding the rest")
            for g in SeqIO.parse(path, 'genbank'):
                genome = g
                break
        self.reference = genome # added to be able to add reference
    def load_features( self, gene_list = None, feature = 'CDS', nb = None, cod_tab = 11 ):
        '''
        This function will take all the genes in gene list or produced by parsing snippy core file and parse them out of the genome information.
        '''
        if gene_list != None:
            try:
                fn = open(gene_list, 'r')
                feature_set = set([line.strip() for line in fn])
                fn.close()
                # take a random sample if necessary
                if nb != None:
                    feature_set = random.sample(feature_set, nb)
            except IOError:
                print("Could not open file {}".format(feature_list))
        else:
            # if this function is run after parse snippy_core
            # otherwise, it will make an empty list. This will be a
            # a problem below
            feature_set = self.snippy_list

        # create a string for features to allow for some fuzzy pattern matching
        # of possible genes
        #feature_string = "|".join(feature_set).lower()
        feature_string = [f.lower() for f in feature_set]
        features_found = 0
        features_ok = 0
        ref_fasta = "ref.fa"
        ref_seq = ""
        #import pdb; pdb.set_trace()
        n_cds = 0
        n_gene = 0
        n_other = 0
        for feat in genome.features:
            try:
                if feat.qualifiers[feature][0].lower() in feature_string:
                    if feat.type == "CDS":
                        n_cds = n_cds + 1
                    elif feat.type == 'gene':
                        n_gene = n_gene + 1
                    else:
                        n_other = n_other + 1
            except:
                pass
            if feat.type == 'CDS':
                try:
                    # here we have to find based on some fuzzy set of possible
                    # genes. So, we use a string of features, and try to use
                    # the gene as the pattern
                    #import pdb; pdb.set_trace()
                    pat = feat.qualifiers[feature][0].lower()
                    #if re.search(feature_string, pat):
                    if pat in feature_string:
                        features_found += 1
                        try:
                            # test if it is a complete CDS
                            tmp = feat.extract(genome.seq).translate(cds=True, table = cod_tab)
                            self.features[feat.qualifiers[feature][0]] = feat
                            features_ok += 1
                            ref_seq += str(feat.extract(genome.seq))
                        except:
                            print("Incomplete CDS {}".format(feat.qualifiers[feature][0]))
                            print("Feature's length divided by 3 had remainder = {}".format(len(feat) % 3))
                            tmp = feat.extract(genome.seq)
                            print("Features first codon was {}".format(tmp[0:3]))
                            print("Features last codon was {}".format(tmp[-3:]))
                            print(str(tmp))
                            pass
                except:
                    pass
        print n_cds, n_gene, n_other
        print("Found {} features".format(features_found))
        print("of which {} were ok".format(features_ok))
        for k in self.features:
            print("Found gene {}".format(k))
        fn = open(ref_fasta, 'w')
        fn.write(">ref_seq\n")
        fn.write(ref_seq + "\n")
        fn.close()
    def parse_snippycore(self, path, corefn = "core.tab", nb = None, sample = 'random'):
        snippy_corefn = corefn
        path_tofile = os.path.join(path, snippy_corefn)
        if os.path.isfile(path_tofile):
            fn = open(path_tofile, 'r')
        else:
            raise IOError("Could not find file {}".format(path_tofile))
        snps = [line.strip().split("\t") for line in fn]
        header = snps[0]
        # assuming that the following headers are present
        chrom = header.index('CHR')
        pos = header.index('POS')
        ref = header.index('Reference')
        locus_tag = header.index('LOCUS_TAG')
        gene = header.index('GENE')
        product = header.index('PRODUCT')
        # keep only SNPs in coding regions
        snps = [s for s in snps[1:] if len(s) == len(header) and s[locus_tag] != '']

        # indexing the variable sites among the isolates
        ix_variable = [len(set([a for n,a in enumerate(snp_list) \
                if n not in [chrom, pos, ref, locus_tag, gene, product]])) > 1 \
            for snp_list in snps]
        # pull out variable sites among the isolates
        cds_var = list(set([s[locus_tag] for n,s in enumerate(snps) if ix_variable[n]]))
        # for k in sorted(count_snps_cds, key=count_snps_cds.get, reverse = False):
        #     print("{}, {}".format(k, count_snps_cds[k]))

        # give some choices about how to output coding regions
        # 1. random, meaning take a random subset of CDS that are variable
        # 2. top, order the CDS by number of variable sites, and take the top N entries
        # 3. None of the above, meaning give all the CDS **NOT RECOMMENDED** as
        #   the input file would be too long
        if sample == 'random':
            print("{}, {}, {}".format(len(cds_var), nb, len(cds_var) > nb))
            self.snippy_list = random.sample(cds_var, nb)
        elif sample == 'top':
            count_snps_cds = {}
            for cds in cds_var:
                try:
                    count_snps_cds[cds] += 1
                except:
                    count_snps_cds[cds] = 1
            self.snippy_list = sorted(count_snps_cds, key=count_snps_cds.get, reverse = True)[0:nb]
        else:
            self.snippy_list = cds_var
        print("Found {} variable genes, and picked {}.".format(len(cds_var), len(self.snippy_list)))
    def index_locations( self, feature_type = 'CDS' ):
        '''
        Function will index sites by codon position.
        There will be 5 categories:
            first_codon <- first codon position
            second_codon <- second codon position
            third_codon <- third codon position
            fourth_codon <- position is located in overlapping CDS regions and could be at different codon positions depending on which annotation is considered
            non_coding <- a position not found within an annotated CDS region

        The function will first annotate positions into one of the first four categories, and then figure out which positions are not annotated to put in non_coding.

        The code will cycle through features in self.genome (should allow for multiple genomes i.e., multi-genbank file)
        if feature.type == 'CDS',
        (will need here six lists, one for each category above, and another to keep track of all annotated positions (union of first-fourth codon positions))
        then figure out strand
            if positive, take location start = start and end = end ( should check what it is doing because Python is 0-base)
            else, take location start = end, end = start
        for each codon position, check
        if position is already in annotated_positions list
            check if already in fourth_codon,
                if not in fourth_codon but already in annotated list,
                    add to fourth_codon list, and deleted from appropriate list
                if presend in fourth_codon, do nothing
        else
            add to appropriate first-third_codon list, and to annotated_list

        for non_coding, create iterator with positions 1 to length of genome.
        iterate through returning only positions not in annotated_positions
        feature_start = genome.feature.location.start
        feature_end = genome.feature.location.end + 1 (we add the +1 because BioPython codes the positions as starting at 0, but we to have things starting at 1)
        THESE POSITIONS ARE THE SAME IF FEATURE IS IN POSITIVE OR NEGATIVE STRAND, BUT CODON POSITIONS ARE REVERSED
        if strand is positive (i.e., not reverse complement), then first, second, and third codon positions can be found thus:
            for codon_pos in 1:3: (codon positions are coded from 1 to 3, we can then add to the codon start position to have a natural, starting at 1, positions)
                codon_positions = range( feature_start + codon_pos, feature_end, 3 )
        else:
            for codon_pos in 1:3: (codon positions are coded from 1 to 3, BUT MEAN 3 TO 1 IN THIS CASE, we can then add to the codon start position to have a natural, starting at 1, positions)
                codon_positions = range( feature_start + codon_pos, feature_end, 3 )

        SOME ASSUMPTIONS:
            1. NO FUZZY POSITIONS --- START AND END OF FEATURES IS CODED WITH CERTAINY
            2. CODING REGIONS ALWAYS START AT POSITION -1 OR 1, SO STRAND IS ALWAYS -1 OR 1

        The naive implementation takes about 60 seconds to index about 100K  sites. Not ideal. It should be faster.
        A Pandas implementation might be faster:
        Here, we would mirror something close to the implementation done in R.
        '''
        # using a pandas approach
        # first step is to create a dataframe which unpacks individual CDS locations, and assigns codon positions
        # here, we assign codon positions on whether the strand is 1 or -1. with one getting counted as [1,2,3], and -1 being counted as [3,2,1] from start
        # at the moment, we assume that all CDS regions have a locus_tag --- that may not always be the case
        # the next step is to filter out duplicated positions, and assign these positions a codon position of 4
        gff_list = [ pd.DataFrame( {'positions': range( feature.location.start + 1, feature.location.end + 1 ),
                                    'codon_pos': ( [1,2,3] if feature.strand == 1 else [3,2,1] ) * ((feature.location.end - feature.location.start + 1)/3),
                                    'locus_tag': feature.qualifiers['locus_tag'][0] } ) for feature in self.reference.features if feature.type in feature_type ]
        gff_df = pd.concat( gff_list, ignore_index = True )
        # count the number of found features
        n_features = len( gff_df.groupby('locus_tag').groups )
        #find duplicated sites (needs at least pandas version 0.17)
        ix_dup=gff_df.duplicated( 'positions', keep = False )
        # make all codon_pos for duplicated sites = 4
        gff_df.loc[ix_dup,'codon_pos'] = 4
        # generate a deduplicated list
        gff_dedup = gff_df[['positions', 'codon_pos']].drop_duplicates('positions')
        # now create all positions dataframe
        all_pos = pd.DataFrame( {'positions': np.arange( 1, len( self.reference ) + 1 )} )
        # add codon positions, and fill nan with 5, and make integer
        # this produces a single data frame with two columns: 'positions', and 'codon_pos'. Positions is counted from 1 to length of chromosome, and 'codon_pos' is one of [1,2,3,4,5] depending on where the base is relative to one or more CDS annotations
        self.indexed_positions = pd.merge( all_pos, gff_dedup, how = 'left', on = 'positions' ).fillna(5).astype(int)
        #import pdb; pdb.set_trace()
        print( "Found {} features that match CDS".format( n_features ) )
        self.summary_index = self.indexed_positions.groupby('codon_pos').size()
        print( "Found {} sites in the first codon position".format( self.summary_index[1]))
        print( "Found {} sites in the second codon position".format( self.summary_index[2]))
        print( "Found {} sites in the third codon position".format( self.summary_index[3]))
        print( "Found {} sites in that are in multiple CDS annotations".format( self.summary_index[4]))
        print( "Found {} sites that are not in any annotation".format( self.summary_index[5]))
        return n_features

class Isolate:
    def __init__(self):
        self.id = ''
        self.seq = ''
        self.genes = {}
    def __str__(self):
        if self.id != '':
            return("Isolate: {} (Total bases: {})".format(self.id, len(self.seq)))
        else:
            return("Object is empty.")
    def __getitem__(self,key):
        return self.genes[key]
    def load_fasta(self, path, isolate_id):
        if not os.path.isfile(path):
            raise IOError("Could not open file {}".format(path))
        try:
            self.seq = SeqIO.read(path, "fasta")
            self.seq.id = isolate_id
            self.id = self.seq.id
        except IOError:
            print "Could not parse file {}".format(path)
            raise
        except ValueError:
            # introduced this to mainly deal with the case when there a FASTA
            # file has a Chromosome and one or more plamids
            # possibly not the most general approach but it solves the issue
            # at hand.
            print("#### WARNING ####")
            print("Found more than one sequence in the FASTA file for {}".format(isolate_id))
            print("Picking the first one, and ignoring the rest.")
            for seq in SeqIO.parse(path, "fasta"):
                self.seq = seq
                self.seq.id = isolate_id
                self.id = self.seq.id
                break
    def load_seqRec(self, record, isolate_id):
        '''
        Load a sequence record. To be used when loading the reference, which
        is first parsed using the Gene class, or from an alignment object.
        '''
        self.id = isolate_id
        self.seq = record
        self.seq.id = isolate_id
    def get_genes(self, genes_obj):
        '''
        Generates a string with the gene_id using the Genes object
        '''
        #generate sequence as a string, and remove the stop codon
        for f in genes_obj.features:
            seq = str(genes_obj.features[f].extract(self.seq.seq))[:-3]
            seq = re.sub("N", "-", seq)
            self.genes[f] = seq
            #print(self.id, f, self.genes[f], len(self.genes[f]))
        return

class Collection:
    def __init__(self):
        self.isolates = collections.OrderedDict({})
    def __getitem__(self, key):
        return self.isolates[key]
    def __str__(self):
        n_isolats = len( self.isolates.keys() )
        if n_isolates  == 0:
            print( 'An empty collection of isolates' )
        else:
            print( 'A collection of {} isolates.'.format( n_isolates ))
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
                    print("Skipping {}".format(d))
                    continue
                print("Trying to load {} at {}".format(d, os.path.join(path, d, default_file)))
                tmp_fn = os.path.join(path, d, default_file)
                self.isolates[d] = Isolate()
                self.isolates[d].load_fasta(tmp_fn, d)
                print("Successfully loaded: ")
                print(self.isolates[d])
            break
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
            aln_open = open( aln_file, 'r' )
        except IOError:
            print( "Could not open file {}.".format( aln_file ))
            raise
        except:
            print( "Something while trying to open file {}.".format( aln_file ))
            raise
        # try to parse alignment file
        try:
            aln = AlignIO.read( aln_open, aln_format )
        except ValueError:
            print( "Opened the file {}, but could not load the alignment. Is the format {} correct?".format( aln_file, aln_format ))
            raise
        except:
            print( "Something happend while trying to parse {}.".format( aln_file ))
            raise
        aln_list = list( aln )
        for s in aln_list:
            print( "Trying to load isolate {} into collection.".format( s.id ))
            self.isolates[s.id] = Isolate()
            self.isolates[s.id].load_seqRec( s, s.id )
        print( "Successfully loaded {} sequences from file {}.".format( len( self.isolates.keys() ), aln_file ) )
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
    def make_nexus( self, outfile, gene_obj ):
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
        #pdb.set_trace()
        partitions = [1, 2, 3, 4, 5]
        pos_index = gene_obj.indexed_positions
        # create a pandas DataFrame of the alignment
        self.aln_pd = pd.DataFrame( np.array([list(self.isolates[rec].seq) for rec in self.isolates], np.character, order="F"), \
                                    index = [self.isolates[rec].id for rec in self.isolates])
        # sort the columns so codon position categories are contiguous
        pos_index = pos_index.sort_values( ['codon_pos', 'positions'])
        new_column_order = pos_index.index
        self.aln_pd = self.aln_pd.iloc[:, new_column_order ]
        # figure out the ranges for each partition
        #pdb.set_trace()
        pos_index = pos_index.reset_index() # this generates a new index for the sorted DataFrame
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
        out += "beastify.py version {}\n".format( VERSION )
        out += "Date: {}\n".format( time.strftime("%d/%m/%Y") )
        out += "Reference: {}\n".format( gene_obj.reference.id )
        out += "]\n\n"
        out+= "begin taxa;\n"
        out += sp + "dimensions ntax={};\n".format(len(self.isolates.keys()))
        out += sp + "taxlabels\n"
        try:
            for i in self.isolates.keys():
                out += i + ":" + self.info[i] + "\n"
        except:
            for i in self.isolates.keys():
                out += i + "\n"
        out += sp + ";\n"
        out += "end;\n\n"
        out += "begin characters;\n"
        out += sp + "dimensions nchar={};\n".format( self.aln_pd.shape[ 1 ])
        out += sp + "format missing=? gap=- datatype=dna;\n"
        out += sp + "gapmode=missing;\n"
        out += sp + "matrix\n"
        #same as above when outputting the tax labels
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
        print("Finished parsing individual genes!")
    def gen_align(self, outfile, gene_obj):
        '''
        iterates over the keys in the self.__dict__ to return
        individual gen alignments, and format the output so it is
        suitable for BEAST.
        '''
        #import pdb; pdb.set_trace()
        concat_seqs = {}
        gene_lens = {}
        for i in self.isolates.keys():
            genome = self.isolates[i]
            genome.get_genes(gene_obj)
            concat_seqs[i] = ""
            for g in sorted(genome.genes.keys()):
                concat_seqs[i] += genome.genes[g]
                try:
                    gene_lens[g]
                    if (gen_lens[g] != len(genome.genes[g])):
                        print("Gene {} has different length in isolate {} than in previously recorded: {}.".format(g, len(genome.genes[g], gen_lens[g])))
                except:
                    gene_lens[g] = len(genome.genes[g])
                    pass
        seq_len = len(concat_seqs[concat_seqs.keys()[0]])
        var_sites = 0
        for n in xrange(seq_len):
            tmp = []
            for i in concat_seqs:
                tmp.append(concat_seqs[i][n])
            tmp = [t for t in tmp if t not in ['-']]
            tmp = set(tmp)
            if len(tmp) > 1:
                var_sites += 1
                # print("Site {} is variable.".format(n))
                # print(tmp)
        print("Total variable sites found: {}.".format(var_sites))
        # outputting nexus file
        sp = "    "
        out = "#NEXUS\n"
        out += "[Data from:\n"
        out += "beastify.py\n"
        out += "]\n\n"
        out+= "begin taxa;\n"
        out += sp + "dimensions ntax={};\n".format(len(self.isolates.keys()))
        out += sp + "taxlabels\n"
        # if self.dates exists, then add it to the id name
        # otherwise, ignore
        # for i in self.dates:
        #     print(i, self.dates[i])
        try:
            for i in self.isolates.keys():
                out += i + ":" + self.info[i] + "\n"
        except:
            for i in self.isolates.keys():
                out += i + "\n"
        out += sp + ";\n"
        out += "end;\n\n"
        out += "begin characters;\n"
        out += sp + "dimensions nchar={};\n".format(len(concat_seqs[concat_seqs.keys()[0]]))
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
        print("Finished parsing individual genes!")
        # creating an nexus file
        return

class TestBeastify(unittest.TestCase):
    '''
    A class for unit testing...
    To run:
    $> python -m unittest beastify.TestBeastify
    '''
    @classmethod
    def setUpClass(cls):
        '''
        Create class Genes to test
        '''
        # elements required for testing the class Gene
        reference = "test_data/test.gbk"
        cls.genes = Genes()
        cls.genes.load_genome( path = reference )
        cls.genes.index_locations()
        # elements required for testing the class Collection
        cls.alignment_multifasta = "test_data/test_multi.fasta"
        cls.isolate_collection = Collection()
        cls.isolate_collection.load_alignment( cls.alignment_multifasta, "fasta" )
        cls.isolate_collection.make_nexus( 'test_data/test.nexus', cls.genes )
        # load test alignment_multifasta
        aln_ordered_open = open( "test_data/test_reorder_multi.fasta", 'r')
        ordered_aln = AlignIO.read( aln_ordered_open, 'fasta' )
        aln_ordered_open.close()
        cls.aln_ordered_df = pd.DataFrame( np.array([list(rec) for rec in ordered_aln], np.character, order="F"), \
                                           index = [rec.id for rec in ordered_aln] )
    def test_1genes_load_genome(self):
        '''
        Tests if Genbank file is loaded correctly
        '''
        self.assertEqual( self.genes.reference.id , 'NC_007795.1' )
    def test_2genes_genome_size(self):
        '''
        Test if Genbank file loaded with appropriate size
        '''
        self.assertEqual( len( self.genes.reference ), 100801 )
    def test_3genes_index_locations(self):
        '''
        Test if genome indexing works
        Should result in the following count for each partition in the test dataset (results based on R implementation)
        First codon position: 29,513
        Second codon position: 29,514
        Third codon position: 29,513
        Fourth codon position: 109
        Unannotated: 12,152
        '''
        self.assertEqual( list( self.genes.summary_index ), [29513, 29514, 29513, 109, 12152] )
    def test_4collection_load_multifasta_alignment(self):
        '''
        Test if multifasta alignment works. Should produce 5 sequences of length
        equal to the reference sequence (100801bp).
        '''
        self.assertEqual( len( self.isolate_collection.isolates ), 5 )
    def test_5collection_make_pandas_align(self):
        '''
        Test if Pandas Alignment DataFrame is successfully made.
        '''
        self.assertEqual( self.isolate_collection.aln_pd.shape, (5, 100801 ) )
    def test_6collection_check_alignment(self):
        '''
        Check if reordered matrix matches the expectation produced with R
        '''
        self.isolate_collection.aln_pd.columns = list( self.aln_ordered_df.columns )
        self.assertTrue( pd_assert_df_equal( self.isolate_collection.aln_pd, self.aln_ordered_df ) == None  )
def run_tests(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    suite = unittest.TestLoader().loadTestsFromTestCase(TestBeastify)
    unittest.TextTestRunner(verbosity=2, buffer = False).run(suite)
    ctx.exit()


@click.command()
@click.option("--nb", \
                help = "Select N random genes from list", \
                default = None,
                type = int)
@click.option("--out", \
                help="Outfile name (default: out.nexus)", \
                default = 'out.nexus')
@click.option("--snippy", \
                help="Use snippy-core, and filter with 'random' or 'top' or 'all' (default: all). If specified 'random' or 'top', then --n must be specified too.", \
                default = "all")
@click.option("--feature", \
                help="Feature name to search in Genbank file (default: CDS)", \
                default = "CDS")
@click.option("--info",
                help = "Path to a tab-delimited file with two or more columns. The first column has the isolate ID, and other columns have dates, location, etc. The information will be added to the isolate ID in the same order as the columns",\
                default = None)
@click.option("--core_filename", \
                help = "If using --snippy, name of filenam containing SNP data (default: core.tab)",\
                default = "core.tab")
@click.option("--seq_filename", \
                help = "The name of the snippy alignment file to search for (default: snps.consensus.subs.fa)", \
                default = "snps.consensus.subs.fa")
@click.option("--gene_list", \
                help = "A list of genes to include in the nexus file.", \
                default = None)
@click.option("--exclude", \
                help = "Comma separated list of isolates to exclude (default: None). Example: Iso22,Iso34", \
                default = "")
@click.option("--inc_ref", \
                help = "Whether to include the reference in the final out file (default: False)", \
                is_flag = True, \
                default = False)
@click.option("--aln_file", \
                help = "A sequence alignment file to give in lieu of folder with snippy output.")
@click.option("--aln_file_format", \
                help = "If providing an alignment file with --aln_file, set the format of the alignment. Any format supported by BioPython:AlignIO could be valid. Default: fasta. Tested: fasta.", \
                default = 'fasta')
@click.option("--test", is_flag=True, default=False, callback=run_tests, expose_value=False, is_eager = True)
@click.argument("reference")
@click.argument("path")
def beastify(reference, path, gene_list, nb, out, snippy, feature, info, core_filename, seq_filename, exclude, inc_ref, test, aln_file, aln_file_format):
    '''
    REFERENCE: a path to reference Genbank file\n
    PATH: a path to a collection of snippy alignment files\n

    By Anders Goncalves da Silva
    '''
    #transforming exclude list to a Python list
    if exclude != None:
        exclude = exclude.split(",")

    #import pdb; pdb.set_trace()
    # start by parsing genbank file, and loading features
    genes = Genes()
    genes.load_genome( path = reference )
    if snippy != None and nb != None:
        genes.parse_snippycore(path = path, \
            corefn = core_filename, \
            nb = nb, \
            sample = snippy)
        genes.load_features( \
                            gene_list = None, \
                            nb = nb, \
                            feature = feature)
    elif snippy.lower() in ['random', 'top'] and nb == None:
        raise ValueError("If specifying --snippy, then --nb must be specified too.")
    elif snippy.lower() == 'all':
        genes.parse_snippycore(path = path, \
                    corefn = core_filename, \
                    nb = nb, \
                    sample = 'all')
        genes.load_features( \
                            gene_list = None, \
                            nb = None, \
                            feature = feature)
    else:
        genes.load_features( \
                            gene_list = gene_list, \
                            nb = nb, \
                            feature = feature)
    collection = Collection()
    if inc_ref:
        collection.load_reference(reference = genes.reference)
    collection.load_isolates(path = path, \
                            seq_file = seq_filename, \
                            ignore = exclude)
    if info != None:
        print('#'*80)
        collection.add_info(info)
    collection.gen_align(out, genes)
    return

if __name__=="__main__":
    beastify()
