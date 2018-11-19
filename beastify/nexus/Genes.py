'''
This is a class for handling genes
'''

from Bio import SeqIO
import pandas as pd
import numpy as np

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
            print(("Could not open file {}".format(path)))
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
                print(("Could not open file {}".format(feature_list)))
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
                            print(("Incomplete CDS {}".format(feat.qualifiers[feature][0])))
                            print(("Feature's length divided by 3 had remainder = {}".format(len(feat) % 3)))
                            tmp = feat.extract(genome.seq)
                            print(("Features first codon was {}".format(tmp[0:3])))
                            print(("Features last codon was {}".format(tmp[-3:])))
                            print((str(tmp)))
                            pass
                except:
                    pass
        print(n_cds, n_gene, n_other)
        print(("Found {} features".format(features_found)))
        print(("of which {} were ok".format(features_ok)))
        for k in self.features:
            print(("Found gene {}".format(k)))
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
            print(("{}, {}, {}".format(len(cds_var), nb, len(cds_var) > nb)))
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
        print(("Found {} variable genes, and picked {}.".format(len(cds_var), len(self.snippy_list))))
    
    def _feature_parsing(self, feature):
        '''
        A function to return a pandas.DataFrame of a feature
        '''
        codon_pos = [1,2,3]
        total_codons = int((feature.location.end - feature.location.start)/3)
        data = {'positions': list(range(feature.location.start + 1, feature.location.end + 1)),
         'codon_pos': (codon_pos if feature.strand == 1 else codon_pos[::-1]) * total_codons,
         'locus_tag': feature.qualifiers['locus_tag'][0]}
        return pd.DataFrame(data)

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
        gff_list = [ self._feature_parsing(feature) for feature in self.reference.features if feature.type in feature_type ]
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
        print(( "Found {} features that match CDS".format( n_features ) ))
        self.summary_index = self.indexed_positions.groupby('codon_pos').size()
        print(( "Found {} sites in the first codon position".format( self.summary_index[1])))
        print(( "Found {} sites in the second codon position".format( self.summary_index[2])))
        print(( "Found {} sites in the third codon position".format( self.summary_index[3])))
        print(( "Found {} sites in that are in multiple CDS annotations".format( self.summary_index[4])))
        print(( "Found {} sites that are not in any annotation".format( self.summary_index[5])))
        return n_features
