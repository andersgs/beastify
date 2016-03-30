'''
beastify.py

Anders Goncalves da Silva
October 2015

Updated on 7 January 2016
'''

import click
from Bio import SeqIO
import os
import random
import re
import sys

class Genes:
    def __init__(self):
        self.features = {}
        self.snippy_list = []
    def load_features(self, path, gene_list = None, feature = 'gene', nb = None, cod_tab = 11):
        # load features
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
        # load features to keep
        #import pdb; pdb.set_trace()
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
class Isolate:
    def __init__(self):
        self.id = ''
        self.seq = ''
        self.genes = {}
#        self.concat = ''
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
        is first parsed using the Gene class
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
        self.isolates = {}
    def __getitem__(self, key):
        return self.isolates[key]
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
                help="Feature name to search in Genbank file (default: gene)", \
                default = "gene")
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
@click.argument("reference")
@click.argument("path")
def beastify(reference, path, gene_list, nb, out, snippy, feature, info, core_filename, seq_filename, exclude, inc_ref):
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
    if snippy != None and nb != None:
        genes.parse_snippycore(path = path, \
            corefn = core_filename, \
            nb = nb, \
            sample = snippy)
        genes.load_features(path = reference, \
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
        genes.load_features(path = reference, \
                            gene_list = None, \
                            nb = None, \
                            feature = feature)
    else:
        genes.load_features(path = reference, \
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
