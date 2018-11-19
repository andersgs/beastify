#!/home/agoncalves/dev/beastify/venv/bin/python
'''
beastify.py

Anders Goncalves da Silva
October 2015

Updated on 7 January 2016
'''

import click
import os
import random
import re
import sys
import unittest
import warnings
from pandas.util.testing import assert_frame_equal as pd_assert_df_equal
import pdb

from beastify.nexus.Genes import Genes
from beastify.nexus.Collection import Collection
from beastify.nexus.tests.Test_Nexus import TestBeastify
from beastify import __VERSION__ as version_string


def run_tests(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    suite = unittest.TestLoader().loadTestsFromTestCase(TestBeastify)
    unittest.TextTestRunner(verbosity=2, buffer=False).run(suite)
    ctx.exit()


@click.command()
# @click.option("--path", \
#                 help = "Path to SNIPPY output folder.",
#                 default = None)
# @click.option("--nb", \
#                 help = "Select N random genes from list", \
#                 default = None,
#                 type = int)
@click.option("--out",
              help="Outfile name (default: out.nexus)",
              default='out.nexus')
# @click.option("--snippy", \
#                 help="Use snippy-core, and filter with 'random' or 'top' or 'all' (default: all). If specified 'random' or 'top', then --n must be specified too.", \
#                 default = "all")
# @click.option("--feature", \
#                 help="Feature name to search in Genbank file (default: CDS)", \
#                 default = "CDS")
@click.option("--info",
              help="Path to a tab-delimited file with two or more columns. The first column has the isolate ID, and other columns have dates, location, etc. The information will be added to the isolate ID in the same order as the columns",
              default=None)
# @click.option("--core_filename", \
#                 help = "If using --snippy, name of filenam containing SNP data (default: core.tab)",\
#                 default = "core.tab")
# @click.option("--seq_filename", \
#                 help = "The name of the snippy alignment file to search for (default: snps.consensus.subs.fa)", \
#                 default = "snps.consensus.subs.fa")
# @click.option("--gene_list", \
#                 help = "A list of genes to include in the nexus file.", \
#                 default = None)
# @click.option("--exclude", \
#                 help = "Comma separated list of isolates to exclude (default: None). Example: Iso22,Iso34", \
#                 default = "")
@click.option("--inc_ref",
              help="Whether to include the reference in the final out file (default: False)",
              is_flag=True,
              default=False)
@click.option("--aln_file",
              help="A sequence alignment file to give in lieu of folder with snippy output.")
@click.option("--aln_file_format",
              help="If providing an alignment file with --aln_file, set the format of the alignment. Any format supported by BioPython:AlignIO could be valid. Default: fasta. Tested: fasta.",
              default='fasta')
@click.option("--subsample",
              help="Subsample X number of bases at random from each partition. default: all bases",
              default=None, type=int)
@click.option("--subsample_seed",
              help="Set the seed when subsampling sites. Default:42",
              default=42, type=int)
@click.option("--parts",
              help="Comma-separated list of partitions to include. default:1,2,3,4,5",
              default='1,2,3,4,5')
@click.option("--test", is_flag=True, default=False, callback=run_tests, expose_value=False, is_eager=True,
              help="Run beastify tests and exit")
@click.option("--mask", help="A BED file indicating regions to mask from the genome", default=None)
@click.version_option(version=version_string, message="beastify v{}".format(version_string))
@click.argument("reference")
def beastify(reference, out, info, inc_ref, aln_file, aln_file_format, subsample, subsample_seed, parts, mask):
    '''
    REFERENCE: a path to reference Genbank file\n

    By Anders Goncalves da Silva
    '''
    # transforming exclude list to a Python list
    # if exclude != None:
    #     exclude = exclude.split(",")
    # start by parsing genbank file, and loading features
    #import pdb; pdb.set_trace()
    if parts == '1,2,3,4,5':
        parts = [1, 2, 3, 4, 5]
    else:
        parts = [int(part.strip()) for part in parts.split(',')]
    genes = Genes()
    genes.load_genome(path=reference)
    genes.index_locations()
    collection = Collection()
    if inc_ref:
        collection.load_reference(reference=genes.reference)

    if aln_file != None:
        collection.load_alignment(aln_file, aln_file_format)
    else:
        raise ValueError("beastify only accepts alignments for the moment.")
    # elif snippy != None and nb != None:
    #     genes.parse_snippycore(path = path, \
    #         corefn = core_filename, \
    #         nb = nb, \
    #         sample = snippy)
    #     genes.load_features( \
    #                         gene_list = None, \
    #                         nb = nb, \
    #                         feature = feature)
    # elif snippy.lower() in ['random', 'top'] and nb == None:
    #     raise ValueError("If specifying --snippy, then --nb must be specified too.")
    # elif snippy.lower() == 'all':
    #     genes.parse_snippycore(path = path, \
    #                 corefn = core_filename, \
    #                 nb = nb, \
    #                 sample = 'all')
    #     genes.load_features( \
    #                         gene_list = None, \
    #                         nb = None, \
    #                         feature = feature)
    # else:
    #     genes.load_features( \
    #                         gene_list = gene_list, \
    #                         nb = nb, \
    #                         feature = feature)
    # collection.load_isolates(path = path, \
    #                         seq_file = seq_filename, \
    #                         ignore = exclude)
    if info != None:
        print(('#'*80))
        collection.add_info(info)
    collection.make_nexus(out, genes, subsample=subsample,
                          subsample_seed=subsample_seed, partitions=parts, mask=mask)
    return


if __name__ == "__main__":
    beastify()
