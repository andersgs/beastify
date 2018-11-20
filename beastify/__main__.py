'''
The main script for beastify
'''

import os
import random
import re
import sys
import unittest

import click

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
@click.option("--out",
              help="Outfile name (default: out.nexus)",
              default='out.nexus')
@click.option("--info",
              help="Path to a tab-delimited file with two or more columns. The first column has the isolate ID, and other columns have dates, location, etc. The information will be added to the isolate ID in the same order as the columns",
              default=None)
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
    if info != None:
        print(('#'*80))
        collection.add_info(info)
    collection.make_nexus(out, genes, subsample=subsample,
                          subsample_seed=subsample_seed, partitions=parts, mask=mask)
    return


if __name__ == "__main__":
    beastify()
