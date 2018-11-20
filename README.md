# beastify: Generate input file for BEAST from whole-genome alignmennt

## Background

Sometimes you want to partiion your alignment in to distinct codon positions
(i.e., 1, 2, and 3), and you want to also model non-coding positions in your
BEAST analysis.

`beastify` does that for you. It will:

1. Figure out all the codon positions in your reference (including overlapping positions)
2. Optionally, label your sequences with any metadata (e.g., dates)
3. Optionally, allows you to remove one or more positions from the alignment
4. Optionally, allows you to mask positions form the alignment
5. Optionally, allows you to sub-sample the alignment (if you want to work on a smaller dataset to test your models before throwing the whole kitchen at BEAST).
6. Output a NEXUS files with the partitions ready for running BEAUTi.

Partitions are labelled:

1. For the first codon position
2. For the second codon position
3. For the third codon position
4. For any overlapping codons (sometimes CDS annotations overlap because sometimes bacterial genes will share codons)
5. If position is not found in a CDS.

## Installation

### Dependencies

- Python >=3.6
- Click
- Pandas
- Numpy
- BioPython

### Using pip

```
pip3 install beastify
```

### Testing your installation

```
beastify --test
```

## Input

1. Genbank reference
2. `snippy` \*.consensus.subs.fa files
3. List of genes to include in the final alignment
4. N (optional) --- random number of genes to select and include

### Command list

```
  --out TEXT                Outfile name (default: out.nexus)
  --info TEXT               Path to a tab-delimited file with two or more
                            columns. The first column has the isolate ID, and
                            other columns have dates, location, etc. The
                            information will be added to the isolate ID in the
                            same order as the columns
  --inc_ref                 Whether to include the reference in the final out
                            file (default: False)
  --aln_file TEXT           A sequence alignment file to give in lieu of
                            folder with snippy output.
  --aln_file_format TEXT    If providing an alignment file with --aln_file,
                            set the format of the alignment. Any format
                            supported by BioPython:AlignIO could be valid.
                            Default: fasta. Tested: fasta.
  --subsample INTEGER       Subsample X number of bases at random from each
                            partition. default: all bases
  --subsample_seed INTEGER  Set the seed when subsampling sites. Default:42
  --parts TEXT              Comma-separated list of partitions to include.
                            default:1,2,3,4,5
  --test                    Run beastify tests and exit
  --mask TEXT               A BED file indicating regions to mask from the
                            genome
  --version                 Show the version and exit.
  --help                    Show this message and exit.
```

## Output

A `nexus` formatted file ready for `beast`.

## Script outline

1. Parse coordinates of genes from Genbank into a `Genes` Class
   - Methods:
     - load_features: a method to load the Genbank features into a
       dictionary. **Method should check that there
       the length is a multiple of 3**, and that the
       **start** and **end** codons are in place. **stop**
       codon should be stripped.
     - parse_snippycore: a method to load the snippy core.tab data
       and identify all variable SNPs among the data that are in
       coding regions --- has options to return a 'random' sample
       of size N genes, 'top' genes with the most SNPs, with the
       N top genes with most SNPs.
   - Data:
     - features: a dictionary with key = genename and value
       set by seqFeature object --- **IF** N is provided, only
       keep a random set of gene*coords of size \_N*
2. Load `snippy` alignment into an `Isolate` Class
   - Methods:
     - load_fasta: will load the sequence into the object
     - cat_genes: given an isolate id, and a genes object,
       return a concatenated sequence (NOT IMPLEMENTED YET)
     - get_gene: return a string with the sequence for the gene specified
       by gene_id using a `Genes` object
     - **str**: print the sequence ID and length, if there is one
     - **getitem**: return the sequence string associated with the key
     - add_dates: the user supplies a table of isolate IDs, and dates in
       a format suitable for BEAST, and the script adds them to the
       identifier
   - Data:
     - seq: A SeqRecord
     - id: The isolate id
     - genes: a dictionary with 'gene_name' as keys and sequence string as
       value
3. `Collection` class to store all the `Isolate` objects
   - Methods:
     - load_isolates: given a list of isolate files, creates
       and stores individual `Isolate` objects for each.
     - gen_align: given a `Genes` object, generate the
       alignment --- uses `cat_genes`
     - **getitem**: given an isolate ID as a key, return the `Isolate`
       object
   - Data:
     - isolates: a dictionary with isolate id as keys and
       `Isolate` objects as values
