# `beastify`: Generate input file for BEAST from whole-genome alignmennt

## Input

1. Genbank reference
2. `snippy` \*.consensus.subs.fa files
3. List of genes to include in the final alignment
4. N (optional) --- random number of genes to select and include

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
            keep a random set of gene_coords of size *N*
2. Load `snippy` alignment into an `Isolate` Class
    - Methods:
        - load_fasta: will load the sequence into the object
        - cat_genes: given an isolate id, and a genes object,
            return a concatenated sequence (NOT IMPLEMENTED YET)
        - get_gene: return a string with the sequence for the gene specified
            by gene_id using a `Genes` object
        - __str__: print the sequence ID and length, if there is one
        - __getitem__: return the sequence string associated with the key
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
        - __getitem__: given an isolate ID as a key, return the `Isolate`
            object
    - Data:
        - isolates: a dictionary with isolate id as keys and
            `Isolate` objects as values
