#!/usr/bin/env Rscript
################################################################################
### snp_sim.R
### A script to simulate SNPs for testing software
###
### Anders Gonçalves da Silva and Torsten Seemann
### 18 April 2016
### version 0.1
################################################################################

################################################################################
### required libraries

suppressPackageStartupMessages( require( ape ) )
suppressPackageStartupMessages( require( phyclust ) )
suppressPackageStartupMessages( require( dplyr ) )
suppressPackageStartupMessages( require( magrittr ) )
suppressPackageStartupMessages( library( magrittr ) )
suppressPackageStartupMessages( require( optparse ) )
suppressPackageStartupMessages( require( stringr ) )

################################################################################
### version number

version_number = 0.1
args = commandArgs(trailingOnly = TRUE)
if( grepl( pattern = '-v', x = args) || grepl( pattern = '--version', x = args)) {
  write( "snp_sim v0.1", file = stderr() )
  q( save = 'no', status = 0, runLast = FALSE )
}

################################################################################
### set arguments

option_list <- list(
  optparse::make_option( c("-o", "--outprefix"), 
                         action = 'store', type = "character", default = 'snp_sim_outfile', 
                         help = 'Prefix for outfiles. [default %default]'),
  optparse::make_option( c( '-n', '--n_samples' ), 
                         action = 'store', type = 'integer', default = 5, 
                         help = 'Number of isolates to simulate. [default %default]' ),
  optparse::make_option( c( '-s', '--seed' ), 
                         type = 'integer', default = 42, 
                         help = 'Set the random seed. [default %default]' ),
  optparse::make_option( c("-v", '--version'), action = 'store_true', default = FALSE, help = 'Display version and exit' )
)
 
parser <- OptionParser(usage = "%prog [options] FASTA", option_list=option_list, 
                       epilogue = "By Anders Goncalves da Silva and Torsten Seemann on 18 April 2016", description = "An R script to produce simulated sequences with SNPs evolved along a tree given an ancestral sequence. Output will be tree file, a multifasta with the mutated sequences, and a VCF file with the introduced SNPs. ")

#arguments <- parse_args(parser, positional_arguments = 1, args = c("test_data/test.fasta"))
arguments <- parse_args(parser, positional_arguments = 1)
opt <- arguments$options
FASTA <- arguments$args

################################################################################
### try to read ancestral sequence file
if( file.access(FASTA) == -1) {
  stop(sprintf("Could not find specified FASTA file ( %s ).", FASTA))
} else {
 ancestral_seq <- ape::read.dna( FASTA, format = 'fasta' )
}
write( sprintf( "Successfully read %s in to memory!\n", FASTA ), stderr() )
#create chrom label for VCF file
chrom = stringr::str_split( labels( ancestral_seq ), pattern = ' +' )[[1]][1]

################################################################################
### load a couple aux functions

#function will identify alternate alleles to the reference
# used in producing the VCF file
find_alternate_allele <- function( site ) {
  ref = site$Ref
  alleles = sort( toupper( unique( unlist( site[ which( site != ref ) ] ) ) ) )
  if ( length( alleles ) > 1) {
    return( stringr::str_c( alleles, collapse = ',' ) )
  } else {
    return( alleles )
  }
}

#transform the genotype base calls, to a coded number of alleles
# used in making VCF file. 
# reference allele gets a 0, alternate alleles get coded values 1 
# or larger in the same order they appear in the ALT column
genotype_calls = function( site ) {
  ref <-  site$Ref
  alleles <-  sort( unique( unlist( site[ which( site != ref ) ] ) ) )
  alleles <-  c(ref, alleles )
  code <-  1:length(alleles) - 1
  names( code ) <-  alleles
  seq_id <-  names( site )
  site <- as.character( site )
  coded_sites <- code[ site ]
  names( coded_sites ) <- seq_id
  coded_sites <- as.data.frame( t( coded_sites ), stringsAsFactors = F )
  coded_sites <- coded_sites %>% dplyr::select( -Ref )
  return( coded_sites )
}


################################################################################
### set file names
out_multifasta_file <- paste( opt$outprefix, '.fasta', sep = '' )
out_tree_file <- paste( opt$outprefix, '.nwk', sep = '' )
out_vcf_file <- paste( opt$outprefix, '.vcf', sep = '' )

################################################################################
### simulate tree and sequences

### simulation parameters
seed_number = opt$seed
n_seq = opt$n_samples
theta = 0.1
kappa = 2
base_freq = c( 0.25, 0.25, 0.25, 0.25 )
n_ancestral_seq = 1

### simulation steps
set.seed( seed_number )
write( sprintf( "Setting the seed to %d!\n", opt$seed ), stderr() )
anc_seq = code2nid( toupper( as.character( ancestral_seq )))
seq_len = length( anc_seq )
write( "Simulating a tree!\n", stderr() )
tree <- ape::read.tree( text = phyclust::ms(nsam = n_seq, nreps = 1, opts = '-T' )[3] )
ape::write.tree(phy =  tree, file = out_tree_file )
write( "Simulating sequences!\n", stderr() )
sim_seqs <- gen.seq.HKY(rooted.tree = tree, 
                        pi = base_freq,
                        rate.scale = theta, 
                        anc.seq = anc_seq, 
                        kappa = kappa, L = seq_len )
write( "Writing simulated sequences!\n", stderr() )
cat(stringr::str_c( sapply( 2:length(sim_seqs), function( i ) {
  tmp <- stringr::str_split(sim_seqs[i], pattern = "\ +")
  stringr::str_c( stringr::str_c( '>', tmp[[1]][1] ), tmp[[1]][2], sep = '\n' )
}), collapse = '\n' ), file = out_multifasta_file )

## add ancestral/reference sequence to file

ancestral_string <- stringr::str_c( stringr::str_c( '\n>', 'Ref' ), 
                                    stringr::str_c( toupper( as.character( ancestral_seq ) ), 
                                                    collapse = '' ), sep = '\n' )

cat( ancestral_string, file = out_multifasta_file, append = T )
################################################################################
### generate the VCF file

write( "Writing VCF file!\n", stderr() )
# upload the simulated sequences
sim_seqs = ape::read.dna( out_multifasta_file, format = 'fasta' )
# find variable sites
var_sites = ape::seg.sites( sim_seqs )

isolate_ids = stringr::str_c( rownames( sim_seqs )[! rownames( sim_seqs ) == 'Ref' ], collapse = '\t' )

header1 = '##fileFormat=VCFv4.0'
header2 = stringr::str_c( '##fileDate=', format(Sys.Date(), '%Y%m%d'), sep = '' )
header3 = stringr::str_c( '##source=snp_sim', version_number, sep = '' )
header4 = stringr::str_c( '##reference=', labels(ancestral_seq))
header5 = '##phasing=complete'
header6 = '##INFO=<ID=SYNTHETIC,Description=\"Alelle simulated with snp_sim.\">'
header7 = '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
header8 = stringr::str_c( '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT', isolate_ids, sep = '\t' )

genotype_df <- as.data.frame( t( as.character( sim_seqs[, var_sites ] ) ), stringsAsFactors = F )

genotype_df <- genotype_df %>%
  purrr::by_row(..f = find_alternate_allele, .collate = 'cols', .to = 'ALT' )

genotypes_code <-genotype_df %>% 
      dplyr::select( -ALT ) %>%
      purrr::by_row(..f = genotype_calls ) %>% 
      dplyr::select( `.out` )
genotypes_code <- dplyr::bind_rows( genotypes_code$.out )
vcf_df <- data.frame( chrom = chrom, 
                      pos = as.numeric( var_sites ), 
                      ID = '.', 
                      ref = toupper( genotype_df$Ref ), 
                      alt = genotype_df$ALT, 
                      QUAL = 99, 
                      filter = 'PASS', 
                      info = 'SYNTHETIC', 
                      format = 'GT', 
                      genotypes_code)

header = stringr::str_c( header1, 
                         header2, 
                         header3, 
                         header4, 
                         header5,
                         header6,
                         header7,
                         header8, 
                         sep = '\n')

genotype_string = stringr::str_c( apply( vcf_df, 1, function( x ) {
  stringr::str_c( x, collapse = '\t' )
}), collapse = '\n' )

cat( header, genotype_string, sep = '\n', file = out_vcf_file )

write( "Work here is done!\n", stderr() )