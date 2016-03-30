################################################################################
### generate some simulated data to test beastify

library( dplyr )
library( phyclust )
library( ape )

################################################################################
### load sequence data

test_fasta <- ape::read.dna( file = "test_data/test.fasta", format = 'fasta' )
test_gff <- readr::read_tsv( file = "test_data/test.gff", comment = "#", col_names = F )
names( test_gff ) <- c( "reference", "source","feature", "start", "end", "score", "strand", "phase", "description" )

# remove first row, which has information about the whole sequence
test_gff <- test_gff %>% 
  dplyr::filter( feature != 'source' )

test_gff %>% dplyr::group_by( feature ) %>% dplyr::summarise( count = n() )
# table of features

test_gff %>% 
  dplyr::group_by( feature ) %>%
  dplyr::mutate( feature_length = end - start + 1 ) %>%
  dplyr::select( description )
  dplyr::summarise( sum( feature_length ) )
# table of total length of features

# intergenic regions
test_gff_cds <- test_gff %>%
  dplyr::filter( feature == 'CDS' ) %>%
  dplyr::select( end, start )

test_gff_cds
intergenetic_regions <- data.frame( start = (test_gff_cds$end[-nrow( test_gff_cds )] + 1), 
                                    end = (test_gff_cds$start[-1] - 1) ) 

intergenetic_regions %>%
  dplyr::mutate( length_intergenic = end - start + 1) %>%
  dplyr::filter( length_intergenic > 0 ) %>%
  dplyr::summarise( sum( length_intergenic ) )

head( test_gff_cds, n = 8 )
################################################################################
### index all positions
n_bases <- length( test_fasta )

make_site_table <- function( positions, n_bases, features ) {
  # partitions are numeric:
  # 1 - 1st codon position
  # 2 - 2nd codon position
  # 3 - 3rd codon position
  # 4 - shared codon positions
  # 5 - other genetic positions
  filtered_sites <- positions %>%
    dplyr::filter( feature %in% features )
  references = unique( positions$reference )
  index_sites <- dplyr::bind_rows(
    mapply(FUN = function(p, n) { data.frame( key = 1:n, chr = p, pos = 1:n ) }, 
           references, n_bases, SIMPLIFY = F)
  )
  
  return( index_sites )
}

head( make_site_table( test_gff, n_bases = n_bases, feature = 'CDS' ) )

unpack_sites <- function( feature_row ) {
  if (feature_row$strand == '+' ) {
    start = feature_row$start + feature_row$phase
    end = feature_row$end
    total_length <- end - start + 1
    increment = 3
  } else {
    start = feature_row$end - feature_row$phase
    end = feature_row$start
    total_length <- start - end + 1
    increment = -3
  }
  description <- strsplit( feature_row$description, split = ';' )[[1]]
  locus_tag_ix <- pmatch( x = 'locus_tag', table = description )
  locus_tag <- gsub(pattern = " ", replacement = '', 
       gsub(pattern = "\"", replacement = '', 
            gsub( pattern = 'locus_tag', 
                  replacement = "", description[locus_tag_ix] ) ) )
  index_sites <- data.frame( pos = start:end, partition = rep( c(1, 2, 3), total_length/3 ), locus_tag = locus_tag, stringsAsFactors = F )
  return( index_sites )
}

tmp <- test_gff %>%
  dplyr::filter( feature == 'CDS' ) %>%
  purrr::by_row( unpack_sites, .collate = "rows" ) %>%
  dplyr::select( pos, partition, locus_tag ) %>%
  dplyr::arrange( pos )

tmp <- tmp %>%
  dplyr::group_by(pos) %>%
  dplyr::mutate( partition = ifelse( length( pos ) > 1, 4, partition ) )

tmp_dist <- tmp %>% 
  dplyr::select( pos, partition ) %>%
  dplyr::distinct()
nrow( tmp )
nrow( tmp_dist )

all_sites = data.frame( reference = "NCXXXX", pos = 1:n_bases )
all_sites <- all_sites %>%
  dplyr::left_join( tmp_dist, c('pos' = 'pos' ) ) %>%
  tidyr::replace_na( list( partition = 5 ) )