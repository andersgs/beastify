import unittest
from ..Genes import Genes
from ..Collection import Collection
from ..Isolate import Isolate
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
