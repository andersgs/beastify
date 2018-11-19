'''
A class to handle individual isolate information
'''

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
            print("Could not parse file {}".format(path))
            raise
        except ValueError:
            # introduced this to mainly deal with the case when there a FASTA
            # file has a Chromosome and one or more plamids
            # possibly not the most general approach but it solves the issue
            # at hand.
            print("#### WARNING ####")
            print(("Found more than one sequence in the FASTA file for {}".format(isolate_id)))
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
