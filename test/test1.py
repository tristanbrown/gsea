import unittest

from gsea import config
from gsea import dataprep
from gsea import analysis

class FileStructureTestCase(unittest.TestCase):
    def setUp(self):
        pass
    
    def test_gep_data_extract(self):
        gep_file = dataprep.IO('leukemia.txt', 'data') 
        gep_data = gep_file.load_array_with_labels(delim='\t')
        # print("Gene Expression Data:")
        # print(gep_data)
        data, genes, phenos = gep_data
        assert len(data) == len(genes)
        assert len(data[0]) == len(phenos)
        assert data[0][0] == -13
        assert phenos[0] == 'ALL'
        assert genes[0] == 'STAT1'
        assert genes[-1] == 'TUBA1'
    
    def test_gene_set_looping(self):
        geneset_file = dataprep.IO('pathways.txt', 'data')
        f = lambda x,y: x
        geneset_data = geneset_file.row_analysis(f, [], delim='\t')
        # print(geneset_data)
        assert len(geneset_data) == 522
        assert len(geneset_data[0]) == 20
        assert len(geneset_data[-1]) == 77
    
    def test_geneset_extract(self):
        geneset_file = dataprep.IO('pathways.txt', 'data')
        data = geneset_file.load_arb_rows(delim='\t', type='str', skipcol=1)
        print(data)
        
    
    def test_full_analysis(self):
        A = analysis.Analysis(config.analysis)
        A.analyzefiles('leukemia.txt', 'pathways.txt', 'test_out.csv',
                        config.path)
