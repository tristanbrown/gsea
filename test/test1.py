import unittest
import numpy as np

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
        data, phenos, genes = gep_data
        assert len(data) == len(genes)
        assert len(data[0]) == len(phenos)
        assert data[0][0] == -13
        assert phenos[0] == 'ALL'
        assert genes[0] == 'STAT1'
        assert genes[-1] == 'TUBA1'
    
    def test_geneset_extract(self):
        geneset_file = dataprep.IO('pathways.txt', 'data')
        labels, data = geneset_file.load_arb_rows(delim='\t', 
                                                    type='str', skipcol=1)
        
        assert len(data) == 522
        
        # print(self.calc_ES(ranked, corr,
                # self.geneset['41bbPathway']))
        # print(self.calc_ES(ranked, corr,
                # self.geneset['CBF_LEUKEMIA_DOWNING_AML']))
        # print(self.calc_ES(ranked, corr,
            # self.geneset['MAP00532_Chondroitin_Heparan_sulfate_biosynthesis']))
        
    
    def test_full_analysis(self):
        A = analysis.Analysis(
            'leukemia.txt', 'pathways.txt', 'data', config.analysis
                )
        sum_hits = np.cumsum(A.hits, axis=-1)
        assert sum_hits[0][-1] == 12
        assert sum_hits[-1][-1] == 37
        
        A.analyzefiles('test_out.csv', 'data')
