import unittest

import gsea

class FileStructureTestCase(unittest.TestCase):
    def setUp(self):
        pass
    
    def test_file_io(self):
        A = gsea.batchanalysis.Analyzer('data/leukemia.txt', 'data/pathways.txt',
                                    'data/test_out.txt', gsea.config.analysis)
        A.writecsv(A.outfile, A.gepdata)
