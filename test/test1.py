import unittest

class FileStructureTestCase(unittest.TestCase):
    def setUp(self):
        pass
    
    def test_file_io(self):
        from gsea import batchanalysis
        A = batchanalysis.Analyzer('data/leukemia.txt', 'data/pathways.txt',
                                    'data/test_out.txt')
        A.writecsv(A.outfile, A.gepdata)
