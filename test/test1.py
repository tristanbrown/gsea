import unittest

class FileStructureTestCase(unittest.TestCase):
    def setUp(self):
        pass
    
    def test_filepath_handling(self):
        from gsea import batchanalysis
        batchanalysis.Analyzer('leukemia.txt', 'pathways.txt', 'test_out.txt')
