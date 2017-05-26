import unittest

from gsea import config
from gsea import dataprep
from gsea import batchanalysis

class FileStructureTestCase(unittest.TestCase):
    def setUp(self):
        pass
    
    def test_file_io(self):
        A = batchanalysis.Analyzer('data/leukemia.txt', 'data/pathways.txt',
                                    'data/test_out.txt', config.analysis)
        A.writecsv(A.outfile, A.gepdata)
