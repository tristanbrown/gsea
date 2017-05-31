import unittest

from gsea import config
from gsea import dataprep
from gsea import batchanalysis

class FileStructureTestCase(unittest.TestCase):
    def setUp(self):
        pass
    
    def test_file_io(self):
        gep_file = dataprep.IO('leukemia.txt', 'data') 
        gep_data = gep_file.load_array_with_labels(delim='\t')
        print("Gene Expression Data:")
        print(gep_data)
