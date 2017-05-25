import unittest
from gsea import batchanalysis

batchanalysis.Analyzer('leukemia.txt', 'pathways.txt', 'test_out.txt')