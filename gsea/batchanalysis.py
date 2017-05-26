"""
This module loops the Gene Set Enrichment Analysis over a file containing 
multiple gene sets, producing an output file tabulating the results. 
"""

import numpy as np
import csv
from gsea.ranker import Ranker

class Analyzer():
    def __init__(self, expression_file, geneset_file, outfile, parameters):
        
        self.outfile = outfile
        
        self.gepdata = self.loadgep(expression_file)
        self.genesets = self.loadgeneset(geneset_file)
        
        self.ranked = Ranker(self.gepdata, parameters['rankby'])
        
    def loadgep(self, filename):
        return np.genfromtxt(filename, delimiter='\t', dtype=None,
                                skip_header=0)
    
    def loadgeneset(self, filename):
        pass
    
    def writecsv(self, filename, data):
        with open(filename, 'w', newline='\n') as f:
            writer = csv.writer(f)
            writer.writerows(data)