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
        
        self.params = parameters
    
    def rank(self):
        self.ranked = Ranker(self.gepdata, self.params['rankby'])
        # self.permuted = Permuter(self.gepdata, self.params['permutations'])
    
    def analyzesets(self):
        
        self.rank()
        self.nes_table = self.loop_genesets(self.genesets)
    
    def loop_genesets(self, genesets):
        return np.array([self.analyzeset(a_set) for a_set in genesets])
        
        
    def analyzeset(self, geneset):
        # es = ES_calc(self.ranked, self.params['p_weight'])
        # es_p = ES_P_calc(self.permuted, self.params['p_weight'])
        # p_stat = P_calc(es, es_p)
        # nes = NES_calc(es, es_p)
        # return (nes, p_stat)
        return None
        
    def loadgep(self, filename):
        return np.genfromtxt(filename, delimiter='\t', dtype=None,
                                skip_header=0)
    
    def loadgeneset(self, filename):
        return []
    
    def writecsv(self, filename, data):
        with open(filename, 'w', newline='\n') as f:
            writer = csv.writer(f)
            writer.writerows(data)
    
    def writeresults(self):
        self.writecsv(self.outfile, self.nes_table)