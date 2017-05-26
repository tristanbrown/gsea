"""
This module loops the Gene Set Enrichment Analysis over a file containing 
multiple gene sets, producing an output file tabulating the results. 
"""

import numpy as np
from gsea.ranker import Ranker

class Analyzer():
    """
    """
    def __init__(self, expression_data, geneset_data, parameters):
        
        self.gepdata = expression_data
        self.genesets = geneset_data
        self.params = parameters
    
    def rank(self):
        self.ranked = Ranker(self.gepdata, self.params['rankby'])
        # self.permuted = Permuter(self.gepdata, self.params['permutations'])
    
    def analyzesets(self):
        
        self.rank()
        return self.loop_genesets(self.genesets)
    
    def loop_genesets(self, genesets):
        return np.array([self.analyzeset(a_set) for a_set in genesets])
        
        
    def analyzeset(self, geneset):
        # es = ES_calc(self.ranked, self.params['p_weight'])
        # es_p = ES_P_calc(self.permuted, self.params['p_weight'])
        # p_stat = P_calc(es, es_p)
        # nes = NES_calc(es, es_p)
        # return (nes, p_stat)
        return None
        