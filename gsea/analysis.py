"""
This module loops the Gene Set Enrichment Analysis over a file containing 
multiple gene sets, producing an output file tabulating the results. 
"""

import numpy as np
from gsea.dataprep import IO
from gsea.gep import Gene_Expression_Profile

class Analysis():
    """
    """
    def __init__(self, params):
        
        self.rankby = params['rankby']
        self.permut = params['permut']
        self.p_weight = params['p_weight']
    
    def analyzefiles(self, gep_name, set_name, out_name, paths):
        """Performs the complete GSEA analysis, given filenames and paths.
        Creates ranked gene expression profiles and uses them for each row in a
        file containing multiple gene sets. 
        """
        
        # Find the input and output file locations. 
    
        gep_file = IO(gep_name, paths['input'])
        geneset_file = IO(set_name, paths['input'])
        out_file = IO(out_name, paths['output'])
        
        # Create ranked and permuted/ranked data.
        
        self.generate_gep(gep_file)
        
        self.gep.metric = self.rankby
        scores = self.gep.score
        ranked = self.gep.ranked
        
        permuted, perm_scores = self.gep.permutations(self.permut)
        
        print(self.gep.genes)
        print(ranked)
        print(scores)
        print(permuted)
        print(perm_scores)
        
        # Loop analysis over gene sets. 
        
        results = geneset_file.row_analysis(self.analyzeset, ranked, permuted,
                                                delim='\t')
        
        # Write the data to an output file. 
        
        out_file.writecsv(results)
        
    def generate_gep(self, file):
        """Creates a Gene Expression Profile object from an IO file object.
        """
        gep_data = file.load_array_with_labels(delim='\t')
        self.gep = Gene_Expression_Profile(*gep_data)
    
    def analyzeset(self, geneset, ranked, permuted):
        # es = ES_calc(self.ranked, self.params['p_weight'])
        # es_p = ES_P_calc(self.permuted, self.params['p_weight'])
        # p_stat = P_calc(es, es_p)
        # nes = NES_calc(es, es_p)
        # return (nes, p_stat)
        return geneset[:2]
        