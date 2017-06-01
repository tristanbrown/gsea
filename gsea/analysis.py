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
        
        gep_data = gep_file.load_array_with_labels(delim='\t')
        gep = Gene_Expression_Profile(*gep_data)
        
        gep.select_metric(self.rankby)
        ranked = gep.rank_genes()
        permuted = gep.permutations(self.permut)
        
        # Loop analysis over gene sets. 
        
        results = geneset_file.row_analysis(self.analyzeset, ranked, permuted,
                                                delim='\t')
        
        # Write the data to an output file. 
        
        out_file.writecsv(results)
        
    def generate_gep(self):
        """Creates a Gene Expression Profile object and uses it to acquire a
        ranked 
        """
        pass
    
    def analyzeset(self, geneset, ranked, permuted):
        # es = ES_calc(self.ranked, self.params['p_weight'])
        # es_p = ES_P_calc(self.permuted, self.params['p_weight'])
        # p_stat = P_calc(es, es_p)
        # nes = NES_calc(es, es_p)
        # return (nes, p_stat)
        return geneset[:2]
        