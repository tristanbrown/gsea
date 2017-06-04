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
        self.gep.num_perm = self.permut
        
        ranked = self.gep.ranked
        corr = self.gep.score
        
        permuted = self.gep.permutations
        permcorr = self.gep.permscores
        
        print(self.gep.genes)
        print(ranked)
        print(corr)
        print(permuted)
        print(permcorr)
        
        self.generate_genesets(geneset_file)
        
        ES = {name : self.calc_ES(ranked, corr, S)
                        for name, S in self.geneset.items()}
        print(ES.values())
        
        # Calculate NES and P-stat
        
        # Write the data to an output file. 
        
        results = ES
        out_file.writecsv(results)
    
    #Convert these generate functions into @property types. 
    def generate_gep(self, file):
        """Creates a Gene Expression Profile object from an IO file object.
        """
        gep_data = file.load_array_with_labels(delim='\t')
        self.gep = Gene_Expression_Profile(*gep_data)
    
    def generate_genesets(self, file):
        """Creates a dictionary of gene sets from an IO file object."""
        self.geneset = file.load_arb_rows(delim='\t', type='str', skipcol=1)
    
    def analyzeset(self, geneset, ranked, permuted):
        # es = ES_calc(self.ranked, self.params['p_weight'])
        # es_p = ES_P_calc(self.permuted, self.params['p_weight'])
        # p_stat = P_calc(es, es_p)
        # nes = NES_calc(es, es_p)
        # return (nes, p_stat)
        return geneset[:2]
    
    def calc_ES(self, ranked, corr, geneset):
        """Takes an array of gene labels from a GEP, the correlation scores by
        which they were ranked, and an independent gene set. Returns the ES, as
        the maximum value of the running sum of the correlation-weighted 
        fraction of ranked genes present in the geneset. 
        """
        N = len(ranked)
        Nh = len(geneset)
        
        hits = np.in1d(ranked, geneset)
        if not np.any(hits):
            P_hit = hits
        else:
            hits_wtd = np.cumsum(abs(hits * corr)**self.p_weight)
            P_hit = hits_wtd / hits_wtd[-1]
            
        P_miss = np.cumsum((1 - hits)/(N - Nh))
        
        return np.amax(P_hit - P_miss)