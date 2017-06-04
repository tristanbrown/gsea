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
    def __init__(self, gep_name, geneset_name, in_path, params):
        
        self.rankby = params['rankby']
        self.permut = params['permut']
        self.p_weight = params['p_weight']
        
        self.gep_file = IO(gep_name, in_path)
        self.gs_file = IO(geneset_name, in_path)
    
    def analyzefiles(self, out_name, out_path):
        """Performs the complete GSEA analysis, given an output file location.
        Creates ranked gene expression profiles and uses them for each row in a
        file containing multiple gene sets. 
        """
        # Find the output file location. 
    
        out_file = IO(out_name, out_path)
        
        # Create ranked and permuted/ranked data.
        
        self.gep.metric = self.rankby
        self.gep.num_perm = self.permut

        print(self.hits)
        print(np.argsort(self.gep.corr))
        print(self.ES)
        print(len(self.ES))
        
        # Calculate NES and P-stat
        
        self.ES_null
        print(self.ES_null)
        
        # Write the data to an output file. 
        
        # results = self.ES
        # out_file.writecsv(results)
    
    @property
    def gep(self):
        """If not found, creates a Gene Expression Profile object from an IO
        file object.
        """
        try:
            return self._gep
        except:
            data, phenos, self._genes =( 
                self.gep_file.load_array_with_labels(delim='\t')
                )
            self._gep = Gene_Expression_Profile(data, phenos)
            return self._gep
    
    @property
    def genes(self):
        """Gives the array of gene labels from the Gene Expression Profile.
        """
        try:
            return self._genes
        except:
            self.gep
            return self._genes
    
    @property
    def genesets(self):
        """If not found, creates a dictionary of gene sets from an IO file
        object."""
        try:
            return self._genesets
        except:
            self._setlabels, self._genesets =(
                self.gs_file.load_arb_rows(delim='\t', type='str', skipcol=1)
                )
            return self._genesets
    
    @property
    def setlabels(self):
        """Returns the labels of the genesets."""
        try:
            return self._setlabels
        except:
            self.genesets
            return self._setlabels
    
    @property
    def ES(self):
        """If not found, runs the calculation to generate an array of
        Enrichment Scores, labeled by gene set title. 
        """
        try:
            return self._ES
        except:
            self._ES = self.calc_ES(self.gep.corr)
            return self._ES
    
    @property
    def hits(self):
        """Gives a 2D boolean array in which each row represents the presence of
        the genes in a particular gene set.
        """
        try:
            return self._hits
        except:
            self._hits =(
                np.array([np.in1d(self.genes, S) for S in self.genesets])
                )
            return self._hits
    
    @property
    def ES_null(self):
        """If not found, generates a 2D array of enrichment scores from the 
        phenotype-permuted ranked genes for the gene sets.
        """
        try:
            return self._ES_null
        except:
            self._ES_null = [self.calc_ES(row)
                for row in self.gep.permcorrs]
            return self._ES_null
    
    def calc_ES(self, corr):
        """Takes an array of correlation scores. Returns the ES, as
        the maximum value of the running sum of the correlation-weighted 
        fraction of ranked genes present in the geneset. 
        """
        N = len(self.genes)
        Nh = np.vstack(np.array([len(S) for S in self.genesets]))
        
        sort = np.argsort(corr)
        hits = self.hits[:,sort]

        hits_wtd = np.cumsum(abs(hits * corr)**self.p_weight, axis=-1)
        P_hit = np.nan_to_num(hits_wtd / np.vstack(hits_wtd[:,-1]))
            
        P_miss = np.cumsum((1 - hits)/(N - Nh), axis=-1)
        
        return np.amax(P_hit - P_miss, axis=-1)