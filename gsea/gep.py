"""
This module contains the class Ranker, which ranks and sorts the genes in a gene
expression profile according to a given algorithm.  
"""

import numpy as np

class Gene_Expression_Profile():
    def __init__(self, data, genes, phenos):
        self.data = data
        self.genes = genes
        self.phenos = phenos
        
        # Setting the default metric:
        self.select_metric('s2n')
        
    def rank_genes(self):
        """Gives the ranked and sorted gene labels according to the initialized
        data."""
        print("Ranking genes.")
        return self.rank_by_metric(self.phenos, self.metric)
    
    def permuted_rank(self):
        """Gives the ranked and sorted gene labels after permuting the
        phenotype classes."""
        print("Permuting and ranking genes.")
        return self.rank_by_metric(np.random.permutation(self.phenos),
                                        self.metric)
    
    def permutations(self, m):
        """Returns an m x n array of m ranked n-length gene arrays, each
        generated from a permutation of the phenotype classes."""
        print("Permuter here.")
        n = len(self.genes)
        permuted_genes = np.tile(self.genes, (m, 1))
        permuted_genes[0][0] = 'mod'
        
        print(permuted_genes)
        return permuted_genes
    
    def rank_by_metric(self, categories, metric):
        """Returns a ranked and sorted array of gene labels according to a
        metric used to score each line of data. The category (phenotype) labels
        must be taken into account in the metrics.
        """
        ranked_genes = np.copy(self.genes)
        
        return ranked_genes
    
    def select_metric(self, label=None):
        """Select a method for assigning a score to 1d arrays of numbers."""
        if label == None:
            raise
        elif label == 's2n':
            self.metric = self.signal2noise
        print("Metric set as %s." % label)
    
    def signal2noise(self, data, categories):
        return 0