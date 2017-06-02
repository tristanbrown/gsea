"""
This module contains the class Ranker, which ranks and sorts the genes in a gene
expression profile according to a given algorithm.  
"""

import numpy as np
import time

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
        return self.rank_by_metric(self.genes, self.phenos)
    
    def permuted_rank(self):
        """Gives the ranked and sorted gene labels after permuting the
        phenotype classes."""
        # print("Permuting and ranking genes.")
        return self.rank_by_metric(self.genes, 
                                        np.random.permutation(self.phenos))
    
    def permutations(self, m):
        """Returns an m x n array of m ranked n-length gene arrays, each
        generated from a permutation of the phenotype classes."""
        print("Permuter here.")
        
        n = len(self.genes)
        maxstr = len(max(self.genes, key=len))
        # permuted_genes = np.tile(self.genes, (m, 1))
        # permuted_genes = np.empty([m, n], dtype=('str', maxstr))
        permuted_genes = []
        # phenos = np.tile(self.phenos, (m, 1))
        print(permuted_genes)
        for row in range(m):
            # np.random.shuffle(phenos[row]) # 2.423s
            # permuted_genes[row] = self.rank_by_metric(self.genes, phenos[row])
            permuted_genes.append(self.permuted_rank())
        
        # print(phenos)
        permuted_genes = np.array(permuted_genes)
        print(permuted_genes)
        
        return permuted_genes
    
    def rank_by_metric(self, unranked, categories):
        """Returns an array of ranking indices for gene labels according to a
        metric used to score each line of data. The category (phenotype) labels
        must be taken into account in the metrics.
        """
        cat1 = np.where(categories == self.phenos[0])
        cat2 = np.where(categories != self.phenos[0])
        
        data1 = self.data[:,cat1[0]]
        data2 = self.data[:,cat2[0]]
        
        scores = self.metric(data1, data2)
        indices = np.argsort(scores)
        
        return unranked[indices]
    
    def select_metric(self, label=None):
        """Select a method for assigning a score to 1d arrays of numbers."""
        if label == None:
            raise
        elif label == 's2n':
            self.metric = self.signal2noise
        print("Metric set as %s." % label)
    
    def signal2noise(self, data1, data2):
        mean1 = np.mean(data1, axis=1)
        mean2 = np.mean(data2, axis=1)
        
        std1 = np.std(data1, axis=1)
        std2 = np.std(data2, axis=1)
        
        return (mean2 - mean1)/(std1 + std2)