"""
This module contains the class Ranker, which ranks and sorts the genes in a gene
expression profile according to a given algorithm.  
"""

import numpy as np
import time

class Gene_Expression_Profile():
    """An object that contains expression profile data, each with a gene and 
    one of two phenotype labels. These gene labels can be ranked through
    metrics generated from statistical comparisons between the data from each
    phenotype. The phenotype labels can be permuted an given number of
    times to re-rank the gene labels and establish a baseline for significance
    of the unpermuted data. 
    """
    def __init__(self, data, genes, phenos):
        self.data = data
        self.genes = genes
        self.phenos = phenos
        
    @property
    def ranked(self):
        """Gives the ranked and sorted gene labels according to the initialized
        data."""
        try:
            return self._ranked
        except:
            print("Ranking genes.")
            self._ranked, self._corr =( 
                self.rank_by_metric(self.genes, self.phenos)
                )
            return self._ranked
    
    @property
    def corr(self):
        """Gives the sorted correlation scores for ranked gene labels."""
        try:
            return self._corr
        except:
            print("Ranking genes.")
            self.ranked
            return self._corr
    
    @property
    def num_perm(self):
        try:
            return self._num_perm
        except:
            self._num_perm = 1000
            return self._num_perm
    
    @num_perm.setter
    def num_perm(self, value):
        self._num_perm = value
    
    @property
    def permutations(self):
        try:
            return self._permutations
        except:
            self._permutations, self._permcorrs =(
                self.permute(self.num_perm)
                )
            return self._permutations
    
    @property
    def permcorrs(self):
        try:
            return self._permcorrs
        except:
            self.permutations
            return self._permcorrs
    
    def permute(self, m):
        """Returns an m x n array of m ranked n-length gene arrays, each
        generated from a permutation of the phenotype classes."""
        print("Permuting phenotypes and ranking genes %s times." % str(m))
        
        n = len(self.genes)
        
        maxstr = len(max(self.genes, key=len))
        genes = np.empty([m, n], dtype=('str', maxstr))
        corrs = np.empty([m, n], dtype='float')
        
        for i in range(m):
            genes[i], corrs[i] = self.permuted_rank()
        
        return (genes, corrs)
        
    def permuted_rank(self):
        """Gives the ranked and sorted gene labels after permuting the
        phenotype classes."""
        # print("Permuting and ranking genes.")
        return self.rank_by_metric(self.genes, 
                                        np.random.permutation(self.phenos))
    
    def rank_by_metric(self, unranked, categories):
        """Returns sorted arrays of gene labels and correlation scores,
        generated according to a ranking metric used to score each
        line of data. The category (phenotype) labels must be taken into account
        in the metrics.
        """
        cat1 = np.where(categories == self.phenos[0])
        cat2 = np.where(categories != self.phenos[0])
        
        data1 = self.data[:,cat1[0]]
        data2 = self.data[:,cat2[0]]
        
        corrs = self.metric(data1, data2)
        indices = np.argsort(corrs)
        
        return (unranked[indices], corrs[indices])
    
    @property
    def metric(self):
        try:
            return self._metric
        except:
            self.metric = 'rat' # Default metric.
            return self._metric
    
    @metric.setter
    def metric(self, label=None):
        """Select a method for assigning a correlation score to 1d arrays of
        numbers."""
        if label == 's2n':
            self._metric = self.signal2noise
        elif label == 'diff':
            self._metric = self.diff_classes
        elif label == 'rat':
            self._metric = self.ratio_classes
        else:
            raise
        print("Metric set as %s." % label)
    
    ##### Ranking Metrics #####
    
    def signal2noise(self, data1, data2):
        mean1 = np.mean(data1, axis=1)
        mean2 = np.mean(data2, axis=1)
        
        std1 = np.std(data1, axis=1)
        std2 = np.std(data2, axis=1)
        
        return (mean2 - mean1)/(std1 + std2)
    
    def diff_classes(self, data1, data2):
        mean1 = np.mean(data1, axis=1)
        mean2 = np.mean(data2, axis=1)
        
        return mean2 - mean1
    
    def ratio_classes(self, data1, data2):
        mean1 = np.mean(data1, axis=1)
        mean2 = np.mean(data2, axis=1)
        
        return mean2/mean1