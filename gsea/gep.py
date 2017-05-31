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
        
    def rank(self, method):
        print("Ranker here.")
        return []
    
    def permute(self, iters):
        print("Permuter here.")
        return []