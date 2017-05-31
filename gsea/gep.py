"""
This module contains the class Ranker, which ranks and sorts the genes in a gene
expression profile according to a given algorithm.  
"""

import numpy as np

class Gene_Expression_Profile():
    def __init__(self, data):
        self.data = data
        
    def rank(self, method):
        print("Ranker here.")