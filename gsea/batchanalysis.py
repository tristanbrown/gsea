"""
This module loops the Gene Set Enrichment Analysis over a file containing 
multiple gene sets, producing an output file tabulating the results. 
"""

import numpy as np

class Analyzer():
    def __init__(self, expression_file, geneset_file, outfile,
                    permut=1000, p_weight=1.0):
        
        print("Parameters loaded:")
        print(expression_file)
        print(geneset_file)
        print(outfile)
        print(permut)
        print(p_weight)
        