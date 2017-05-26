"""
Main routine of GSEA. 
"""

import sys
import os
from . import config
from .batchanalysis import Analyzer

def main(args=None):
    if args is None:
        args = sys.argv[1:]
        
    infile1 = os.path.join(config.path['input'], args[0])
    infile2 = os.path.join(config.path['input'], args[1])
    outfile = os.path.join(config.path['output'], 'results.csv')
    
    print("File found?")
    print(os.path.exists(infile1))
    print(os.path.exists(infile2))
    print(os.path.exists(outfile))
    
    A = Analyzer(infile1, infile2, outfile, config.analysis)

if __name__ == "__main__":
    main()