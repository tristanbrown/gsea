"""
Main routine of GSEA. 
"""

import sys
import config
from dataprep import IO
from batchanalysis import Analyzer

def main(args=None):
    if args is None:
        args = sys.argv[1:]
        
    gep_file = IO(args[0], config.path['input'])
    geneset_file = IO(args[1], config.path['input'])
    out_file = IO('results.csv', config.path['output'])
    
    A = Analyzer(gep_file.loadarray(),
                    geneset_file.loadgeneset(),
                    config.analysis
                )
    results = A.analyzesets()
    out_file.writecsv(results)

if __name__ == "__main__":
    main()