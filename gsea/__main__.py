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
    
    # Acquire the data and check the output path. 
    
    gep_data = IO(args[0], config.path['input']).loadarray()
    geneset_data = IO(args[1], config.path['input']).loadgeneset()
        #Acquire ahead of time and loop over the array, or just acquire one row
        #at a time? 
    out_file = IO('results.csv', config.path['output'])
    
    # Analyze the data.
    
    A = Analyzer(gep_data, geneset_data, config.analysis)
    results = A.analyzesets()
    
    # Write the results.
    
    out_file.writecsv(results)

if __name__ == "__main__":
    main()