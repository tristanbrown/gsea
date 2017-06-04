"""
Main routine of GSEA. 
"""

import sys
import config
from analysis import Analysis

def main(args=None):
    
    # Determine filenames.
    
    if args is None:
        args = sys.argv[1:]
 
    out_name = 'results.csv'
    
    # Analyze the data and write the results.
    
    A = Analysis(args[0], args[1], config.path['input'], config.analysis)
    A.analyzefiles(out_name, config.path['output'])

if __name__ == "__main__":
    main()