"""
This module contains the functions necessary to extract data from text files
and put it into numpy arrays to be usable by other functions. 
"""

import os
import csv
import numpy as np

class IO():
    """An object that points to a desired file location and either extracts data
    from an existing text file as a numpy array, or writes data to a text file. 
    """
    def __init__(self, filename, dir=''):
        print("Extractor here.")
        self.fn = os.path.join(dir, filename)
        print("File found?")
        print(os.path.exists(self.fn))
    
    def loadarray(self):
        return np.genfromtxt(self.fn, delimiter='\t', dtype=None,
                                skip_header=0)
    
    def loadgeneset(self):
        return []
    
    def writecsv(self, data):
        with open(self.fn, 'w', newline='\n') as f:
            writer = csv.writer(f)
            writer.writerows(data)
    