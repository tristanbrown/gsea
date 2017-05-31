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

        if os.path.exists(self.fn):
            print("%s successfully loaded" % self.fn)
        else:
            raise
    
    def load_array(self):
        return np.genfromtxt(self.fn, delimiter='\t', dtype=None,
                                skip_header=0)
    
    def load_array_with_labels(self, delim=',', datatype='int',
                                    coltype='str', rowtype='str'):
        """Used for extracting data from text files with row and column labels.
        Returns a tuple containing 3 arrays: the data, the column labels,
        and the row labels.
        
        The optional arguments allow the user to choose the delimiter, as well
        as specifying the types for the labels and data. 
        """
        with open(self.fn, 'r') as f:
            col_labels = np.array(f.readline().split()[1:], dtype=coltype)
        num_col = len(col_labels)
        row_labels = np.genfromtxt(self.fn, delimiter=delim, dtype=rowtype, 
                                skip_header=1, usecols=0)
        data = np.genfromtxt(self.fn, delimiter=delim, dtype=datatype, 
                                skip_header=1, usecols=range(1, num_col+1))
        
        return (data, col_labels, row_labels)
    
    def loadgeneset(self):
        return []
    
    def writecsv(self, data):
        with open(self.fn, 'w', newline='\n') as f:
            writer = csv.writer(f)
            writer.writerows(data)
    