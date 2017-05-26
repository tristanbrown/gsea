"""
This module contains the functions necessary to extract data from text files
and put it into numpy arrays to be usable by other functions. 
"""

import os
import numpy as np

class Extractor():
    def __init__(self, filename, dir):
        print("Extractor here.")