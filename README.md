# GSEA
Python implementation of the Gene Set Enrichment Analysis method

# Installation and dependencies
Above the directory containing setup.py, install using:
    >>pip install -e gsea

This software has been tested in an Anaconda environment containing the 
following packages and their respective dependencies:
- python (v. 3.6.1)
- numpy (v. 1.12.1)

# Usage
Currently, the package can be run from the root directory as follows:

    >> python gsea <filename1> <filename2>
    
where filename1 refers to the Gene Expression Profile data, and filename2 refers
to the Gene Set data. 
The gsea/config.py file determines the location of the data (must be in the 
python path), the output filename, as well as several analysis parameters. 

config.py is meant to be easily edited by the user. 

The default filename for outputting the results is "results.csv". 

# Testing
Run a testing script (e.g. test1.py) as follows:

    >> python -m unittest test.test1

or all tests via:

    >> python -m unittest discover