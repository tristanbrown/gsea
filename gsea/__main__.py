"""
Main routine of GSEA. 
"""


import sys
import config

def main(args=None):
    if args is None:
        args = sys.argv[1:]
    print("Hello World")
    print(args)
    print(config.get['inputpath'])

if __name__ == "__main__":
    main()