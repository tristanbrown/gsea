"""
This is the config file for GSEA, containing various parameters the user may
wish to modify. 
"""

path = dict(
    input = 'data',
    output = 'data',

    )

analysis = dict(
    permutations = 1000,
    p_weight = 1.0,
    rankby = 's2n'
    )