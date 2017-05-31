"""
This is the config file for GSEA, containing various parameters the user may
wish to modify. 
"""

path = dict(
    input = 'data',
    output = 'data',

    )

analysis = dict(
    rankby = 's2n',
    permut = 1000,
    p_weight = 1.0
    )