'''
Owner: Qian Yang 06/20/2016
Last edit: Enze Chen 10/22/2018

This script is a wrapper used to extract molecular dynamics data.
It is capable of handling multiple files at once. Place it in the
same folder as the 'log_read' files.
    - XXX.80 are examples that you can run on.
'''

import log_read_md as readMols
import log_read_reacts as readReacts

maxframenum = 47000   # varies with file. Just set to the maximum

molfiles = []
reactfiles = []
sets = ['80']       # could be useful if you want to process multiple files at once

for i in sets:
    print(i)    
    molfiles.append('mine.out')
    # reactfiles.append('mine.out')
    
# [moleculedict, newmoldict, bonded_reaction, bond_dist] = readMols.getMols(molfiles, maxframenum)
readMols.getMols(molfiles, maxframenum)
# readReacts.getReacts(reactfiles, maxframenum, moleculedict, newmoldict, bonded_reaction, bond_dist)
