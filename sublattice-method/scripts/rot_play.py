import numpy as np
import sys
sys.path.insert(0, '/home/david/Work/DFT-Surface-Science/sublattice-method/funcs')
from funcs import rot
for i in range(2):
    for j in range(2):
        for k in range(2):
            print([i,j,k],[90,0,0],np.round(rot([[i,j,k]],90,0,0),2))
            print([i,j,k],[0,90,0],np.round(rot([[i,j,k]],0,90,0),2))
            print([i,j,k],[0,0,90],np.round(rot([[i,j,k]],0,0,90),2))
