"""
Compute the number of moulins our domain should have using
observation-constrained density
"""

import numpy as np
import scipy.stats

# Domain
surf = lambda x: 6*(np.sqrt(x + 5e3) - np.sqrt(5e3))

# PDF: empirically derived from satellite data
pdf = lambda x: 42.12111317*scipy.stats.norm.pdf(x, loc=1138.25114462 - 400, scale=280.12484981)


dy = 25e3
dx = 10
x = np.arange(0, 100e3, dx)
z = surf(x)

n_moulins = 0
for (i, xi) in enumerate(x):
    zmean = z[i]
    n_moulins += pdf(zmean)*dx*dy/1e6

print(n_moulins)
