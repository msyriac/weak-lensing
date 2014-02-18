#!/usr/bin/env python

from tabulate_like import *

G=Gridder(200)

G.load_many('grids/grid')

import pylab

if True:
    for i,gf in enumerate(G.grid):
        pylab.figure()
        pylab.subplot(1,2,1)
        pylab.imshow(gf[0],extent=(-1,1,-1,1))
        pylab.colorbar()
        pylab.subplot(1,2,2)
        pylab.imshow(gf[0]-G.grid[0][0],extent=(-1,1,-1,1))
        pylab.colorbar()
        pylab.savefig('fig%i.pdf'%(i))

d1=(G.grid[3][0]-G.grid[4][0])/0.02
#d2=(G.grid[3][0]+G.grid[4][0]-2*G.grid[0][0])/0.02
d12=(G.grid[1][0]-G.grid[0][0])/0.001

pylab.figure()

pylab.subplot(2,2,1)
pylab.imshow(d1,extent=(-1,1,-1,1))
pylab.colorbar()
pylab.subplot(2,2,2)
pylab.imshow(d12,extent=(-1,1,-1,1))
pylab.colorbar()

pylab.subplot(2,2,3)
pylab.imshow(d1-d12,extent=(-1,1,-1,1))
pylab.colorbar()
pylab.subplot(2,2,4)
pylab.imshow(d1/d12,extent=(-1,1,-1,1),vmin=0.5, vmax=1.5)
pylab.colorbar()


pylab.savefig('figd.pdf')
