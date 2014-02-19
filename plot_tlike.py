#!/usr/bin/env python

from tabulate_like import *
import pylab, cPickle, scipy
print "starting"
G=cPickle.load(open('grids/combined.pickle'))




if False:
    for i,gf in enumerate(G.grid):
        pylab.figure()
        pylab.subplot(1,2,1)
        pylab.imshow(gf[0],extent=(-1,1,-1,1))
        pylab.colorbar()
        pylab.subplot(1,2,2)
        pylab.imshow(gf[0]-G.grid[0][0],extent=(-1,1,-1,1))
        pylab.colorbar()
        pylab.savefig('fig%i.pdf'%(i))

diffB=(G.grid[3][0]-G.grid[0][0])
#d2=(G.grid[3][0]+G.grid[4][0]-2*G.grid[0][0])/0.02
diffA=(G.grid[2][0]-G.grid[0][0])
s=0.005
d1=(8*diffA-diffB)/(6*s)
d3=(diffB-2*diffA)/(6*s**3)

pylab.figure()

pylab.subplot(2,2,1)
pylab.imshow(d1,extent=(-1,1,-1,1))
pylab.colorbar()
pylab.subplot(2,2,2)
pylab.imshow(d3,extent=(-1,1,-1,1))
pylab.colorbar()

pylab.subplot(2,2,3)
pylab.imshow(diffB,extent=(-1,1,-1,1))
pylab.colorbar()
pylab.subplot(2,2,4)
pylab.imshow(diffA,extent=(-1,1,-1,1))
pylab.colorbar()

pylab.savefig('figd.pdf')

d1=diffB=(G.grid[3][0]-G.grid[4][0])/(0.02)

pylab.figure()
 
pylab.subplot(2,2,1)
pylab.imshow(d1,extent=(-1,1,-1,1))
pylab.colorbar()
pylab.subplot(2,2,2)
pylab.imshow((d1[::-1,:])/2.0,extent=(-1,1,-1,1))
pylab.colorbar()
pylab.subplot(2,2,3)
pylab.imshow((d1[:,::-1]),extent=(-1,1,-1,1))
pylab.colorbar()
pylab.subplot(2,2,4)
pylab.imshow(df1,extent=(-1,1,-1,1))
pylab.colorbar()
pylab.savefig('figd2.pdf')



pylab.figure()
pylab.subplot(2,2,1)
pylab.imshow(df1a/Pa,extent=(-1,1,-1,1))
pylab.colorbar()

pylab.subplot(2,2,2)
pylab.imshow(df1b/Pb,extent=(-1,1,-1,1))
pylab.colorbar()

pylab.subplot(2,2,3)
pylab.imshow(df1a*df1b/(Pa*Pb)*P,extent=(-1,1,-1,1))
pylab.colorbar()


pylab.subplot(2,2,4)
pylab.imshow(df1*df1/P,extent=(-1,1,-1,1))
pylab.colorbar()
pylab.savefig('figQPx.pdf')

pylab.savefig('figQPx.pdf')


