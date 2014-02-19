#!/usr/bin/env python

from tabulate_like import *
import pylab, cPickle, scipy
print "starting"
G=cPickle.load(open('grids/combined.pickle'))



print G.shearlist
for i,gf in enumerate(G.grid):
    print i,gf[0].sum(), gf[1].sum()
    s=gf[0].sum()
    gf[0]/=s
    gf[1]/=s

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
df1a=(d1+d1[:,::-1])/2.0
df1b=(-d1[::-1,:]-d1[::-1,::-1])/2.0
df1=(df1a+df1b)/2.0
pylab.imshow(df1,extent=(-1,1,-1,1))
pylab.colorbar()
pylab.savefig('figd2.pdf')

P=G.grid[0][0]
Pa=(P+P[:,::-1])/2.0
Pb=(P[::-1,:]+P[::-1,::-1])/2.0
P=(Pa+Pb)/2.0


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

F1=df1a*df1b/(Pa*Pb)*P
F1[scipy.isnan(F1)]=0.0
print F1
print F1.sum()
for i,row in enumerate(F1):
    print i, row,row.sum()
F2=df1*df1/P


I=TabLike(G.grid[0][0], d1, G.N)
cPickle.dump(I,open('grids/tablike.pickle','w'))
