#!/usr/bin/env python

from tabulate_like import *
import sys,scipy
print "starting"

G=Gridder(200)

##option 1
root='grids'
pi=3 ## index of plus st
mi=4 #3 index of minus st
zi=0 #index of 0 shear
st=0.01

#option2
root='grids2'
pi=4
mi=1
zi=0
st=0.01


G.load_many(root+'/grid')

#normalize

print G.shearlist
for i,gf in enumerate(G.grid):
    print i,gf[0].sum(), gf[1].sum()
    s=gf[0].sum()
    gf[0]/=s
    gf[1]/=s

G.dump_to_file(root+'/combined.pickle')

# get fishers, etc.
d1=(G.grid[pi][0]-G.grid[mi][0])/(2*st)

df1a=(d1+d1[:,::-1])/2.0
df1b=(-d1[::-1,:]-d1[::-1,::-1])/2.0
df1=(df1a+df1b)/2.0

P=G.grid[zi][0]
Pa=(P+P[:,::-1])/2.0
Pb=(P[::-1,:]+P[::-1,::-1])/2.0
P=(Pa+Pb)/2.0


F1=df1a*df1b/(Pa*Pb)*P
F1[scipy.isnan(F1)]=0.0
F1[scipy.isinf(F1)]=0.0

F2=df1*df1/P
F2[scipy.isnan(F2)]=0.0
F2[scipy.isinf(F2)]=0.0

## two ways of getting fisher.
print F1.sum(), F2.sum()



F3=df1*transpose(df1)/P
F3[scipy.isnan(F3)]=0.0
F3[scipy.isinf(F3)]=0.0
##off-diag term
print F3.sum()


Fisher=array([[F1.sum(),0],[0,F1.sum()]])
L_i=[df1/P, transpose(df1/P)]
L_ij=[[df1a*df1b/(Pa*Pb),df1a*transpose(df1b)/(Pa*transpose(Pb))],
      [df1a*transpose(df1b)/(Pa*transpose(Pb)),transpose(df1a*df1b/(Pa*Pb))]]

I=TabLike(G.N,P,L_i, L_ij, Fisher)
cPickle.dump(I,open(root+'/tablike.pickle','w'),-1)
print "done"
