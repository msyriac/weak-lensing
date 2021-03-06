#!/usr/bin/env python

from tabulate_like import *
import sys,scipy
import numpy as np

print "starting"

G=Gridder(200)

##option 1
#root='/astro/u/anze/work/mat/wl/'+'grids2'
#root='grids'
#pi=3 ## index of plus st
#mi=4 #3 index of minus st
#zi=0 #index of 0 shear
#st=0.01

#option2
root='/astro/u/anze/work/mat/wl/'+'grids2'
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

#G.dump_to_file(root+'/combined.pickle')

# get fishers, etc.
d1=(G.grid[pi][0]-G.grid[mi][0])/(2*st)

df1a=(d1+d1[:,::-1])/2.0
df1b=(-d1[::-1,:]-d1[::-1,::-1])/2.0
df1=(df1a+df1b)/2.0

P=G.grid[zi][0]
Pa=(P+P[:,::-1])/2.0
Pb=(P[::-1,:]+P[::-1,::-1])/2.0
P=(Pa+Pb)/2.0
#print P
#sys.exit()


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

for v in L_i:
    v[scipy.isnan(v)]=0.0
    v[scipy.isinf(v)]=0.0

for vl in L_ij:
    for v in vl:
        v[scipy.isnan(v)]=0.0
        v[scipy.isinf(v)]=0.0


#N=G+H
#K1 = \int L,, L, L,
#K1_1111 = L11 * L1 * L1

b=[0,1]
K1 = [[[[np.zeros((200,200)) for u in b] for v in b] for w in b] for x in b]
K2 = [[[[np.zeros((200,200)) for u in b] for v in b] for w in b] for x in b]
for i,j,k,l in [(u,v,w,x) for u in b for v in b for w in b for x in b]:
    print "Finding Ks for index ", i, j, k, l
    K1[i][j][k][l] = (L_ij[i][j]*L_i[k]*L_i[l]*P).sum()
    K2[i][j][k][l] = (L_ij[i][j]*L_ij[k][l]*P).sum()

print K1
print K2

#print Fisher.shape, L_i[0].shape, L_ij[0][0].shape

I=TabLike(G.N,P,L_i, L_ij, Fisher)
#cPickle.dump(I,open(root+'/tablike.pickle','w'),-1)
print "done"
