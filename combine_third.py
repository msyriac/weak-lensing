#!/usr/bin/env python

from tabulate_like import *
import sys,scipy
import pylab
print "starting"

G=Gridder(200)
root='grids4'
G.load_many(root+'/grid')
print G.shearlist
#normalize
print G.shearlist
for i,gf in enumerate(G.grid):
    print i,gf[0].sum(), gf[1].sum()
    s=gf[0].sum()
    gf[0]/=s
    gf[1]/=s
    if False:
        pylab.figure()
        pylab.imshow(gf[0],extent=(-1,1,-1,1))
        pylab.colorbar()
        pylab.savefig('pic%i.pdf'%i)
        pylab.figure()
        pylab.imshow(gf[1],extent=(-1,1,-1,1))
        pylab.colorbar()
        pylab.savefig('picP%i.pdf'%i)
    
grid_0_0=G.grid[0][0]
grid_1_0=G.grid[1][0]
grid_2_0=G.grid[2][0]
grid_1_1=G.grid[3][0]
grid_2_1=G.grid[4][0]
grid_m1_0=grid_1_0[::-1,:]
grid_m2_0=grid_2_0[::-1,:]

grid_m1_1=grid_1_1[::-1,:]
grid_1_m1=grid_1_1[:,::-1]
grid_m1_m1=grid_1_1[::-1,::-1]
st=0.05

if False:
    for i,g in enumerate([grid_1_1,grid_m1_1, grid_1_m1, grid_m1_m1]):
            pylab.figure()
            pylab.imshow(g,extent=(-1,1,-1,1))
            pylab.colorbar()
            pylab.savefig('dbg%i.pdf'%i)
    

dersc=[[ 0,       0,    1.,       0,      0],
      [1./12,	-2./3,	0,	2./3,	-1./12],
      [-1./12,	4./3,	-5./2,	4./3,	-1./12],
      [-1./2,	1.,	0,	-1.,	1./2],
      [1.,	-4.,	6.,	-4.,	1. ]]

pvec=[grid_m2_0, grid_m1_0, grid_0_0, grid_1_0, grid_2_0]
ders=[]
for do,l in enumerate(dersc):
    pic=reduce(lambda x,y:x+y, map(lambda x,y:x*y,pvec,l))/st**do
    ders.append(pic)
    pylab.figure()
    pylab.imshow(pic,extent=(-1,1,-1,1))
    pylab.colorbar()
    pylab.savefig('ders%i.pdf'%do)


crossder=(grid_1_1+grid_m1_m1-grid_1_m1-grid_m1_1)/(4*st**2)
pylab.figure()
pylab.imshow(crossder,extent=(-1,1,-1,1))
pylab.colorbar()
pylab.savefig('crossder.pdf')


def fixar(ar):
    print ar
    ar[scipy.isnan(ar)]=0.0
    ar[scipy.isinf(ar)]=0.0
    return ar

P=ders[0]
P=0.25*(P+P[::-1,:]+P[:,::-1]+P[::-1,::-1])
L_i=[fixar(ders[1]/ders[0]), fixar(transpose(ders[1]/ders[0]))]
L_ij=[[fixar(ders[2]/ders[0]-L_i[0]**2), fixar(crossder/ders[0]-L_i[0]*L_i[1])],
      [fixar(crossder/ders[0]-L_i[0]*L_i[1]), fixar(transpose(ders[2]/ders[0]))-L_i[1]**2]]


F=zeros((2,2))
for i in range(2):
    for j in range(2):
        print "zero=",((L_i[i]*L_i[j]+L_ij[i][j])*P).sum()
        Fisher=((L_i[i]*L_i[j])*P).sum()
        print "Fisher=",Fisher
        F[i,j]=Fisher
        
print F
## now let's fix F according to what we think it should be
Fd=(F[0,0]+F[1,1])/2
F=array([[Fd,0],[0,Fd]])
FI=array([[1/Fd,0],[0,1/Fd]])
print F

# b=[0,1]
# K1 = zeros((2,2,2,2))
# K2 = zeros((2,2,2,2))
# for i,j,k,l in [(u,v,w,x) for u in b for v in b for w in b for x in b]:
#     print "Finding Ks for index ", i, j, k, l
#     K1[i,j,k,l] = (L_ij[i][j]*L_i[k]*L_i[l]*P).sum()
#     K2[i,j,k,l] = (L_ij[i][j]*L_ij[k][l]*P).sum()

# print K1
# print K2


## Let's do one dimension, can rotate later, use letter Q for derivative
Q_i=L_i[0]
Q_ii=L_ij[0][0]
Q_iii=fixar(ders[3]/ders[0]-3*ders[2]*ders[1]/ders[0]**2 + 2*ders[1]**3/ders[0]**3)

QFi=1/(Q_i*Q_i*P).sum()
print QFi,FI[0,0]

b3a=QFi*(Q_i*(Q_iii+3*Q_ii*Q_i+Q_i**3)*P).sum()/6.0

print "bias b3a=",b3a

I=TabLike(G.N,P,L_i, L_ij, F, K1, K2, b3)
cPickle.dump(I,open(root+'/tablike.pickle','w'),-1)
print "done"

