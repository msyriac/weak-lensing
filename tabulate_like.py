#!/usr/bin/env python

import toy
import time, glob
import cPickle
from scipy import *
from math import *
import sys

class Gridder:
    def __init__ (self,N):
        self.T=toy.ToyGenerator(0,0)
        self.shearlist=[(0.0,0.0),(0.05, 0.00), (0.1, 0.0),  (0.05,0.05),
           (0.1,0.05) ]

        self.grid=[]
        self.N=N
        self.Ns=len(self.shearlist)
        for i in range(self.Ns):
            self.grid.append((zeros((N,N)),zeros((N,N))))

    def add_samples(self,N):

        for i in range(N):
            lst=self.T.generate_multishear(self.shearlist)
            for e12,grd in zip(lst, self.grid):
                if e12!=None:
                    e1,e2=e12

                    i1=int((1+e1)/2.*self.N)
                    i2=int((1+e2)/2.*self.N)
                    grd[0][i1,i2]+=1

                    esq=e1**2+e2**2
                    theta=atan2(e2,e1)+pi
                    i1=int(esq*self.N)
                    i2=int(theta/(2*pi)*self.N)
                    try:

                        grd[1][i1,i2]+=1
                    except:
                        print i1,i2,e1,e2,esq
            
    def dump_to_file(self,fn):
        cPickle.dump (self,open(fn,'w'),-1)

    def load_many(self,root):
        fnl=glob.glob(root+"_*.pickle")
        fnl.sort()
        for i,fn in enumerate(fnl):
            print "Loading ",fn
            if i==0:
                cself=cPickle.load(open(fn))
                # convert tuples to list:
                ng=[]
                for g in cself.grid:
                    ng.append(list(g))
                cself.grid=ng
                print len(ng)
            else:
                tmp=cPickle.load(open(fn))
                for grd,grd2 in zip(cself.grid,tmp.grid):
                    grd[0]+=grd2[0]
                    grd[1]+=grd2[1]
                print cself.grid[0][1]
        self.grid=cself.grid
        self.N=cself.N
        self.Ns=cself.Ns
        self.shearlist=cself.shearlist


class TabLike:
    def __init__(self, N,P, L_i, L_ij, Fisher, K1, K2, b3):
        self.N=N
        self.P=P
        self.L_i=L_i
        self.L_ij=L_ij
        self.Fisher=Fisher
        self.K1=K1
        self.K2=K2
        self.b3=b3
        ## we assume Fisher is diagonal here
        self.IFisher=array([[1/Fisher[0,0], 0],[0,1/Fisher[1,1]]])

    
    def P(self,e1m,e2m):
        i1=int((1+e1m)/2.*self.N)
        i2=int((1+e2m)/2.*self.N)
        P=self.P[i,i2]
        
    def FD(self,e1m,e2m):
        i1=int((1+e1m)/2.*self.N)
        i2=int((1+e2m)/2.*self.N)
        return array([self.L_i[0][i1,i2], self.L_i[1][i1,i2]])

    def SD(self,e1m,e2m):
        i1=int((1+e1m)/2.*self.N)
        i2=int((1+e2m)/2.*self.N)
        i00=self.L_ij[0][0][i1,i2]
        i10=self.L_ij[1][0][i1,i2]
        i11=self.L_ij[1][1][i1,i2]
        return array([[i00,i10],[i10,i11]])

def main():
    G=Gridder(200)
    Ni=200
    N=400000
    for i in range(Ni):
        print i,'/',Ni
        G.add_samples(N)
        G.dump_to_file('grids4/grid_%s.pickle'%sys.argv[1])



if __name__=="__main__":
    main()
