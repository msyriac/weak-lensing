#!/usr/bin/env python
from scipy import *
import scipy.linalg as la
import random

class Pairs:
    def __init__(self, sigmad=0.02):
        cross11=0.1
        cross12=0.2
        cross22=0.05
        sigmad2=sigmad**2
        self.cov = zeros((4,4))
        for i  in range(4):
            self.cov[i,i]=sigmad2
        self.cov[2,0]=sigmad2*cross11
        self.cov[3,1]=sigmad2*cross22
        self.cov[3,0]=sigmad2*cross12
        self.cov[2,1]=sigmad2*cross12

        ## and the other ones
        self.cov[0,2]=sigmad2*cross11
        self.cov[1,3]=sigmad2*cross22
        self.cov[0,3]=sigmad2*cross12
        self.cov[1,2]=sigmad2*cross12
        
        self.chol=la.cholesky(self.cov)
        print self.cov


    def draw (self):
        rv=[random.gauss(0,1.),random.gauss(0,1.),random.gauss(0,1.),random.gauss(0,1.)]
        return dot(rv,self.chol)



if __name__=="__main__":
    p=Pairs()
    sw=0
    sv=zeros(4)
    svv=zeros((4,4))
    for c in arange(100000):
        v=p.draw()
        sv+=v
        svv+=array([[x*y for x in v] for y in v])
        sw+=1

    print sv/sw
    print svv/sw
        
        
