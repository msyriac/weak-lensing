#!/usr/bin/env python

import glob
import csv
import numpy
from math import *
import matplotlib.pyplot as p
import scipy.linalg as l
from scipy import *
from numba import autojit


#import generate_pairs as gen

filesG = glob.glob("datan/g*.csv")
filesH = glob.glob("datan/h*.csv")

CinvG=numpy.zeros((2,2))
CinvH=numpy.zeros((2,2))
Cinv_sumG=numpy.zeros((2,2))
Cinv_sumH=numpy.zeros((2,2))
GHsum=numpy.zeros((4,1))
GHsumw=numpy.zeros((4,1))

Q=numpy.zeros((4,4))
Q_sumG=numpy.zeros((2,1))
Q_sumH=numpy.zeros((2,1))

#p=gen.Pairs()


i=0
j=0

freq=1000

x=[]
w=[]


def getMats(valsG):
    PG=valsG[0]
    QG=numpy.array([[valsG[1]],
                    [valsG[2]]])
    RG=numpy.array( [ [valsG[3], valsG[4]],
                      [valsG[4], valsG[5]] ] )
    
    Q_normG=(QG/PG)
    
    CinvG=((numpy.dot(Q_normG,Q_normG.transpose()))-(RG/PG))
    return Q_normG, CinvG


def process(valsG, valsH, Cinv_sumG, Cinv_sumH, x):
    
    Q_normG,CinvG=getMats(valsG)
    Cinv_sumG+=CinvG
    Q_normH,CinvH=getMats(valsH)
    Cinv_sumH+=CinvH

    x.append( (Q_normG.flat, Q_normH.flat) )
    


for fname in filesG[:][:]:
    fG = open(fname, 'r')
    fH = open(fname.replace('g','h'),'r')
              
    readerG = csv.reader(fG,delimiter=',')
    readerH = csv.reader(fH,delimiter=',')
    for rowG, rowH in zip(readerG,readerH):
        i+=1
        if (i % freq)==0: print i
        valsG = [float(r) for r in (rv for rv in rowG if rv!='')]
        valsH = [float(r) for r in (rv for rv in rowH if rv!='')]
        
        process(valsG, valsH,Cinv_sumG, Cinv_sumH, x)

        
    fG.close()
    fH.close()
        


Cinv=Cinv_sumG+Cinv_sumH
CC=l.inv(Cinv)*2*i
Cd=(CC[0,0]+CC[1,1])/2.0
CC=array([[Cd,0],[0,Cd]])

print Cinv_sumG
print Cinv_sumH
print CC


Nch=30

Q=[numpy.zeros((2,2)) for a in range(Nch)]
sw=zeros(Nch)
i=0

for QG,QH in x:
    j=i%Nch
    M=numpy.outer(dot(CC,QG),dot(CC,QH))
    Q[j]+=M
    sw[j]+=1
    i+=1

S=zeros((2,2))
SS=zeros((2,2))

for j in range(Nch):
    Q[j]/=sw[j]
    print Q[j]
    S+=Q[j]
    SS+=Q[j]*Q[j]
    

S/=Nch
SS/=Nch
SS-=S*S
print '--------'
print S
print sqrt(SS)/sqrt(Nch)
print CC[0,0]/sqrt(i)





