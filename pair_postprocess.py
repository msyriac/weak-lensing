#!/usr/bin/env python

import glob
import csv
import numpy
from math import *
import matplotlib.pyplot as p
import scipy.linalg as l
from scipy import *

#import generate_pairs as gen

filesG = glob.glob("data/pairs_1/g*.csv")
filesH = glob.glob("data/pairs_1/h*.csv")

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

        PG=valsG[0]
        QG=numpy.array([[valsG[1]],
                        [valsG[2]]])
        RG=numpy.array( [ [valsG[3], valsG[4]],
                          [valsG[4], valsG[5]] ] )

        Q_normG=(QG/PG)
        #Q_sumG+=Q_normG ####
        CinvG=((numpy.dot(Q_normG,Q_normG.transpose()))-(RG/PG))
        #Cinv_sumG+=CinvG ###
        CmG=numpy.linalg.inv(CinvG)
        infGp=numpy.dot(CmG,Q_normG)



        PH=valsH[0]
        QH=numpy.array([[valsH[1]],
                        [valsH[2]]])
        RH=numpy.array( [ [valsH[3], valsH[4]],
                          [valsH[4], valsH[5]] ] )

        Q_normH=(QH/PH)
        #Q_sumH+=Q_normH ####
        CinvH=((numpy.dot(Q_normH,Q_normH.transpose()))-(RH/PH))
        #Cinv_sumH+=CinvH ####
        CmH=numpy.linalg.inv(CinvH)
        infHp=numpy.dot(CmH,Q_normH)


        #print infHp,CmH
        #print infGp,CmG

        DG=CmG[0,0]*CmG[1,1]-CmG[0,1]**2
        # TrG=CmG[0,0]+CmG[1,1]
        # E1G=TrG/2.0+sqrt(TrG**2/4-DG)
        # E2G=TrG/2.0-sqrt(TrG**2/4-DG)
        # GpD=(E1G>0) and (E2G>0)

        DH=CmH[0,0]*CmH[1,1]-CmH[0,1]**2
        # TrH=CmH[0,0]+CmH[1,1]
        # E1H=TrH/2.0+sqrt(TrH**2/4-DH)
        # E2H=TrH/2.0-sqrt(TrH**2/4-DH)
        # HpD=(E1H>0) and (E2H>0)

        if (any(abs(infGp)>1)) or (any(abs(infHp)>1)):
            continue
        
        
        weG=abs(1/(DG))
        vecG=numpy.array(infGp.flat)
        weH=abs(1/(DH))
        vecH=numpy.array(infHp.flat)
        x.append( (weG,vecG,weH,vecH) )
        

        
    fG.close()
    fH.close()
        

Nch=12

Q=[numpy.zeros((2,2)) for a in range(Nch)]
sw=zeros(Nch)
i=0

for weG,vecG,weH,vecH in x:
    j=i%Nch
    M=numpy.outer(vecG,vecH)
    Q[j]+=M*weG*weH
    sw[j]+=weG*weH
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

print S
print sqrt(SS)/sqrt(Nch)





