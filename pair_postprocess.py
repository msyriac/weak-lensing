#!/usr/bin/env python

import glob
import csv
import numpy
from math import *
import matplotlib.pyplot as p

filesG = glob.glob("data/g*.csv")
filesH = glob.glob("data/h*.csv")

CinvG=numpy.zeros((2,2))
CinvH=numpy.zeros((2,2))
GHsum=numpy.zeros((4,1))
Q=numpy.zeros((4,4))
i=0

freq=1000

x=[]

for fname in filesG:
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
        CinvG=((numpy.dot(Q_normG,Q_normG.transpose()))-(RG/PG))
        CmG=numpy.linalg.inv(CinvG)
        infGp=numpy.dot(CmG,Q_normG)



        PH=valsH[0]
        QH=numpy.array([[valsH[1]],
                        [valsH[2]]])
        RH=numpy.array( [ [valsH[3], valsH[4]],
                          [valsH[4], valsH[5]] ] )

        Q_normH=(QH/PH)
        CinvH=((numpy.dot(Q_normH,Q_normH.transpose()))-(RH/PH))
        CmH=numpy.linalg.inv(CinvH)
        infHp=numpy.dot(CmH,Q_normH)

        x.append( numpy.array( [[infGp.flat[0]],
                                [infGp.flat[1]],
                                [infHp.flat[0]],
                                [infHp.flat[1]] ] )  )
        GHsum+=x[-1]

GHavg=GHsum/i
for g in x:
    Q+=(1/(i-1.)) * numpy.dot(g-GHavg,(g-GHavg).transpose())  


#print Q


    
print i, " lines read."




fG.close()
fH.close()
