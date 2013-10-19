#!/usr/bin/env python

import glob
import csv
import numpy
from math import *
import matplotlib.pyplot as p
#import generate_pairs as gen

filesG = glob.glob("data/g*.csv")
filesH = glob.glob("data/h*.csv")

CinvG=numpy.zeros((2,2))
CinvH=numpy.zeros((2,2))
Cinv_sumG=numpy.zeros((2,2))
Cinv_sumH=numpy.zeros((2,2))
GHsum=numpy.zeros((4,1))
Q=numpy.zeros((4,4))
Q_sumG=numpy.zeros((2,1))
Q_sumH=numpy.zeros((2,1))

#p=gen.Pairs()


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

        x.append( numpy.array( [[infGp.flat[0]],
                                [infGp.flat[1]],
                                [infHp.flat[0]],
                                [infHp.flat[1]] ] )  )
        
        #x.append(p.draw())
        #x[-1].resize(4,1)

        GHsum+=x[-1]
        
        


#CmsumG=numpy.linalg.inv(Cinv_sumG)
#CmsumH=numpy.linalg.inv(Cinv_sumH)
#infGsump=numpy.dot(CmsumG,Q_sumG)
#infHsump=numpy.dot(CmsumH,Q_sumH)
#GHavg=numpy.array([infGsump.flat[0],infGsump.flat[1],infHsump.flat[0],infHsump.flat[1]])
#print GHavg

GHavg=GHsum/i
for g in x:
    Q+=numpy.dot(g-GHavg,(g-GHavg).transpose())
Q/=(i-1.) 


print Q




    
print i, " lines read."




fG.close()
fH.close()



