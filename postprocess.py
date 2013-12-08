#!/usr/bin/env python

import glob
import csv
import numpy
from math import *
import matplotlib.pyplot as p
import scipy.linalg  as l

files = glob.glob("run_working_1/*.csv")

Cinv=numpy.zeros((2,2))
Qsum=0.
appG=[-0.01,0.02]

x=[]
y=[]
y2=[]

i=0
k=0

freq=1000

ped=numpy.zeros(2)
pedsw=0

QNl=[]

for fname in files:
    f = open(fname, 'r')
    reader = csv.reader(f,delimiter=',')
    for row in reader:
        vals = [float(r) for r in (rv for rv in row if rv!='')]
        if len(vals)!=6:
            print "Invalid line, ", row
            print fname
            k+=1
        else:
            i+=1
            P=vals[0]
            Q=numpy.array([[vals[1]],
                          [vals[2]]])
            R=numpy.array( [ [vals[3], vals[4]],
                             [vals[4], vals[5]] ] )

            Q_norm=(Q/P)
            QNl.append(Q_norm.flat)
            Qsum+=Q_norm
            CinvC=((numpy.dot(Q_norm,Q_norm.transpose()))-(R/P))
            Cinv+=CinvC
            if False:
                detCinv=CinvC[0,0]*CinvC[1,1]-CinvC[0,1]**2
                CC=l.inv(CinvC)
                infGX=numpy.dot(CC,Q_norm.flat)
                we=detCinv
                ped+=infGX*we
                pedsw+=we
                
                #we=abs(detCinv)
                #if (detCinv>0):
                #    ped+=infGX*we
                #    pedsw+=we
                ##else:
                #    ped-=infGX*we
                #    pedsw+=we
            
            

            if (i % freq)==0:
                print i
Cm=numpy.linalg.inv(Cinv)
infGp=numpy.dot(Cm,Qsum)
infG=[infGp.flat[0],infGp.flat[1]]
err=[sqrt((Cm[0,0])),sqrt((Cm[1,1]))]
sigma = [(infG[0]-appG[0])/err[0],(infG[1]-appG[1])/err[1]]
#err=[0.,0.]
#sigma=[0.,0.]
bias=[infG[0]-appG[0],infG[1]-appG[1]]
biasp=[(bias[0]*100/appG[0]),(bias[1]*100/appG[1])]
x.append(i)
y.append(biasp[0])
y2.append(biasp[1])

print i, " valid lines read."
print k, " invalid lines skipped."
print i
print infG[0], err[0],  sigma[0],  biasp[0]
print infG[1], err[1],  sigma[1],  biasp[1]
print Cm

Cmo=Cm*len(QNl)
Cmo[1,0]=0
Cmo[0,1]=0
Cmo[1,1]=(Cmo[0,0]+Cmo[1,1])/2
Cmo[0,0]=Cmo[1,1]

for Qn in QNl:
    ped+=numpy.dot(Cmo,Qn)
    pedsw+=1

print Cmo/len(QNl)

ped/=pedsw
print ped

if (freq!=i):
    p.clf()
    p.plot(x,y)
    p.plot(x,y2)
    p.show()
f.close()
