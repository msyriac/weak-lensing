#!/usr/bin/env python

import glob, sys
import csv
import numpy
from math import *
import matplotlib.pyplot as p
import scipy.linalg  as l
import pickle

#files = glob.glob("data/g_zero_hightol/*.csv")
files = glob.glob("data/run_working_5/*.csv")

Cinv=numpy.zeros((2,2))
F=numpy.zeros((2,2))
Qsum=0.
appG=[0.06,-0.07]

x=[]
y=[]
y2=[]
i=0
k=0

freq=1000

ped=numpy.zeros(2)
pedsw=0

QNl=[]

for fname in files[:]:
    f = open(fname, 'r')
    reader = csv.reader(f,delimiter=',')
    for row in reader:
        try:
            vals = [float(r) for r in (rv for rv in row if rv!='')]
        except:
            continue
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
            F+=(numpy.dot(Q_norm,Q_norm.transpose()))

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

#F/=i
#print F
#Fr=numpy.zeros((2,2))
#Fr[0,0]=9.782
#Fr[1,1]=9.782
#Fr[0,0]=(F[0,0]+F[1,1])/2.
#Fr[1,1]=(F[0,0]+F[1,1])/2.
#print Fr
#print numpy.dot(numpy.linalg.inv(Fr),Qsum/i)

F=pickle.load(open("F.pickle",'r'))
#F[0,0]=(F[0,0]+F[1,1])/2.
#F[1,1]=F[0,0]
#F[0,1]=0.
#F[1,0]=0.
print F

infGp=numpy.dot(numpy.linalg.inv(F),Qsum/i)
print infGp
#sys.exit()
Cm=numpy.linalg.inv(Cinv)
infGpalt=numpy.dot(Cm,Qsum)
infG=[infGp.flat[0],infGp.flat[1]]
infGalt=[infGpalt.flat[0],infGpalt.flat[1]]
#err=[sqrt((Cm[0,0])),sqrt((Cm[1,1]))]
#sigma = [(infG[0]-appG[0])/err[0],(infG[1]-appG[1])/err[1]]
#err=[0.,0.]
#sigma=[0.,0.]
bias=[infG[0]-appG[0],infG[1]-appG[1]]
biasalt=[infGalt[0]-appG[0],infGalt[1]-appG[1]]
print "Bias with F estimator is ", [-b for b in bias]
print "Bias is otherwise ", [-b for b in biasalt]
pickle.dump([-b for b in bias],open("g06m05bias.pickle",'w'))
sys.exit()

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
