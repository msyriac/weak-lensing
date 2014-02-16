#!/usr/bin/env python

import glob
import csv
import numpy
from math import *
import matplotlib.pyplot as p
import scipy.linalg as l
from scipy import *

import cPickle as pickle
import sys


#import generate_pairs as gen

root='/astro/u/anze/work/mat/wl/'

filesG = glob.glob(root+"datan/g*.csv")
filesH = glob.glob(root+"datan/h*.csv")


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

freq=10000

x=[]
w=[]

Pickled=True ###################### change to True after first run on new data to save time
PickleIt=True #change to false if you don't want to pickle a first run (>1.5x speedup)
pickleroot='' #change if you want to pickle a different set!

if (not Pickled):
    print "WARNING: Pickling will take place since it is enabled and you have indicated this is a first run. Pickling results in a >1.5x slowdown on first run but a significant speedup in subsequent runs, so if you're really in a hurry now, open the source and disable pickling."
    print "Reading for first time..."
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
    
            CinvG=((numpy.dot(Q_normG,Q_normG.transpose()))-(RG/PG))
            Cinv_sumG+=CinvG 
    
            PH=valsH[0]
            QH=numpy.array([[valsH[1]],
                            [valsH[2]]])
            RH=numpy.array( [ [valsH[3], valsH[4]],
                              [valsH[4], valsH[5]] ] )
    
            Q_normH=(QH/PH)
            CinvH=((numpy.dot(Q_normH,Q_normH.transpose()))-(RH/PH))
            Cinv_sumH+=CinvH
            x.append( (Q_normG, Q_normH) ) #removed .flat
            
    
            
        fG.close()
        fH.close()
            
    if PickleIt:
        print "Pickling. This will take some time..."
        pickle.dump(x,open(pickleroot+'x.pickle','wb'))
        pickle.dump(Cinv_sumG,open(pickleroot+'Cinv_sumG.pickle','wb'))
        pickle.dump(Cinv_sumH,open(pickleroot+'Cinv_sumH.pickle','wb'))
        pickle.dump(i,open(pickleroot+'i.pickle','wb'))


if Pickled:
    print "Unpickling. This will take some time... but not as much as if you hadn't pickled! (If you've generated new data files, Ctrl+C and change Pickled to False in source.)"
    x=pickle.load(open(pickleroot+'x.pickle','rb'))
    Cinv_sumG=pickle.load(open(pickleroot+'Cinv_sumG.pickle','rb'))
    Cinv_sumH=pickle.load(open(pickleroot+'Cinv_sumH.pickle','rb'))
    i=pickle.load(open(pickleroot+'i.pickle','rb'))

Cinv=Cinv_sumG+Cinv_sumH
CC=l.inv(Cinv)*2*i ######## should this be i or 2i?
Cd=(CC[0,0]+CC[1,1])/2.0
CC=array([[Cd,0],[0,Cd]])

##
Cg=l.inv(Cinv_sumG)*i
Ch=l.inv(Cinv_sumH)*i
##

print Cinv_sumG
print Cinv_sumH
print CC

Full=True

Nch=30

if Full:
    Q=[numpy.zeros((4,4)) for a in range(Nch)]
else:
    Q=[numpy.zeros((2,2)) for a in range(Nch)]
sw=zeros(Nch)
i=0

for QG,QH in x:
    j=i%Nch
    if Full:
        #gest=dot(CC,QG)
        #hest=dot(CC,QH)
        gest=dot(Cg,QG)
        hest=dot(Ch,QH)
        joint=numpy.array([[gest[0]],
                           [gest[1]],
                           [hest[0]],
                           [hest[1]]])

        M=numpy.outer(joint,joint.transpose())

    else:
        M=numpy.outer(dot(CC,QG),dot(CC,QH))

    Q[j]+=M
    sw[j]+=1
    i+=1


if Full:
    S=zeros((4,4))
    SS=zeros((4,4))
else:
    S=zeros((2,2))
    SS=zeros((2,2))

for j in range(Nch):
    Q[j]/=(sw[j]-1.) # was sw[j] !!!
    print Q[j]
    S+=Q[j]
    SS+=Q[j]*Q[j]
    

S/=Nch
SS/=Nch
SS-=S*S
print '--------'
print S
print sqrt(SS)/sqrt(Nch-1.) # was Nch !!!
print CC[0,0]/sqrt(i)







