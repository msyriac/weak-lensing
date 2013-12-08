#!/usr/bin/env python

import glob
import csv
import numpy
from math import *
import matplotlib.pyplot as p
import scipy.linalg as l
import generate_pairs as gen ###
import matplotlib.pyplot as pl

def LnLikelihood(Qlist,Nlist,C):
    LnLsum=0
    for Q,N in zip(Qlist,Nlist):
        Ninv=numpy.linalg.inv(N)
        NinvC=numpy.dot(Ninv,C)
        F=-0.5*NinvC.trace()
        S=0.25*numpy.dot(NinvC,NinvC).trace()
        QC=numpy.dot(Q,C)
        T=0.5*numpy.dot(QC,Q.transpose()) #check order!
        T2=-0.5*numpy.dot(numpy.dot(QC,NinvC),Q.transpose()) #check order!
        LnL=  F + S + T + T2
        LnLsum+=LnL
    return LnLsum[0][0]



    
root="/astro/u/anze/work/mat/wl/datan"
#root="data/pairs_1"

filesG = glob.glob(root+"/g*.csv")
filesH = glob.glob(root+"/h*.csv")

CinvG=numpy.zeros((2,2))
CinvH=numpy.zeros((2,2))
Cinv_sumG=numpy.zeros((2,2))
Cinv_sumH=numpy.zeros((2,2))
GHsum=numpy.zeros((4,1))
Q=numpy.zeros((4,4))
Q_sumG=numpy.zeros((2,1))
Q_sumH=numpy.zeros((2,1))

#p=gen.Pairs() ###


i=0

freq=1000

x=[]
Q=[]
N=[]

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
        #infGp=numpy.dot(CmG,Q_normG)



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
        #infHp=numpy.dot(CmH,Q_normH)



        Q.append(numpy.array([[Q_normG.flat[0],
                               Q_normG.flat[1],
                               Q_normH.flat[0],
                               Q_normH.flat[1]]]))
        N.append(l.block_diag(CmG,CmH))
        
        #print Q,N

        #raise SystemExit(0)
        #x.append( numpy.array( [[infGp.flat[0]],
        #[infGp.flat[1]],
        #                        [infHp.flat[0]],
        #                        [infHp.flat[1]] ] )  )
        
        #x.append(p.draw()) ###
        #x[-1].resize(4,1) ###

        #GHsum+=x[-1]
        
        





print i, " lines read."

p=gen.Pairs()
trueC=p.cov
xtrue=trueC[0,3]
npt=5
step=xtrue*0.1
x=numpy.arange(xtrue-step*npt/2,xtrue+step*npt/2,step)
y=[]
for xi in x:
    trueC[0,3]=xi
    y.append(LnLikelihood(Q,N,trueC))


print x,y
pl.plot(x.tolist(),y)
pl.show()


fG.close()
fH.close()





