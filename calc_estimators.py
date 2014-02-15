#!/usr/bin/env python

import glob, sys
import csv
import numpy
from math import *
import matplotlib.pyplot as p
import scipy.linalg  as linalg
import pickle

from itertools import product

# Tell the program where you have stored P Q R values in csv
# Each row contains the respective derivative (zeroth, first, second)
# of the full marginalized likelihood evaluated for ZERO shear
# and a value of e_m drawn from that likelihood for ZERO shear using toy.py. 
# (In practice, draw e_s from intrinsic ellipticity prior distribution,
# combine with (zero) shear and add gaussian error to it.)

# Which means that if you have some quantity, F, G, H that you want to calculate
# that is an expectation value (at zero shear) for some product of derivatives
# of likelihoods calculated at zero shear, just find that quantity for each
# for each set of P,Q,R values and take an average (sample mean).

#files = glob.glob("data/g_zero_hightol/*.csv")
#files1 = glob.glob("data/g_zero_good/*.csv")
#files2 = glob.glob("data/g_zero_good2/*.csv")
#files = files1 #+ files2
files = glob.glob("data/*.csv")

F=numpy.zeros((2,2))
Qsum=0.
#appG=[0.0,0.0]

i=0
k=0

freq=10000
#Gijk
G=numpy.zeros((2,2,2))
#(Gijk + Hijk)
N=numpy.zeros((2,2,2))

b=[0,1]
s=[b,b,b]
Nch=30

RP=numpy.zeros((2,2)) 

Fc=[numpy.zeros((2,2)) for a in range(Nch)]
Rc=[numpy.zeros((2,2)) for a in range(Nch)]
Gc=[numpy.zeros((2,2,2)) for a in range(Nch)]
Nc=[numpy.zeros((2,2,2)) for a in range(Nch)]
sw=numpy.zeros(Nch)


for fname in files:
    f = open(fname, 'r')
    reader = csv.reader(f,delimiter=',')
    for row in reader:
        try:
            vals = [float(r) for r in (rv for rv in row if rv!='')]
        except:
            print ":("
            continue
        if len(vals)!=6:
            print ":( Invalid line, ", row
            print fname
            k+=1
        else:
            P=vals[0]
            Q=numpy.array([[vals[1]],
                          [vals[2]]])
            R=numpy.array( [ [vals[3], vals[4]],
                             [vals[4], vals[5]] ] )

            Q_norm=(Q/P)
            #Qsum+=Q_norm

            Find=numpy.outer(Q_norm,Q_norm)-(R/P)
            #F += Find #(numpy.dot(Q_norm,Q_norm.transpose()))


            #(G+H)
            P2=P*P
            P3=P2*P
            j=i%Nch
            
            for l,m,n in product(*s):
                Nind = Q[l]*R[m,n]/P2
                Gind = Q[l]*Q[m]*Q[n]/P3
                #N[l,m,n] += Nind
                #G[l,m,n] += Gind
                Nc[j][l,m,n]+=Nind
                Gc[j][l,m,n]+=Gind
            

            RP+=(R/P)
            Fc[j]+=Find
            Rc[j]+=(R/P)
            sw[j]+=1

            i+=1

            if (i % freq)==0:
                print i

#N/=i
#G/=i
#F/=i
#print F
#Fr=numpy.zeros((2,2))
#Fr[0,0]=(F[0,0]+F[1,1])/2.
#Fr[1,1]=(F[0,0]+F[1,1])/2.

print RP/i
#print Fr
#print numpy.dot(numpy.linalg.inv(Fr),Qsum/i)


print i, " valid lines read."
print k, " invalid lines skipped."

f.close()


i=0


S=numpy.zeros((2,2))
SS=numpy.zeros((2,2))

RS=numpy.zeros((2,2))
RSS=numpy.zeros((2,2))

GS=numpy.zeros((2,2,2))
GSS=numpy.zeros((2,2,2))

NS=numpy.zeros((2,2,2))
NSS=numpy.zeros((2,2,2))

for j in range(Nch):
    Fc[j]/=(sw[j]-1.) # was sw[j] !!!
    Rc[j]/=(sw[j]-1.) # was sw[j] !!!
    for l,m,n in product(*s):
        Nc[j][l,m,n]/=(sw[j]-1.)
        Gc[j][l,m,n]/=(sw[j]-1.)
        NS[l,m,n]+=Nc[j][l,m,n]
        NSS[l,m,n]+=Nc[j][l,m,n]*Nc[j][l,m,n]
        GS[l,m,n]+=Gc[j][l,m,n]
        GSS[l,m,n]+=Gc[j][l,m,n]*Gc[j][l,m,n]


    S+=Fc[j]
    SS+=Fc[j]*Fc[j]
    RS+=Rc[j]
    RSS+=Rc[j]*Rc[j]

S/=Nch
SS/=Nch
SS-=S*S

RS/=Nch
RSS/=Nch
RSS-=RS*RS

for l,m,n in product(*s):
    GS[l,m,n]/=Nch
    GSS[l,m,n]/=Nch
    GSS[l,m,n]-=GS[l,m,n]*GS[l,m,n]

    NS[l,m,n]/=Nch
    NSS[l,m,n]/=Nch
    NSS[l,m,n]-=NS[l,m,n]*NS[l,m,n]


print '--------'
print "F matrix"
print S
pickle.dump(S,open('F.pickle','w'))
print "F matrix stddev"
print numpy.sqrt(SS)/sqrt(Nch-1.) # was Nch !!!

print "G matrix"
print GS
pickle.dump(GS,open('G.pickle','w'))
print "G matrix stddev"
print numpy.sqrt(GSS)/sqrt(Nch-1.)


print "N matrix"
print NS
pickle.dump(NS,open('N.pickle','w'))
print "N matrix stddev"
print numpy.sqrt(NSS)/sqrt(Nch-1.)

print "R/P matrix"
print RS
pickle.dump(RS,open('RP.pickle','w'))
print "R/P matrix stddev"
print numpy.sqrt(RSS)/sqrt(Nch-1.)

#print numpy.dot(linalg.inv(S),RS)
