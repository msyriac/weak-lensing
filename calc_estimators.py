#!/usr/bin/env python

import glob, sys
import csv
import numpy
from math import *
import matplotlib.pyplot as p
import scipy.linalg  as linalg
import pickle
from accumulator import *
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
files = glob.glob("data/g_zero_good1/*.csv")

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

Fc=Accumulator("Fisher matrix acc.",Nch)
Oc=Accumulator("Other quantity", Nch)
Oc2=Accumulator("Other quantity2", Nch)

Rc=Accumulator("Rc matrix acc.",Nch)
Gc=Accumulator("G matrix acc.",Nch)
Nc=Accumulator("N matrix acc.",Nch)

Nind=zeros((2,2,2))
Gind=zeros((2,2,2))

## this is the average of R/P for uncorrected quantity
Rcorrector=-array([[0.697824,0],[0,0.697824]])

for fname in files[:]:
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
            ##correction due to shitty normalization
            R+=Rcorrector*P

            L_i=(Q/P)
            L_iL_iT=numpy.outer(L_i,L_i)
            L_ij=(R/P)-L_iL_iT

            Find=-L_ij
            Fc.accumulate(-L_ij)
            #Oc.accumulate(R/P)
            #Oc2.accumulate(Q/P)
            Rc.accumulate(R/P)

            #if False
            if True:
                P2=P*P
                P3=P2*P
                for l,m,n in product(*s):
                    Nind[l,m,n] = Q[l]*R[m,n]/P2
                    Gind[l,m,n] = Q[l]*Q[m]*Q[n]/P3
                Nc.accumulate(Nind)
                Gc.accumulate(Gind)

            i+=1
            if (i % freq)==0:
                print i

print i, " valid lines read."
print k, " invalid lines skipped."

f.close()

for mat in [Fc,Oc,Oc2, Rc,Gc,Nc]:
    mat.print_stat()
