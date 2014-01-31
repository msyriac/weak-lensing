import random
from math import *
from scipy.integrate import dblquad
import matplotlib.pyplot as pl
import numpy as np

import mcint
import random
import sys, getopt
#from numba import autojit
import time

def sampler():
    while True:
        y1=2.*random.random()-1.
        x1=2.*random.random()-1.
        #y2=2.*random.random()-1.
        #x2=2.*random.random()-1.
        if ((x1**2.+y1**2. <= 1)):# and (x2**2.+y2**2. <= 1)):
            yield x1,y1 #,x2,y2)


def MonteCarloEval(g1,g2,es1,es2,ids,nsamples):
    #result, error=mcint.integrate(TestIntegrand,sampler(),measure=pi,n=1000000)
    global sm, es1_,es2_,g1_,g2_

    g1_=g1
    g2_=g2
    es1_=es1
    es2_=es2

    if ids=='11':
        result, error=mcint.integrate(mfisherIntegrand11,sampler(),measure=pi,n=nsamples)
    elif ids=='22':
        result, error=mcint.integrate(mfisherIntegrand22,sampler(),measure=pi,n=nsamples)
    elif ids=='12':
        result, error=mcint.integrate(mfisherIntegrand12,sampler(),measure=pi,n=nsamples)
    else:
        print "Invalid integrand ID. Exiting."
        sys.exit()
    #result, error=mcint.integrate(TestIntegrand,sampler(),measure=1.0,n=1000)

    #print result, error
    return result,error

def TestIntegrand((x1,y1)):
    if ((x1**2.+y1**2.)<=1.0):
        return pi #(x1**2.+y1**2.)
    else:
        return 0.

def main(argv):
    global sm, es1_,es2_,g1_,g2_

    g1_=0.00
    g2_=0.00
    sm=0.05

    ntheta=50
    nsamples = 50000
    i=20
    v=1


    try:
        opts, args = getopt.getopt(argv,"n:1:2:e:t:i:v:")
    except getopt.GetoptError:
        print "I don't understand the arguments you passed. Run with -h to see available options."
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'python cluster_fisher.py -n <number of MC samples=50000> -1 <shear component 1=0.0> -2 <shear component 2=0.00> -e <sigma_e=0.05> -t <number of theta samples=50> -i <file index=20, E_S=i*0.005> -v <1 for F11, 2 for F22, 3 for F12>'
            sys.exit()
        elif opt == "-n":
            nsamples=int(arg)
        elif opt == "-1":
            g1_=float(arg)
        elif opt == "-2":
            g2_=float(arg)
        elif opt == "-e":
            sm=float(arg)
        elif opt == "-t":
            ntheta=int(arg)
        elif opt == "-i":
            i=int(arg)
        elif opt == "-v":
            v=int(arg)



    fH=open('cluster_fisher_opts.csv', 'w')
    row=[nsamples,ntheta,g1_,g2_,i]
    fH.write(','.join(str(j) for j in row) + '\n')
    fH.close()

    #i goes from 0 to 80 so that E goes from 0 to 0.4
    E=i*0.005

    if v==1:
        ids='11'
        pref='va/'
    elif v==2:
        ids='22'
        pref='vb/'
    elif v==3:
        ids='12'
        pref='vc/'
    else:
        print "Invalid v option."
        sys.exit()

    F,eF=avgTheta(E,ntheta,ids,nsamples)
    #F22,eF22=avgTheta(E,ntheta,'22',nsamples)
    #F12,eF12=avgTheta(E,ntheta,'12',nsamples)
    row=[F,eF]

    fH=open('data_fisher/'+pref+'i'+str(i)+'.csv', 'w')
    fH.write(','.join(str(j) for j in row) + '\n')
    fH.close()



def avgTheta(E,n,ids,nsamples):
    Integ=0.
    err=0.
    for i in range(0,n):

        theta=random.uniform(0,2.*pi);

        es1_=E*cos(theta)
        es2_=E*sin(theta)
        val,erro=MonteCarloEval(g1_,g2_,es1_,es2_,ids,nsamples)
        Integ+=val
        err+=erro

    return Integ/n,err/n

def initFisher(s):
    global sm
    sm=s

def LogLike(es1,es2,g1,g2,em1,em2):
    global sm
    Ans = -(es1**2. + es2**2. - 2.*es1*g1 + g1**2. - 2.*es2*g2 + g2**2. + 2.*em2*(es2**2.*g2 + (1. + es1**2. - 2.*es1*g1)*g2 + es2*(-1. + g1**2. - g2**2.)) + 2.*em1*(es1**2.*g1 + g1*(1. + es2**2. - 2.*es2*g2) + es1*(-1. - g1**2. + g2**2.)) + em1**2.*(1. - 2.*es1*g1 - 2.*es2*g2 + es1**2.*(g1**2. + g2**2.) + es2**2.*(g1**2. + g2**2.)) + em2**2.*(1. - 2.*es1*g1 - 2.*es2*g2 + es1**2.*(g1**2. + g2**2.) + es2**2.*(g1**2. + g2**2.)))/ (2.*(1. - 2.*es1*g1 - 2.*es2*g2 + es1**2.*(g1**2. + g2**2.) + es2**2.*(g1**2. + g2**2.))*sm**2.)
    return Ans

def DLogLikeg1g2(es1,es2,g1,g2,em1,em2):
    global sm
    Ans=-((((-1 + es1**2 - 2*em2*es2 + es2**2 + 2*em2*(es1**2 + es2**2)*g2)*(g1 + es1**2*g1 + es2*g1*(es2 - 2*g2) + es1*(-1 - g1**2 + g2**2)) +em1*(es1**4*(g1 - g2)*(g1 + g2) + 2*es1**3*g1*(-1 + 2*g2**2) +(1 + es2*(g1 - g2))*(1 + es2**2 - 2*es2*g2)*(-1 + es2*(g1 + g2)) + 2*es1*g1*(1 - 2*es2*g2 + es2**2*(-1 + 2*g2**2)) +es1**2*(1 + g1**2*(-1 + 2*es2*(es2 - g2)) + g2*(-3*g2 + 2*es2*(1 - es2*g2 + g2**2)))))*(-((-1 - 2*em1*es1 + es1**2 + es2**2 + 2*em1*(es1**2 + es2**2)*g1)*(es2*(-1 + g1**2) + (1 + es1**2 + es2**2 - 2*es1*g1)*g2 - es2*g2**2)) +em2*(1 + es1**4*(g1 - g2)*(g1 + g2) - 2*es1**3*(g1 + g1**3 - g1*g2**2) -2*es1*g1*(2 - 2*es2*g2 + es2**2*(1 + g1**2 - g2**2)) +es1**2*(1 + g1**2*(5 + 2*es2*(es2 - 2*g2)) - g2*(g2 + 2*es2*(-1 + es2*g2))) +es2*(-2*g2 + es2*(-1 + g1**2*(3 + es2**2 - 4*es2*g2) + g2*(g2 + es2*(2 - es2*g2)))))))/((1 - 2*es1*g1 - 2*es2*g2 + es1**2*(g1**2 + g2**2) + es2**2*(g1**2 + g2**2))**4*sm**4))
    return Ans

def DLogLikeg1g1(es1,es2,g1,g2,em1,em2):
    global sm
    Ans=((-1. + es1**2. - 2.*em2*es2 + es2**2. + 2.*em2*(es1**2. + es2**2.)*g2)*(g1 + es1**2.*g1 + es2*g1*(es2 - 2.*g2) + es1*(-1. - g1**2 + g2**2)) +em1*(es1**4.*(g1 - g2)*(g1 + g2) + 2.*es1**3.*g1*(-1. + 2.*g2**2.) +(1. + es2*(g1 - g2))*(1. + es2**2. - 2.*es2*g2)*(-1. + es2*(g1 + g2)) + 2.*es1*g1*(1. - 2.*es2*g2 + es2**2.*(-1. + 2.*g2**2.)) +es1**2.*(1. + g1**2.*(-1. + 2.*es2*(es2 - g2)) + g2*(-3.*g2 + 2.*es2*(1. - es2*g2 + g2**2.)))))**2./((1. - 2.*es1*g1 - 2.*es2*g2 + es1**2.*(g1**2. + g2**2.) + es2**2.*(g1**2. + g2**2.))**4.*sm**4.)
    return Ans


def DLogLikeg2g2(es1,es2,g1,g2,em1,em2):
    global sm
    Ans=((-1 - 2*em1*es1 + es1**2 + es2**2 + 2*em1*(es1**2 + es2**2)*g1)*(es2*(-1 + g1**2) + (1 + es1**2 + es2**2 - 2*es1*g1)*g2 - es2*g2**2) +em2*(-1 + es1**4*(-g1**2 + g2**2) + 2*es1**3*(g1 + g1**3 - g1*g2**2) + 2*es1*g1*(2 - 2*es2*g2 + es2**2*(1 + g1**2 - g2**2)) +es1**2*(-1 + g1**2*(-5 - 2*es2**2 + 4*es2*g2) + g2*(g2 + 2*es2*(-1 + es2*g2))) +es2*(2*g2 - es2*(-1 + g1**2*(3 + es2**2 - 4*es2*g2) + g2*(g2 + es2*(2 - es2*g2))))))**2/((1 - 2*es1*g1 - 2*es2*g2 + es1**2*(g1**2 + g2**2) + es2**2*(g1**2 + g2**2))**4*sm**4)
    return Ans

def DDLogLikeg1g1(es1,es2,g1,g2,em1,em2):
    global sm

    Ans=(8*(-es1 + es1**2*g1 + es2**2*g1)*(em2**2*es1**2*g1 + (1 + em2*es2)**2*g1 + em1**2*(-es1 + es1**2*g1 + es2**2*g1) -es1*(1 + em2**2 + 2*em2*g2) + em1*(1 + es1**2 + es2**2 - 2*es1*g1 - 2*es2*g2))* (1 - 2*es1*g1 - 2*es2*g2 + es1**2*(g1**2 + g2**2) + es2**2*(g1**2 + g2**2)) - 2*(1 - 2*em1*es1 + 2*em2*es2 + em1**2*(es1**2 + es2**2) + em2**2*(es1**2 + es2**2))*(1 - 2*es1*g1 - 2*es2*g2 + es1**2*(g1**2 + g2**2) + es2**2*(g1**2 + g2**2))**2 -8*(-es1 + es1**2*g1 + es2**2*g1)**2*(es1**2 + es2**2 - 2*es1*g1 + g1**2 - 2*es2*g2 + g2**2 +2*em2*(es2**2*g2 + (1 + es1**2 - 2*es1*g1)*g2 + es2*(-1 + g1**2 - g2**2)) + 2*em1*(es1**2*g1 + g1*(1 + es2**2 - 2*es2*g2) + es1*(-1 - g1**2 + g2**2)) + em1**2*(1 - 2*es1*g1 - 2*es2*g2 + es1**2*(g1**2 + g2**2) + es2**2*(g1**2 + g2**2)) + em2**2*(1 - 2*es1*g1 - 2*es2*g2 + es1**2*(g1**2 + g2**2) + es2**2*(g1**2 + g2**2))) + 2*(es1**2 + es2**2)*(1 - 2*es1*g1 - 2*es2*g2 + es1**2*(g1**2 + g2**2) + es2**2*(g1**2 + g2**2))* (es1**2 + es2**2 - 2*es1*g1 + g1**2 - 2*es2*g2 + g2**2 + 2*em2*(es2**2*g2 + (1 + es1**2 - 2*es1*g1)*g2 + es2*(-1 + g1**2 - g2**2)) + 2*em1*(es1**2*g1 + g1*(1 + es2**2 - 2*es2*g2) + es1*(-1 - g1**2 + g2**2)) + em1**2*(1 - 2*es1*g1 - 2*es2*g2 + es1**2*(g1**2 + g2**2) + es2**2*(g1**2 + g2**2)) + em2**2*(1 - 2*es1*g1 - 2*es2*g2 + es1**2*(g1**2 + g2**2) + es2**2*(g1**2 + g2**2))))/ (2.*(1 - 2*es1*g1 - 2*es2*g2 + es1**2*(g1**2 + g2**2) + es2**2*(g1**2 + g2**2))**3*sm**2)
    return -Ans


def fisherIntegrand12(em1,em2):
    global es1_,es2_,g1_,g2_
    return (DLogLikeg1g2(es1_,es2_,g1_,g2_,em1,em2)*exp(LogLike(es1_,es2_,g1_,g2_,em1,em2)))
def mfisherIntegrand12((em1,em2)):
    global es1_,es2_,g1_,g2_
    return (DLogLikeg1g2(es1_,es2_,g1_,g2_,em1,em2)*exp(LogLike(es1_,es2_,g1_,g2_,em1,em2)))
def fisherIntegrand11(em1,em2):
    global es1_,es2_,g1_,g2_
    return (DLogLikeg1g1(es1_,es2_,g1_,g2_,em1,em2)*exp(LogLike(es1_,es2_,g1_,g2_,em1,em2)))
def mfisherIntegrand11((em1,em2)):
    global es1_,es2_,g1_,g2_
    return (DLogLikeg1g1(es1_,es2_,g1_,g2_,em1,em2)*exp(LogLike(es1_,es2_,g1_,g2_,em1,em2)))
def mfisherIntegrand22((em1,em2)):
    global es1_,es2_,g1_,g2_
    return (DLogLikeg2g2(es1_,es2_,g1_,g2_,em1,em2)*exp(LogLike(es1_,es2_,g1_,g2_,em1,em2)))
def fisherIntegrand22(em1,em2):
    global es1_,es2_,g1_,g2_
    return (DLogLikeg2g2(es1_,es2_,g1_,g2_,em1,em2)*exp(LogLike(es1_,es2_,g1_,g2_,em1,em2)))
def AltfisherIntegrand11(em1,em2):
    global es1_,es2_,g1_,g2_
    return ((DDLogLikeg1g1(es1_,es2_,g1_,g2_,em1,em2)-DLogLikeg1g1(es1_,es2_,g1_,g2_,em1,em2))*exp(LogLike(es1_,es2_,g1_,g2_,em1,em2)))
def mAltfisherIntegrand11((em1,em2)):
    global es1_,es2_,g1_,g2_
    return ((DDLogLikeg1g1(es1_,es2_,g1_,g2_,em1,em2)-DLogLikeg1g1(es1_,es2_,g1_,g2_,em1,em2))*exp(LogLike(es1_,es2_,g1_,g2_,em1,em2)))
#def mfisherIntegrand11((em1,em2)):
#    global es1_,es2_,g1_,g2_
#    return (-1.*DDLogLikeg1g1(es1_,es2_,g1_,g2_,em1,em2)*exp(LogLike(es1_,es2_,g1_,g2_,em1,em2)))
def dfisherIntegrand11(em1,em2):
    global es1_,es2_,g1_,g2_
    return (-1.*DDLogLikeg1g1(es1_,es2_,g1_,g2_,em1,em2)*exp(LogLike(es1_,es2_,g1_,g2_,em1,em2)))


if __name__=="__main__":
    main(sys.argv[1:])
