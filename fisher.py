import random
from math import *
from scipy.integrate import dblquad
import matplotlib.pyplot as pl
import numpy as np

import mcint
import random
import sys
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

def main():
    global sm, es1_,es2_,g1_,g2_
    tol=1.49e-2

    g1_=0.00
    g2_=0.00
    sm=0.05 #

    #MonteCarloEval(g1_,g2_,0.3,-0.3)
    #sys.exit()

    Erange=np.arange(0.0,0.4,0.01)
    lF11=[]
    lF12=[]
    lF22=[]
    leF11=[]
    leF12=[]
    leF22=[]
    i=1
    l=len(Erange)
    start=time.time()

    nsamples = 50000
    pl.clf()

    for E in Erange:

        #E=0.3

        F11,eF11=avgTheta(E,50,'11',nsamples)
        F22,eF22=avgTheta(E,50,'22',nsamples)
        #F11, eF11=dblquad(dfisherIntegrand11,-1.0, 1.0, lambda x:-sqrt(1.-x*x), lambda x:sqrt(1.-x*x),epsabs=tol, epsrel=tol)
        #F11, eF11 = MonteCarloEval(g1_,g2_,es1_,es2_)
        #F11, eF11=dblquad(fisherIntegrand11,-1.0, 1.0, lambda x:-sqrt(1.-x*x), lambda x:sqrt(1.-x*x),epsabs=tol, epsrel=tol)
        #F12, eF12=dblquad(fisherIntegrand12,-1.0, 1.0, lambda x:-sqrt(1.-x*x), lambda x:sqrt(1.-x*x),epsabs=tol, epsrel=tol)
        #F22, eF22=dblquad(fisherIntegrand22,-1.0, 1.0, lambda x:-sqrt(1.-x*x), lambda x:sqrt(1.-x*x),epsabs=tol, epsrel=tol)
        lF11.append(F11)
        lF22.append(F22)
        #lF12.append(F12)
        leF11.append(eF11)
        leF22.append(eF22)
        #leF12.append(eF12)
        elapsed=time.time()-start
        avgtime=elapsed/i
        timeleft=avgtime*(l-i)
        print "Done ", i*100/l, " %"
        print "Time left is ", timeleft/60., " minutes."
        #print F11,F12,F22
        i+=1

    #pl.plot(Erange,lF11)
    # pl.plot(Erange,lF12)
    # pl.plot(Erange,lF22)
    pl.errorbar(Erange,lF11,yerr=leF11)
    #pl.errorbar(Erange,lF12,yerr=leF12)
    pl.errorbar(Erange,lF22,yerr=leF22)
    pl.savefig('plto.png')
    #pl.show()
    #print LogLike(es1_,es2_,g1_,g2_,em1,em2)
    #print DLogLikeg1g2(es1_,es2_,g1_,g2_,em1,em2)
    #print fisherIntegrand(em1,em2)

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
    main()

