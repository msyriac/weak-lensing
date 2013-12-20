import random
from math import *
from scipy.integrate import dblquad 
import matplotlib.pyplot as pl
import numpy as np


def main():
    global sm, es1_,es2_,g1_,g2_
    tol=1.49e-06

    g1_=-0.01
    g2_=0.02
    sm=0.05

    Erange=np.arange(0.0,0.4,0.01)
    lF11=[]
    lF12=[]
    lF22=[]
    i=1
    l=len(Erange)
    for E in Erange:
        
        #E=0.3
        theta=random.uniform(0,2.*pi);

        es1_=E*cos(theta)
        es2_=E*sin(theta)

        F11=dblquad(fisherIntegrand11,-1.0, 1.0, lambda x:-sqrt(1.-x*x), lambda x:sqrt(1.-x*x),epsabs=tol, epsrel=tol)[0]
        F12=dblquad(fisherIntegrand12,-1.0, 1.0, lambda x:-sqrt(1.-x*x), lambda x:sqrt(1.-x*x),epsabs=tol, epsrel=tol)[0]
        F22=dblquad(fisherIntegrand22,-1.0, 1.0, lambda x:-sqrt(1.-x*x), lambda x:sqrt(1.-x*x),epsabs=tol, epsrel=tol)[0]
        lF11.append(F11)
        lF22.append(F22)
        lF12.append(F12)
        print "Done ", i*100/l, " %"
        #print F11,F12,F22
        i+=1

    pl.plot(Erange,lF11)
    pl.plot(Erange,lF12)
    pl.plot(Erange,lF22)
    pl.show()
    #print LogLike(es1_,es2_,g1_,g2_,em1,em2)
    #print DLogLikeg1g2(es1_,es2_,g1_,g2_,em1,em2)
    #print fisherIntegrand(em1,em2)




def LogLike(es1,es2,g1,g2,em1,em2):
    global sm
    Ans = -(es1**2 + es2**2 - 2*es1*g1 + g1**2 - 2*es2*g2 + g2**2 + 2*em2*(es2**2*g2 + (1 + es1**2 - 2*es1*g1)*g2 + es2*(-1 + g1**2 - g2**2)) + 2*em1*(es1**2*g1 + g1*(1 + es2**2 - 2*es2*g2) + es1*(-1 - g1**2 + g2**2)) + em1**2*(1 - 2*es1*g1 - 2*es2*g2 + es1**2*(g1**2 + g2**2) + es2**2*(g1**2 + g2**2)) + em2**2*(1 - 2*es1*g1 - 2*es2*g2 + es1**2*(g1**2 + g2**2) + es2**2*(g1**2 + g2**2)))/ (2.*(1 - 2*es1*g1 - 2*es2*g2 + es1**2*(g1**2 + g2**2) + es2**2*(g1**2 + g2**2))*sm**2)
    return Ans

def DLogLikeg1g2(es1,es2,g1,g2,em1,em2):
    global sm
    Ans=-((((-1 + es1**2 - 2*em2*es2 + es2**2 + 2*em2*(es1**2 + es2**2)*g2)*(g1 + es1**2*g1 + es2*g1*(es2 - 2*g2) + es1*(-1 - g1**2 + g2**2)) +em1*(es1**4*(g1 - g2)*(g1 + g2) + 2*es1**3*g1*(-1 + 2*g2**2) +(1 + es2*(g1 - g2))*(1 + es2**2 - 2*es2*g2)*(-1 + es2*(g1 + g2)) + 2*es1*g1*(1 - 2*es2*g2 + es2**2*(-1 + 2*g2**2)) +es1**2*(1 + g1**2*(-1 + 2*es2*(es2 - g2)) + g2*(-3*g2 + 2*es2*(1 - es2*g2 + g2**2)))))*(-((-1 - 2*em1*es1 + es1**2 + es2**2 + 2*em1*(es1**2 + es2**2)*g1)*(es2*(-1 + g1**2) + (1 + es1**2 + es2**2 - 2*es1*g1)*g2 - es2*g2**2)) +em2*(1 + es1**4*(g1 - g2)*(g1 + g2) - 2*es1**3*(g1 + g1**3 - g1*g2**2) -2*es1*g1*(2 - 2*es2*g2 + es2**2*(1 + g1**2 - g2**2)) +es1**2*(1 + g1**2*(5 + 2*es2*(es2 - 2*g2)) - g2*(g2 + 2*es2*(-1 + es2*g2))) +es2*(-2*g2 + es2*(-1 + g1**2*(3 + es2**2 - 4*es2*g2) + g2*(g2 + es2*(2 - es2*g2)))))))/((1 - 2*es1*g1 - 2*es2*g2 + es1**2*(g1**2 + g2**2) + es2**2*(g1**2 + g2**2))**4*sm**4))
    return Ans

def DLogLikeg1g1(es1,es2,g1,g2,em1,em2):
    global sm
    Ans=((-1 + es1**2 - 2*em2*es2 + es2**2 + 2*em2*(es1**2 + es2**2)*g2)*(g1 + es1**2*g1 + es2*g1*(es2 - 2*g2) + es1*(-1 - g1**2 + g2**2)) +em1*(es1**4*(g1 - g2)*(g1 + g2) + 2*es1**3*g1*(-1 + 2*g2**2) +(1 + es2*(g1 - g2))*(1 + es2**2 - 2*es2*g2)*(-1 + es2*(g1 + g2)) + 2*es1*g1*(1 - 2*es2*g2 + es2**2*(-1 + 2*g2**2)) +es1**2*(1 + g1**2*(-1 + 2*es2*(es2 - g2)) + g2*(-3*g2 + 2*es2*(1 - es2*g2 + g2**2)))))**2/((1 - 2*es1*g1 - 2*es2*g2 + es1**2*(g1**2 + g2**2) + es2**2*(g1**2 + g2**2))**4*sm**4)
    return Ans


def DLogLikeg2g2(es1,es2,g1,g2,em1,em2):
    global sm
    Ans=((-1 - 2*em1*es1 + es1**2 + es2**2 + 2*em1*(es1**2 + es2**2)*g1)*(es2*(-1 + g1**2) + (1 + es1**2 + es2**2 - 2*es1*g1)*g2 - es2*g2**2) +em2*(-1 + es1**4*(-g1**2 + g2**2) + 2*es1**3*(g1 + g1**3 - g1*g2**2) + 2*es1*g1*(2 - 2*es2*g2 + es2**2*(1 + g1**2 - g2**2)) +es1**2*(-1 + g1**2*(-5 - 2*es2**2 + 4*es2*g2) + g2*(g2 + 2*es2*(-1 + es2*g2))) +es2*(2*g2 - es2*(-1 + g1**2*(3 + es2**2 - 4*es2*g2) + g2*(g2 + es2*(2 - es2*g2))))))**2/((1 - 2*es1*g1 - 2*es2*g2 + es1**2*(g1**2 + g2**2) + es2**2*(g1**2 + g2**2))**4*sm**4)
    return Ans



def fisherIntegrand12(em1,em2):
    global es1_,es2_,g1_,g2_
    return (DLogLikeg1g2(es1_,es2_,g1_,g2_,em1,em2)*exp(LogLike(es1_,es2_,g1_,g2_,em1,em2)))
def fisherIntegrand11(em1,em2):
    global es1_,es2_,g1_,g2_
    return (DLogLikeg1g1(es1_,es2_,g1_,g2_,em1,em2)*exp(LogLike(es1_,es2_,g1_,g2_,em1,em2)))
def fisherIntegrand22(em1,em2):
    global es1_,es2_,g1_,g2_
    return (DLogLikeg2g2(es1_,es2_,g1_,g2_,em1,em2)*exp(LogLike(es1_,es2_,g1_,g2_,em1,em2)))

if __name__=="__main__":
    main()

