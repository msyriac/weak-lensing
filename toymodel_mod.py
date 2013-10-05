#!/usr/bin/env python
from scipy import *
from  scipy.integrate import quad
from scipy.integrate import dblquad 
from scipy.interpolate import InterpolatedUnivariateSpline
import cmath
import numpy
import matplotlib
matplotlib.use('Agg')
import pylab, math,random, sys, time, getopt
import matplotlib.pyplot as p
import csv


#Shear e with g to return eo1 and eo2
def shear (e,g):
    ee=e[0]+1j*e[1]
    gg=g[0]+1j*g[1]
    sh = (ee-gg)/(1.-conjugate(gg)*ee)
    return [sh.real,sh.imag]


    
#Return prior as a function of |es|
def Eprior(e2):
    e=math.sqrt(e2)
    return e*((1-e2)**2.)*exp(-e2/(sp2twice))

def fEprior(es0,es1):
    return Eprior((es0**2.+es1**2.))

class ToyGenerator:
    def __init__ (self, g1, g2):
        self.g1=g1
        self.g2=g2
        
        ## generate lookup table so that we can sample from 
        ## prob given in eq 18 of Bernstein and Armstrong:

        x=linspace(0,1.0,100)
        y=array([quad(Eprior,0,emax)[0] for emax in x])
        y/=y[-1]


        self.genI=InterpolatedUnivariateSpline(y,x) # This is to invert the CDF right? Just making sure. - Mat
        #print self.genI
        #xx=arange(0,1,0.0001)
        #pylab.plot(y,x,'ro')
        #pylab.plot(xx,map(self.genI,xx),'r-')
        #pylab.show()


    def generateE(self):
        return self.genI(random.uniform(0.,1.0))


    def generate (self):
        E=self.generateE();
        theta=random.uniform(0,2.*math.pi);
        e1=E*cos(theta)
        e2=E*sin(theta)
        #print e1,e2

        ##intrinsic ellipticty now shear
        e1,e2=shear([e1,e2],[self.g1,self.g2])
        
        ### now add error
        ok=False
        while not ok:

            e1m=e1+random.gauss(0,sigma_e)
            e2m=e2+random.gauss(0,sigma_e)
        
            ee=sqrt(e1**2.+e2**2.)
            if (ee<1.): ## requier they fall onto <1
                ok=True
            else:
                print "Retrying em"
            
        return e1m,e2m
        


#def PriorFirstDer1(es1,es2):
    #-(es1*(-1 + es1**2 + es2**2)**2*(es1**4 + es2**4 + sp**2 + es1**2*(-1 + 2*es2**2 - 5*sp**2) - es2**2*(1 + 5*sp**2)))s

    #def PriorFirstDer2(es1,es2):
    
def LikeZeroDer(es1,es2):
    eo=[es1,es2]
    #eo=shear([es1,es2],gg_)
    dlta2=(em1-eo[0])**2.+(em2-eo[1])**2.
    return exp(-dlta2/(2.*sm2))




def LikeFirstDer1(es1,es2):

    dlta2=(em1-es1)**2.+(em2-es2)**2.
    es12=es1**2.

    es22=es2**2.

    return exp(-dlta2/(2.*sm2))*(-em1+em1*es12+es1-es12*es1+2.*em2*es1*es2-em1*es22-es1*es22)/sm2

def LikeFirstDer2(es1,es2):

    dlta2=(em1-es1)**2.+(em2-es2)**2.
    es12=es1**2.
    es22=es2**2.
    return exp(-dlta2/(2.*sm2))*(-em2+em2*es22+es2-es22*es2+2.*em1*es2*es1-em2*es12-es2*es12)/sm2

def LikeSecondDer11(es1,es2):

    dlta2=(em1-es1)**2.+(em2-es2)**2.
    es12=es1**2.
    es22=es2**2.
    return exp(-dlta2/(2.*sm2))*(((em1 - es1)*(-1. + es12) + 2.*em2*es1*es2 - (em1 + es1)*es22)**2. + (2.*em1*es1*(-1. + es12 - 3.*es22) - (-1. + 3.*es12 - es22)*(-1. + es12 - 2.*em2*es2 + es22))*sm2)/sm4
    

def LikeSecondDer22(es1,es2):


    dlta2=(em1-es1)**2.+(em2-es2)**2.
    es12=es1**2.
    es14=es12**2.

    es22=es2**2.
    return exp(-dlta2/(2.*sm2))*((em2*(1. + es12 - es22) + es2*(-1. - 2.*em1*es1 + es12 + es22))**2. + (es14 - 2.*es12*es2*(3.*em2 + es2) + (1. + 2.*em2*es2 - 3.*es22)*(-1. + es22) - 2.*em1*(es1 + es12*es1 - 3.*es1*es22))*sm2)/sm4


def LikeSecondDer12(es1,es2):
    dlta2=(em1-es1)**2.+(em2-es2)**2.
    es12=es1**2.
    es14=es12**2.

    es22=es2**2.
    es24=es22**2.

    return exp(-dlta2/(2.*sm2))*(2.*em12*es1*es2*(-1. + es12 - es22) + es1*(2.*em22*es2*(-1. - es12 + es22) + es2*(-1. + es12 + es22)*(-1. + es12 + es22 - 4.*sm2) + em2*(-1. + 4.*es22 + (es12 - 3.*es22)*(es12 + es22 - 2.*sm2))) + em1*(-(em2*(-1. + es14 - 6.*es12*es22 + es24)) + es2*(-1 - 3.*es14 + es24 - 2.*es22*sm2 + es12*(4. - 2.*es22 + 6.*sm2))))/sm4
    

'''
def getDerivative (d,em,gg):
    global gg_
    gg_=gg
    if d==0:
        return dblquad (IntegrandP,-1.0, 1.0, lambda x:-sqrt(1-x*x), lambda x:sqrt(1-x*x),epsabs=tol, epsrel=tol)[0]
    if d==1:
        fder1=dblquad (IntegrandQ1,-1.0, 1.0, lambda x:-sqrt(1-x*x), lambda x:sqrt(1-x*x),epsabs=tol, epsrel=tol)[0]
        fder2=dblquad (IntegrandQ2,-1.0, 1.0, lambda x:-sqrt(1-x*x), lambda x:sqrt(1-x*x),epsabs=tol, epsrel=tol)[0]
        return [fder1,fder2]
    if d==2:
        sder11=dblquad (IntegrandR11,-1.0, 1.0, lambda x:-sqrt(1-x*x), lambda x:sqrt(1-x*x),epsabs=tol, epsrel=tol)[0]
        sder12=dblquad (IntegrandR12,-1.0, 1.0, lambda x:-sqrt(1-x*x), lambda x:sqrt(1-x*x),epsabs=tol, epsrel=tol)[0]
        sder22=dblquad (IntegrandR22,-1.0, 1.0, lambda x:-sqrt(1-x*x), lambda x:sqrt(1-x*x),epsabs=tol, epsrel=tol)[0]
        return [sder11,sder12,sder22]
'''



def IntegrandP(es0,es1):
    return fEprior(es0,es1)*LikeZeroDer(es0,es1)


def IntegrandQ1(es0,es1):
    return fEprior(es0,es1)*LikeFirstDer1(es0,es1)

def IntegrandQ2(es0,es1):
    return fEprior(es0,es1)*LikeFirstDer2(es0,es1)

def IntegrandR11(es0,es1):
    return fEprior(es0,es1)*LikeSecondDer11(es0,es1)

def IntegrandR12(es0,es1):
    return fEprior(es0,es1)*LikeSecondDer12(es0,es1)

def IntegrandR22(es0,es1):
    return fEprior(es0,es1)*LikeSecondDer22(es0,es1)


def main(argv):

    n=50
    tollr=1.49e-03
    gg1=-0.01
    gg2=0.02
    sige=0.05
    sigp=0.3
    i=1
    
    try:
        opts, args = getopt.getopt(argv,"n:1:2:e:s:t:i:")
    except getopt.GetoptError:
        print "I don't understand the arguments you passed. Run with -h to see available options."
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'python toymodel_mod.py -n <number of galaxies> -1 <shear component 1> -2 <shear component 2> -e <sigma_e> -s <sigma_pr> -t <integrator tolerance> -i <file index>'
            sys.exit()
        elif opt == "-n":
            n=int(arg)
        elif opt == "-1":
            gg1=float(arg)
        elif opt == "-2":
            gg2=float(arg)
        elif opt == "-e":
            sige=float(arg)
        elif opt == "-s":
            sigp=float(arg)
        elif opt == "-t":
            tollr=float(arg)
        elif opt == "-i":
            ind=int(arg)

    makeDerivs(n,gg1,gg2,sige,sigp,tollr,ind)


def makeDerivs(n,gg1,gg2,sige,sigp,tollr,ind):
    global em1, em2, sm2, sm4, em12, em22, sp2twice, tol,sig_pr, sigma_e

    sig_pr=sigp
    sigma_e=sige #0.05


    tol=tollr #1.49e-03 #integrator tolerance

    # some global stuff to speed up function calls by integrator
    sp2twice=2*sig_pr**2.
    sm2=sigma_e**2.
    sm4=sm2**2.
    ##
    
    appG=[gg1,gg2]#[-0.01,0.02]

    print "Applied shear is ",appG
    
    A=ToyGenerator(appG[0],appG[1])
    eps=0.001
    #random.seed(12296)#3)
    random.seed(time.time())
    
    Cinv=numpy.zeros((2,2))
    Qsum=0


    started=time.time()
    fileappend=str(started)


    tottime=0
    k=0


    f = open('data/'+fileappend+'i'+str(ind)+'.csv', 'w')
    #writer = csv.writer(f,delimiter=' ')


    for i in range(n):
        
        em=A.generate()
        
        #Stuff for global use to speed up function calls from integrator
        em1=em[0]
        em2=em[1]
        em12=em1**2.
        em22=em2**2.
        ##


        P=dblquad (IntegrandP,-1.0, 1.0, lambda x:-sqrt(1.-x*x), lambda x:sqrt(1.-x*x),epsabs=tol, epsrel=tol)[0]
        Q1=dblquad (IntegrandQ1,-1.0, 1.0, lambda x:-sqrt(1.-x*x), lambda x:sqrt(1.-x*x),epsabs=tol, epsrel=tol)[0]
        Q2=dblquad (IntegrandQ2,-1.0, 1.0, lambda x:-sqrt(1.-x*x), lambda x:sqrt(1.-x*x),epsabs=tol, epsrel=tol)[0]
        R11=dblquad (IntegrandR11,-1.0, 1.0, lambda x:-sqrt(1.-x*x), lambda x:sqrt(1.-x*x),epsabs=tol, epsrel=tol)[0]
        R12=dblquad (IntegrandR12,-1.0, 1.0, lambda x:-sqrt(1.-x*x), lambda x:sqrt(1.-x*x),epsabs=tol, epsrel=tol)[0]
        R22=dblquad (IntegrandR22,-1.0, 1.0, lambda x:-sqrt(1.-x*x), lambda x:sqrt(1.-x*x),epsabs=tol, epsrel=tol)[0]

        
        #writer.writerow([P,Q1,Q2,R11,R12,R22])
        row=[P,Q1,Q2,R11,R12,R22]
        f.write(','.join(str(j) for j in row) + '\n')

        '''
        #Code to check correctness of analytical derivatives
        #Zeroth derivative
    
        P0=getDerivative(0,em,[0.,0.])

        #Finite Differences
        Forward1=getDerivative(0,em,[eps,0.])
        Backward1=getDerivative(0,em,[-eps,0.])
        Forward2=getDerivative(0,em,[0.0,eps])
        Backward2=getDerivative(0,em,[0.0,-eps])
        For1For2=getDerivative(0,em,[eps,eps])
        Bak1Bak2=getDerivative(0,em,[-eps,-eps])

        #Central finite difference for first derivative
        D1N1=(Forward1-Backward1)/(2*eps)
        D1N2=(Forward2-Backward2)/(2*eps)

        #Central finite differences for second derivative
        D2N11=(Forward1-2*P0+Backward1)/eps**2
        D2N12=(For1For2-Forward1-Forward2+2*P0-Backward1-Backward2+Bak1Bak2)/(2*eps**2)
        #(getDerivative(0,em,[eps,eps]) - getDerivative(0,em,[eps,-eps]) - getDerivative(0,em,[-eps,eps]) + getDerivative(0,em,[-eps,-eps]) )/(4*eps**2) # more efficient exists
        D2N22=(Forward2-2*P0+Backward2)/eps**2

        #Analytical first and second derivatives
        D1=getDerivative(1,em,[0.,0.])
        D2=getDerivative(2,em,[0.,0.])


        print "Input em is ", em
        print "Zeroth derivative is ",P0
        print "Analytic  first derivative is ", D1
        print "Numerical first derivative is " ,[D1N1,D1N2]
        print "Percentage difference is ",[(D1N1-D1[0])*100/D1[0], (D1N2-D1[1])*100/D1[1]], "%"
        print "Analytic  second derivative is ", D2
        print "Numerical second derivative is " ,[D2N11,D2N12,D2N22]
        print "Percentage difference is ",[(D2N11-D2[0])*100/D2[0], (D2N12-D2[1])*100/D2[1], (D2N22-D2[2])*100/D2[2]], "%"
        '''

    #print "Average time = ", (time.time()-started)/n, " seconds."
    f.close()


if __name__=="__main__":
    main(sys.argv[1:])
