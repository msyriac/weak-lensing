#!/usr/bin/env python
from scipy import *
from  scipy.integrate import quad
from scipy.integrate import dblquad 
from scipy.interpolate import InterpolatedUnivariateSpline
import cmath
import numpy
import matplotlib
matplotlib.use('Agg')
import pylab, math,random, sys, time
import matplotlib.pyplot as p

sig_pr=0.3
sigma_e=0.05

#Shear e with g to return eo1 and eo2
def shear (e,g):
    ee=e[0]+1j*e[1]
    gg=g[0]+1j*g[1]
    sh = (ee-gg)/(1-conjugate(gg)*ee)
    return [sh.real,sh.imag]


    
#Return prior as a function of |es|
def Eprior(e2):
    e=math.sqrt(e2)
    return e*((1-e2)**2)*exp(-e2/(sp2twice))

def fEprior(es0,es1):
    return Eprior((es0**2+es1**2))

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
        theta=random.uniform(0,2*math.pi);
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
        
            ee=sqrt(e1**2+e2**2)
            if (ee<1): ## requier they fall onto <1
                ok=True
            else:
                print "Retrying em"
            
        return e1m,e2m
        


def LikeZeroDer(es1,es2): #,ggg_):
    #if (ggg_!=[0.,0.]): print ggg_
    #eo=shear([es1,es2],ggg_)
    eo=[es1,es2] #####
    dlta2=(em1-eo[0])**2+(em2-eo[1])**2
    #print es1, es2
    #print em1, em2, eo[0], eo[1]
    #print sm2
    #print dlta2
    return exp(-dlta2/(2*sm2))




def LikeFirstDer1(es1,es2):

    dlta2=(em1-es1)**2+(em2-es2)**2

    
    #em12=em1**2
    #em14=em1s**2

    #em22=em2**2
    #em24=em2s**2

    es12=es1**2
    #es14=es1s**2

    es22=es2**2
    #es24=es2s**2

    return exp(-dlta2/(2*sm2))*(-em1+em1*es12+es1-es12*es1+2*em2*es1*es2-em1*es22-es1*es22)/sm2

def LikeFirstDer2(es1,es2):

    dlta2=(em1-es1)**2+(em2-es2)**2

    
    #em12=em1**2
    #em14=em1s**2

    #em22=em2**2
    #em24=em2s**2

    es12=es1**2
    #es14=es1s**2

    es22=es2**2
    #es24=es2s**2



    
    return exp(-dlta2/(2*sm2))*(-em2+em2*es22+es2-es22*es2+2*em1*es2*es1-em2*es12-es2*es12)/sm2

def LikeSecondDer11(es1,es2):

    dlta2=(em1-es1)**2+(em2-es2)**2

    
    #em12=em1**2
    #em14=em1s**2

    #em22=em2**2
    #em24=em2s**2

    es12=es1**2
    #es14=es1s**2

    es22=es2**2
    #es24=es2s**2

    #return exp(-dlta*dlta/(2*sm**2))*(((em1 - es1)*(-1 + es1**2) + 2*em2*es1*es2 - (em1 + es1)*es2**2)**2 + (2*em1*es1*(-1 + es1**2 - 3*es2**2) - (-1 + 3*es1**2 - es2**2)*(-1 + es1**2 - 2*em2*es2 + es2**2))*sm**2)/sm**4
    

    return exp(-dlta2/(2*sm2))*(((em1 - es1)*(-1 + es12) + 2*em2*es1*es2 - (em1 + es1)*es22)**2 + (2*em1*es1*(-1 + es12 - 3*es22) - (-1 + 3*es12 - es22)*(-1 + es12 - 2*em2*es2 + es22))*sm2)/sm4
    

def LikeSecondDer22(es1,es2):


    dlta2=(em1-es1)**2+(em2-es2)**2

    
    #em12=em1**2
    #em14=em1s**2

    #em22=em2**2
    #em24=em2s**2

    es12=es1**2
    es14=es12**2

    es22=es2**2
    #es24=es2s**2

    #return exp(-dlta*dlta/(2*sm**2))*((em2*(1 + es1**2 - es2**2) + es2*(-1 - 2*em1*es1 + es1**2 + es2**2))**2 + (es1**4 - 2*es1**2*es2*(3*em2 + es2) + (1 + 2*em2*es2 - 3*es2**2)*(-1 + es2**2) - 2*em1*(es1 + es1**3 - 3*es1*es2**2))*sm**2)/sm**4
    return exp(-dlta2/(2*sm2))*((em2*(1 + es12 - es22) + es2*(-1 - 2*em1*es1 + es12 + es22))**2 + (es14 - 2*es12*es2*(3*em2 + es2) + (1 + 2*em2*es2 - 3*es22)*(-1 + es22) - 2*em1*(es1 + es12*es1 - 3*es1*es22))*sm2)/sm4


def LikeSecondDer12(es1,es2):
    dlta2=(em1-es1)**2+(em2-es2)**2

    

    es12=es1**2
    es14=es12**2

    es22=es2**2
    es24=es22**2


    #return exp(-dlta*dlta/(2*sm**2))*(2*em1**2*es1*es2*(-1 + es1**2 - es2**2) + es1*(2*em2**2*es2*(-1 - es1**2 + es2**2) + es2*(-1 + es1**2 + es2**2)*(-1 + es1**2 + es2**2 - 4*sm**2) + em2*(-1 + 4*es2**2 + (es1**2 - 3*es2**2)*(es1**2 + es2**2 - 2*sm**2))) + em1*(-(em2*(-1 + es1**4 - 6*es1**2*es2**2 + es2**4)) + es2*(-1 - 3*es1**4 + es2**4 - 2*es2**2*sm**2 + es1**2*(4 - 2*es2**2 + 6*sm**2))))/sm**4
    return exp(-dlta2/(2*sm2))*(2*em12*es1*es2*(-1 + es12 - es22) + es1*(2*em22*es2*(-1 - es12 + es22) + es2*(-1 + es12 + es22)*(-1 + es12 + es22 - 4*sm2) + em2*(-1 + 4*es22 + (es12 - 3*es22)*(es12 + es22 - 2*sm2))) + em1*(-(em2*(-1 + es14 - 6*es12*es22 + es24)) + es2*(-1 - 3*es14 + es24 - 2*es22*sm2 + es12*(4 - 2*es22 + 6*sm2))))/sm4
    


def getDerivative (d,em,gg):
    global em_
    #global gg_
    em_=em
    #gg_=gg
    if d==0:
        #return dblquad (zeroDer,-1.0, 1.0, lambda x:-sqrt(1-x*x), lambda x:sqrt(1-x*x),epsabs=tol, epsrel=tol)[0]
        return dblquad (IntegrandP,-1.0, 1.0, lambda x:-sqrt(1-x*x), lambda x:sqrt(1-x*x),epsabs=tol, epsrel=tol)[0]
    if d==1:
        #fder1=dblquad (firstDer1,-1.0, 1.0, lambda x:-sqrt(1-x*x), lambda x:sqrt(1-x*x),epsabs=tol, epsrel=tol)[0]
        #fder2=dblquad (firstDer2,-1.0, 1.0, lambda x:-sqrt(1-x*x), lambda x:sqrt(1-x*x),epsabs=tol, epsrel=tol)[0]
        fder1=dblquad (IntegrandQ1,-1.0, 1.0, lambda x:-sqrt(1-x*x), lambda x:sqrt(1-x*x),epsabs=tol, epsrel=tol)[0]
        fder2=dblquad (IntegrandQ2,-1.0, 1.0, lambda x:-sqrt(1-x*x), lambda x:sqrt(1-x*x),epsabs=tol, epsrel=tol)[0]
        return [fder1,fder2]
    if d==2:
        #sder11=dblquad (secondDer11,-1.0, 1.0, lambda x:-sqrt(1-x*x), lambda x:sqrt(1-x*x),epsabs=tol, epsrel=tol)[0]
        #sder12=dblquad (secondDer12,-1.0, 1.0, lambda x:-sqrt(1-x*x), lambda x:sqrt(1-x*x),epsabs=tol, epsrel=tol)[0]
        #sder22=dblquad (secondDer22,-1.0, 1.0, lambda x:-sqrt(1-x*x), lambda x:sqrt(1-x*x),epsabs=tol, epsrel=tol)[0]
        sder11=dblquad (IntegrandR11,-1.0, 1.0, lambda x:-sqrt(1-x*x), lambda x:sqrt(1-x*x),epsabs=tol, epsrel=tol)[0]
        sder12=dblquad (IntegrandR12,-1.0, 1.0, lambda x:-sqrt(1-x*x), lambda x:sqrt(1-x*x),epsabs=tol, epsrel=tol)[0]
        sder22=dblquad (IntegrandR22,-1.0, 1.0, lambda x:-sqrt(1-x*x), lambda x:sqrt(1-x*x),epsabs=tol, epsrel=tol)[0]
        return [sder11,sder12,sder22]




def IntegrandP(es0,es1):
    #ggg_=[0.,0.] ####
    return fEprior(es0,es1)*LikeZeroDer(es0,es1) #,[0.,0.])

#def derivP(es0,es1,gg):
    #ggg_=gg
    #print ggg_
    #return fEprior(es0,es1)*LikeZeroDer(es0,es1,gg)

def IntegrandQ1(es0,es1):

    ####
    #P0=derivP(es0,es1,[0.,0.])

    #Finite Differences
    #Forward1=derivP(es0,es1,[eps,0.])
    #Backward1=derivP(es0,es1,[-eps,0.])
    #Central finite difference for first derivative
    #return ((Forward1-Backward1)/(2*eps))

    #Central finite differences for second derivative

  
    #####    
    
    
    return fEprior(es0,es1)*LikeFirstDer1(es0,es1)

def IntegrandQ2(es0,es1):
    #P0=derivP(es0,es1,[0.,0.])

    #Forward2=derivP(es0,es1,[0.0,eps])
    #Backward2=derivP(es0,es1,[0.0,-eps])
    #return ((Forward2-Backward2)/(2*eps))

    
    return fEprior(es0,es1)*LikeFirstDer2(es0,es1)

def IntegrandR11(es0,es1):

    #P0=derivP(es0,es1,[0.,0.])

    #Finite Differences
    #Forward1=derivP(es0,es1,[eps,0.])
    #Backward1=derivP(es0,es1,[-eps,0.])
    #return ((Forward1-2*P0+Backward1)/eps**2)
     
    return fEprior(es0,es1)*LikeSecondDer11(es0,es1)

def IntegrandR12(es0,es1):
    #P0=derivP(es0,es1,[0.,0.])

    #Finite Differences
    #Forward2=derivP(es0,es1,[0.0,eps])
    #Backward2=derivP(es0,es1,[0.0,-eps])
    #return ((Forward2-2*P0+Backward2)/eps**2)

    
    return fEprior(es0,es1)*LikeSecondDer12(es0,es1)

def IntegrandR22(es0,es1):
    #P0=derivP(es0,es1,[0.,0.])

    #Forward1=derivP(es0,es1,[eps,0.])
    #Backward1=derivP(es0,es1,[-eps,0.])
    #Forward2=derivP(es0,es1,[0.0,eps])
    #Backward2=derivP(es0,es1,[0.0,-eps])
    #For1For2=derivP(es0,es1,[eps,eps])
    #Bak1Bak2=derivP(es0,es1,[-eps,-eps])
    #return ((For1For2-Forward1-Forward2+2*P0-Backward1-Backward2+Bak1Bak2)/(2*eps**2))

    return fEprior(es0,es1)*LikeSecondDer22(es0,es1)


if __name__=="__main__":
    global em1
    global em2

    global sm2
    global sm4

    global em12
    global em22

    global sp2twice
    global tol
    tol=1.49e-01

    #global eps
    #global ggg_
    #ggg_=[0.,0.]

    sp2twice=2*sig_pr**2

    

    sm2=sigma_e**2
    sm4=sm2**2
    
    appG=numpy.matrix([[0.01],[0.02]])

    print appG
    
    A=ToyGenerator(appG.flat[0],appG.flat[1])
    #eps=0.001
    random.seed(1241)#3)
    #P=[]
    #Q=[]
    #R=[]
    #Q_norm=[]

    Cinv=numpy.zeros((2,2))
    Qsum=0
    biasp1=[]
    biasp2=[]
    x=[]

    started=time.time()

    freq=5
    tottime=0
    k=0


    n=10000
    #print LikeZeroDer(0.1,-0.1)
    #print gg_
    
    for i in range(n):
        
        em=A.generate()
        em1=em[0]
        em2=em[1]

        em12=em1**2
        #em14=em1s**2

        em22=em2**2
        #em24=em2s**2

        
        P=(dblquad (IntegrandP,-1.0, 1.0, lambda x:-sqrt(1-x*x), lambda x:sqrt(1-x*x),epsabs=tol, epsrel=tol)[0])
        Q=(numpy.matrix( [ [dblquad (IntegrandQ1,-1.0, 1.0, lambda x:-sqrt(1-x*x), lambda x:sqrt(1-x*x),epsabs=tol, epsrel=tol)[0]],
                            [dblquad (IntegrandQ2,-1.0, 1.0, lambda x:-sqrt(1-x*x), lambda x:sqrt(1-x*x),epsabs=tol, epsrel=tol)[0]]
                            ] ))
        Rcross=dblquad (IntegrandR12,-1.0, 1.0, lambda x:-sqrt(1-x*x), lambda x:sqrt(1-x*x),epsabs=tol, epsrel=tol)[0]
        
        R=(numpy.matrix( [ [dblquad (IntegrandR11,-1.0, 1.0, lambda x:-sqrt(1-x*x), lambda x:sqrt(1-x*x),epsabs=tol, epsrel=tol)[0], Rcross],
                             [Rcross, dblquad (IntegrandR22,-1.0, 1.0, lambda x:-sqrt(1-x*x), lambda x:sqrt(1-x*x),epsabs=tol, epsrel=tol)[0]] ]))


        Q_norm=(Q/P)
        Qsum+=Q_norm
        Cinv+= Q_norm*(Q_norm.transpose())-(R/P)

        if (((i+1) % freq)==0):
            now=time.time()
            elapsed=now-started
            started=now
            tottime+=elapsed
            k+=1
            avgtime=(tottime/k)/freq
            print i+1, " out of ", n, " galaxies done."
             
            print "Average time for 1 galaxy is ", avgtime, " seconds."

            predict=(n-i-1)*avgtime
            print "Predicted time of ", predict/60., " minutes left."
            
            
            x.append(i+1)
            infG=numpy.linalg.inv(Cinv)*Qsum
            bias=infG-appG
            biasp1.append(bias.flat[0]*100/appG.flat[0])
            biasp2.append(bias.flat[1]*100/appG.flat[1])
            print infG
            print biasp1[-1]
            print biasp2[-1]
            p.clf()
            p.plot(x,biasp1)
            p.plot(x,biasp2)
            p.savefig('bias.png')

        '''    
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

        
    print infG
    print biasp1[-1]
    print biasp2[-1]
    p.clf()
    p.plot(x,biasp1)
    p.plot(x,biasp2)
    p.savefig('bias.png')
        
    p.show()


    
    



        





