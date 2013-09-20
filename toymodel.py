#!/usr/bin/env python
from scipy import *
from  scipy.integrate import quad
from scipy.integrate import dblquad 
from scipy.interpolate import InterpolatedUnivariateSpline
import pylab, math,random
import cmath

sig_pr=0.3
sigma_e=0.05

def shear (e,g):
    ee=e[0]+1j*e[1]
    gg=g[0]+1j*g[1]
    sh = (ee-gg)/(1-conjugate(gg)*ee)
    return [sh.real,sh.imag]


    

def Eprior(e):
    return e*(1-e*e)**2*exp(-e**2/(2*sig_pr**2))


class ToyGenerator:
    def __init__ (self, g1, g2): # Changed sig_i to sig_pr
        self.g1=g1
        self.g2=g2
        
        ## generate lookup table so that we can sample from 
        ## prob given in eq 18 of Bernstein and Armstrong:

        x=linspace(0,1.0,100)
        #y=array([quad(lambda x:x*(1-x*x)**x*exp(-x**2/(2*0.3**2)),0,a)[0] for a in x]) # Typo? Should be x^2 not x^x, also changed 0.3 to sig_pr -Mat
        y=array([quad(Eprior,0,emax)[0] for emax in x])
        y/=y[-1]
        #print x
        #print y


        self.genI=InterpolatedUnivariateSpline(y,x) # This is to invert the CDF right? Just making sure. - Mat
        print self.genI
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

        ##intrinsic ellipticty now shear
        ### MAT CHECK THE SHIT BELOW # Corrected to agree with Seitz and Schneider 1997 Eq 3.2 -Mat
        # This is cut off at second order. Try using full expression if encountering problems later.

        

        e1,e2=shear([e1,e2],[self.g1,self.g2])
        
        ### now add error
        ok=False
        while not ok:

            e1m=e1+random.gauss(0,sigma_e)
            e2m=e2+random.gauss(0,sigma_e)
        
            ee=sqrt(e1**2+e2**2)
            if (ee<1): ## requier they fall onto <1
                ok=True
            
        return e1m,e2m
        


def zeroDer(es0,es1):
    es=sqrt(es0**2+es1**2)
    eo=shear([es0,es1],gg_)
    dlta=sqrt((em_[0]-eo[0])**2+(em_[1]-eo[1])**2)
    return Eprior(es)*exp(-dlta*dlta/(2*sigma_e**2))




def firstDer1(es0,es1):
    es=sqrt(es0**2+es1**2)
    eo=[es0,es1] #expanding around g=0
    dlta=sqrt((em_[0]-eo[0])**2+(em_[1]-eo[1])**2)
    Gprime= exp(-dlta*dlta/(2*sigma_e**2))*(-1.*dlta/sigma_e**2)

    den = dlta
    top = -1*(es0**3+em_[0]-es0**2*em_[0]+es1**2*em_[0]+es0*(-1+es1**2-2*es1*em_[1]))
    return Eprior(es)*Gprime*top/den

def firstDer2(es0,es1):
    es=sqrt(es0**2+es1**2)
    eo=[es0,es1] #expanding around g=0
    dlta=sqrt((em_[0]-eo[0])**2+(em_[1]-eo[1])**2)
    Gprime= exp(-dlta*dlta/(2*sigma_e**2))*(-1.*dlta/sigma_e**2)

    den = dlta
    top = -1*(es1**3+em_[1]-es1**2*em_[1]+es0**2*em_[1]+es1*(-1+es0**2-2*es0*em_[0]))
    return Eprior(es)*Gprime*top/den



def getDerivative (d,em,gg):
    global em_
    global gg_
    em_=em
    gg_=gg
    if d==0:
        return dblquad (zeroDer,-1.0, 1.0, lambda x:-sqrt(1-x*x), lambda x:sqrt(1-x*x),epsabs=1.49e-02, epsrel=1.49e-02)[0]
    if d==1:
        fder1=dblquad (firstDer1,-1.0, 1.0, lambda x:-sqrt(1-x*x), lambda x:sqrt(1-x*x),epsabs=1.49e-02, epsrel=1.49e-02)[0]
        fder2=dblquad (firstDer2,-1.0, 1.0, lambda x:-sqrt(1-x*x), lambda x:sqrt(1-x*x),epsabs=1.49e-02, epsrel=1.49e-02)[0]
        return [fder1,fder2]


if __name__=="__main__":
    A=ToyGenerator(0.01,0.02)
    eps=0.001
    random.seed(3)
    for x in range(3):
        em=A.generate()

        P0=getDerivative(0,em,[0.,0.])
        #D1N1=(getDerivative(0,em,[eps,0.])-P0)/eps
        #D1N2=(getDerivative(0,em,[0.0,eps])-P0)/eps
        D1=getDerivative(1,em,[0.,0.])
        
        print em,P0, D1



        

