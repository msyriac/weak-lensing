#!/usr/bin/env python
from scipy import *
from  scipy.integrate import quad
from scipy.interpolate import InterpolatedUnivariateSpline
import pylab, math,random

class ToyGenerator:
    def __init__ (self, g1, g2, sig_i):
        self.g1=g1
        self.g2=g2
        
        ## generate lookup table so that we can sample from 
        ## prob given in eq 18 of Bernstein and Armstrong:

        x=linspace(0,1.0,100)
        y=array([quad(lambda x:x*(1-x*x)**x*exp(-x**2/(2*0.3**2)),0,a)[0] for a in x])
        y/=y[-1]
        print x
        print y


        self.genI=InterpolatedUnivariateSpline(y,x)
        xx=arange(0,1,0.0001)
        #pylab.plot(y,x,'ro')
        #pylab.plot(xx,map(self.genI,xx),'r-')
        #pylab.show()


    def generateE(self):
        return self.genI(random.uniform(0.,1.0))


    def generate (self, sigma):
        E=self.generateE();
        theta=random.uniform(0,2*math.pi);
        e1=E*cos(theta)
        e2=E*sin(theta)

        ##intrinsic ellipticty now shear
        ### MAT CHECK THE SHIT BELOW
        corrfact=1/(1-(e1*self.g1+e2*self.g2))
        e1=(e1+self.g1)*corrfact
        e2=(e2+self.g2)*corrfact
        
        ### now add error
        e1+=random.gauss(0,sigma)
        e2+=random.gauss(0,sigma)
        
        ee=sqrt(e1**2+e2**2)
        if (ee>1): ## requier they fall onto <1
            e1/=ee
            e2/=ee
            
        return e1,e2
        
        

if __name__=="__main__":
    A=ToyGenerator(0.01,0.02,0.3)
    
    for x in range(10):
        print A.generate(0.05)
