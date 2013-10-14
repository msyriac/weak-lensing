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

import generate_pairs as gen
import galfunctions as fn #module containg all the analytical priors, likelihoods, derivatives


#Shear e with g to return eo1 and eo2
def shear (e,g):
    ee=e[0]+1j*e[1]
    gg=g[0]+1j*g[1]
    sh = (ee-gg)/(1.-conjugate(gg)*ee)
    return [sh.real,sh.imag]


class ToyGenerator:
    def __init__ (self, g1, g2):
        self.g1=g1
        self.g2=g2
        
        ## generate lookup table so that we can sample from 
        ## prob given in eq 18 of Bernstein and Armstrong:

        x=linspace(0,1.0,100)
        y=array([quad(fn.Eprior,0,emax)[0] for emax in x])
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

            e1m=e1+random.gauss(0,fn.sm)
            e2m=e2+random.gauss(0,fn.sm)
        
            ee=sqrt(e1**2.+e2**2.)
            if (ee<1.): ## requier they fall onto <1
                ok=True
            else:
                print "Retrying em"
            
        return e1m,e2m
        

def main(argv):

    n=3
    tollr=1.49e-06
    gg1=-0.01
    gg2=0.02
    sige=0.05#0.05
    sigp=0.3
    ind=1
    doPairs=False
    showTime=False
    
    try:
        opts, args = getopt.getopt(argv,"dhpn:1:2:e:s:t:i:")
    except getopt.GetoptError:
        print "I don't understand the arguments you passed. Run with -h to see available options."
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'python toy.py -n <number of galaxies=3> -1 <shear component 1=-0.01> -2 <shear component 2=0.02> -e <sigma_e=0.05> -s <sigma_pr=0.3> -t <integrator tolerance=1.49e-06> -i <file index=1> \n Add -p if you want to do g,h pairs instead \n Add -d if you want to display average time per galaxy'
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
        elif opt == "-p":
            doPairs=True
        elif opt == "-d":
            showTime=True



    started=time.time()
    fileappend=str(started).translate(None,'.')

    fG = open('data/'+"g"+fileappend+'i'+str(ind)+'.csv', 'w')
    
    fn.initGalFuncs(sigp,sige)

    if doPairs:
        fH = open('data/'+"h"+fileappend+'i'+str(ind)+'.csv', 'w')
        p=gen.Pairs()

        for i in range(n):
            random.seed(started*ind+i)
            [gg1,gg2,hh1,hh2]=p.draw()
            A=ToyGenerator(gg1,gg2)
            emA=A.generate()
            fn.updateEs(emA[0],emA[1])
            row=getDerivs(tollr)
            fG.write(','.join(str(j) for j in row) + '\n')

            random.seed(started*ind+i+1)
            B=ToyGenerator(hh1,hh2)
            emB=B.generate()
            fn.updateEs(emB[0],emB[1])
            row=getDerivs(tollr)
            fH.write(','.join(str(j) for j in row) + '\n')
        if showTime:
            ended=time.time()
            print "Average time per galaxy pair was ", (ended-started)/n, " seconds."
        fH.close()

    else:
        random.seed(started*ind)
        A=ToyGenerator(gg1,gg2)
        for i in range(n):
            em=A.generate()
            fn.updateEs(em[0],em[1])
            row=getDerivs(tollr)
            fG.write(','.join(str(j) for j in row) + '\n')
        if showTime:
            ended=time.time()
            print "Average time per galaxy was ", (ended-started)/n, " seconds."

    fG.close()    


def getDerivs(tol):


    P=dblquad (fn.IntegrandP,-1.0, 1.0, lambda x:-sqrt(1.-x*x), lambda x:sqrt(1.-x*x),epsabs=tol, epsrel=tol)[0]
    Q1=dblquad (fn.IntegrandQ1,-1.0, 1.0, lambda x:-sqrt(1.-x*x), lambda x:sqrt(1.-x*x),epsabs=tol, epsrel=tol)[0]
    Q2=dblquad (fn.IntegrandQ2,-1.0, 1.0, lambda x:-sqrt(1.-x*x), lambda x:sqrt(1.-x*x),epsabs=tol, epsrel=tol)[0]
    R11=dblquad (fn.IntegrandR11,-1.0, 1.0, lambda x:-sqrt(1.-x*x), lambda x:sqrt(1.-x*x),epsabs=tol, epsrel=tol)[0]
    R12=dblquad (fn.IntegrandR12,-1.0, 1.0, lambda x:-sqrt(1.-x*x), lambda x:sqrt(1.-x*x),epsabs=tol, epsrel=tol)[0]
    R22=dblquad (fn.IntegrandR22,-1.0, 1.0, lambda x:-sqrt(1.-x*x), lambda x:sqrt(1.-x*x),epsabs=tol, epsrel=tol)[0]
        
    row=[P,Q1,Q2,R11,R12,R22]

    return row


    

if __name__=="__main__":
    main(sys.argv[1:])
