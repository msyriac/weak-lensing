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
try:
    from numba import autojit
except:
    pass

#Shear e with g to return eo1 and eo2
##@autojit
def shear (e,g):
    ee=e[0]+1j*e[1]
    gg=g[0]+1j*g[1]
    sh = (ee-gg)/(1.-conjugate(gg)*ee)
    return [sh.real,sh.imag]


class ToyGenerator:
    def __init__ (self, g1, g2,sigma_pr=0.3, sm=0.05):
        self.g1=g1
        self.g2=g2
        self.do_shear=not ((g1==0) and (g2==0))
        self.sp2twice=2*sigma_pr**2
        self.sm=sm

        ## generate lookup table so that we can sample from 
        ## prob given in eq 18 of Bernstein and Armstrong:

        x=linspace(0,1.0,100)
        y=array([quad(self.Eprior,0,emax)[0] for emax in x])
        y/=y[-1]
        
        self.genI=InterpolatedUnivariateSpline(y,x) # This is to invert the CDF right? Just making sure. - Mat
        #print self.genI
        #xx=arange(0,1,0.0001)
        #pylab.plot(y,x,'ro')
        #pylab.plot(xx,map(self.genI,xx),'r-')
        #pylab.show()

    ##@autojit
    def Eprior(self,e):
        e2=e**2
        return e*((1-e2)**2.)*exp(-e2/(self.sp2twice))

    ##@autojit
    def generateE(self):
        return self.genI(random.uniform(0.,1.0))


    ##@autojit
    def generate (self):
        E=self.generateE();
        theta=random.uniform(0,2.*math.pi);
        e1=E*cos(theta)
        e2=E*sin(theta)
        #print e1,e2

        ##intrinsic ellipticty now shear
        if (self.do_shear):
            e1,e2=shear([e1,e2],[self.g1,self.g2])
        
        ### now add error
        ok=False
        while not ok:

            e1m=e1+random.gauss(0,self.sm)
            e2m=e2+random.gauss(0,self.sm)
        
            ee=sqrt(e1m**2.+e2m**2.)
            if (ee<1.): ## requier they fall onto <1
                ok=True
            else:
                pass
                #print "Retrying em"
            
        return e1m,e2m

    ##@autojit
    def generate_multishear (self, shearlist):
        E=self.generateE();
        theta=random.uniform(0,2*math.pi);
        e1g=E*cos(theta)
        e2g=E*sin(theta)
        e12g=[e1g,e2g]
        e1m,e2m=0.0,0.0
        toret=[]
        err1=random.gauss(0,self.sm)
        err2=random.gauss(0,self.sm)
        for g1,g2 in shearlist:
            #g1,g2=g12
            e1,e2=shear(e12g,[g1,g2])
            e1m=e1+err1
            e2m=e2+err2
            ee=sqrt(e1m**2.+e2m**2.)
            if (ee<1.): ## requier they fall onto <1
                toret.append((e1m,e2m))                
            else:
                toret.append(None)
                    
        return toret



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
        #random.seed(started*ind)
        #[gg1,gg2,hh1,hh2]=p.draw()
        #[gg1,gg2,hh1,hh2]=[ 0.01017254,  0.02438026,  0.01457512, -0.00170212]
        #[gg1,gg2,hh1,hh2]=[0.003,0.004,0.01,0.01]
        #fini=open('data/gval.ini', 'w')
        #fini.write(','.join(str(j) for j in [gg1,gg2,hh1,hh2]) + '\n')
        #fini.close()

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

    #start=time.time()
    P,Perr=dblquad (fn.IntegrandP,-1.0, 1.0, lambda x:-sqrt(1.-x*x), lambda x:sqrt(1.-x*x),epsabs=0., epsrel=tol)
    Q1,Q1err=dblquad (fn.IntegrandQ1,-1.0, 1.0, lambda x:-sqrt(1.-x*x), lambda x:sqrt(1.-x*x),epsabs=0., epsrel=tol)
    Q2,Q2err=dblquad (fn.IntegrandQ2,-1.0, 1.0, lambda x:-sqrt(1.-x*x), lambda x:sqrt(1.-x*x),epsabs=0., epsrel=tol)
    R11,R11err=dblquad (fn.IntegrandR11,-1.0, 1.0, lambda x:-sqrt(1.-x*x), lambda x:sqrt(1.-x*x),epsabs=0., epsrel=tol)
    R12,R12err=dblquad (fn.IntegrandR12,-1.0, 1.0, lambda x:-sqrt(1.-x*x), lambda x:sqrt(1.-x*x),epsabs=0., epsrel=tol)
    R22,R22err=dblquad (fn.IntegrandR22,-1.0, 1.0, lambda x:-sqrt(1.-x*x), lambda x:sqrt(1.-x*x),epsabs=0., epsrel=tol)
    #end=time.time()

    #print abs(Perr)/P, abs(Q1err)/Q1, abs(Q2err)/Q2
    #print abs(R11err)/R11, abs(R12err)/R12, abs(R22err)/R22
    #print "Time taken is ", end-start, " seconds"
        
    row=[P,Q1,Q2,R11,R12,R22]

    return row


    

if __name__=="__main__":
    main(sys.argv[1:])
