#!/usr/bin/env python
from random import *
from scipy import *
from math import *

class toytoy:
    ## toy toy is a toy,
    ## where we draw from a gaussain with mean mu, 
    ## but clipped at -1..1
    ### PQR expressions from maple

    def __init__ (self, mu):
        self.mu=mu
        
    def sample(self):
        while True:
            g=gauss(self.mu,1)
            if (abs(g)<1):
                return g
        

    def PQR_unnorm(self, x):
        P=exp(-(x-self.mu)**2/2)
        Q=(x-self.mu)*P
        R=(x*x-1-2*x*self.mu+self.mu**2)
        return P,Q,R

    def PQR(self, x):
        mu=self.mu
        norm=0.5*2**(0.5)*pi**(0.5)*(erf(0.5*2**(0.5)+0.5*mu*2**(0.5))-erf(-0.5*2**(0.5)+0.5*mu*2**(0.5)))
        P=exp(-(x-self.mu)**2/2)/norm
        Q=exp(-0.5*(x-mu)**2)*2**(0.5)*(pi**(0.5)*x*erf(0.5*2**(0.5)+0.5*mu*2**(0.5))-pi**(0.5)*x*erf(-0.5*2**(0.5)+0.5*mu*2**(0.5))-pi**(0.5)*mu*erf(0.5*2**(0.5)+0.5*mu*2**(0.5))+pi**(0.5)*mu*erf(-0.5*2**(0.5)+0.5*mu*2**(0.5))-exp(-0.5*(1+mu)**2)*2**(0.5)+exp(-0.5*(-1+mu)**2)*2**(0.5))/pi/(erf(0.5*2**(0.5)+0.5*mu*2**(0.5))-erf(-0.5*2**(0.5)+0.5*mu*2**(0.5)))**2
        R=exp(-0.5*(x-mu)**2)*2**(0.5)*(2*pi*erf(0.5*2**(0.5)+0.5*mu*2**(0.5))*erf(-0.5*2**(0.5)+0.5*mu*2**(0.5))+pi*x**2*erf(0.5*2**(0.5)+0.5*mu*2**(0.5))**2+pi*x**2*erf(-0.5*2**(0.5)+0.5*mu*2**(0.5))**2+pi*mu**2*erf(0.5*2**(0.5)+0.5*mu*2**(0.5))**2-pi*erf(0.5*2**(0.5)+0.5*mu*2**(0.5))**2-pi*erf(-0.5*2**(0.5)+0.5*mu*2**(0.5))**2-2*pi*x**2*erf(0.5*2**(0.5)+0.5*mu*2**(0.5))*erf(-0.5*2**(0.5)+0.5*mu*2**(0.5))-2*pi*x*mu*erf(0.5*2**(0.5)+0.5*mu*2**(0.5))**2-2*pi*x*mu*erf(-0.5*2**(0.5)+0.5*mu*2**(0.5))**2-2*pi*mu**2*erf(0.5*2**(0.5)+0.5*mu*2**(0.5))*erf(-0.5*2**(0.5)+0.5*mu*2**(0.5))+2**(0.5)*pi**(0.5)*exp(-0.5*(1+mu)**2)*erf(0.5*2**(0.5)+0.5*mu*2**(0.5))-2**(0.5)*pi**(0.5)*exp(-0.5*(1+mu)**2)*erf(-0.5*2**(0.5)+0.5*mu*2**(0.5))+2**(0.5)*pi**(0.5)*exp(-0.5*(-1+mu)**2)*erf(0.5*2**(0.5)+0.5*mu*2**(0.5))-2**(0.5)*pi**(0.5)*exp(-0.5*(-1+mu)**2)*erf(-0.5*2**(0.5)+0.5*mu*2**(0.5))+4*pi*x*mu*erf(0.5*2**(0.5)+0.5*mu*2**(0.5))*erf(-0.5*2**(0.5)+0.5*mu*2**(0.5))+pi*mu**2*erf(-0.5*2**(0.5)+0.5*mu*2**(0.5))**2-2*2**(0.5)*pi**(0.5)*x*exp(-0.5*(1+mu)**2)*erf(0.5*2**(0.5)+0.5*mu*2**(0.5))+2*2**(0.5)*pi**(0.5)*x*exp(-0.5*(1+mu)**2)*erf(-0.5*2**(0.5)+0.5*mu*2**(0.5))+2*2**(0.5)*pi**(0.5)*x*exp(-0.5*(-1+mu)**2)*erf(0.5*2**(0.5)+0.5*mu*2**(0.5))-2*2**(0.5)*pi**(0.5)*x*exp(-0.5*(-1+mu)**2)*erf(-0.5*2**(0.5)+0.5*mu*2**(0.5))+3*2**(0.5)*pi**(0.5)*exp(-0.5*(1+mu)**2)*mu*erf(0.5*2**(0.5)+0.5*mu*2**(0.5))-3*2**(0.5)*pi**(0.5)*exp(-0.5*(1+mu)**2)*mu*erf(-0.5*2**(0.5)+0.5*mu*2**(0.5))-3*2**(0.5)*pi**(0.5)*exp(-0.5*(-1+mu)**2)*mu*erf(0.5*2**(0.5)+0.5*mu*2**(0.5))+3*2**(0.5)*pi**(0.5)*exp(-0.5*(-1+mu)**2)*mu*erf(-0.5*2**(0.5)+0.5*mu*2**(0.5))+4*exp(-(1+mu)**2)+4*exp(-(-1+mu)**2)-8*exp(-1-mu**2))/pi**(3/2)/(erf(0.5*2**(0.5)+0.5*mu*2**(0.5))-erf(-0.5*2**(0.5)+0.5*mu*2**(0.5)))**3
        return P,Q,R



from accumulator import *

QPx=Accumulator("Q/P unnormed")
RPx=Accumulator("R/P unnormed")

QP=Accumulator("Q/P")
RP=Accumulator("R/P")

FishX=Accumulator("FisherX")
Fish=Accumulator("Fisher")

N=100000
T=toytoy(0.0)

for cc in range(N):
    v=T.sample()

    P,Q,R=T.PQR_unnorm(v)
    QPx.accumulate(Q/P)
    RPx.accumulate(R/P)
    FishX.accumulate(Q**2/P**2)

    P,Q,R=T.PQR(v)
    QP.accumulate(Q/P)
    RP.accumulate(R/P)
    Fish.accumulate(Q**2/P**2)
    
for x in [QPx,RPx,QP,RP,Fish,FishX]:
    x.print_stat()


print '------------------------------------'
FishI=1./Fish.x
FishIx=1./FishX.x
## Now let's try estimator
Tt=toytoy(0.1)
E=Accumulator("Estimate,true=0.1")
Ex=Accumulator("Estimate from unnormed,true=0.1")

FD=Accumulator("First loglike der")
SD=Accumulator("Second loglike der")

FDx=Accumulator("First loglike der (unnorm)")
SDx=Accumulator("Second loglike der (unnorm)")

for cc in range(N):
    v=Tt.sample()

    P,Q,R=T.PQR_unnorm(v)
    Ex.accumulate(FishIx*(Q/P))
    FDx.accumulate(Q/P)
    SDx.accumulate(R/P-Q**2/P**2)

    P,Q,R=T.PQR(v)
    E.accumulate(FishI*(Q/P))
    FD.accumulate(Q/P)
    SD.accumulate(R/P-Q**2/P**2)


E.print_stat()
Ex.print_stat()
    
print '--------'
for x in [FD,SD,FDx,SDx]:
    x.get_stat()

print "Bernstein like estimator"
print -FD.mean/SD.mean
print -FDx.mean/SDx.mean
