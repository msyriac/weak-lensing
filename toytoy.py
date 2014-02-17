#!/usr/bin/env python
from random import *
from scipy import *
from math import *

class toytoy:
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

N=100000
T=toytoy(0.4)
for cc in range(N):
    v=T.sample()

    P,Q,R=T.PQR_unnorm(v)
    QPx.accumulate(Q/P)
    RPx.accumulate(R/P)

    P,Q,R=T.PQR(v)
    QP.accumulate(Q/P)
    RP.accumulate(R/P)


for x in [QPx,RPx,QP,RP]:
    x.print_stat()

