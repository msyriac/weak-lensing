#!/usr/bin/env python
from toytoy import *
from scipy.integrate import quad

T=toytoy(0.0)
# Basic identities

#print quad(T.P, -1,1)

def avg (func,T=T):
    v,e=quad(lambda x:T.P(x)*func(x),-1,1)
    #print e
    return v

print avg(lambda x:1) # int P
print avg(lambda x:T.L_ii(x)+T.L_i(x)**2) #<L_ii+L_i L_j>
Fi=1/avg(lambda x:T.L_i(x)**2)

b2=Fi*avg(lambda x:T.L_ii(x)*T.L_i(x))/2.0
b3=Fi*avg(lambda x:T.L_ii(x)*(T.L_ii(x)+T.L_i(x)**2))/6.0
b3a=Fi*avg(lambda x:T.L_i(x)*(T.L_iii(x)+3*T.L_ii(x)*T.L_i(x)+T.L_i(x)**3))/6.0




#print T.L_iii(0.3),'A', T.L_i(0.99), Fi
#stop()

for s in arange(0.001,1.01,0.01):
    E=avg (T.L_i,T=toytoy(s))*Fi
    print s, E, E+s**3*b3a, "EE"

