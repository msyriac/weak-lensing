#!/usr/bin/env python
from toytoy import *
from scipy.integrate import quad
import scipy.linalg as la

T=toytoy(0.0)
# Basic identities

#print quad(T.P, -1,1)

def avg (func,Tx=T):
    v,e=quad(lambda x:Tx.P(x)*func(x),-1,1)
    #print e
    return v

print avg(lambda x:1) # int P
print avg(lambda x:T.L_ii(x)+T.L_i(x)**2) #<L_ii+L_i L_j>
Fi=1/avg(lambda x:T.L_i(x)**2)

b2=Fi*avg(lambda x:T.L_ii(x)*T.L_i(x))/2.0
b3=Fi*avg(lambda x:T.L_ii(x)*(T.L_ii(x)+T.L_i(x)**2))/6.0
b3a=Fi*avg(lambda x:T.L_i(x)*(T.L_iii(x)+3*T.L_ii(x)*T.L_i(x)+T.L_i(x)**3))/6.0

U=[1,2,3]
U[0]=lambda x:T.L_i(x)
U[1]=lambda x:T.L_ii(x)+T.L_i(x)**2
U[2]=lambda x:T.L_iii(x)+3*T.L_ii(x)*T.L_i(x)+T.L_i(x)**3



#print U[2](0.3), T.L_ii(0.3), T.L_i(0.3)
#stop()
#print avg(U[0]), avg(U[1]), avg(U[2])
#stop()


# for i in [0,1,2]:
#     for j in [0,1,2]:
#         print avg(lambda x:U[i](x)*U[j](x)),
#     print
# stop()


M=zeros((2,2))
for i in [0,2]:
    for j in [0,2]:
        M[i/2,j/2]=avg(lambda x:U[i](x)*U[j](x))

M[:,1]/=6
print M
#stop()
print 1/Fi, b3a/Fi
#stop()

MI=la.inv(M)
print MI,1/M
#stop()

WW=dot(MI,array([1,0]))
print WW, Fi

#print U[0](0.1), U[2](1), T.L_i(1), T.L_ii(1), T.L_iii(1)

#print T.L_iii(0.3),'A', T.L_i(0.99), Fi
#stop()

for s in arange(0.001,1.01,0.01):

    #print s, 
    #print avg(U[2], T=toytoy(s)), M[1,0]*s, M[1,1]*s**3
    #continue


    E=avg (T.L_i,Tx=toytoy(s))*Fi
    E2=avg(U[0],Tx=toytoy(s))*MI[0,0]+avg(U[2],toytoy(s))*MI[0,1]

    print s, E, E+s**3*b3a, "EE",E2
    #print s, E, E+s**3*b3a, "EE",((avg(U[0], toytoy(s)))- M[0,0]*s)/(M[0,1]*s**3),"EF"
    #print s, E, E+s**3*b3a, "EE",((avg(U[2], toytoy(s)))- M[1,0]*s)/(M[1,1]*s**3),"EH"

