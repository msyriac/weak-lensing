#!/usr/bin/env python
#!/usr/bin/env python

from tabulate_like import *
import pylab, cPickle, toy
from accumulator import *
import scipy.linalg as la

print "starting"
I=cPickle.load(open('grids2/tablike.pickle'))


N=10000
for shear in [0.02]:#0.005,0.01,0.02,0.03,0.04, 0.05,0.07, 0.1, 0.2]:
    T=toy.ToyGenerator(shear,0.00)
    A=Accumulator("estimator")
    FDA=Accumulator("Bern 1st d")
    SDA=Accumulator("Bern 2nd d")
    terr=0.0001#shear*0.01
    cc=0
    while True:
        cc+=1
        for i in range(N):
            e1m,e2m=T.generate()
            fd=I.FD(e1m,e2m)
            sd=I.SD(e1m,e2m)
            E=dot(I.IFisher,fd)
            A.accumulate(E)
            FDA.accumulate(fd)
            SDA.accumulate(sd)
        A.get_stat()
        FDA.get_stat()
        SDA.get_stat()
        print shear, A.mean, A.err, ([shear,0]-A.mean)/A.err, cc*N
        print dot(la.inv(SDA.mean),FDA.mean)

        if all(A.err<terr):
            break
#A.print_stat()

