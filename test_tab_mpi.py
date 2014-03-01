#!/usr/bin/env python
#!/usr/bin/env python

from tabulate_like import *
import pylab, cPickle, toy
from accumulator import *
import scipy.linalg as la
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
if (rank==0):
    print "starting ", rank, size
    f=open('tabres.txt','w')

I=cPickle.load(open('grids4/tablike.pickle'))


N=30000

MM=array(I.M13I)
print I.M13I
OO=la.inv(MM)


for shear in [0.1]:#arange(0.005, 0.2, 0.005):
    T=toy.ToyGenerator(shear,0.00)
    A=Accumulator("estimator")
    A2=Accumulator("3rd order estimator")
    FDA=Accumulator("Bern 1st d")
    SDA=Accumulator("Bern 2nd d")
    terr=shear*0.001
    cc=0
    while True:
        cc+=1
        for i in range(N):
            e1m,e2m=T.generate()
            fd=I.FD(e1m,e2m)
            sd=I.SD(e1m,e2m)
            E=dot(I.IFisher,fd)

            e1mx=sqrt(e1m**2+e2m**2)
            e2mx=0
            E1=I.FD(e1mx, e2mx)[0]
            E3=I.E3D(e1mx, e2mx)
            #E13=E1*I.M13I[0,0]+E3*I.M13I[0,1]
            E13=E3/OO[0,1]
            Efx=[e1m,e2m]/e1mx*E13

            A.accumulate(E)
            A2.accumulate(Efx)
            FDA.accumulate(fd)
            SDA.accumulate(sd)
        A.get_stat(mpireduce=comm, tagofs=10)
        A2.get_stat(mpireduce=comm, tagofs=40)
        FDA.get_stat(mpireduce=comm, tagofs=20)
        SDA.get_stat(mpireduce=comm, tagofs=30)
        if (rank==0):
            print shear, A.mean, A.err, A2.mean, A2.err, ([shear,0]-A.mean)/A.err, cc*N
            ber=dot(-la.inv(SDA.mean),FDA.mean)
            print ber
        if all(A.err<terr):
            break
    if (rank==0):
        f.write ("%g %g %g %g %g %g %g %g %g %g %g \n"%(shear, A.mean[0], A.mean[1], A.err[0], A.err[1], ber[0], ber[1], SDA.mean[0,0], SDA.mean[1,1],A2.mean[0], A2.err[0]))

#A.print_stat()
if (rank==0):
    f.close()
