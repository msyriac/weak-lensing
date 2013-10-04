#!/usr/bin/env python

import glob
import csv
import numpy
import math

files = glob.glob("data/*.csv")

Cinv=numpy.zeros((2,2))
Qsum=0.
appG=[-0.01,0.02]

i=0
k=0
for fname in files:
    f = open(fname, 'r')
    reader = csv.reader(f,delimiter=',')
    for row in reader:
        vals = [float(r) for r in (rv for rv in row if rv!='')]
        if len(vals)!=6:
            print "Invalid line, ", vals
            print fname
            k+=1
        else:
            i+=1
            P=vals[0]
            Q=numpy.array([[vals[1]],
                          [vals[2]]])
            R=numpy.array( [ [vals[3], vals[4]],
                             [vals[4], vals[5]] ] )

            Q_norm=(Q/P)
            Qsum+=Q_norm
            Cinv+=((numpy.dot(Q_norm,Q_norm.transpose()))-(R/P))
        


print i, " valid lines read."
print k, " invalid lines skipped."
Cm=numpy.linalg.inv(Cinv)
infGp=numpy.dot(Cm,Qsum)
infG=[infGp.flat[0],infGp.flat[1]]
err=[sqrt(-Cm[0,0]),sqrt(-Cm[1,1])]
sigma = [(infG[0]-appG[0])/err[0],(infG[1]-appG[1])/err[1]]
#err=[0.,0.]
#sigma=[0.,0.]
bias=[infG[0]-appG[0],infG[1]-appG[1]]
biasp=[(bias[0]*100/appG[0]),(bias[1]*100/appG[1])]
print str((i, infG[0], err[0],  sigma[0],  biasp[0],infG[1], err[1],  sigma[1],  biasp[1])).strip('()')
            
f.close()
