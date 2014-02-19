#!/usr/bin/env python
#!/usr/bin/env python

from tabulate_like import *
import pylab, cPickle
print "starting"
I=cPickle.load(open('grids/tablike.pickle'))

print I.P.sum()
print (I.P*I.Q**2).sum()

#print I.Fisher()
