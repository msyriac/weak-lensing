#!/usr/bin/env python

from tabulate_like import *
import sys
print "starting"
G=Gridder(200)
G.load_many(sys.argv[1]+'/grid')
G.dump_to_file(sys.argv[1]+'/combined.pickle')

