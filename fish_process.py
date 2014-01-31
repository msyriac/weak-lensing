
import glob
import csv
from math import *
import matplotlib.pyplot as pl

import re
def get_trailing_number(s):
    m = re.search(r'\d+$', s)
    return int(m.group()) if m else None

pref=['va','vb','vc']

pl.clf()

for p in pref:
    files = glob.glob("data_fisher/"+p+"/*.csv")

    E=[]
    F11=[]
    eF11=[]

    for fname in files:
        f = open(fname, 'r')
        reader = csv.reader(f,delimiter=',')
        i = get_trailing_number(fname[:(len(fname)-4)])
        for row in reader:
            vals = [float(r) for r in (rv for rv in row if rv!='')]
        E.append(i*0.005)
        F11.append(vals[0])
        eF11.append(vals[1])

    sortedF11=[x for (y,x) in sorted(zip(E,F11))]
    sortedeF11=[x for (y,x) in sorted(zip(E,eF11))]
    E=sorted(E)
    pl.errorbar(E,sortedF11,yerr=eF11)

pl.savefig('fisher.png')#show()
