
import glob
import csv
from math import *
import matplotlib.pyplot as p

import re
def get_trailing_number(s):
    m = re.search(r'\d+$', s)
    return int(m.group()) if m else None



files = glob.glob("data_fisher/*.csv")

E=[]
F11=[]
eF11=[]
F22=[]
eF22=[]
F12=[]
eF12=[]

for fname in files:
    f = open(fname, 'r')
    reader = csv.reader(f,delimiter=',')
    i = get_trailing_number(fname[:(len(fname)-4)])
    for row in reader:
        vals = [float(r) for r in (rv for rv in row if rv!='')]
    E.append(i*0.005)
    F11.append(vals[0])
    eF11.append(vals[1])
    F22.append(vals[2])
    eF22.append(vals[3])
    F12.append(vals[4])
    eF12.append(vals[5])

sortedF11=[x for (y,x) in sorted(zip(E,F11))]
sortedF22=[x for (y,x) in sorted(zip(E,F22))]
sortedF12=[x for (y,x) in sorted(zip(E,F12))]
sortedeF11=[x for (y,x) in sorted(zip(E,eF11))]
sortedeF22=[x for (y,x) in sorted(zip(E,eF22))]
sortedeF12=[x for (y,x) in sorted(zip(E,eF12))]
E=sorted(E)
p.clf()
p.errorbar(E,sortedF11,yerr=eF11)
p.errorbar(E,sortedF22,yerr=eF22)
p.errorbar(E,sortedF12,yerr=eF12)
p.show()
