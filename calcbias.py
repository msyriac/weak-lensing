import pickle

g1=0.06
g2=-0.07

g=[g1,g2]

b=[0,1]

G=pickle.load(open("G.pickle",'r'))
N=pickle.load(open("N.pickle",'r'))
F=pickle.load(open("F.pickle",'r'))

print G
print N
print F

bias=[0.,0.]

for i in b:
    bias[i]=0.
    for j in b:
        for k in b:
            for l in b:
                bias[i] += F[i,j]*N[j,k,l]*g[k]*g[l]
          

print bias
