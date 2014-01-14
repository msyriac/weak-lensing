import galfunctions as gal
import fisher as f
import toy as t
from math import *
from scipy.integrate import dblquad
import matplotlib.pyplot as plt
import numpy as np
import random
from matplotlib import rc
import pickle

rc('font',**{'family':'serif','serif':['Times'], 'size':'18'})
rc('text', usetex=True)   

def main():
    global g1_,g2_,em1_,em2_

    tol=1.49e-06
    sige=0.05
    sigp=0.3

    g=np.arange(0.0,0.3,0.01)#[0.01, 0.02, 0.03, 0.04, 0.05]
    
    gal.initGalFuncs(sigp,sige)
    f.initFisher(sige)

    #emA=[0.1,0.1]
    emlist=[0.05,0.1,0.3]
    thetag=random.uniform(0,2.*pi);

    plt.clf()
    ax = plt.subplot(1,1,1)

    for em in emlist:
        theta=random.uniform(0,2.*pi);
        emA=[em*cos(theta),em*sin(theta)]

        Plist=[]
        Perrlist=[]
    	for gg in g:
    	    #A=t.ToyGenerator(gg1,gg1)
    	    #emA=A.generate()
            g1_=gg*cos(thetag)
            g2_=gg*sin(thetag)

    	    em1_=emA[0]
    	    em2_=emA[1]
    	    P, Perr = dblquad (Integrand,-1.0, 1.0, lambda x:-sqrt(1.-x*x), lambda x:sqrt(1.-x*x),epsabs=tol, epsrel=tol)
    	    Plist.append(P)
    	    Perrlist.append(Perr)
    	    print P
    	
    	
        #plt.errorbar(g,Plist,yerr=Perrlist, label="$|\\mathbf{e}_{m}|$ = "+str(em) )
    	ax.plot(g,Plist,label="$|\\mathbf{e}_{m}|$ = "+str(em) )
    	
    ax.set_xlabel('Shear $|\\mathbf{g}|$')
    ax.set_ylabel('Likelihood $\\mathcal{L}(\\mathbf{e}_{m}|\\mathbf{g})$')
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels,loc=4)
    #plt.gcf().subplots_adjust(left=0.15)
    plt.savefig('liket.png')
    pickle.dump((Plist,Perrlist,g),open('plot.pickle','wb'))
    

def Integrand(es1,es2):
    global g1_,g2_,em1_,em2_
    loglike=f.LogLike(es1,es2,g1_,g2_,em1_,em2_)
    prior=gal.fEprior(es1,es2)
    return exp(loglike)*prior


if __name__=="__main__":
    main()
