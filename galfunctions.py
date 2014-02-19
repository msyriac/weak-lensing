from math import *
try:
    from numba import autojit
except:
    pass

def initGalFuncs(spr,se):
    global sp, sm, sm2, sm4, sp2,sp4, sp2twice

    sp=spr
    sp2=sp**2
    sp4=sp2**2
    sp2twice=2*sp2
    sm=se
    sm2=sm**2
    sm4=sm2**2

def updateEs(e1,e2):
    global em1, em2, em12, em22
    em1=e1
    em2=e2
    em12=em1**2
    em22=em2**2
    
def PriorFirstDer1(es1,es2):
    es12=es1**2
    es22=es2**2
    dlta2=es12+es22
    
    return (-(es1*(-1 + es12 + es22)**2*(es12**2 + es22**2 + sp2 + es12*(-1 + 2*es22 - 5*sp2) - es22*(1 + 5*sp2))))*exp(-dlta2/(sp2twice))/(sqrt(dlta2)*sp2)


def PriorFirstDer2(es1,es2):
    es12=es1**2
    es22=es2**2
    dlta2=es12+es22

    return (-(es2*(-1 + es12 + es22)**2*(es12**2 + es22**2 + sp2 + es12*(-1 + 2*es22 - 5*sp2) - es22*(1 + 5*sp2))))*exp(-dlta2/(sp2twice))/(sqrt(dlta2)*sp2)

def PriorSecondDer11(es1,es2):
    es12=es1**2
    es22=es2**2
    dlta2=es12+es22
        
    return (es1**10 + es1**8*(-2 + 4*es22 - 13*sp2) + es22*(1 + es22)*sp2*(es2**4 + sp2 - es22*(1 + 5*sp2)) + es1**6*(1 + 6*es2**4 + 16*sp2 + 30*sp4 - 2*es22*(3 + 19*sp2)) + es12*es22*(es2**6 - 2*es2**4*(1 + 5*sp2) - 2*sp2*(2 + 9*sp2) + es22*(1 + 16*sp2 + 20*sp4)) + es1**4*(4*es2**6 - 6*es2**4*(1 + 6*sp2) - sp2*(3 + 14*sp2) + es22*(2 + 32*sp2 + 55*sp4)))*exp(-dlta2/(sp2twice))/(((es12 + es22)**1.5*sp4)/(-1 + es12 + es22)**2)

def PriorSecondDer22(es1, es2):
    dlta2=es1**2+es2**2

    return (es1**8*(es2**2 + sp2) + es1**6*(4*es2**4 - 5*sp**4 - 2*es2**2*(1 + 5*sp2)) + es1**4*(6*es2**6 - sp2*(1 + 4*sp2) - 6*es2**4*(1 + 6*sp2) + es2**2*(1 + 16*sp2 + 20*sp**4)) + es2**4*(es2**6 - es2**4*(2 + 13*sp2) - sp2*(3 + 14*sp2) + es2**2*(1 + 16*sp2 + 30*sp**4)) + es1**2*(4*es2**8 + sp**4 - 2*es2**2*sp2*(2 + 9*sp2) - 2*es2**6*(3 + 19*sp2) + es2**4*(2 + 32*sp2 + 55*sp**4)))*exp(-dlta2/(2*sp2))/(((es1**2 + es2**2)**1.5*sp**4)/(-1 + es1**2 + es2**2)**2)

def PriorSecondDer12(es1,es2):
    dlta2=es1**2+es2**2

    return (es1*es2*(es1**8 + es2**8 - sp**4 + 2*es1**6*(-1 + 2*es2**2 - 7*sp2) - 2*es2**6*(1 + 7*sp2) - 2*es2**2*(sp2 + 5*sp**4) + es2**4*(1 + 16*sp2 + 35*sp**4) + es1**4*(1 + 6*es2**4 + 16*sp2 + 35*sp**4 - 6*es2**2*(1 + 7*sp2)) + 2*es1**2*(2*es2**6 - sp2*(1 + 5*sp2) - 3*es2**4*(1 + 7*sp2) + es2**2*(1 + 16*sp2 + 35*sp**4))))*exp(-dlta2/(2*sp2))/(((es1**2 + es2**2)**1.5*sp**4)/(-1 + es1**2 + es2**2)**2)

#@autojit    
def LikeZeroDer(es1,es2):
    #eo=[es1,es2]
    #eo=shear([es1,es2],gg_)
    dlta2=(em1-es1)**2.+(em2-es2)**2.
    return exp(-dlta2/(2.*sm2))



#@autojit
def LikeFirstDer1(es1,es2):

    dlta2=(em1-es1)**2.+(em2-es2)**2.
    es12=es1**2.

    es22=es2**2.

    return exp(-dlta2/(2.*sm2))*(-em1+em1*es12+es1-es12*es1+2.*em2*es1*es2-em1*es22-es1*es22)/sm2

#@autojit
def LikeFirstDer2(es1,es2):

    dlta2=(em1-es1)**2.+(em2-es2)**2.
    es12=es1**2.
    es22=es2**2.
    return exp(-dlta2/(2.*sm2))*(-em2+em2*es22+es2-es22*es2+2.*em1*es2*es1-em2*es12-es2*es12)/sm2

#@autojit
def LikeSecondDer11(es1,es2):

    dlta2=(em1-es1)**2.+(em2-es2)**2.
    es12=es1**2.
    es22=es2**2.
    return exp(-dlta2/(2.*sm2))*(((em1 - es1)*(-1. + es12) + 2.*em2*es1*es2 - (em1 + es1)*es22)**2. + (2.*em1*es1*(-1. + es12 - 3.*es22) - (-1. + 3.*es12 - es22)*(-1. + es12 - 2.*em2*es2 + es22))*sm2)/sm4
    

#@autojit
def LikeSecondDer22(es1,es2):


    dlta2=(em1-es1)**2.+(em2-es2)**2.
    es12=es1**2.
    es14=es12**2.

    es22=es2**2.
    return exp(-dlta2/(2.*sm2))*((em2*(1. + es12 - es22) + es2*(-1. - 2.*em1*es1 + es12 + es22))**2. + (es14 - 2.*es12*es2*(3.*em2 + es2) + (1. + 2.*em2*es2 - 3.*es22)*(-1. + es22) - 2.*em1*(es1 + es12*es1 - 3.*es1*es22))*sm2)/sm4


#@autojit
def LikeSecondDer12(es1,es2):
    dlta2=(em1-es1)**2.+(em2-es2)**2.
    es12=es1**2.
    es14=es12**2.

    es22=es2**2.
    es24=es22**2.

    return exp(-dlta2/(2.*sm2))*(2.*em12*es1*es2*(-1. + es12 - es22) + es1*(2.*em22*es2*(-1. - es12 + es22) + es2*(-1. + es12 + es22)*(-1. + es12 + es22 - 4.*sm2) + em2*(-1. + 4.*es22 + (es12 - 3.*es22)*(es12 + es22 - 2.*sm2))) + em1*(-(em2*(-1. + es14 - 6.*es12*es22 + es24)) + es2*(-1 - 3.*es14 + es24 - 2.*es22*sm2 + es12*(4. - 2.*es22 + 6.*sm2))))/sm4
    

def Gaussn(es1,es2):
    return exp((-es1**2-es2**2)/(2*sp**2))

def dGaussn1(es1,es2):
    return exp((-es1**2-es2**2)/(2*sp**2))*(-es1/sp**2)


def dGaussn2(es1,es2):
    return exp((-es1**2-es2**2)/(2*sp**2))*(-es2/sp**2)


def dGaussn11(es1,es2):
    return exp((-es1**2-es2**2)/(2*sp**2))*((-sp**2 + es1**2)/sp**4)

def dGaussn22(es1,es2):
    return exp((-es1**2-es2**2)/(2*sp**2))*((-sp**2 + es2**2)/sp**4)


def dGaussn12(es1,es2):
    return exp((-es1**2-es2**2)/(2*sp**2))*((es1*es2)/sp**4)

#@autojit
def IntegrandP(es0,es1):
    return fEprior(es0,es1)*LikeZeroDer(es0,es1)


#@autojit
def IntegrandQ1(es0,es1):
    return fEprior(es0,es1)*LikeFirstDer1(es0,es1)

#@autojit
def IntegrandQ2(es0,es1):
    return fEprior(es0,es1)*LikeFirstDer2(es0,es1)

#@autojit
def IntegrandR11(es0,es1):
    return fEprior(es0,es1)*LikeSecondDer11(es0,es1)

#@autojit
def IntegrandR12(es0,es1):
    return fEprior(es0,es1)*LikeSecondDer12(es0,es1)

#@autojit
def IntegrandR22(es0,es1):
    return fEprior(es0,es1)*LikeSecondDer22(es0,es1)



#Return prior as a function of |es|
#@autojit
def Eprior(e):
    e2=e**2
    return e*((1-e2)**2.)*exp(-e2/(sp2twice))

#@autojit
def fEprior(es1,es2):
    
    return Eprior(sqrt(es1**2.+es2**2.))
