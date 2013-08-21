#!/usr/bin/env python
import galsim
import matplotlib.pyplot as plt
import numpy as np
from math import *
#import timeit
import sys, getopt
import cmath



def storeXY(nx,ny):
    global XX, XY, YY

    XX=np.zeros((nx,ny))
    XY=np.zeros((nx,ny))
    YY=np.zeros((nx,ny))


    for i in range(0,nx):
        x=0.5+i-(nx)/2.0
        for j in range(0,ny):
            y=0.5+j-(ny)/2.0
            XX[i,j]=x*x
            XY[i,j]=x*y
            YY[i,j]=y*y
            
            
    

def getQuad(img, nx, ny): # returns the (unweighted) quadrupole matrix of an (nx x ny) img 
    global XX, XY, YY
    
    quad=np.array([[0.0,0.0],[0.0,0.0]])

    quad[1][0]=(img*XY).sum()
    quad[0][0]=(img*YY).sum() ##### WARNING: X AND Y HAVE BEEN SWAPPED TO
    quad[1][1]=(img*XX).sum() ##### ACCOUNT FOR NUMPY BEING (Y,X)
    quad[0][1]=quad[1][0]
    
    return quad
                
def polE(quad): # returns the KSB "polarization parameters" defined in KSB Eq 3.2
    e=np.array([[0.0],[0.0]])
    q1=quad[0][0]-quad[1][1]
    q2=2.0*quad[1][0]
    T=(quad[0][0]+quad[1][1]) #/2.
    #T+=2.*sqrt(np.linalg.det(quad))/2.
    e[0]=q1/T
    e[1]=q2/T
    return e



def main(argv):

    #DEFAULTS
    n=1000
    p=20
    g1=0.02
    g2=-0.02
    el=0.3
    show_last_gal=False
    plot_progression=False
    update_freq=100
    gal_sigma=2.0
    make_own=False
    M=1
    
    try:
        opts, args = getopt.getopt(argv,"oPIhn:p:1:2:e:s:u:M:")
    except getopt.GetoptError:
        print "My programmer is stupid, so I don't understand the arguments you passed. Run with -h to see available options."
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'python ksb.py -n <number of galaxies> -p <stamp pixels> -1 <shear component 1> -2 <shear component 2> -e <ellipticity> -s <galaxy sigma in pixels> -M <number of rotations>'
            print 'Add -P if you want to see a plot of estimated shear as a function of n, -I if you want to see an image of the last galaxy -- don\'t use both together. -o if you want to make your own ellipse'
            sys.exit()
        elif opt == "-n":
            n=int(arg)
        elif opt == "-p":
            p=int(arg)
        elif opt == "-1":
            g1=float(arg)
        elif opt == "-2":
            g2=float(arg)
        elif opt == "-e":
            el=float(arg)
        elif opt == "-s":
            gal_sigma=float(arg)
        elif opt == "-u":
            update_freq=int(arg)
        elif opt == "-P":
            plot_progression=True
        elif opt == "-I":
            show_last_gal=True
        elif opt == "-M":
            M=int(arg)
        elif opt == "-o":
            make_own=True

    estimate(n,p,g1,g2,el,gal_sigma,show_last_gal,plot_progression,update_freq,M,make_own)


def estimate(n,p,g1,g2,el,s,show_last_gal,plot_progression,update_freq,M,make_own):


    shear=[]
    shear.append(g1) # lensing shear to apply
    shear.append(g2)
    stamp_xsize = p
    stamp_ysize = p
    random_seed = 1738401
    pixel_scale = 1.0
    sky_level = 1.e6
    gal_sigma=s
    gal_flux=1.0
    ###########
    
    print "\nPixel width and height is", p
    print "Applied shear is", shear
    print "Galaxy sigma in pixels", s
    print "Ellipticity is", el
    print "Number of galaxies is", n, '\n'


    #gi=(sqrt(1+el)-sqrt(1-el))/(sqrt(1+el)+sqrt(1-el)) #shear due to intrinsic ellipticity

 
    ppy=[]
    ppy2=[]
    


    storeXY(stamp_xsize,stamp_ysize)
    #storeXY(3,3)
    #print XX
    #print YY
    #print XY

    
    #start = timeit.default_timer()
    for t in range(0,M):

        phi = t * 2. * pi  / M
        #print "Rotated "+ str(phi) + " radians"
        py=[]
        py2=[]
        py3=[]
        py4=[]
        px=[]


        first_in_pair=True
        img = galsim.ImageD(stamp_xsize, stamp_ysize,scale=pixel_scale)
        profile = galsim.Gaussian(sigma=gal_sigma, flux=gal_flux)


        for i in range(0,n):

            ud = galsim.UniformDeviate(random_seed+i)
            if first_in_pair:
                theta = ud() * 2. * pi * galsim.radians
                first_in_pair = False
            else:
                theta += 90 * galsim.degrees
                first_in_pair = True


    
        
            gal=profile.createSheared(e=el, beta=theta)
    
            s1=(shear[0]*cos(phi)-shear[1]*sin(phi))
            s2=(shear[0]*sin(phi)+shear[1]*cos(phi))
            g=sqrt(s1**2.+s2**2.)
                
            gal.applyLensing(g1=s1,g2=s2,mu=1./(1.-2.*(g**2.)))
    
    
            #psf = galsim.Gaussian(sigma=1.0)
            pix = galsim.Pixel(pixel_scale)
            final=galsim.Convolve([gal,pix])
            final.draw(img)

            if make_own==False:
                imga=img.array
            else:
                imga=drawEllipse(s1,s2,stamp_xsize,stamp_ysize,gal_sigma)

            q = getQuad(imga,*np.shape(imga))
            E=polE(q)
            
    
    
            if ((i % update_freq) is 0) and (i>0):
                print "Done " + str(i) + " galaxies..."
        

            py.append(0.5*E[0])
            py2.append(0.5*E[1])
            if plot_progression:
                px.append(i+1)
                py3.append(np.array(py).mean())
                py4.append(np.array(py2).mean())

        ppy.append(py)
        ppy2.append(py2)
    
        
    #stop = timeit.default_timer()
    #print "\nAverage time spent on one galaxy: ", (stop-start)/(i+1)
   
    if plot_progression:
        plt.clf()
        ax = plt.subplot(1,1,1)

        ax.set_xlabel("Number of galaxies")
        ax.set_ylabel("Inferred shear")
        ax.plot(px,py3,'ro')
        ax.plot(px,py4,'bx')
        ax.axhline(y=shear[0])
        ax.axhline(y=shear[1])
    


    print "\nApplied shear is " + str(shear)


    nnpy=[]
    nnpy2=[]
    

    # THIS SNIPPET UNDOES THE ROTATION ON THE SHEAR
    t=0
    for lpy, lpy2 in zip(ppy, ppy2):
        phi = t * 2. * pi  / M
        npy=[]
        npy2=[]
    
        for s1, s2 in zip(lpy, lpy2):
            sh1=(s1*cos(-phi)-s2*sin(-phi)) # NOTE THE NEGATIVE PHI
            sh2=(s1*sin(-phi)+s2*cos(-phi))
            npy.append(sh1)
            npy2.append(sh2)
    
        nnpy.append(npy)
        nnpy2.append(npy2)
        
        t=t+1            
    ##########################################
    
    shea1=[]
    shea2=[]
    k=1
    for npy, npy2 in zip(nnpy,nnpy2):
        #print "Rotation ", k
        k=k+1
        for i, pp in enumerate([npy,npy2]):
            ## we are doing this clever rotation
            ## so we need to reduce by 2
            pp=[ (pp[2*c]+pp[2*c+1])/2.0 for c in range(len(pp)/2)]
            p=np.array(pp)
            inf=p.mean()
            sh=shear[i]
            err=p.std()/sqrt(len(p))
        
            #print "Inferred shear"+str(i+1)+" is", inf, "+/-", err, "off by ", (inf-sh)*100/sh, " %", (inf-sh)/err, " sigma"
            if i%2==0:
                shea1.append(inf)
            else:
                shea2.append(inf)

    print "Final mean shear1 is ", np.array(shea1).mean()
    print "Final mean shear2 is ", np.array(shea2).mean()




    if show_last_gal:
        plt.imshow(img.array, interpolation='nearest')
    if show_last_gal or plot_progression:
        plt.show()


def save_gals(n, stamp_xsize, stamp_ysize, gal_sigma, sh1, sh2):
    pixel_scale=1.0
    img = galsim.ImageD(stamp_xsize, stamp_ysize,scale=pixel_scale)

    for i in range(0,n):

        theta = i * 2. * pi  / n
        

        gal = galsim.Gaussian(sigma=gal_sigma, flux=1.0)
    
        
        #gal=profile.createSheared(e=el, beta=theta)
    
        s1=(sh1*cos(theta)-sh2*sin(theta))
        s2=(sh1*sin(theta)+sh2*cos(theta))
        g=sqrt(s1**2.+s2**2.)
        
        gal.applyLensing(g1=s1,g2=s2,mu=1./(1.-2.*(g**2.)))


        pix = galsim.Pixel(pixel_scale)
        final=galsim.Convolve([gal,pix])
        final.draw(img)
        plt.clf()
        plt.imshow(img.array, interpolation='nearest')
        plt.savefig('gal'+str(i)+'.png')
        print "Saved galaxy " + str(i)


def drawEllipse(s1,s2, nx, ny,sigma):
    z=s1*2+s2*2*1j
    phi=-cmath.phase(z)/2.
    r=abs(z)
    boa=sqrt((1.-r)/(1.+r))
    img=np.zeros((nx,ny))
    a=sigma
    b=boa*a
    for i in range(0,nx):
        xp=0.5+i-(nx)/2.0        
        for j in range(0,ny):
            yp=0.5+j-(ny)/2.0
            x=cos(phi)*xp-sin(phi)*yp
            y=sin(phi)*xp+cos(phi)*yp
            img[j,i]=exp(-((x*x)/(a*a))-((y*y)/(b*b)))
    return img

    
if __name__=="__main__":
    main(sys.argv[1:])
    #save_gals(20,20,20,2.0,0.3, -0.1)
