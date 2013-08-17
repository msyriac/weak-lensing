
import galsim
import matplotlib.pyplot as plt
import numpy as np
from math import *
#import timeit
import sys, getopt




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
    
    try:
        opts, args = getopt.getopt(argv,"PIhn:p:1:2:e:s:u:")
    except getopt.GetoptError:
        print "My programmer is stupid, so I don't understand the arguments you passed. Run with -h to see available options."
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'python ksb.py -n <number of galaxies> -p <stamp pixels> -1 <shear component 1> -2 <shear component 2> -e <ellipticity> -s <galaxy sigma in pixels>'
            print 'Add -P if you want to see a plot of estimated shear as a function of n, -I if you want to see an image of the last galaxy -- don\'t use both together'
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

    estimate(n,p,g1,g2,el,gal_sigma,show_last_gal,plot_progression,update_freq)


def estimate(n,p,g1,g2,el,s,show_last_gal,plot_progression,update_freq):


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

    px=[]
    py=[]
    py2=[]
    py3=[]
    py4=[]

    
    img = galsim.ImageD(stamp_xsize, stamp_ysize,scale=pixel_scale)
    profile = galsim.Gaussian(sigma=gal_sigma, flux=gal_flux)

    first_in_pair=True

    storeXY(stamp_xsize,stamp_ysize)
    #storeXY(3,3)
    #print XX
    #print YY
    #print XY

    
    #start = timeit.default_timer()
    for i in range(0,n):

        ud = galsim.UniformDeviate(random_seed+i)
        if first_in_pair:
            theta = ud() * 2. * pi * galsim.radians
            first_in_pair = False
        else:
            theta += 90 * galsim.degrees
            first_in_pair = True


    
        
        gal=profile.createSheared(e=el, beta=theta)
    
    
        g=sqrt(shear[0]**2.+shear[1]**2.)
        gal.applyLensing(g1=shear[0],g2=shear[1],mu=1./(1.-2.*(g**2.)))
    
    
        #psf = galsim.Gaussian(sigma=1.0)
        pix = galsim.Pixel(pixel_scale)
        final=galsim.Convolve([gal,pix])
        final.draw(img)

        
        #myarray=np.array([[0.,1.,0.],[1.,1.,1.],[0.,1.,0.]])
        #print img.array
        q = getQuad(img.array,*np.shape(img.array))
        E=polE(q)
        #print img.array
        #print q
        #print E
    
    
        if ((i % update_freq) is 0) and (i>0):
            print "Done " + str(i) + " galaxies..."
        

        py.append(0.5*E[0])
        py2.append(0.5*E[1])
        if plot_progression:
            px.append(i+1)
            py3.append(np.array(py).mean())
            py4.append(np.array(py2).mean())
    
        
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

    p=np.array(py)
    inf=p.mean()
    print "Inferred shear1 is", inf, "+/-", p.std()/sqrt(len(p)), "off by ", (inf-shear[0])*100/shear[0], " %"
    p=np.array(py2)
    inf=p.mean()
    print "Inferred shear2 is", inf, "+/-", p.std()/sqrt(len(p)), "off by ", (inf-shear[1])*100/shear[1], " %\n"

    if show_last_gal:
        plt.imshow(img.array, interpolation='nearest')
    if show_last_gal or plot_progression:
        plt.show()


if __name__=="__main__":
    main(sys.argv[1:])
