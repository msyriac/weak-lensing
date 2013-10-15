#!/usr/bin/env python
import galsim
import pylab
from scipy import *

gal_sigma=2.0
gal_flux=1.0
el=0
shear=[0.02,0.02]
theta = 123 * galsim.degrees
g=sqrt(shear[0]**2.+shear[1]**2.)
p=40
stamp_xsize = p
stamp_ysize = p
random_seed = 1738401
pixel_scale = 1.0

profile = galsim.Gaussian(sigma=gal_sigma, flux=gal_flux)
profile2 = galsim.Gaussian(sigma=gal_sigma, flux=gal_flux)
gal=profile.createSheared(e=el, beta=theta)
gal2=profile2

img = galsim.ImageD(stamp_xsize, stamp_ysize,scale=pixel_scale)
img2 = galsim.ImageD(stamp_xsize, stamp_ysize,scale=pixel_scale)
gal.applyLensing(g1=shear[0],g2=shear[1],mu=1./(1.-2.*(g**2.)))
gal2.applyLensing(g1=shear[0],g2=shear[1],mu=1./(1.-2.*(g**2.)))
gal.draw(img)
gal2.draw(img2)

pylab.imshow(img.array-img2.array)
pylab.colorbar()
print img.array-img2.array
pylab.show()



