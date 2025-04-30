"""
TNH, 04/2025

Script with function that finds v_min and v_max velocity ranges for Galactic emission for a given set of Galactic l and b coordinates.
Also has a helper function which returns a (not very precise) l and b value from the name of the field. (For more precise values, can code up a way to access cube header instead but this is better in case cube is not accessible).

Credit: adapted from Betsey Adams' deviation velocity code, based on Chapter 2.1 of High-Velocity Clouds by van Woerden, Wakker, Schwarz and de Boer (2005) (https://www.cita.utoronto.ca/~amarchal/pdf/High-Velocity-Clouds.pdf)

Usage:
1. from get_exp_gal_em import get_lb, get_gal_vel
2. l, b = get_lb(<field_name>)
    e.g., l, b = get_lb('S1044+5550')
3. get_gal_vel(l, b)
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, Column
from astropy.io import ascii
from astropy import coordinates as coord
from astropy import units as u
from astropy.coordinates import SkyCoord

def get_lb(field):
    #extracting rough RA, DEC from field name and creating astropy SkyCoord
     if len(field.split('+'))>1:
         field_ra = field.split('+')[0]
         field_dec = field.split('+')[1]
         field_coord = SkyCoord('{hour}h{minute}m00.0s'.format(hour=field_ra[1:3],minute=field_ra[3:5]),
                                '+{d}d{minute}m00.0s'.format(d=field_dec[0:2],minute=field_dec[2:4]))
     else:
         field_ra = field.split('-')[0]
         field_dec = (field.split('-')[1])
         field_coord = SkyCoord('{hour}h{minute}m00.0s'.format(hour=field_ra[1:3],minute=field_ra[3:5]),
                                '-{d}d{minute}m00.0s'.format(d=field_dec[0:2],minute=field_dec[2:4]))
     #transforming to Galactic coordinates
     field_gal_coord = field_coord.transform_to('galactic')
     
     #returning Galactic coordinates
     return field_gal_coord.l, field_gal_coord.b

def get_gal_vel(l, b):
    #checking if l, b are array of values or single value
    #in either case, putting them into numpy arrays
    if hasattr(l, "__len__") & hasattr(b, "__len__"):
        larr = np.array(l.value)
        barr = np.array(b.value)
    else:
        larr = np.array([l.value])
        barr = np.array([b.value])

    #initializing parameters which will be used in the Galactic emission model
    r0 = 8.5 #distance of the Sun to the Galactic Center (kpc)
    rmax=26. #diameter of Milky Way disk (kpc)
    rs=0.5 #galactocentric radius cutoff (kpc) for solid-body rotation
    vr=220. #constant velocity (km/s) assumed for rotation outside R_s
    z1=1 # (thickness of MW disk for R<R0)/2 (kpc)
    z2=3 # (thickness of MW disk st R=3R0)/2 (kpc)
    
    #define arrays for the distance along los plus radius and z at that distance
    #use an arbitrary range for distance to start; will have to calculate max distance and cutoff in next step
    darr=np.arange(0,20,0.1)
    rarr = r0 * np.sqrt((np.cos(barr*np.pi/180.)**2) * ((darr/r0)**2) - (2*np.cos(barr*np.pi/180.)*np.cos(larr*np.pi/180.)*(darr/r0)) +1 )
    zarr=darr*np.sin(barr*np.pi/180.)
    
    #now finding the index to truncate arrays
    #this happens if r>rmax or abs(z)>zmax
    #first calculate zmax based on rarr
    #need to account for warp for points exterior to sun
    rext = (rarr>=rmax)
    zmax=np.full(len(darr),z1)
    zmax[rext] = z1+((z2-z1)*(((rarr[rext]/r0)-1)**2)/4)
    
    #truncating radius array at the distance where the sightline leaves the disk (r>rmax or abs(z)>zmax)
    dmax_mask = ((rarr <= rmax)&(abs(zarr)<=zmax))
    rvals = rarr[dmax_mask]
    
    #now calculating the velocity array based on this
    #assuming we are probing flat part of rotation curve/line of sight will never cross r=0.5
    varr = ((r0/rvals) -1)*vr*np.sin(larr*np.pi/180.)*np.cos(barr*np.pi/180.)
    
    #returning min and max of these velocities as the velocity range for Galactic emission
    return min(varr), max(varr)
    
    
