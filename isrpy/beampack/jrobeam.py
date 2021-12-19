"""
#   jrobeam.py
#
#    module for calculating geometry parameters and magnetic aspect
#    angle of radar targets monitored by jro
#
#    use aspect_elaz or aspect_txty to calculate aspect angles of targets
#    specified by (el,az) or (tx,ty) angles
#
#  Created by Erhan Kudeki on 11/29/08.
#  Copyright (c) 2008 ECE, UIUC. All rights reserved.
#
"""

from .beam import *
import numpy as np



def aspect_txty(year,rr,tx,ty):
    """
    # returns magnetic aspect angle and geocentric coordinates of a target tracked by jro at
    # range rr (km)
    # tx along jro building
    # ty into the building
    """

    tz=sqrt(1-tx**2.-ty**2.)
    xyz=xyz0+rr*(tx*ux+ty*uy+tz*uo)            #geocentric coordinates of target

    [r,lat,lon,aspect]=aspect_angle(year,xyz)
    [dec,ha]=xyz2dec_ha(xyz-xyz0)
    return r,lon,lat,dec,ha,aspect

def aspect_elaz(year,rr,el,az):
    """
    # returns magnetic aspect angle and geocentric coordinates of a target tracked by jro at
    # range       rr (km)
    # elevation   el (rad above local tangent plane to ellipsoid)
    # azimuth     az (rad east of local north)
    """

    tx=cos(el)*sin(az)                    # direction cosines wrt east and north
    ty=cos(el)*cos(az)
    tz=sin(el)
    xyz=xyz0+rr*(tx*east0+ty*north0+tz*zenith0)        #geocentric coordinates of target

    [r,lat,lon,aspect]=aspect_angle(year,xyz)
    [dec,ha]=xyz2dec_ha(xyz-xyz0)
    return r,lon,lat,dec,ha,aspect

def cosBs(year,rr,el,az):
    # decomposes the radial unit vector to the target to direction cosines of magnetic North, East, and Up

    tx=cos(el)*sin(az)                              # direction cosines wrt east and north
    ty=cos(el)*cos(az)
    tz=sin(el)
    xyz=xyz0+rr*(tx*east0+ty*north0+tz*zenith0)     # target vector
    r=sqrt(dot(xyz,xyz))
    lat,lon,h=xyz2llh(xyz[0],xyz[1],xyz[2])         # target lat, lon, height

    radial=xyz/r; # unit vector to target
    p=sqrt(xyz[0]**2+xyz[1]**2)
    east=array([-xyz[1],xyz[0],0])/p # unit vector to east from target
    north=-cross(east,radial) # unit vector to north from target
    rr_=xyz-xyz0 # vector from radar to target
    rr_u=rr_/sqrt(dot(rr_,rr_)) # unit vector from radar to target

    [bX,bY,bZ,bB]=igrf.igrf_B(year,r-a_igrf,lon/deg,lat/deg)
    bfield=array([bX,bY,bZ])
    B=bX*north+bY*east-bZ*radial # magnetic field vector B
    bn=B/sqrt(dot(B,B)) # "magnetic north" unit vector since B points by definition in "magnetic north" direction
    be=cross(bn,radial)
    be=be/sqrt(dot(be,be)) # magnetic east unit vector
    bu=cross(be,bn) # magnetic up unit vector

    cosBn=dot(bn,rr_u) # magnetic north direction-cosine of rr_u
    aspect_angle=arccos(cosBn)
    cosBe=dot(be,rr_u) # magnetic east direction-cosine of rr_u
    cosBu=dot(bu,rr_u) # magnetic up direction-cosine of rr_u

    """
    uLOS=cosBe*U(h)+cosBn*V(h)+cosBu*W(h) ... LOS wind model in terms of wind components to calculate and direction cosines
    """

    return r,lat,lon,h,xyz,B,aspect,cosBn,cosBe,cosBu

# --------------------------------------------------------------
#from pylab import *  # bad idea... name space clobbered
#from pyigrf import igrf

# ------------ jro radar specifications -------------------------
jrospecs = RadarSpecs(
    lat0 = -11.947917 * deg,   # geodetic, the usual map or GPS latitude
    lon0 = -76.872306 * deg,   # east of Greenwich
    h0 = 0.463,                # local height above reference ellipsoid
    dec = -12.88 * deg,        # antenna declination
    ha = -(4+37./60.) * deg    # hour angle on-axis direction at JRO
    )

# Make this variables accesible to the module for backwards compatibility

lat0 = jrospecs.lat0    # geodetic, the usual map or GPS latitude
lon0 = jrospecs.lon0    # east of Greenwich
h0   = jrospecs.h0      # local height above reference ellipsoid
dec = jrospecs.dec      # antenna declination
ha  = jrospecs.ha       # hour angle on-axis direction at JRO

n0 = jrospecs.n0
x0 = jrospecs.x0
y0 = jrospecs.y0
z0 = jrospecs.z0
xyz0 = jrospecs.xyz0

# unit vectors from jro
east0 = jrospecs.east0
zenith0 = jrospecs.zenith0
north0 = jrospecs.north0

# radar methods from beam.py
xyz2dec_ha = jrospecs.xyz2dec_ha
aspect_angle = jrospecs.aspect_angle

# orthonormal basis vectors including the jro on-axis direction
uo = jrospecs.uo        # on axis

ux = np.cross(zenith0, uo)
ux = ux / np.sqrt(np.dot(ux, ux))  # along the building to the right
uy = np.cross(uo, ux)              # away from the building into the valley
