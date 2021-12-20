#
#
#    module for calculating geometry parameters and magnetic aspect
#    angle of radar targets monitored by ALTAIR
#
#    use aspect_elaz or aspect_txty to calculate aspect angles of targets
#    specified by (el,az) or (tx,ty) angles
#
#  Created by Erhan Kudeki on 11/29/08.
#  Copyright (c) 2008 ECE, UIUC. All rights reserved.
#


# --------------------------------------------------------------
from .beam import *
import numpy as np

# ------------ ALTAIR radar specifications -------------------------
altair_specs = RadarSpecs(
    lat0 = 9.39541666667 * deg, # this is geodetic, the usual map/GPS latitude
    lon0 = 167.479333333 * deg, # east of Greenwich
    h0 = 0.012 # local height above reference ellipsoid
    )

lat0 = altair_specs.lat0    # geodetic, the usual map or GPS latitude
lon0 = altair_specs.lon0    # east of Greenwich
h0   = altair_specs.h0      # local height above reference ellipsoid

n0 = altair_specs.n0
x0 = altair_specs.x0
y0 = altair_specs.y0
z0 = altair_specs.z0
xyz0 = altair_specs.xyz0

# unit vectors from jro
east0 = altair_specs.east0
zenith0 = altair_specs.zenith0
north0 = altair_specs.north0

# radar methods from beam.py
xyz2dec_ha = altair_specs.xyz2dec_ha
dec_ha2el_az = altair_specs.dec_ha2el_az
aspect_angle = altair_specs.aspect_angle
aspect_elaz = altair_specs.aspect_elaz
cosBs = altair_specs.cosBs
# ------------ ALTAIR specifications -------------------------

n0=a_WGS/sqrt(1-flatness*(2-flatness)*sin(lat0)**2.)
x0=(n0+h0)*cos(lat0)*cos(lon0)            # cartesian geocentric coordinates wrt Greenwich
y0=(n0+h0)*cos(lat0)*sin(lon0)
z0=(n0*(1-eccentricity**2)+h0)*sin(lat0)
xyz0=array([x0,y0,z0])
xy0=array([x0,y0])
r0=sqrt(dot(xyz0,xyz0))
p0=sqrt(dot(xy0,xy0))

# unit vectors from ALTAIR
east0=array([-y0,x0,0])/p0                # zenith and north directions wrt local ellipsoid
zenith0=array([cos(lat0)*cos(lon0),cos(lat0)*sin(lon0),sin(lat0)])
north0=cross(zenith0,east0)
