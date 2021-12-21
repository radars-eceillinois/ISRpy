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

from . import beamtools
import numpy as np # needed for JRO specific calculations below

# ------------ jro radar specifications -------------------------
deg = beamtools.deg             # to express angles in degree values
lat0 = -11.947917 * deg    # geodetic, the usual map or GPS latitude
lon0 = -76.872306 * deg    # east of Greenwich
h0   = 0.463      # local height above reference ellipsoid

jromodel = beamtools.TargetGeometry(lat0, lon0, h0)

# Radar location in ECEF coordinates
xyz0 = jromodel.xyz0

# unit vectors
east0 = jromodel.east0
zenith0 = jromodel.zenith0
north0 = jromodel.north0

# radar methods from beamtools
xyz2dec_ha = jromodel.xyz2dec_ha
dec_ha2el_az = jromodel.dec_ha2el_az
aspect_angle = jromodel.aspect_angle
aspect_elaz = jromodel.aspect_elaz
cosBs = jromodel.cosBs

# orthonormal basis vectors including the jro on-axis direction
dec = -12.88 * deg        # antenna declination
ha = -(4+37./60.) * deg    # hour angle on-axis direction at JRO

uo = np.array([np.cos(dec) * np.cos(ha/4. + lon0),
               np.cos(dec) * np.sin(ha/4. + lon0),
               np.sin(dec)])    # on axis

ux = np.cross(zenith0, uo)
ux = ux / np.sqrt(np.dot(ux, ux))  # along the building to the right
uy = np.cross(uo, ux)              # away from the building into the valley

def aspect_txty(year, rr, tx, ty):
    """Returns magnetic aspect angle and geocentric coordinates of a target tracked by jro at
    range rr (km)
    tx along jro building
    ty into the building
    """

    tz = np.sqrt(1 - tx ** 2. - ty ** 2.)
    xyz = xyz0 + rr * (tx * ux + ty * uy + tz * uo)
    #geocentric coordinates of target

    [r, lat, lon, aspect] = aspect_angle(year, xyz)
    [dec, ha] = xyz2dec_ha(xyz - xyz0)

    return r,lon,lat,dec,ha,aspect
