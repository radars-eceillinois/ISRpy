#
#
#    module for calculating geometry parameters and magnetic aspect
#    angle of radar targets monitored by SANYA
#
#    use aspect_elaz or aspect_txty to calculate aspect angles of targets
#    specified by (el,az) or (tx,ty) angles
#
#  Created by Erhan Kudeki on 11/29/08.
#  Copyright (c) 2008 ECE, UIUC. All rights reserved.
#

from . import beamtools

# ------------ SANYA specifications -------------------------
deg = beamtools.deg             # to express angles in degree values
lat0 = 18.34 * deg  #this is geodetic, the usual map or GPS latitude
lon0 = 109.62 * deg #east of Greenwich
h0 = 0.                        # local height above reference ellipsoid

sanya_model = beamtools.TargetGeometry(lat0, lon0, h0)

# Radar location in ECEF coordinates
xyz0 = sanya_model.xyz0

# unit vectors
east0 = sanya_model.east0
zenith0 = sanya_model.zenith0
north0 = sanya_model.north0

# radar methods from beamtools
xyz2dec_ha = sanya_model.xyz2dec_ha
dec_ha2el_az = sanya_model.dec_ha2el_az
aspect_angle = sanya_model.aspect_angle
aspect_elaz = sanya_model.aspect_elaz
cosBs = sanya_model.cosBs
