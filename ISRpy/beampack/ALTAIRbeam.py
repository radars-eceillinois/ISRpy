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

from . import beamtools

# ------------ ALTAIR radar specifications -------------------------
deg = beamtools.deg             # to express angles in degree values
lat0 = 9.39541666667 * deg # this is geodetic, the usual map/GPS latitude
lon0 = 167.479333333 * deg # east of Greenwich
h0 = 0.012 # local height above reference ellipsoid

altair_model = beamtools.TargetGeometry(lat0, lon0, h0)

# Radar location in ECEF coordinates
xyz0 = altair_model.xyz0

# unit vectors
east0 = altair_model.east0
zenith0 = altair_model.zenith0
north0 = altair_model.north0

# radar methods from beamtools
xyz2dec_ha = altair_model.xyz2dec_ha
dec_ha2el_az = altair_model.dec_ha2el_az
aspect_angle = altair_model.aspect_angle
aspect_elaz = altair_model.aspect_elaz
cosBs = altair_model.cosBs
