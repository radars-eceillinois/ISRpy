#
#
#    module for calculating geometry parameters and magnetic aspect
#    angle of radar targets monitored by Arecibo Observatory
#
#    use aspect_elaz or aspect_txty to calculate aspect angles of targets
#    specified by (el,az) or (tx,ty) angles
#
#  Created by Erhan Kudeki on 11/29/08.
#  Copyright (c) 2008 ECE, UIUC. All rights reserved.
#

from . import beamtools

# ------------ AO specifications -------------------------
deg = beamtools.deg             # to express angles in degree values
lat0 = 18.3464 * deg #this is geodetic, the usual map or GPS latitude
lon0 = -66.7528 * deg #east of Greenwich
h0 = 0                # local height above reference ellipsoid

ao_model = beamtools.TargetGeometry(lat0, lon0, h0)

# Radar location in ECEF coordinates
xyz0 = ao_model.xyz0

# unit vectors
east0 = ao_model.east0
zenith0 = ao_model.zenith0
north0 = ao_model.north0

# radar methods from beamtools
xyz2dec_ha = ao_model.xyz2dec_ha
dec_ha2el_az = ao_model.dec_ha2el_az
aspect_angle = ao_model.aspect_angle
aspect_elaz = ao_model.aspect_elaz
cosBs = ao_model.cosBs
