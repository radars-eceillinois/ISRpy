#
#
#    module for calculating geometry parameters and magnetic aspect
#    angle of radar targets monitored by EISCAT Svalbard Radar (esr)
#
#    use aspect_elaz or aspect_txty to calculate aspect angles of targets
#    specified by (el,az) or (tx,ty) angles
#
#  Created by Erhan Kudeki on 11/29/08.
#  Copyright (c) 2008 ECE, UIUC. All rights reserved.
#

from . import beamtools

# ------------ ESR specifications -------------------------
deg = beamtools.deg             # to express angles in degree values
lat0 = (78 + 11 / 60. + 14.9 / 3600.) * deg
lon0 = (16 + 8 / 60. + 25.0 / 3600.) * deg
h0 = 82. / 1000.                            # local height above reference ellipsoid
esr_model = beamtools.TargetGeometry(lat0, lon0, h0)

# Radar location in ECEF coordinates
xyz0 = esr_model.xyz0

# unit vectors
east0 = esr_model.east0
zenith0 = esr_model.zenith0
north0 = esr_model.north0

# radar methods from beamtools
xyz2dec_ha = esr_model.xyz2dec_ha
dec_ha2el_az = esr_model.dec_ha2el_az
aspect_angle = esr_model.aspect_angle
aspect_elaz = esr_model.aspect_elaz
cosBs = esr_model.cosBs
