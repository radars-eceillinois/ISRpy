#
#
#    module for calculating geometry parameters and magnetic aspect
#    angle of radar targets monitored by EISCAT Tromso UHF (esr)
#
#    use aspect_elaz or aspect_txty to calculate aspect angles of targets
#    specified by (el,az) or (tx,ty) angles
#
#  Created by Erhan Kudeki on 11/29/08.
#  Copyright (c) 2008 ECE, UIUC. All rights reserved.
#

from . import beamtools

# ------------ Tromso UHF specifications -------------------------
deg = beamtools.deg             # to express angles in degree values
lat0 = (69 + 35/60. + 11.12 / 3600.) * deg
lon0 = (19 + 13/60. + 35.15 / 3600.) * deg
h0   = 70. / 1000.                            # local height above reference ellipsoid

tromso_model = beamtools.TargetGeometry(lat0, lon0, h0)

# Radar location in ECEF coordinates
xyz0 = tromso_model.xyz0

# unit vectors
east0 = tromso_model.east0
zenith0 = tromso_model.zenith0
north0 = tromso_model.north0

# radar methods from beamtools
xyz2dec_ha = tromso_model.xyz2dec_ha
dec_ha2el_az = tromso_model.dec_ha2el_az
aspect_angle = tromso_model.aspect_angle
aspect_elaz = tromso_model.aspect_elaz
cosBs = tromso_model.cosBs
