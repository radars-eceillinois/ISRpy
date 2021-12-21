#
#
#    module for calculating geometry parameters and magnetic aspect
#    angle of radar targets monitored by laPlata AMISR
#
#    use aspect_elaz or aspect_txty to calculate aspect angles of targets
#    specified by (el,az) or (tx,ty) angles
#
#  Created by Erhan Kudeki on 11/29/08.
#  Copyright (c) 2008 ECE, UIUC. All rights reserved.
#

from . import beamtools

# ------------ laPlata AMISR specifications -------------------------
deg = beamtools.deg             # to express angles in degree values
lat0 = -(33 + 51/60.+57.6/3600.) * deg
lon0 = -(58 + 8 / 60.+11.04/3600.) * deg
h0=0.                            # local height above reference ellipsoid

laplata_model = beamtools.TargetGeometry(lat0, lon0, h0)

# Radar location in ECEF coordinates
xyz0 = laplata_model.xyz0

# unit vectors
east0 = laplata_model.east0
zenith0 = laplata_model.zenith0
north0 = laplata_model.north0

# radar methods from beamtools
xyz2dec_ha = laplata_model.xyz2dec_ha
dec_ha2el_az = laplata_model.dec_ha2el_az
aspect_angle = laplata_model.aspect_angle
aspect_elaz = laplata_model.aspect_elaz
cosBs = laplata_model.cosBs

