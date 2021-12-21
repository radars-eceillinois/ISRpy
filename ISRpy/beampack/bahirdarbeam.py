#
#
#    module for calculating geometry parameters and magnetic aspect
#    angle of radar targets monitored by BAHIRDAR
#
#    use aspect_elaz or aspect_txty to calculate aspect angles of targets
#    specified by (el,az) or (tx,ty) angles
#
#  Created by Erhan Kudeki on 11/29/08.
#  Copyright (c) 2008 ECE, UIUC. All rights reserved.
#

from . import beamtools

# ------------ BAHIRDAR radar specifications -------------------------
deg = beamtools.deg             # to express angles in degree values
lat0 = (11. + 34 / 60. + 18.34 / 3600) * deg  # geodetic, the usual map or GPS latitude
lon0 = (37. + 23. / 60 + 39.80 / 3600) * deg  # east of Greenwich
h0 = (5879.0 * 0.3048) / 1000.                # local height above reference ellipsoid
bahirdar_model = beamtools.TargetGeometry(lat0, lon0, h0)

# Radar location in ECEF coordinates
xyz0 = bahirdar_model.xyz0

# unit vectors
east0 = bahirdar_model.east0
zenith0 = bahirdar_model.zenith0
north0 = bahirdar_model.north0

# radar methods from beamtools
xyz2dec_ha = bahirdar_model.xyz2dec_ha
dec_ha2el_az = bahirdar_model.dec_ha2el_az
aspect_angle = bahirdar_model.aspect_angle
aspect_elaz = bahirdar_model.aspect_elaz
cosBs = bahirdar_model.cosBs
