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

# ------------ BAHIRDAR radar specifications -------------------------
bahirdar_specs = RadarSpecs(
    lat0 = (11. + 34 / 60. + 18.34 / 3600) * deg  # geodetic, the usual map or GPS latitude
    lon0 = (37. + 23. / 60 + 39.80 / 3600) * deg  # east of Greenwich
    h0 = (5879.0 * 0.3048) / 1000.                # local height above reference ellipsoid
    )

lat0 = bahirdar_specs.lat0    # geodetic, the usual map or GPS latitude
lon0 = bahirdar_specs.lon0    # east of Greenwich
h0   = bahirdar_specs.h0      # local height above reference ellipsoid

n0 = bahirdar_specs.n0
x0 = bahirdar_specs.x0
y0 = bahirdar_specs.y0
z0 = bahirdar_specs.z0
xyz0 = bahirdar_specs.xyz0

# unit vectors from jro
east0 = bahirdar_specs.east0
zenith0 = bahirdar_specs.zenith0
north0 = bahirdar_specs.north0

# radar methods from beam.py
xyz2dec_ha = bahirdar_specs.xyz2dec_ha
dec_ha2el_az = bahirdar_specs.dec_ha2el_az
aspect_angle = bahirdar_specs.aspect_angle
aspect_elaz = bahirdar_specs.aspect_elaz
cosBs = bahirdar_specs.cosBs
