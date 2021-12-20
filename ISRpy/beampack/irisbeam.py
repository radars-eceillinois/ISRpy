#
#
#    module for calculating geometry parameters and magnetic aspect
#    angle of radar targets monitored by IRIS
#
#    use aspect_elaz or aspect_txty to calculate aspect angles of targets
#    specified by (el,az) or (tx,ty) angles
#
#  Created by Erhan Kudeki on 11/29/08.
#  Copyright (c) 2008 ECE, UIUC. All rights reserved.
#


# ------------ IRIS radar specifications -------------------------
from .beam import *
import numpy as np

# ------------ iris radar specifications -------------------------
iris_specs = RadarSpecs(
    lat0 = (40 + 10 / 60. + 0.61 / 3600.) * deg                    # geodetic, the usual map or GPS latitude
    lon0 = -(88 + 9 / 60. + 30.95 / 3600.) * deg                    # east of Greenwich
    h0 = 221. / 1000.                            # local height above reference ellipsoid
    )

# Make this variables accesible to the module for backwards compatibility

lat0 = iris_specs.lat0    # geodetic, the usual map or GPS latitude
lon0 = iris_specs.lon0    # east of Greenwich
h0   = iris_specs.h0      # local height above reference ellipsoid

n0 = iris_specs.n0
x0 = iris_specs.x0
y0 = iris_specs.y0
z0 = iris_specs.z0
xyz0 = iris_specs.xyz0

# unit vectors
east0 = iris_specs.east0
zenith0 = iris_specs.zenith0
north0 = iris_specs.north0

# radar methods from beam.py
xyz2dec_ha = iris_specs.xyz2dec_ha
dec_ha2el_az = iris_specs.dec_ha2el_az
aspect_angle = iris_specs.aspect_angle
aspect_elaz = iris_specs.aspect_elaz
cosBs = iris_specs.cosBs

