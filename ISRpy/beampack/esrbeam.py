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

from .beam import *
import numpy as np

# ------------ ESR specifications -------------------------
esr_specs = RadarSpecs(
    lat0 = (78 + 11 / 60. + 14.9 / 3600.) * deg
    lon0 = (16 + 8 / 60. + 25.0 / 3600.) * deg
    h0 = 82. / 1000.                            # local height above reference ellipsoid
    )

# Make this variables accesible to the module for backwards compatibility

lat0 = esr_specs.lat0    # geodetic, the usual map or GPS latitude
lon0 = esr_specs.lon0    # east of Greenwich
h0   = esr_specs.h0      # local height above reference ellipsoid

n0 = esr_specs.n0
x0 = esr_specs.x0
y0 = esr_specs.y0
z0 = esr_specs.z0
xyz0 = esr_specs.xyz0

# unit vectors
east0 = esr_specs.east0
zenith0 = esr_specs.zenith0
north0 = esr_specs.north0

# radar methods from beam.py
xyz2dec_ha = esr_specs.xyz2dec_ha
dec_ha2el_az = esr_specs.dec_ha2el_az
aspect_angle = esr_specs.aspect_angle
aspect_elaz = esr_specs.aspect_elaz
cosBs = esr_specs.cosBs
