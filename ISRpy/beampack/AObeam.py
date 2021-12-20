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

from .beam import *
import numpy as np


# ------------ AO specifications -------------------------
ao_specs = RadarSpecs(
    lat0 = 18.3464 * deg, #this is geodetic, the usual map or GPS latitude
    lon0 = -66.7528 * deg, #east of Greenwich
    h0 = 0,                # local height above reference ellipsoid
    )

# Make this variables accesible to the module for backwards compatibility

lat0 = ao_specs.lat0    # geodetic, the usual map or GPS latitude
lon0 = ao_specs.lon0    # east of Greenwich
h0   = ao_specs.h0      # local height above reference ellipsoid


n0 = ao_specs.n0
x0 = ao_specs.x0
y0 = ao_specs.y0
z0 = ao_specs.z0
xyz0 = ao_specs.xyz0

# unit vectors from AO
east0 = ao_specs.east0
zenith0 = ao_specs.zenith0
north0 = ao_specs.north0

# radar methods from beam.py
xyz2dec_ha = ao_specs.xyz2dec_ha
dec_ha2el_az = ao_specs.dec_ha2el_az
aspect_angle = ao_specs.aspect_angle
aspect_elaz = ao_specs.aspect_elaz
cosBs = ao_specs.cosBs
