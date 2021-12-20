#
#
#    module for calculating geometry parameters and magnetic aspect
#    angle of radar targets monitored by SANYA
#
#    use aspect_elaz or aspect_txty to calculate aspect angles of targets
#    specified by (el,az) or (tx,ty) angles
#
#  Created by Erhan Kudeki on 11/29/08.
#  Copyright (c) 2008 ECE, UIUC. All rights reserved.
#


from .beam import *
import numpy as np

# ------------ SANYA specifications -------------------------
sanya_specs = RadarSpecs(
    lat0 = 18.34 * deg,  #this is geodetic, the usual map or GPS latitude
    lon0 = 109.62 * deg, #east of Greenwich
    h0 = 0.                        # local height above reference ellipsoid
    )

# Make this variables accesible to the module for backwards compatibility

lat0 = sanya_specs.lat0    # geodetic, the usual map or GPS latitude
lon0 = sanya_specs.lon0    # east of Greenwich
h0   = sanya_specs.h0      # local height above reference ellipsoid

n0 = sanya_specs.n0
x0 = sanya_specs.x0
y0 = sanya_specs.y0
z0 = sanya_specs.z0
xyz0 = sanya_specs.xyz0

# unit vectors
east0 = sanya_specs.east0
zenith0 = sanya_specs.zenith0
north0 = sanya_specs.north0

# radar methods from beam.py
xyz2dec_ha = sanya_specs.xyz2dec_ha
dec_ha2el_az = sanya_specs.dec_ha2el_az
aspect_angle = sanya_specs.aspect_angle
aspect_elaz = sanya_specs.aspect_elaz
cosBs = sanya_specs.cosBs
