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



# ------------ KWAJ specifications -------------------------

from .beam import *
import numpy as np

# ------------ jro radar specifications -------------------------
kwaj_specs = RadarSpecs(
    lat0 = 9.3935605 * deg, #this is geodetic, the usual map or GPS latitude
    lon0 = 167.4763514 * deg, #east of Greenwich
    h0 = 0,                # local height above reference ellipsoid
    )
#lat0=9.3935605*deg                        # geodetic, the usual map or GPS latitude
#lon0=167.4763514*deg
# ----IRIS@Roi----
#lat0=(9+23/60.+53.818/3600.)*deg
#lon0=(167+28/60.+9.123/3600.)*deg
#h0=54.21
# ---IRIS@Roi///GPS Shack ---
#lat0=(9+23/60.+53.822/3600.)*deg
#lon0=(167+28/60.+9.128/3600.)*deg
#h0=53.91/1000.

# Make this variables accesible to the module for backwards compatibility

lat0 = kwaj_specs.lat0    # geodetic, the usual map or GPS latitude
lon0 = kwaj_specs.lon0    # east of Greenwich
h0   = kwaj_specs.h0      # local height above reference ellipsoid

n0 = kwaj_specs.n0
x0 = kwaj_specs.x0
y0 = kwaj_specs.y0
z0 = kwaj_specs.z0
xyz0 = kwaj_specs.xyz0

# unit vectors
east0 = kwaj_specs.east0
zenith0 = kwaj_specs.zenith0
north0 = kwaj_specs.north0

# radar methods from beam.py
xyz2dec_ha = kwaj_specs.xyz2dec_ha
dec_ha2el_az = kwaj_specs.dec_ha2el_az
aspect_angle = kwaj_specs.aspect_angle
aspect_elaz = kwaj_specs.aspect_elaz
cosBs = kwaj_specs.cosBs

