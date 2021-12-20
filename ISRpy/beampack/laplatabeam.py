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

# ------------ laPlata AMISR specifications -------------------------
laplata_specs = RadarSpecs(
    lat0 = -(33 + 51/60.+57.6/3600.) * deg,
    lon0 = -(58 + 8 / 60.+11.04/3600.) * deg,
    h0=0.                            # local height above reference ellipsoid
    )

# Make this variables accesible to the module for backwards compatibility

lat0 = laplata_specs.lat0    # geodetic, the usual map or GPS latitude
lon0 = laplata_specs.lon0    # east of Greenwich
h0   = laplata_specs.h0      # local height above reference ellipsoid

n0 = laplata_specs.n0
x0 = laplata_specs.x0
y0 = laplata_specs.y0
z0 = laplata_specs.z0
xyz0 = laplata_specs.xyz0

# unit vectors
east0 = laplata_specs.east0
zenith0 = laplata_specs.zenith0
north0 = laplata_specs.north0

# radar methods from beam.py
xyz2dec_ha = laplata_specs.xyz2dec_ha
dec_ha2el_az = laplata_specs.dec_ha2el_az
aspect_angle = laplata_specs.aspect_angle
aspect_elaz = laplata_specs.aspect_elaz
cosBs = laplata_specs.cosBs

