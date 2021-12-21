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

from . import beamtools

# ------------ iris radar specifications -------------------------
deg = beamtools.deg             # to express angles in degree values
lat0 = (40 + 10 / 60. + 0.61 / 3600.) * deg                    # geodetic, the usual map or GPS latitude
lon0 = -(88 + 9 / 60. + 30.95 / 3600.) * deg                    # east of Greenwich
h0 = 221. / 1000.                            # local height above reference ellipsoid

iris_model = beamtools.TargetGeometry(lat0, lon0, h0)

# Radar location in ECEF coordinates
xyz0 = iris_model.xyz0

# unit vectors
east0 = iris_model.east0
zenith0 = iris_model.zenith0
north0 = iris_model.north0

# radar methods from beamtools
xyz2dec_ha = iris_model.xyz2dec_ha
dec_ha2el_az = iris_model.dec_ha2el_az
aspect_angle = iris_model.aspect_angle
aspect_elaz = iris_model.aspect_elaz
cosBs = iris_model.cosBs

