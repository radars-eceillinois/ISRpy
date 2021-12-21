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

from . import beamtools

# ------------ KWAJ specifications -------------------------
deg = beamtools.deg             # to express angles in degree values
lat0 = 9.3935605 * deg #this is geodetic, the usual map or GPS latitude
lon0 = 167.4763514 * deg #east of Greenwich
h0 = 0                # local height above reference ellipsoid

kwaj_model = beamtools.TargetGeometry(lat0, lon0, h0)

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

# Radar location in ECEF coordinates
xyz0 = kwaj_model.xyz0

# unit vectors
east0 = kwaj_model.east0
zenith0 = kwaj_model.zenith0
north0 = kwaj_model.north0

# radar methods from beam.py
xyz2dec_ha = kwaj_model.xyz2dec_ha
dec_ha2el_az = kwaj_model.dec_ha2el_az
aspect_angle = kwaj_model.aspect_angle
aspect_elaz = kwaj_model.aspect_elaz
cosBs = kwaj_model.cosBs

