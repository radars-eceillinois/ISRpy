#
#
#    module for calculating geometry parameters and magnetic aspect
#    angle of radar targets monitored by PFISR
#
#    use aspect_elaz or aspect_txty to calculate aspect angles of targets
#    specified by (el,az) or (tx,ty) angles
#
#  Created by Pablo Reyes on 1/4/2022.
#

from . import beamtools

# ------------ PFISR specifications -------------------------
deg = beamtools.deg             # to express angles in degree values
lat0 = 65.12992 * deg #this is geodetic, the usual map or GPS latitude
lon0 = -147.47104 * deg #east of Greenwich
h0 = 0.213                # local height in km above reference ellipsoid

pfisr_model = beamtools.TargetGeometry(lat0, lon0, h0)

# Radar location in ECEF coordinates
xyz0 = pfisr_model.xyz0

# unit vectors
east0 = pfisr_model.east0
zenith0 = pfisr_model.zenith0
north0 = pfisr_model.north0

# radar methods from beamtools.py
#xyz2dec_ha = pfisr_model.xyz2dec_ha
#dec_ha2el_az = pfisr_model.dec_ha2el_az
aspect_angle = pfisr_model.aspect_angle
aspect_elaz = pfisr_model.aspect_elaz
cosBs = pfisr_model.cosBs
elaz2xyz = pfisr_model.elaz2xyz
xyz2elaz = pfisr_model.xyz2elaz
