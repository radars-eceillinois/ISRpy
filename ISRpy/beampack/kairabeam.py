# -*- coding: utf-8 -*-
#
#
#    module for calculating geometry parameters and magnetic aspect
#    angle of radar targets monitored by KAIRA VHF at Kilpisjärvi - Finland.
#
#    use aspect_elaz or aspect_txty to calculate aspect angles of targets
#    specified by (el,az) or (tx,ty) angles
#
#  Created by Erhan Kudeki on 11/29/08.
#  Copyright (c) 2008 ECE, UIUC. All rights reserved.
#

from . import beamtools

# ------------ KAIRA specifications -------------------------
deg = beamtools.deg             # to express angles in degree values
lat0 = 69.0707445 * deg
lon0 = 20.7620758 * deg
h0 = 493. / 1000.                            # local height above reference ellipsoid

kaira_model = beamtools.TargetGeometry(lat0, lon0, h0)

# Radar location in ECEF coordinates
xyz0 = kaira_model.xyz0

# unit vectors
east0 = kaira_model.east0
zenith0 = kaira_model.zenith0
north0 = kaira_model.north0

# radar methods from beamtools
xyz2dec_ha = kaira_model.xyz2dec_ha
dec_ha2el_az = kaira_model.dec_ha2el_az
aspect_angle = kaira_model.aspect_angle
aspect_elaz = kaira_model.aspect_elaz
cosBs = kaira_model.cosBs
