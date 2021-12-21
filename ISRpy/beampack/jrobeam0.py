"""
#   jrobeam.py
#
#	module for calculating geometry parameters and magnetic aspect
#	angle of radar targets monitored by jro
#
#	use aspect_elaz or aspect_txty to calculate aspect angles of targets
#	specified by (el,az) or (tx,ty) angles
#
#  Created by Erhan Kudeki on 11/29/08.
#  Copyright (c) 2008 ECE, UIUC. All rights reserved.
#
"""

def llh2xyz(latg,lon,h):
	""" returns geocentric xyz coordinates in km of a target with
	#	latitude       latg (rad) --- geodetic
	#	longitude      lon (rad)
	#	height         h (km above local ellipsoid)
	"""

	n=a_WGS/sqrt(1-flatness*(2-flatness)*sin(latg)**2.)
	x=(n+h)*cos(latg)*cos(lon)      # cartesian geocentric coordinates wrt Greenwich
	y=(n+h)*cos(latg)*sin(lon)
	z=(n*(1-eccentricity**2.)+h)*sin(latg)
	return x,y,z

def xyz2llh(x,y,z):
	"""
	# returns longitude 'lon', geodetic latitude 'lat', and height 'h'
	# of position (x,y,z) defined in geocentric coordinate system
	"""

	p=sqrt(x**2.+y**2.)
	lon=arctan2(y,x)
	lat=arctan2(z,p)
	latp=lat
	for i in range(10):
		n=a_WGS/sqrt(1-flatness*(2-flatness)*sin(latp)**2.)
		h=p/cos(latp)-n
		lat=arctan(z/(p*(1-n*eccentricity**2./(n+h))))
		if abs(lat-latp)<3*eps:
			n=a_WGS/sqrt(1-flatness*(2-flatness)*sin(lat)**2.)
			h=p/cos(lat)-n
			break
		latp=lat
	return lat,lon,h

def dec_ha2el_az(dec,ha):
	"""
	# returns elevation and azimuth angles of a radar beam with respect to local tangent plane.
	# the beam is specified by:
	#		declination dec (deg)
	#		hour angle  ha  (min)
	# with respect to radar location at longitude lon0 and height h0
	# above reference ellipsiod at geodetic latitude lat0
	"""

	lat=dec*deg	                        	# on celestial sphere
	lon=2*pi*(ha/(24*60))
	lon=lon+lon0					# on celestial sphere
	vec=array([cos(lat)*cos(lon),cos(lat)*sin(lon),sin(lat)])
	hor=vec-dot(vec,zenith0)*zenith0
	hor=hor/sqrt(dot(hor,hor))
	el=arccos(dot(hor,vec))/deg
	north=dot(hor,north0)
	east=dot(hor,east0)
	az=arctan2(east,north)/deg
	return el,az

def xyz2dec_ha(vec):
	"""
	# declination and hour angle in target direction used to describe radar beam
	# direction at JRO, corresponding to latitude and relative longitude of the
	# beam-spot on the celestial sphere, corresponds to rr->\infty, in which case:
	"""
	vec=vec/sqrt(dot(vec,vec))
	p=sqrt(vec[0]**2+vec[1]**2)
	dec=arctan2(vec[2],p)/deg                               # in degrees
	ha=(arctan2(vec[1],vec[0])-lon0)*(24/(2.*pi))*60		# in minutes
	return dec,ha

def aspect_angle(year,xyz):
	"""
	# returns the magnetic aspect angle (rad) of a target with
	# geocentric vector xyz defined in geocentric coordinates
	"""

	r=sqrt(dot(xyz,xyz))
	p=sqrt(xyz[0]**2+xyz[1]**2)
	lat=arctan2(xyz[2],p)
	lon=arctan2(xyz[1],xyz[0])
	radial=xyz/r;						# directions from target
	east=array([-xyz[1],xyz[0],0])/p
	north=-cross(east,radial)
	rr=xyz-xyz0
	u_rr=rr/sqrt(dot(rr,rr))			# unit vector from radar to target

	[bX,bY,bZ,bB]=igrf.igrf_B(year,r-a_igrf,lon/deg,lat/deg)
	bfield=array([bX,bY,bZ])
	B=bX*north+bY*east-bZ*radial
	u_B=B/sqrt(dot(B,B))
	aspect=arccos(dot(u_B,u_rr))
	return r,lat,lon,aspect,B

def aspect_txty(year,rr,tx,ty):
	"""
    # returns magnetic aspect angle and geocentric coordinates of a target tracked by jro at
    # range rr (km)
	# tx along jro building
	# ty into the building
	"""

	tz=sqrt(1-tx**2.-ty**2.)
	xyz=xyz0+rr*(tx*ux+ty*uy+tz*uo)			#geocentric coordinates of target

	[r,lat,lon,aspect,B]=aspect_angle(year,xyz)
	[dec,ha]=xyz2dec_ha(xyz-xyz0)
	return r,lon,lat,dec,ha,aspect,B
def aspect_elaz(year,rr,el,az):
	"""
    # returns magnetic aspect angle and geocentric coordinates of a target tracked by jro at
    # range       rr (km)
    # elevation   el (rad above local tangent plane to ellipsoid)
    # azimuth     az (rad east of local north)
	"""

	tx=cos(el)*sin(az)					# direction cosines wrt east and north
	ty=cos(el)*cos(az)
	tz=sin(el)
	xyz=xyz0+rr*(tx*east0+ty*north0+tz*zenith0)		#geocentric coordinates of target

	[r,lat,lon,aspect]=aspect_angle(year,xyz)
	[dec,ha]=xyz2dec_ha(xyz-xyz0)
	return r,lon,lat,dec,ha,aspect

def cosBs(year,rr,el,az):
    # decomposes the radial unit vector to the target to direction cosines of magnetic North, East, and Up

    tx=cos(el)*sin(az)                              # direction cosines wrt east and north
    ty=cos(el)*cos(az)
    tz=sin(el)
    xyz=xyz0+rr*(tx*east0+ty*north0+tz*zenith0)     # target vector
    r=sqrt(dot(xyz,xyz))
    lat,lon,h=xyz2llh(xyz[0],xyz[1],xyz[2])         # target lat, lon, height

    radial=xyz/r; # unit vector to target
    p=sqrt(xyz[0]**2+xyz[1]**2)
    east=array([-xyz[1],xyz[0],0])/p # unit vector to east from target
    north=-cross(east,radial) # unit vector to north from target
    rr_=xyz-xyz0 # vector from radar to target
    rr_u=rr_/sqrt(dot(rr_,rr_)) # unit vector from radar to target

    [bX,bY,bZ,bB]=igrf.igrf_B(year,r-a_igrf,lon/deg,lat/deg)
    bfield=array([bX,bY,bZ])
    B=bX*north+bY*east-bZ*radial # magnetic field vector B
    bn=B/sqrt(dot(B,B)) # "magnetic north" unit vector since B points by definition in "magnetic north" direction
    be=cross(bn,radial)
    be=be/sqrt(dot(be,be)) # magnetic east unit vector
    bu=cross(be,bn) # magnetic up unit vector

    cosBn=dot(bn,rr_u) # magnetic north direction-cosine of rr_u
    aspect_angle=arccos(cosBn)
    cosBe=dot(be,rr_u) # magnetic east direction-cosine of rr_u
    cosBu=dot(bu,rr_u) # magnetic up direction-cosine of rr_u

    """
    uLOS=cosBe*U(h)+cosBn*V(h)+cosBu*W(h) ... LOS wind model in terms of wind components to calculate and direction cosines
    """

    return r,lat,lon,h,xyz,B,aspect,cosBn,cosBe,cosBu

# --------------------------------------------------------------
from pylab import *
from ..pyigrf import pyigrf
igrf = pyigrf() # loading the latest IGRF model

eps=finfo(float).eps					# float resolution
deg=pi/180.								# to express angles in degree values
a_igrf=6371.2							# mean earth radius (km)

a_WGS=6378.137							# equatorial radius WGS 84
flatness=1/298.257
b_WGS=a_WGS*(1-flatness)				# WGS polar radius
eccentricity=sqrt(a_WGS**2-b_WGS**2)/a_WGS

# ------------ jro radar specifications -------------------------
lat0=-11.947917*deg						# geodetic, the usual map or GPS latitude
lon0=-76.872306*deg						# east of Greenwich
h0=0.463								# local height above reference ellipsoid
n0=a_WGS/sqrt(1-flatness*(2-flatness)*sin(lat0)**2.)
x0=(n0+h0)*cos(lat0)*cos(lon0)			# cartesian geocentric coordinates wrt Greenwich
y0=(n0+h0)*cos(lat0)*sin(lon0)
z0=(n0*(1-eccentricity**2)+h0)*sin(lat0)
xyz0=array([x0,y0,z0])
xy0=array([x0,y0])
r0=sqrt(dot(xyz0,xyz0))
p0=sqrt(dot(xy0,xy0))

# unit vectors from jro
east0=array([-y0,x0,0])/p0				# zenith and north directions wrt local ellipsoid
zenith0=array([cos(lat0)*cos(lon0),cos(lat0)*sin(lon0),sin(lat0)])
north0=cross(zenith0,east0)

# orthonormal basis vectors including the jro on-axis direction
dec=-12.88*deg
ha=-(4+37./60.)*deg						# on-axis direction at JRO
uo=array([cos(dec)*cos(ha/4.+lon0),cos(dec)*sin(ha/4.+lon0),sin(dec)])	# on axis
ux=cross(zenith0,uo)
ux=ux/sqrt(dot(ux,ux))					# along the building to the right
uy=cross(uo,ux)							# away from the building into the valley
