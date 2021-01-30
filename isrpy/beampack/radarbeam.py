#
#  radarbeam.py
#
#	module for calculating geometry parameters and magnetic aspect
#	angle of radar targets monitored by any radar
#
#	use aspect_elaz or aspect_txty to calculate aspect angles of targets
#	specified by (el,az) or (tx,ty) angles
#
#  Created by Erhan Kudeki on 11/29/08 as jrobeam.py
#  Copyright (c) 2008 ECE, UIUC. All rights reserved.
#  history
#  - Aug29,2013 by P. Reyes
#    -Generate a module that accepts the lon,lat,h coordinates for the location
#    of any radar.
#    -flattening has been changed from 1/298.257	to 1./298.257223563
#    using the WGS84 reference in:
#    http://earth-info.nga.mil/GandG/publications/tr8350.2/wgs84fin.pdf
#    - A new routine called enu2xyz to move a point from xr,yr,zr to some
#      direction east, north, up

def llh2xyz(latg,lon,h):
    # returns geocentric xyz coordinates (ECEF) in km of a target with
    # latitude       latg (rad) --- geodetic
    # longitude      lon (rad)
    # height         h (km above local ellipsoid)
    n=a_WGS / np.sqrt(1.-flatness*(2.-flatness) * np.sin(latg)**2.)
    # cartesian geocentric coordinates wrt Greenwich
    x=(n+h)*np.cos(latg)*np.cos(lon)
    y=(n+h)*np.cos(latg)*np.sin(lon)
    z=(n*(1.-eccentricity**2.)+h)*np.sin(latg)
    return x,y,z

def xyz2llh(x,y,z):
    # returns longitude 'lon', geodetic latitude 'lat', and height 'h'
    # of position (x,y,z) defined in geocentric coordinate system (ECEF)

    # on Oct23,2013 by P. Reyes, adding the .all() in order to support
    # arrays

    p=np.sqrt(x**2.+y**2.)
    lon=np.arctan2(y,x)
    lat=np.arctan2(z,p)
    latp=lat.copy()
    for i in range(10):
        n=a_WGS/np.sqrt(1.-flatness*(2-flatness)*np.sin(latp)**2.)
        h=p/np.cos(latp)-n
        lat=np.arctan(z/(p*(1.-n*eccentricity**2./(n+h))))
        if (abs(lat-latp)<3.*eps).all():
            n=a_WGS/np.sqrt(1.-flatness*(2.-flatness)*np.sin(lat)**2.)
            h=p/np.cos(lat)-n
            break
        latp=lat.copy()
    return lat,lon,h

def enu2xyz(xr,yr,zr,east,north,up):
    # moves a point from xr,yr,zr to x,y,z by moving into the direction
    # specified by east,north,up (enu) coordinates in km
    latg,lon,h = xyz2llh(xr,yr,zr)
    A = np.array([[-np.sin(lon),-np.sin(latg)*np.cos(lon),np.cos(latg)*np.cos(lon)],
                  [ np.cos(lon),-np.sin(latg)*np.sin(lon),np.cos(latg)*np.sin(lon)],
                  [           0 , np.cos(latg)             ,np.sin(latg)]])
    x,y,z = np.dot(A,np.array([east,north,up]))+np.array([xr,yr,zr])
    return x,y,z

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

    [bX,bY,bZ,bB]=igrf_B(year,r-a_igrf,lon/deg,lat/deg)
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
import numpy as np
import pyigrf

eps=np.finfo(float).eps         # float resolution
deg=np.pi/180.                  # to express angles in degree values
a_igrf=6371.2                   # mean earth radius (km)

# WGS84 constants
# reference:
# http://earth-info.nga.mil/GandG/publications/tr8350.2/wgs84fin.pdf

a_WGS=6378.137         # equatorial radius WGS 84 (semi-major axis) in km
#flatness=1/298.257
flatness = 1./298.257223563  # flatenning
b_WGS=a_WGS*(1.-flatness)    # WGS polar radius (semi-minor axis) in km
eccentricity=np.sqrt(a_WGS**2-b_WGS**2)/a_WGS

# ------------ radar specifications -------------------------
class radarspecs:
    """Will contain radar coordinates and coordinate conversions
    saved locations:
    JRO : lat: -11.947917 , lon: -76.872306, h0: 0.463 km
    JRO_GE : as zoom in with GoogleEarth to the center of the antenna.
    IRIS@ROI
    ALTAIR
    IRIS@URBANA
    """
    def __init__(self,lat0=None,lon0=None,h0=None,location=None):
        if location!=None:
            if location.upper() == "JRO":
                # geodetic, the usual map or GPS latitude
                self.lat0 = -11.947917 * deg
                self.lon0 = -76.872306 * deg
                self.h0 = 0.463 # local height in km above reference ellipsoid
            elif location.upper() == "JRO_GE":
                # gusing google earth to the center of the Antenna
                # -11.9514944444 = -(11.+57./60.+5.38/3600.) # 11deg57'5.38"S
                self.lat0 = -11.9514944444 * deg
                # -76.8743916667#-(76.+52./60.+27.81/3600.) #  76deg52'27.81"W
                self.lon0 = -76.8743916667 * deg
                self.h0 = 0.463 # local height in km above reference ellipsoid
            elif location.upper() == "IRIS@ROI":
                #  9.39794444444 = (9.+23./60.+52.6/3600.) # 9deg23'52.60"N
                self.lat0 = 9.39794444444 * deg
                #  167.469166667 = (167.+28./60.+9./3600.) # 167deg28'9.00"E
                self.lon0 = 167.469166667 * deg
                self.h0 = 0.012
            elif location.upper() == "ALTAIR":
                #  9.39794444444 = (9.+23./60.+43.5/3600.) # 9deg23'43.50"N
                self.lat0 = 9.39541666667 * deg
                #  167.469166667 = (167.+28./60.+45.6/3600.) # 167deg28'45.60"E
                self.lon0 = 167.479333333 * deg
                self.h0 = 0.012
            elif location.upper() == "IRIS@URBANA":
                # 40.16683888888889 = (40.+10./60.+0.62/3600.) #40deg10'0.62"N
                self.lat0 = 40.16683888888889 * deg
                #-88.1586 = -(88.+9./60.+30.96/3600.) #88deg9'30.96"W
                self.lon0 = (360. -88.1586) * deg
                self.h0 = 0.221
        elif lat0==None or lon0==None or h0==None:
            # By default: JRO center of antenna with google earth
            # -11.9514944444 = -(11.+57./60.+5.38/3600.) # 11deg57'5.38"S
            self.lat0 = -11.9514944444 * deg
            # -76.8743916667#-(76.+52./60.+27.81/3600.) #  76deg52'27.81"W
            self.lon0 = -76.8743916667 * deg
            self.h0 = 0.463 # local height in km above reference ellipsoid
        else:
            self.lat0 = lat0 * deg
            self.lon0 = lon0* deg
            self.h0 = h0 # local height in km above reference ellipsoid
        x0,y0,z0 = llh2xyz(self.lat0,self.lon0,self.h0)
        self.xyz0 = np.array([x0,y0,z0])
        xy0 = np.array([x0,y0])
        p0 = np.sqrt(np.dot(xy0,xy0))

        # unit vectors from jro
        self.east0 = np.array([-y0,x0,0])/p0
        # zenith and north directions wrt local ellipsoid
        self.zenith0 = np.array([np.cos(self.lat0) * np.cos(self.lon0),
                            np.cos(self.lat0) * np.sin(self.lon0),
                            np.sin(self.lat0)])
        self.north0 = np.cross(self.zenith0,self.east0)

        # orthonormal basis vectors including the jro on-axis direction
        dec=-12.88*deg
        ha=-(4.+37./60.)*deg # on-axis direction at JRO
        self.uo = np.array([np.cos(dec) * np.cos(ha/4. + self.lon0), # on axis
                    np.cos(dec) * np.sin(ha/4. + self.lon0), np.sin(dec)])
        self.ux = np.cross(self.zenith0,self.uo)
        # along the building to the right
        self.ux = self.ux / np.sqrt(np.dot(self.ux,self.ux))
        # away from the building into the valley
        self.uy = np.cross(self.uo,self.ux)

    def locations(self):
        return ["JRO","JRO_GE","IRIS@ROI","ALTAIR","IRIS@URBANA"]

    def dec_ha2el_az(dec,ha):
        # returns elevation and azimuth angles of a radar beam
        # with respect to local tangent plane.
        # the beam is specified by:
        #		declination dec (deg)
        #		hour angle  ha  (min)
        # with respect to radar location at longitude lon0 and height h0
        # above reference ellipsiod at geodetic latitude lat0

        lat=dec*deg                 # on celestial sphere
        lon=2.*pi*(ha/(24.*60.))
        lon=lon+lon0                # on celestial sphere
        vec=array([cos(lat)*cos(lon),cos(lat)*sin(lon),sin(lat)])
        hor=vec-dot(vec,zenith0)*zenith0
        hor=hor/sqrt(dot(hor,hor))
        el=arccos(dot(hor,vec))/deg
        north=dot(hor,north0)
        east=dot(hor,east0)
        az=arctan2(east,north)/deg
        return el,az

    def xyz2dec_ha(self,vec):
        # declination and hour angle in target direction used to describe radar
        # beam direction at JRO, corresponding to latitude and relative
        # longitude of the beam-spot on the celestial sphere, corresponds to
        # rr->\infty, in which case:
        vec = vec/np.sqrt(np.dot(vec,vec))
        p = np.sqrt(vec[0]**2.+vec[1]**2.)
        dec = np.arctan2(vec[2],p)/deg                           # in degrees
        ha = (np.arctan2(vec[1],vec[0]) - self.lon0)*(24./(2.*np.pi))*60.    # in minutes
        return dec,ha

    def aspect_angle(self,year,xyz):
        # returns the magnetic aspect angle (rad) of a target with
        # geocentric vector xyz defined in geocentric coordinates

        r   = np.sqrt(np.dot(xyz,xyz))
        p   = np.sqrt(xyz[0]**2. + xyz[1]**2.)
        lat = np.arctan2(xyz[2],p)
        lon = np.arctan2(xyz[1],xyz[0])
        radial = xyz/r;			   # directions from target
        east   = np.array([-xyz[1],xyz[0],0.])/p
        north  = -np.cross(east,radial)
        rr = xyz - self.xyz0
        u_rr = rr / np.sqrt(np.dot(rr,rr))   # unit vector from radar to target

        [bX,bY,bZ,bB] = pyigrf.igrf_B(year, r - a_igrf, lon/deg, lat/deg)
        bfield = np.array([bX,bY,bZ])
        B = bX*north + bY*east - bZ*radial
        u_B = B / np.sqrt(np.dot(B,B))
        aspect = np.arccos(np.dot(u_B, u_rr))
        return r,lat,lon,aspect

    def aspect_txty(self,year,rr,tx,ty):
        # returns magnetic aspect angle and geocentric coordinates of a target
        # tracked by jro at
        # range rr (km)
        # tx along jro building
        # ty into the building

        tz = np.sqrt(1.-tx**2.-ty**2.)
        #geocentric coordinates of target
        xyz = self.xyz0 + rr*(tx*self.ux + ty*self.uy + tz*self.uo)

        [r,lat,lon,aspect] = self.aspect_angle(year,xyz)
        [dec,ha] = self.xyz2dec_ha(xyz - self.xyz0)
        return r,lon,lat,dec,ha,aspect

    def aspect_elaz(self,year,rr,el,az):
        # returns magnetic aspect angle and geocentric coordinates of a target
        # tracked by jro at
        # range       rr (km)
        # elevation   el (rad above local tangent plane to ellipsoid)
        # azimuth     az (rad east of local north)

        tx = np.cos(el) * np.sin(az)		# direction cosines wrt east and north
        ty = np.cos(el) * np.cos(az)
        tz = np.sin(el)
        #geocentric coordinates of target :
        xyz = self.xyz0 + rr*(tx * self.east0 + ty*self.north0+tz*self.zenith0)

        [r,lat,lon,aspect] = self.aspect_angle(year,xyz)
        [dec,ha] = xyz2dec_ha(xyz - self.xyz0)
        return r,lon,lat,dec,ha,aspect
