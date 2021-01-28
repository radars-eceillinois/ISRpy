"""igrfb.py

Created by Erhan Kudeki on 11/29/08.
Copyright (c) 2008 ECE, UIUC. All rights reserved.

    [X,Y,Z,B]=igrf_B(year,ht,lon,lat,ORD)
    returns X Y Z components of geomagnetic field based on igrf-11 model
    and B=sqrt(X**2+Y**2+Z**2), X=north, Y=east, Z=down (nT),
    1900.<year<2015., ht(km above Earth radius a),
    lon(deg, east>0), lat(deg, geocentric, north>0)
    note: geodetic coordinates should be translated to geocentric
    before calling this function.
    based on an earler MATLAB code by Erhan Kudeki, March 2004.
    Ref: The Earth's magnetic field, Merril & McElhinny, Academic

history:
-igrf11.py version based on igrf11coeffs.txt file is created by EK on 4/29/11
-input parameter ORD with a default of 10 added by EK on 12/28/11
-sum over g and h coefficients defined for indices m and n running from
-0 to ORD --- max possible value for ORD in IGRF is 13.
-igrfb.py version by P. Reyes on 10/19/12 prepared for any # of coefficients.
 and allows for years larger than 5 years after last epoch.
"""

import pylab as py
import numpy as np
import os
this_file_folder = os.path.split(os.path.abspath(__file__))[0]
# read the information from the file
txtlines = open(this_file_folder+'/igrf11coeffs.txt','r').read().split('\n')
for line in reversed(txtlines): # start from the bottom to get largest n
    if len(line) < 3: continue # If line is too small skip
    max_n = int(line.split()[1]) # getting the largest n (13 in igrf11)
    break
for line in txtlines:
    if len(line) < 3: continue # If line is too small skip
    if line[0:2] in ['g ', 'h ']: # reading the coefficients
        n = int(line.split()[1])
        m = int(line.split()[2])
        if line[0] == 'g':
            gdat[:,m,n] = py.array(line.split()[3:], dtype=float)
        elif line[0] == 'h':
            hdat[:,m,n] = py.array(line.split()[3:], dtype=float)
    elif line[0:3] == 'g/h': #reading the epochs
        epoch = py.array(line.split()[3:-1],dtype=float) # read the epochs
        gdat = py.zeros([epoch.size+1,max_n + 1, max_n + 1],float) #SV+1
        hdat = py.zeros([epoch.size+1,max_n + 1, max_n + 1],float) #SV+1

# ------ declare and initialize fixed parameters for all epochs ---------
a=6371.2                    # igrf earth radius
[m,n]=py.mgrid[0:max_n + 1,0:max_n + 1]  # set up 14X14(IGRF11) meshgrid

from scipy.misc import factorial
"""
 build up the "schmidt" coefficients !!! careful with this definition
"""
schmidt=py.sqrt(2*factorial(n-m)/factorial(n+m))*(-1)**m
schmidt[0,:]=1.

def igrf_B(year,ht,lon,lat,ORD=max_n):
    """
    [X,Y,Z,B]=igrf_B(year,ht,lon,lat,ORD)
    returns X Y Z components of geomagnetic field based on igrf-11 model
    and B=sqrt(X**2+Y**2+Z**2), X=north, Y=east, Z=down (nT),
    1900.<year<2015., ht(km above Earth radius a),
    lon(deg, east>0), lat(deg, geocentric, north>0)
    note: geodetic coordinates should be translated to geocentric
    before calling this function.
    """
    base=py.find(year >= epoch).max()     # base epoch year index
    y0 = epoch[base]      # starting year
    if year >= epoch[-1]: # If year larger than last epoch, use SV
        gSV = gdat[-1][:][:] # last value is SV
        hSV = hdat[-1][:][:] # last value is SV
    elif year < epoch[-1]: # otherwise linear interpolation
        y1 = epoch[base+1]      # ending year
        gSV = (gdat[base+1][:][:] - gdat[base][:][:])/(y1-y0)
        hSV = (hdat[base+1][:][:] - hdat[base][:][:])/(y1-y0)
    g = gdat[base][:][:] + (year-y0)*gSV
    h = hdat[base][:][:] + (year-y0)*hSV

    phi=lon*py.pi/180.    # set phi=longitude dependence - co-sinusoids
    cp=py.cos(m*phi)
    sp=py.sin(m*phi)
    az=g*cp+h*sp
    az_phi=m*(-g*sp+h*cp)

    r=a+ht        # set geocentric altitude dependence
    amp=a*((a/r)**(n+1))
    amp_r=-(n+1)*amp/r                # r derivative of amp

    from scipy.special import lpmn
    theta=(90.-lat)*py.pi/180.    # set theta=colatitude dependence
    ct=py.cos(theta)
    st=py.sqrt(1.-ct**2.)
    [lPmn,lPmn_der]=lpmn(max_n,max_n,ct)    # assoc legendre and derivative
    lPmn=lPmn*schmidt    # schmidt normalization
    lPmn_theta=-st*lPmn_der*schmidt

    #Z=py.sum((amp_r*lPmn*az)[0:ORD+1,0:ORD+1])       # get field components (nT)
    #Y=-py.sum((amp*lPmn*az_phi)[0:ORD+1,0:ORD+1])/(r*st)
    #X=py.sum((amp*lPmn_theta*az)[0:ORD+1,0:ORD+1])/r
    Z=py.sum((amp_r*lPmn*az))       # get field components (nT)
    Y=-py.sum((amp*lPmn*az_phi))/(r*st)
    X=py.sum((amp*lPmn_theta*az))/r
    B=py.sqrt(X**2.+Y**2.+Z**2.)

    return X,Y,Z,B

def igrf_B_V(year,ht,lon,lat,ORD=max_n):
    """
    [X,Y,Z,B]=igrf_B(year,ht,lon,lat,ORD)
    returns X Y Z components of geomagnetic field based on igrf-11 model
    and B=sqrt(X**2+Y**2+Z**2), X=north, Y=east, Z=down (nT),
    1900.<year<2015., ht(km above Earth radius a),
    lon(deg, east>0), lat(deg, geocentric, north>0)
    note: geodetic coordinates should be translated to geocentric
    before calling this function.
    """
    base=py.find(year >= epoch).max()     # base epoch year index
    y0 = epoch[base]      # starting year
    if year >= epoch[-1]: # If year larger than last epoch, use SV
        gSV = gdat[-1][:][:] # last value is SV
        hSV = hdat[-1][:][:] # last value is SV
    elif year < epoch[-1]: # otherwise linear interpolation
        y1 = epoch[base+1]      # ending year
        gSV = (gdat[base+1][:][:] - gdat[base][:][:])/(y1-y0)
        hSV = (hdat[base+1][:][:] - hdat[base][:][:])/(y1-y0)
    g = gdat[base][:][:] + (year-y0)*gSV
    h = hdat[base][:][:] + (year-y0)*hSV

    phi=np.array(lon*py.pi/180.)  # set phi=longitude dependence - co-sinusoids
    phi.shape = phi.shape + (1,1) # # to prepare for multidimensional inputs
    cp=py.cos(m*phi)
    sp=py.sin(m*phi)
    az=g*cp+h*sp
    az_phi=m*(-g*sp+h*cp)

    r = np.array(a + ht)        # set geocentric altitude dependence
    r.shape = r.shape + (1,1)  # to prepare for multidimensional inputs

    amp=a*((a/r)**(n+1))
    amp_r=-(n+1)*amp/r                # r derivative of amp

    from scipy.special import lpmn
    theta=np.array((90. - lat)*py.pi/180.)    # set theta=colatitude dependence

    ct=np.cos(theta)
    st=np.sqrt(1.-ct**2.)

    arr_shape = ct.shape

    ct = ct.ravel()
    st = st.ravel()
    lPmn_arr = np.empty(ct.shape + schmidt.shape)
    lPmn_theta_arr = np.empty(ct.shape + schmidt.shape)

    for i in range(ct.shape[0]):
        [lPmn,lPmn_der]=lpmn(max_n,max_n,ct[i])    # assoc legendre and derivative
        lPmn_arr[i] = lPmn*schmidt    # schmidt normalization
        lPmn_theta_arr[i] = -st[i]*lPmn_der*schmidt
    lPmn_arr.shape = arr_shape + schmidt.shape
    lPmn_theta_arr.shape = arr_shape + schmidt.shape
    st.shape = arr_shape
    r.shape = arr_shape

    Z=py.sum(amp_r*lPmn_arr*az, axis=(-2,-1))       # get field components (nT)
    Y=-py.sum(amp*lPmn_arr*az_phi, axis=(-2,-1))/(r*st)
    X=py.sum(amp*lPmn_theta_arr*az, axis=(-2,-1))/r
    B=py.sqrt(X**2.+Y**2.+Z**2.)

    return X,Y,Z,B
