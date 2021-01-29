"""igrf.py

Created by Erhan Kudeki on 11/29/08.
Copyright (c) 2008 ECE, UIUC. All rights reserved.

	[X,Y,Z,B]=igrf_B(year,ht,lon,lat,ORD)
	returns X Y Z components of geomagnetic field based on igrf2000 model
	and B=sqrt(X**2+Y**2+Z**2), X=north, Y=east, Z=down (nT),
	1900.<year<2010., ht(km above Earth radius a),
	lon(deg, east>0), lat(deg, geocentric, north>0)
	note: geodetic coordinates should be translated to geocentric
	before calling this function.
	based on an earler MATLAB code by Erhan Kudeki, March 2004.
	Ref: The Earth's magnetic field, Merril & McElhinny, Academic

igrf11.py version based on igrf11coeffs.txt file is created by EK on 4/29/11

input parameter ORD with a default of 10 added by EK on 12/28/11
sum over g and h coefficients defined for indices m and n running from
0 to ORD --- max possible value for ORD in IGRF is 13.
"""

def igrf_B(year,ht,lon,lat,ORD=10):
	base=int((year-1900.)/5)			# base model
	shift=year-(5*base+1900)		# shift years
	g=(1-shift/5.)*gdat[base][:][:]+(shift/5.)*gdat[base+1][:][:]
	h=(1-shift/5.)*hdat[base][:][:]+(shift/5.)*hdat[base+1][:][:]

	phi=lon*pi/180.	# set phi=longitude dependence - co-sinusoids
	cp=cos(m*phi)
	sp=sin(m*phi)
	az=g*cp+h*sp
	az_phi=m*(-g*sp+h*cp)

	r=a+ht		# set geocentric altitude dependence
	amp=a*((a/r)**(n+1))
	amp_r=-(n+1)*amp/r				# r derivative of amp

	from scipy.special import lpmn
	theta=(90.-lat)*pi/180.	# set theta=colatitude dependence
	ct=cos(theta)
	st=sqrt(1.-ct**2.)
	[lPmn,lPmn_der]=lpmn(13,13,ct)	# assoc legendre and derivative
	lPmn=lPmn*schmidt	# schmidt normalization
	lPmn_theta=-st*lPmn_der*schmidt

#	figure()
#	imshow(schmidt,interpolation='nearest')
#	colorbar()
#	title("schmidt_mn")

	Z=sum((amp_r*lPmn*az)[0:ORD+1,0:ORD+1])	   # get field components (nT)
	Y=-sum((amp*lPmn*az_phi)[0:ORD+1,0:ORD+1])/(r*st)
	X=sum((amp*lPmn_theta*az)[0:ORD+1,0:ORD+1])/r
	B=sqrt(X**2.+Y**2.+Z**2.)

	return X,Y,Z,B

"""
--------------- read in igrf schmidt coefficients g and h -------------
the following section of the module establishes the namespace and variables
needed by the module function igrf_B() upon importing the module as in:
import igrf11 as grf
-----------------------------------------------------------------------
adjust the following code after new a coefficients file replaces igrf11
-----------------------------------------------------------------------
"""
from pylab import *
from string import atoi, atof

import os
pwd=os.path.dirname(__file__)
pwd=pwd.encode('utf-8')
if pwd=='':
	file=open('igrf11coeffsORIG.txt')	# open the coefficients textfile
else:
	file=open(pwd+'/igrf11coeffsORIG.txt')	# open the coefficients textfile

line=file.readline()			# read the first line into string
while line[0] in ["#", ' ']:
	line = file.readline()
keys=line.split()
epochs=size(keys[3:])			# number of epochs in current model

gdat=zeros([epochs,14,14],float)	# 3D arrays are declared in [0,13]
hdat=zeros([epochs,14,14],float)

for line in file:			# read the data lines one by one
	#print line
	dat=line.split()
	if size(dat[3:])<epochs-2:	# only the 10X10 part retained in IGRF10
		break                   # but has no effect in IGRF11
	coeff=dat[0]
	n=atoi(dat[1])
	m=atoi(dat[2])
	ep=0
	for item in dat[3:]:
		if coeff=="g":
			gdat[ep][m][n]=atof(item)
		else:
			hdat[ep][m][n]=atof(item)
		ep+=1
file.close()		       # close the coefficients file

for n in range(1,14):	       # extrapolate over next epoch g/hdat are in 0-13
	for m in range(n+1):
		gdat[epochs-1][m][n]=gdat[epochs-2][m][n]+5.*gdat[epochs-1][m][n]
		hdat[epochs-1][m][n]=hdat[epochs-2][m][n]+5.*hdat[epochs-1][m][n]

# ------ declare and initialize fixed parameters for all epochs ---------
a=6371.2					# igrf earth radius
[m,n]=mgrid[0:14,0:14]				# set up 14X14 meshgrid

from scipy.misc import factorial

"""
 build up the "schmidt" coefficients !!! careful with this definition
"""
schmidt=sqrt(2*factorial(n-m)/factorial(n+m))*(-1)**m
schmidt[0,:]=1.
