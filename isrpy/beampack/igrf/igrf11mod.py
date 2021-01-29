#
#  igrf11.py
#  
#
#  Created by Erhan Kudeki on 11/29/08.
#  Copyright (c) 2008 ECE, UIUC. All rights reserved.
#
#	[X,Y,Z,B]=igrf_B(year,ht,lon,lat)
#	returns X Y Z components of geomagnetic field based on igrf2000 model
#	and B=sqrt(X**2+Y**2+Z**2), X=north, Y=east, Z=down (nT),
#	1900.<year<2010., ht(km above Earth radius a), 
#	lon(deg, east>0), lat(deg, geocentric, north>0)
#	note: geodetic coordinates should be translated to geocentric
#	before calling this function.
#	based on an earler MATLAB code by Erhan Kudeki, March 2004.
#	Ref: The Earth's magnetic field, Merril & McElhinny, Academic
#
#  igrf11.py version based on igrf11coeffsORIG.txt file is created by EK on 4/29/11
#

def igrf_B(year,ht,lon,lat):
	base=int((year-1900)/5)			# base model 
	shift=year-(5*base+1900)		# shift years
	g=(1-shift/5.)*gdat[base][:][:]+(shift/5.)*gdat[base+1][:][:]
	h=(1-shift/5.)*hdat[base][:][:]+(shift/5.)*hdat[base+1][:][:]
			
	phi=lon*pi/180.					# set phi=longitude dependence - co-sinusoids
	cp=cos(m*phi)
	sp=sin(m*phi)
	az=g*cp+h*sp
	az_phi=m*(-g*sp+h*cp)
	
	r=a+ht							# set geocentric altitude dependence
	amp=a*((a/r)**(n+1))	
	amp_r=-(n+1)*amp/r				# r derivative of amp
	
	from scipy.special import lpmn
	theta=(90.-lat)*pi/180.			# set theta=colatitude dependence
	ct=cos(theta)
	st=sqrt(1.-ct**2.)
	[lPmn,lPmn_der]=lpmn(13,13,ct)	# associated legendre function and derivative
	lPmn=lPmn*schmidt				# schmidt normalization
	lPmn_theta=-st*lPmn_der*schmidt 

#	figure()
#	imshow(schmidt,interpolation='nearest')
#	colorbar()
#	title("schmidt_mn")
		
	Z=sum(amp_r*lPmn*az)			# get field components (nT)
	Y=-sum(amp*lPmn*az_phi)/(r*st)
	X=sum(amp*lPmn_theta*az)/r		
	B=sqrt(X**2.+Y**2.+Z**2.)
	
	return X,Y,Z,B
		
# --------------- read in igrf schmidt coefficients g and h -------------
# -----------------------------------------------------------------------
# adjust the following code after new a coefficients file replaces igrf10	

from pylab import *
from string import atoi, atof
from pdb import set_trace

import os
#pwd = os.path.dirname(__file__)
pwd = os.path.split(os.path.abspath(__file__))[0]
file=open(pwd+'/igrf11coeffsORIG.txt')		# open the coefficients textfile
line=file.readline()				# read the first line into string

while line[0] in ["#", ' '] :
        line = file.readline()
keys = line.split()
epochs = size(keys[3:])                         # number of epochs in current model

gdat = zeros([epochs,14,14],float)              # 3D Coefficient arrays are declared
hdat = zeros([epochs,14,14],float)

for line in file:					# read the data lines one by one
	print line

	dat=line.split()
	if size(dat[3:])<epochs-2:		# only the 10X10 part retained
		break

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

file.close()						# close the coefficients file

for n in range(1,14):				# extrapolate over the next epoch 
	for m in range(n+1):
		gdat[epochs-1][m][n]=gdat[epochs-2][m][n]+5.*gdat[epochs-1][m][n]
		hdat[epochs-1][m][n]=hdat[epochs-2][m][n]+5.*hdat[epochs-1][m][n]

# ------ declare and initialize fixed parameters for all epochs ---------
a=6371.2							# igrf earth radius
[m,n]=mgrid[0:14,0:14]				# set up 14X14 meshgrid

set_trace()

from scipy.misc.common import factorial

schmidt=sqrt(2*factorial(n-m)/factorial(n+m))*(-1)**m
schmidt[0,:]=1.		# build up the "schmidt" coefficients !!! careful with this definition
#
#fig=figure()
#imshow(schmidt,interpolation='nearest')
#colorbar()
#title("schmidt_mn")
#fig.show()
