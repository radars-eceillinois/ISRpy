"""

Created by Erhan Kudeki on 11/29/08.
Copyright (c) 2008 ECE, UIUC. All rights reserved.

    [Bn,Be,Bd,B]=igrf_B(year,ht,lon,lat,ORD)
    returns X Y Z components of geomagnetic field based on igrf-x model
    and B=sqrt(Bn**2 + Be**2 + Bd**2), Bn=north, Be=east, Bd=down (nT),
    1900.<year<2025., ht(km above Earth radius a), (a=6371.2  km)
    lon(deg, east>0), lat(deg, geocentric, north>0)
    note: geodetic coordinates should be translated to geocentric
    before calling this function.
    based on an earler MATLAB code by Erhan Kudeki, March 2004.
    Ref: The Earth's magnetic field, Merril & McElhinny, Academic

history:
-igrf11.py version based on igrf11coeffs.txt file is created by EK on 4/29/11
-sum over g and h coefficients defined for indices m and n running from
-igrfb.py now searches for the latest file with pattern igrf??coeffs.txt
 date: 02/03/2020
-P. Reyes on 2/9/2020 has converted igrf into a class.
 01/30/2021 - P. Reyes has simplified pyigrf package by having the routines
     directly into the __init__.py file and not using a class. Reversing
     the change from 2/9/2020
 02/19/2021 - P. Reyes has created a branch that brings back the routines
     inside a class. Advantages are, among others, to load a different set
     of coefficients from a file.

"""

import numpy as np
import os

__version__ = "0.0.4"

__IGRF_MODELS__ = [
    # Model, Release Year, Main Field 1, MF2, Secular Variation 1, SV2, file
    ['IGRF-13', 2020, 1900., 2020., 2020., 2025., 'igrf13coeffs.txt'],
    ['IGRF-12', 2015, 1900., 2015., 2015., 2020., 'igrf12coeffs.txt'],
    ['IGRF-11', 2010, 1900., 2010., 2010., 2015., 'igrf11coeffs.txt'],
    ]
__LATEST_IGRF__ = "IGRF-13"


class pyigrf:
    a_igrf = 6371.2 #IGRF Earth's radius
    _avail_models = dict(
        [[x[0], {'release year':x[1],
                 'main field':x[2:4],
                 'secular variation':x[4:6],
                 'coeffs file':x[6]}
         ] for x in __IGRF_MODELS__
        ])
    def __init__(self, model_name = None, coeff_file=None,
            verbose=0):
        """Holds an igrf instance with an  igrf model loaded.

        available method:
        igrf_B
        Usage:
        [Bn,Be,Bd,B] = igrf_B(year, ht, lon, lat)

        Example:
        --------
        >>> import ISRpy
        >>> [Bn,Be,Bd,B] = ISRpy.igrf.igrf_B(year, ht, lon, lat)

        Example: Use a different IGRF model
        -----------------------------------
        To start a new instance with t different available model:

        >>> import ISRpy
        >>> igrf12 = ISRpy.pyigrf("IGRF-12", verbose=1)
        Reading IGRF-12 coefficients
        max_n is 13
        Last Epoch year is: 2015.0
        secular variation: 2015-20
        >>> [Bn,Be,Bd,B] = igrf12.igrf_B(year, ht, lon, lat)
        """
        self.verbose = verbose
        if type(model_name) == type(None):
            self.model_name = __LATEST_IGRF__
        else:
            self.model_name = model_name.upper()
        if not type(coeff_file) == type(None):
            # Coefficients will be obtained from coeff_file
            self.coeff_file = coeff_file
            self.model_name = "fromfile"
            if self.verbose >= 1:
                print(f"IGRF coefficients will be read from file:{coeff_file}")
        else:
            if self.verbose >= 1:
                print(f"Reading {self.model_name} coefficients ")
            assert self.model_name in self._avail_models.keys(), \
                    f"{self.model_name} not one of "\
                    f"{list(self._avail_models.keys())}"
            this_file_folder = os.path.split(os.path.abspath(__file__))[0]
            self.coeff_file = os.path.join(this_file_folder,'igrfdata',
                self._avail_models[self.model_name]['coeffs file'])

        self._read_coeff_file(self.coeff_file)
        self._get_m_n_schmidt()

    def display_available_models(self):
        print("Available IGRF models:",list(self._avail_models.keys()))


    def igrf_B(self,year,ht,lon,lat):
        """
        [Bn,Be,Bd,B] = igrf_B(year, ht, lon, lat)
        returns Bn Be Bd components of geomagnetic field based on igrf-X model
        and B=sqrt(Bn**2 + Be**2 + Bd**2), Bn=north, Be=east, Bd=down (nT),
        1900.<year<max_year.,
        ht: (km above Earth radius a),
        lon: (deg, east>0),
        lat: (deg, geocentric, north>0)
             note: geodetic coordinates should be translated to geocentric
                   before calling this function.
        """
        from scipy.special import lpmn

        epoch = self.epoch
        gdat = self._gdat
        hdat = self._hdat
        m = self._m
        n = self._n
        a = self.a_igrf
        max_n = self.max_n
        schmidt = self._schmidt

        base = (np.nonzero(year >= epoch)[0]).max()     # base epoch year index
        y0 = epoch[base]      # starting year
        if year >= epoch[-1]: # If year larger than last epoch, use SV
            gSV = gdat[-1,:,:] # last value is SV
            hSV = hdat[-1,:,:] # last value is SV
        elif year < epoch[-1]: # otherwise linear interpolation
            y1 = epoch[base+1]      # ending year
            gSV = (gdat[base+1,:,:] - gdat[base,:,:]) / (y1-y0)
            hSV = (hdat[base+1,:,:] - hdat[base,:,:]) / (y1-y0)
        g = gdat[base,:,:] + (year-y0) * gSV
        h = hdat[base,:,:] + (year-y0) * hSV

        phi = lon*np.pi/180.    # set phi=longitude dependence - co-sinusoids
        cp  = np.cos(m * phi)
        sp  = np.sin(m * phi)
        az  = g * cp + h * sp
        az_phi = m * (-g * sp + h * cp)

        r = a + ht        # set geocentric altitude dependence
        amp   = a * ((a / r) ** (n + 1))
        amp_r = -(n + 1) * amp / r                # r derivative of amp


        theta = (90. - lat) * np.pi / 180.    # set theta=colatitude dependence
        ct = np.cos(theta)
        st = np.sqrt(1. - ct ** 2.)
        [lPmn,lPmn_der] = lpmn(max_n, max_n, ct)    # assoc legendre and derivative
        lPmn = lPmn * schmidt    # schmidt normalization
        lPmn_theta = -st * lPmn_der * schmidt

        Z = np.sum((amp_r * lPmn * az))       # get field components (nT)
        Y = -np.sum((amp * lPmn * az_phi)) / (r * st)
        X = np.sum((amp * lPmn_theta * az)) / r
        B = np.sqrt(X ** 2. + Y ** 2. + Z ** 2.)

        return X,Y,Z,B

    def _read_coeff_file(self, coeff_file):
        """Reads the IGRF coefficients from coeff_file."""
        try:
            with open(coeff_file,'r') as fp:
                txtlines = fp.read().splitlines()
        except:
            raise IOError("Problems reading coefficients file:\n%s"%coeff_file)

        for line in reversed(txtlines): #start from bottom to get largest n
            if len(line) < 3:
                continue # If line is too small skip
            max_n = int(line.split()[1]) # getting largest n (13 in igrf11)
            if self.verbose:
                print("max_n is",max_n)
            break
        for line in txtlines:
            if len(line) < 3: continue # If line is too small skip
            if line[0:2] in ['g ', 'h ']: # reading the coefficients
                n = int(line.split()[1])
                m = int(line.split()[2])
                if line[0] == 'g':
                    gdat[:,m,n] = np.array(line.split()[3:], dtype=float)
                elif line[0] == 'h':
                    hdat[:,m,n] = np.array(line.split()[3:], dtype=float)
            elif line[0:3] == 'g/h': #reading the epochs
                all_epochs = line.split()[3:]
                secular_variation = all_epochs[-1]
                epoch = np.array(all_epochs[:-1],dtype=float) # read the epochs
                gdat = np.zeros([epoch.size+1,max_n + 1, max_n + 1],float) #SV+1
                hdat = np.zeros([epoch.size+1,max_n + 1, max_n + 1],float) #SV+1

        if self.verbose:
            print("Last Epoch year is:",epoch[-1])
            print("secular variation:",secular_variation)

        self.max_n = max_n
        self._gdat = gdat
        self._hdat = hdat
        self.secular_variation = secular_variation
        self.epoch = epoch


    def _get_m_n_schmidt(self, max_n=None):
        """
         Builds up the "schmidt" coefficients !!! careful with this definition
         Schmidt quasi-normalized associated Legendre functions of degree n
         and order m. Thebault et al. 2015
        """
        # ------ declare and initialize fixed parameters for all epochs ---------
        if not type(max_n) == type(None):
            self.max_n = max_n
        [m, n]=np.mgrid[0:self.max_n + 1, 0:self.max_n + 1]  # set up 14X14(IGRF11) meshgrid

        from scipy.special import factorial as factorial
        schmidt = np.sqrt(2 * factorial(n - m) / factorial(n + m)) * (-1) ** m
        schmidt[0,:] = 1.
        self._m = m
        self._n = n
        self._schmidt = schmidt

