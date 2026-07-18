# ISRpy

Utility and application modules for incoherent scatter radar (ISR) and coherent scatter radar (CSR) 
data processing and modeling.

ISR's are radar systems with large power-aperture products being used in ionospheric probing to 
detect weak returns scattered from thermal level fluctuations encountered in ionospheres in 
thermal equilibrium. CSR's are smaller radar systems used in atmospheric and ionospheric 
research which can only detect stronger returns from the scattering medium caused by larger
amplitude fluctuations in the medium amplified by instabilities and turbulence.

Try "ISRpy. TAB" in iPython or Jupyter for methods and submodules...

-readpack and beampack provides data read and target positioning
routines for a collection of ISR's and CSR's.

-pyigrf is a package for geomagnetic vector field computataions based on the
IGRF model.

-kmltools provides routines for placing radar beam and data maps in Google Earth.

-enoise is spectral noise level estimator based on Hildebrand-Sekhon algorithm.

## IGRF examples

    import ISRpy
    import datetime
    datetime0 = datetime.datetime(2026,7,1,0,0,0)
    # convert datetime into float year
    fyear = ISRpy.beampack.beamtools.datetime2year(datetime0)
    print(fyear)

2026.495890410959

    # obtain IGRF components at lat,lon = 37.4538,-122.1822
    # see level ht = 0
    lat,lon = 37.4538,-122.1822
    ht = 0
    [Bn,Be,Bd,B] = ISRpy.igrf.igrf_B(fyear, ht, lon, lat)
    Bn,Be,Bd,B

(np.float64(22228.509544177166),  
np.float64(5080.003723167496),  
np.float64(41438.41308344268),  
np.float64(47297.51741113908))  

    # Display available IGRF models:
    ISRpy.igrf.display_available_models()

Available IGRF models: ['IGRF-14', 'IGRF-13', 'IGRF-12', 'IGRF-11', 'IGRF-10', 'IGRF-9']  

    # Display the different epochs from a legacy IGRF model:
    igrf12 = ISRpy.pyigrf("IGRF-12", verbose=1)
    print(igrf12.epoch)

Reading IGRF-12 coefficients  
max_n is 13  
Last Epoch year is: 2015.0  
secular variation: 2015-20  
[1900. 1905. 1910. 1915. 1920. 1925. 1930. 1935. 1940. 1945. 1950. 1955.  
 1960. 1965. 1970. 1975. 1980. 1985. 1990. 1995. 2000. 2005. 2010. 2015.]  
 
