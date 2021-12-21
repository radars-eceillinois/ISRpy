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
