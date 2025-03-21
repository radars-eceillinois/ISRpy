"""
Utility and application modules for ISR (and CSR) data processing and modeling

where:

    ISR=Incoherent Scatter Radar
    CSR=smaller ionospheric/atmospheric radars than ISR's that can only detect
stronger coherent scatter than minimal level ionospheric incoherent scatter


Try "ISRpy. TAB" in ipython for module methods and submodules...

-readpack and beampack provides data read and target positioning
routines for a collection of ISR's and CSR's.

-pyigrf is a packet for geomagnetic vector field computataions based on the
IGRF model.

-kmltools provides routines for placing radar beam and data maps in Google Earth.

-enoise is spectral noise level estimator based on Hildebrand-Sekhon algorithm.

Release 0.0.20
- Include PFISR beam
Release 0.0.21 (1/4/2022)
- Adding elaz2xyz definition to find a target with range, elevation, and azimuth
Release 0.0.22 (1/4/2022)
- making elaz2xyz available to all radars
Release 0.0.28 (3/19/2025)
- updating igrf to igrf-14 (2025)
"""

from . import readpack
from . import beampack
from . import enoise
from .pyigrf import pyigrf

# from . import kmltools # needs to be fixed
__version__ = '0.0.28'

igrf = pyigrf() # instantiating with the latest IGRF coefficients

__all__=['readpack','beampack','kmltools','enoise', 'igrf']
