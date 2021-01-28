"""igrf magnetic field model

igrf.py 	---implements igrf10 model based on igrf10coeff.txt valid until 2010
igrf11.py       ---implements igrf11 model based on igrf11coeff.txt valid until 2015

usage:
import igrf.igrf or igrf.igrf11 as grf 
grf.igrf_B(2011.5,150.,0.,0.)


"""

__all__=['igrf','igrf11']
