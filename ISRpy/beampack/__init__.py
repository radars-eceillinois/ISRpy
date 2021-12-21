"""
For soft-target radars including ISR's the target is effectively the content
and Bragg scale spatial Fourier component of the content density of the radar
scattering volume determined by the radar antenna beam and transmitted pulse
length and receiver filter.

beampack provides routines for describing
(a) the location of the radar scattering volume with respect to the locations
of a collection ISR and also in an Earth centered coordinate system,
(b) geomagnetic field vectors and magnetic aspect angles within radar beams.
"""

from . import   AObeam,\
                bahirdarbeam,\
                esrbeam,\
                irisbeam,\
                jrobeam,\
                jrobeam0,\
                kairabeam,\
                kwajbeam,\
                laplatabeam,\
                sanyabeam,\
                tromsobeam


__all__=['jrobeam','kwajbeam','AObeam','sanyabeam','irisbeam','kairabeam',
'tromsobeam','esrbeam','laplatabeam','bahirdarbeam']
