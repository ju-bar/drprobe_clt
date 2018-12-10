# Dr. Probe command-line tools (drprobe_clt)


## Authors and Copyright

   Juri Barthel
   Forschungszentrum Jülich GmbH, 52425 Jülich, Germany
   Copyright (c) 2008 - 2018 - Forschungszentrum Jülich GmbH
   Published under the GNU General Public License, version 3,
   see <http://www.gnu.org/licenses/> and LICENSE!


## CELSLC

CELSLC is a program to calculate object transmission functions to be used
as phase gratings in a multislice algorithm for electron diffraction
calculations. The calculations require an atomic structure model as input,
including the definition of a calculation box, atomic coordinates,
thermal vibration parameters, and partial occupancy factors. Further
parameters concern numerical sampling and the probing electron energy. The
output produced can be used as input of the program MSA.


## MSA

MSA is a program to calculate the diffraction of beam of probing electrons
through a crystal. The crystal data is input in form of phase gratings or
projected scattering potentials as calculated by the program CELSLC. Further
parameters concern the probe forming, sample thickness, scan settings etc.
Output are electron wave functions or STEM images.


## WAVIMG

WAVIMG is a program used for the calculation of high-reslolution TEM images
from an input electron wave function. 


## Development

The programs are written in Fortran 90 code for Intel Fortran compilers.

### TODOs

* publish on GitHub
* keep up with new code from [JMultiSlice](https://github.com/ju-bar/JMultiSliceLib)
  - investigate if partial C++ code can be used in this project
* replace FFTPACK routines by fftwf routines in MSA and test for increased speed
* add a wrapper tool for applying external transition potentials for EELS and EDX

