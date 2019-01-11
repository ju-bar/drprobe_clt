# Dr. Probe command-line tools (drprobe_clt)


## Authors and Copyright

   Juri Barthel, 
   Forschungszentrum Jülich GmbH, 52425 Jülich, Germany

   Copyright (c) 2008 - 2019 - Forschungszentrum Jülich GmbH
   
   Published under the GNU General Public License, version 3,
   see <http://www.gnu.org/licenses/> and LICENSE!


## CELSLC

[CELSLC](http://www.er-c.org/barthel/drprobe/celslc.html)
is a program to calculate object transmission functions to be used
as phase gratings in a multislice algorithm for electron diffraction
calculations. The calculations require an atomic structure model as input,
including the definition of a calculation box, atomic coordinates,
thermal vibration parameters, and partial occupancy factors. Further
parameters concern numerical sampling and the probing electron energy. The
output produced can be used as input of the program MSA.


## MSA

[MSA](http://www.er-c.org/barthel/drprobe/msa.html)
is a program to calculate the diffraction of beam of probing electrons
through a crystal. The crystal data is input in form of phase gratings or
projected scattering potentials as calculated by the program CELSLC. Further
parameters concern the probe forming, sample thickness, scan settings etc.
Output are electron wave functions or STEM images.


## WAVIMG

[WAVIMG](http://www.er-c.org/barthel/drprobe/wavimg.html)
is a program used for the calculation of high-reslolution TEM images
from an input electron wave function.


## Documentation

Documentation and a few examples can be found on the
[Dr. Probe website](http://www.er-c.org/barthel/drprobe/). In addition, each
tool has its own "howto" text file. These files are used as the primary source
of documentation.


## Development

The programs are written in Fortran 90 code for Intel Fortran compilers.

The code of the program MSA links to "libfftwf-3.3.lib"
MSA uses data output by FFTW and is in no form based on work represented by
the FFTW project. Source code and library binary code of FFTW are available
from http://www.fftw.org/ (accessed April 2018).

### TODOs

* keep up with new code from [JMultiSlice](https://github.com/ju-bar/JMultiSliceLib)
  - investigate if partial C++ code can be used in this project
* add notes on experimental and undocumented features to the howto files.
* add a wrapper tool for applying external transition potentials for EELS and EDX

