# Dr. Probe command-line tools (drprobe_clt)


## Authors and Copyright

   Juri Barthel, 
   Forschungszentrum J�lich GmbH, 52425 J�lich, Germany

   Copyright (c) 2008 - 2023 - Forschungszentrum J�lich GmbH
   
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

## Compiling

### Windows

Solution and project files are included that should work with Visual Studio 2017+ 
and the Intel Parallel Studio or the Intel oneAPI.

The code is using the dfti module (mkl_dfti.f90) of the Intel MKL. This module
must be compiled, which is not the case by default. Add the source file manually
to the CELSLC project once and compile it with a release configuration. Then
copy the mkl_dft*.mod files to the lp64 sub-folder of the MKL include folder
(ilp64 for 64-bit integer code, default is 32-bit integers). Once the mod files
are there, and you always use the same MKL version, mkl_dfti can be removed again
from the project and CELSLC and MSA should both compile and link. This manual
procedure is needed, since the MKL code folder can be different on different
systems.

### Linux

Check the makefiles for correct include folders. The setup given should work in
most cases.

## Testing

Example input data is provided in the subfolder "test". The following calls 
are simple test cases. Please, adopt folders and file names according to
your local situation.

### STEM image simulations

* HAADF, ABF and BF thickness series simulation with subsequent source convolutions:
```
celslc -cif test/STO_001_4x4.cif -nx 625 -ny 625 -nz 2 -nv 50 -ht 300 -slc test/sto_001_4x4_300keV_fl50 -fl
msa -prm test/msa-1.prm -out test/img.dat /3dout
msa -prm test/msa-1.prm -in test/img_HAADF.dat -out test/img_HAADF_psc.dat /3dout
msa -prm test/msa-1.prm -in test/img_ABF.dat -out test/img_ABF_psc.dat /3dout
msa -prm test/msa-1.prm -in test/img_BF.dat -out test/img_BF_psc.dat /3dout
```

* HAADF thickness series simulation with subsequent source convolutions
including separation of elastic and thermal-diffuse scattering with
5 processes running in parallel, each solving a set of different scan lines:
```
celslc -cif test/STO_001_4x4.cif -nx 625 -ny 625 -nz 2 -nv 50 -ht 300 -slc test/sto_001_4x4_300keV_fl50 -fl
start msa -prm test/msa-1-avg.prm -out test/img.dat /3dout -py 0 -ly 3 /silavwaveft /verbose
start msa -prm test/msa-1-avg.prm -out test/img.dat /3dout -py 4 -ly 7 /silavwaveft /verbose
start msa -prm test/msa-1-avg.prm -out test/img.dat /3dout -py 8 -ly 11 /silavwaveft /verbose
start msa -prm test/msa-1-avg.prm -out test/img.dat /3dout -py 12 -ly 15 /silavwaveft /verbose
start msa -prm test/msa-1-avg.prm -out test/img.dat /3dout -py 16 -ly 19 /silavwaveft /verbose
msa -prm test/msa-1.prm -in test/img_HAADF_tot.dat -out test/img_HAADF_tot_psc.dat /3dout
msa -prm test/msa-1.prm -in test/img_HAADF_ela.dat -out test/img_HAADF_ela_psc.dat /3dout
msa -prm test/msa-1.prm -in test/img_HAADF_tds.dat -out test/img_HAADF_tds_psc.dat /3dout
```

### HR-TEM image simulations

* NCSI image:
```
celslc -cif test/STO_001_4x4.cif -nx 256 -ny 256 -nz 2 -ht 300 -slc test/sto_001_4x4_300keV_dwf -dwf -abf 0.07
msa -prm test/msa-2.prm -out test/img.dat /ctem
wavimg -prm test/wavimg-2.prm -out test/img_ctem.dat
```

## Development

The programs are written in Fortran 90 code for Intel Fortran compilers. CELSLC
and MSA are linked against the Intel Math Kernel Libraries (MKL).

Current versions:
* CELSLC: 1.1
* MSA: 1.1
* WAVIMG: 1.0

### TODOs

* add slice loading on demand to MSA
* add a wrapper tool for applying external transition potentials for EELS and EDX
* add partial import of wavefunctions and images from 3d data sets with
  MSA and WAVIMG.

