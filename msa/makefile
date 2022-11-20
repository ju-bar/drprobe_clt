# MAKEFILE for compiling and linking msa
# Defining variables
# - choose the compiler, ifort is recommended
f90c = ifort
# - set compiler options
cflags := -nologo -O3 -recursive -132 -assume byterecl -threads -static-intel -fpp
#cflags := $(cflags) -Ddouble_precision
#   add -Ddouble_precision to compile for double precision calculation, default is single precision without this option
incl = ../common 
# mkl path must be adjusted to the local mkl installation
ifneq ("$(wildcard ${MKLROOT})", "")
    mklpath = ${MKLROOT}
else
    mklpath = /usr/local/intel/mkl
endif
mklincl = $(mklpath)/include
# mkl lib linking must be adjusted to the local mkl installation
# - when linking 64-bit MKL with ifort on OSX:
#lpar = $(mklpath)/lib/libmkl_intel_lp64.a $(mklpath)/lib/libmkl_sequential.a $(mklpath)/lib/libmkl_core.a -lpthread -lm -ldl
# - when linking 64-bit MKL with ifort on Linux:
lpar = -Wl,--start-group $(mklpath)/lib/intel64/libmkl_intel_lp64.a $(mklpath)/lib/intel64/libmkl_sequential.a $(mklpath)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl
# - when linking FFTW
# lpar = /usr/local/lib/libfftw3f.a
# removed switches
#	 -libs:static ???
# Start of the makefile
msa: precision.o emsdata.o MultiSlice.o STEMfunctions.o msaparams.o msasub.o plasmon.o random.o msa.o
	$(f90c) -o msa -nologo precision.o emsdata.o MultiSlice.o STEMfunctions.o msaparams.o msasub.o plasmon.o random.o msa.o $(lpar)
precision.mod: precision.o ../common/precision.f90
	$(f90c) -c $(cflags) -I $(incl) ../common/precision.f90
precision.o: ../common/precision.f90
	$(f90c) -c $(cflags) -I $(incl) ../common/precision.f90
emsdata.mod: emsdata.o ../common/emsdata.f90
	$(f90c) -c $(cflags) -I $(incl) ../common/emsdata.f90
emsdata.o: ../common/emsdata.f90
	$(f90c) -c $(cflags) -I $(incl) ../common/emsdata.f90
MultiSlice.mod: MultiSlice.o MultiSlice.f90
	$(f90c) -c $(cflags) -I $(incl) MultiSlice.f90
MultiSlice.o: plasmon.mod mkl_dfti.mod MultiSlice.f90
	$(f90c) -c $(cflags) -I $(incl) MultiSlice.f90
STEMfunctions.mod: STEMfunctions.o STEMfunctions.f90
	$(f90c) -c $(cflags) -I $(incl) STEMfunctions.f90
STEMfunctions.o: STEMfunctions.f90
	$(f90c) -c $(cflags) -I $(incl) STEMfunctions.f90
msaparams.mod: STEMfunctions.mod MultiSlice.mod msaparams.o msaparams.f90
	$(f90c) -c $(cflags) -I $(incl) msaparams.f90
msaparams.o: STEMfunctions.mod MultiSlice.mod msaparams.f90
	$(f90c) -c $(cflags) -I $(incl) msaparams.f90
msa.o: STEMfunctions.mod MultiSlice.mod msaparams.mod msa.f90
	$(f90c) -c $(cflags) -I $(incl) msa.f90
msasub.o: STEMfunctions.mod MultiSlice.mod emsdata.mod msaparams.mod plasmon.mod msasub.f90
	$(f90c) -c $(cflags) -I $(incl) msasub.f90
plasmon.mod: plasmon.o plasmon.f90
	$(f90c) -c $(cflags) -I $(incl) plasmon.f90
plasmon.o: plasmon.f90
	$(f90c) -c $(cflags) -I $(incl) plasmon.f90
random.o: ../common/random.f90
	$(f90c) -c $(cflags) -I $(incl) ../common/random.f90
mkl_dfti.mod: mkl_dfti.o $(mklpath)/include/mkl_dfti.f90
	$(f90c) -c $(cflags) -I $(incl) -I $(mklincl) $(mklpath)/include/mkl_dfti.f90
mkl_dfti.o: $(mklpath)/include/mkl_dfti.f90
	$(f90c) -c $(cflags) -I $(incl) -I $(mklincl) $(mklpath)/include/mkl_dfti.f90
#
.PHONY: clean
clean:
	rm *.mod *.o msa
# End of the makefile
# J. Barthel, ju.barthel@fz-juelich.de, 2019
# Forschungszentrum Juelich GmbH, 52425 Juelich, Germany
# published under the GNU General Public License, Version 3
