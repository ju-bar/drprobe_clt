# MAKEFILE for compiling and linking msa
# Defining variables
f90c = ifort
cpar = -nologo -O3 -recursive -132 -assume byterecl -threads -static-intel -I../common
lpar = /usr/local/lib/libfftw3f.a
# removed switches
#	 -libs:static ???
# Start of the makefile
msa: emsdata.o MultiSlice.o STEMfunctions.o msaparams.o msasub.o random.o msa.o
	$(f90c) -o msa -nologo emsdata.o MultiSlice.o STEMfunctions.o msaparams.o msasub.o random.o msa.o $(lpar)
emsdata.mod: emsdata.o ../common/emsdata.f90
	$(f90c) -c $(cpar) ../common/emsdata.f90
emsdata.o: ../common/emsdata.f90
	$(f90c) -c $(cpar) ../common/emsdata.f90
MultiSlice.mod: MultiSlice.o MultiSlice.f90
	$(f90c) -c $(cpar) MultiSlice.f90
MultiSlice.o: MultiSlice.f90
	$(f90c) -c $(cpar) MultiSlice.f90
STEMfunctions.mod: STEMfunctions.o STEMfunctions.f90
	$(f90c) -c $(cpar) STEMfunctions.f90
STEMfunctions.o: STEMfunctions.f90
	$(f90c) -c $(cpar) STEMfunctions.f90
msaparams.mod: STEMfunctions.mod MultiSlice.mod msaparams.o msaparams.f90
	$(f90c) -c $(cpar) msaparams.f90
msaparams.o: STEMfunctions.mod MultiSlice.mod msaparams.f90
	$(f90c) -c $(cpar) msaparams.f90
msa.o: STEMfunctions.mod MultiSlice.mod msaparams.mod msa.f90
	$(f90c) -c $(cpar) msa.f90
msasub.o: STEMfunctions.mod MultiSlice.mod emsdata.mod msaparams.mod msasub.f90
	$(f90c) -c $(cpar) msasub.f90
random.o: ../common/random.f90
	$(f90c) -c $(cpar) ../common/random.f90
#
.PHONY: clean
clean:
	rm *.mod *.o msa
# End of the makefile
# J. Barthel, ju.barthel@fz-juelich.de, 2019
# Forschungszentrum Juelich GmbH, 52425 Juelich, Germany
# published under the GNU General Public License, Version 3