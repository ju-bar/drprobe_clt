# MAKEFILE for compiling and linking msa
# Defining variables
f90comp = ifort
switch = -nologo -O3 -recursive -132 -assume byterecl -threads -static-intel
# removed switches
#	 -libs:static ???
# Start of the makefile
msa: emsdata.o MultiSlice.o STEMfunctions.o msaparams.o FFTs_x.o msasub.o random.o msa.o
	$(f90comp) -o msa -nologo emsdata.o MultiSlice.o STEMfunctions.o msaparams.o FFTs_x.o msasub.o random.o msa.o
emsdata.mod: emsdata.o ../common/emsdata.f90
	$(f90comp) -c $(switch) ../common/emsdata.f90
emsdata.o: ../common/emsdata.f90
	$(f90comp) -c $(switch) ../common/emsdata.f90
MultiSlice.mod: MultiSlice.o MultiSlice.f90
	$(f90comp) -c $(switch) MultiSlice.f90
MultiSlice.o: MultiSlice.f90
	$(f90comp) -c $(switch) MultiSlice.f90
STEMfunctions.mod: STEMfunctions.o STEMfunctions.f90
	$(f90comp) -c $(switch) STEMfunctions.f90
STEMfunctions.o: STEMfunctions.f90
	$(f90comp) -c $(switch) STEMfunctions.f90
msaparams.mod: STEMfunctions.mod MultiSlice.mod msaparams.o msaparams.f90
	$(f90comp) -c $(switch) msaparams.f90
msaparams.o: STEMfunctions.mod MultiSlice.mod msaparams.f90
	$(f90comp) -c $(switch) msaparams.f90
msa.o: STEMfunctions.mod MultiSlice.mod msaparams.mod msa.f90
	$(f90comp) -c $(switch) msa.f90
FFTs_x.o: ../common/FFTs_x.f
	$(f90comp) -c $(switch) ../common/FFTs_x.f
msasub.o: STEMfunctions.mod MultiSlice.mod emsdata.mod msaparams.mod msasub.f90
	$(f90comp) -c $(switch) msasub.f90
random.o: ../common/random.f90
	$(f90comp) -c $(switch) ../common/random.f90
clean:
	rm *.mod *.o msa
# End of the makefile
# J. Barthel, ju.barthel@fz-juelich.de, 2018
# Forschungszentrum Juelich GmbH, 52425 Juelich, Germany
# published under the GNU General Public License, Version 3