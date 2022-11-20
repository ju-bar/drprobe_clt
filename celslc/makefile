# MAKEFILE for compiling and linking celslc
# Defining variables
f90comp = ifort
switch3 = -nologo -O3 -recursive -132 -assume byterecl -threads -static-intel
switch2 = -nologo -O2 -recursive -132 -assume byterecl -threads -static-intel
switch1 = -nologo -O1 -recursive -132 -assume byterecl -threads -static-intel
switch0 = -nologo -O0 -recursive -132 -assume byterecl -threads -static-intel
# mkl path must be adjusted to the local mkl installation
ifneq ("$(wildcard ${MKLROOT})", "")
    mklpath = ${MKLROOT}
else
    mklpath = /usr/local/intel/mkl
endif
mklpar = -I$(mklpath)/include
# mkl lib linking must be adjusted to the local mkl installation
# - when linking 64-bit MKL with ifort on OSX:
#lpar = $(mklpath)/lib/libmkl_intel_lp64.a $(mklpath)/lib/libmkl_sequential.a $(mklpath)/lib/libmkl_core.a -lpthread -lm -ldl
# - when linking 64-bit MKL with ifort on Linux:
lpar = -Wl,--start-group $(mklpath)/lib/intel64/libmkl_intel_lp64.a $(mklpath)/lib/intel64/libmkl_sequential.a $(mklpath)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl
# removed switches
#	 -libs:static ???
# Start of the makefile
celslc: binio2.o spacegroups.o symops.o cifio.o emsdata.o AberrationFunctions.o mkl_dfti.o CellSlicer.o 3dpot.o fitfeprm.o fscatab.o fftmkl.o wekoscatt.o wakiscatt.o integration.o random.o celslcsubs.o celslcprm.o celslc.o
	$(f90comp) -o celslc -nologo binio2.o spacegroups.o symops.o cifio.o emsdata.o mkl_dfti.o AberrationFunctions.o CellSlicer.o 3dpot.o fitfeprm.o fscatab.o fftmkl.o wekoscatt.o wakiscatt.o integration.o random.o celslcsubs.o celslcprm.o celslc.o $(lpar)
# Modules
celslcprm.mod: celslcprm.o celslcprm.f90
	$(f90comp) -c $(switch3) celslcprm.f90
celslcprm.o: celslcprm.f90
	$(f90comp) -c $(switch3) celslcprm.f90
emsdata.mod: emsdata.o ../common/emsdata.f90
	$(f90comp) -c $(switch3) ../common/emsdata.f90
emsdata.o: ../common/emsdata.f90
	$(f90comp) -c $(switch3) ../common/emsdata.f90
AberrationFunctions.mod: AberrationFunctions.o ../common/AberrationFunctions.f90
	$(f90comp) -c $(switch3) ../common/AberrationFunctions.f90
AberrationFunctions.o: ../common/AberrationFunctions.f90
	$(f90comp) -c $(switch3) ../common/AberrationFunctions.f90
mkl_dfti.mod: mkl_dfti.o $(mklpath)/include/mkl_dfti.f90
	$(f90comp) -c $(switch3) $(mklpar) $(mklpath)/include/mkl_dfti.f90
mkl_dfti.o: $(mklpath)/include/mkl_dfti.f90
	$(f90comp) -c $(switch3) $(mklpar) $(mklpath)/include/mkl_dfti.f90
CellSlicer.mod: fscatab.o CellSlicer.o ../common/CellSlicer.f90
	$(f90comp) -c $(switch3) ../common/CellSlicer.f90
CellSlicer.o: fscatab.o ../common/CellSlicer.f90
	$(f90comp) -c $(switch3) ../common/CellSlicer.f90
m3dpot.mod: 3dpot.o 3dpot.f90
	$(f90comp) -c $(switch3) $(mklpar) 3dpot.f90
3dpot.o: 3dpot.f90
	$(f90comp) -c $(switch3) $(mklpar) 3dpot.f90
fitfeprm.mod: fitfeprm.o ../common/fitfeprm.f90
	$(f90comp) -c $(switch3) ../common/fitfeprm.f90
fitfeprm.o: ../common/fitfeprm.f90
	$(f90comp) -c $(switch3) ../common/fitfeprm.f90
fscatab.mod: fitfeprm.o fscatab.o fscatab.f90
	$(f90comp) -c $(switch3) fscatab.f90
fscatab.o: fitfeprm.o fscatab.f90
	$(f90comp) -c $(switch3) fscatab.f90
cifio.mod: spacegroups.o symops.o binio2.o cifio.o
	$(f90comp) -c $(switch3) ../common/cifio.f90
cifio.o: spacegroups.o symops.o binio2.o ../common/cifio.f90
	$(f90comp) -c $(switch3) ../common/cifio.f90
symops.mod: spacegroups.o symops.o
	$(f90comp) -c $(switch3) ../common/symops.f90
symops.o: spacegroups.o ../common/symops.f90
	$(f90comp) -c $(switch3) ../common/symops.f90
spacegroups.mod: spacegroups.o
	$(f90comp) -c $(switch0) ../common/spacegroups.f90
spacegroups.o: ../common/spacegroups.f90
	$(f90comp) -c $(switch0) ../common/spacegroups.f90
# Routine sources
binio2.o: ../common/binio2.f90
	$(f90comp) -c $(switch3) ../common/binio2.f90
fftmkl.o: ../common/fftmkl.f90
	$(f90comp) -c $(switch3) $(mklpar) ../common/fftmkl.f90
wekoscatt.o: ../common/wekoscatt.f90
	$(f90comp) -c $(switch3) ../common/wekoscatt.f90
wakiscatt.o: ../common/wakiscatt.f90
	$(f90comp) -c $(switch3) ../common/wakiscatt.f90
integration.o: ../common/integration.f90
	$(f90comp) -c $(switch3) ../common/integration.f90
random.o: ../common/random.f90
	$(f90comp) -c $(switch3) ../common/random.f90
celslcsubs.o: CellSlicer.mod emsdata.mod celslcprm.mod celslcsubs.f90
	$(f90comp) -c $(switch3) celslcsubs.f90
# Main
celslc.o: spacegroups.mod symops.mod cifio.mod emsdata.mod CellSlicer.mod fscatab.mod m3dpot.mod celslcprm.mod celslc.f90
	$(f90comp) -c $(switch3) celslc.f90
clean:
	rm *.mod *.o celslc
# End of the makefile
# J. Barthel, ju.barthel@fz-juelich.de, 2018
# Forschungszentrum Juelich GmbH, 52425 Juelich, Germany
# published under the GNU General Public License, Version 3
