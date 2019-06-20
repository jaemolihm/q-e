# Makefile for EPW

include ../../make.inc
include make.libs

#
# use recursive definitions in Makefile to modify the PATH for modules
#

MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../EPW/src \
          $(MOD_FLAG)../../PW/src \
          $(MOD_FLAG)../../PHonon/PH \
          $(MOD_FLAG)../../LR_Modules \
          $(MOD_FLAG)../../wannier90-3.0.0/src/obj \
          $(MOD_FLAG)../../dft-d3

#
# The above is in alphabetical order except when order matter during compilation
#


# setphases.o \ : SP: We keep it in case its usefull in the future but now depreciated.

default : run_wigner.x


PHOBJS = ../../PHonon/PH/libph.a
PWOBJS = ../../PW/src/libpw.a
W90LIB = ../../wannier90-3.0.0/libwannier.a 
LRMODS = ../../LR_Modules/liblrmod.a
PWOBJS = ../../PW/src/libpw.a 
EPWOBJS = ../../EPW/src/libepw.a 
QEMODS = ../../Modules/libqemod.a ../../KS_Solvers/libks_solvers.a \
         ../../FFTXlib/libqefft.a ../../dft-d3/libdftd3qe.a
LIBOBJS =../../LAXlib/libqela.a ../../UtilXlib/libutil.a ../../clib/clib.a ../../iotk/src/libiotk.a 

run_wigner.x : run_wigner.o $(EPWOBJS) $(PHOBJS) $(LRMODS) $(PWOBJS) $(W90LIB) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
	run_wigner.o  $(EPWOBJS) $(PHOBJS) $(LRMODS) $(W90LIB) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)

clean :
	- /bin/rm -f  *.o *~ *.d *.mod *.i

include make.depend
