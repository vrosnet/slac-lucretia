# 
# Makefile for construction of Lucretia mex files
#  Written for use on a Linux OS (EL5 & SL5 32 and 64-bit supported)
#  - edit MATPATH variable to link to the location of your
#    Matlab  installation
#  - edit MATEXT variable to provide the mex file extension for your
#    environment, either "mexglx" for 32-bit or "mexa64" for 64-bit
#  - Ensure the mex build environment is set up correctly for
#    your environment. From within Matlab: (>> mex -setup)
#  - type "make clean" to cleanup any previous build attempts
#  - type "make" to build mex executable files
#  - type "make install" to move mex files to their correct
#    installation directories
#
# Auth: GW: Aug 4 2011
# Mod: 
#
#===================================================================

# Link to matlab installation directory
 MATPATH = /Applications/MATLAB_R2011b.app
# Name of mex file extension for this architecture (mexglx or mexa64)
 MEXEXT = mexmaci64

  MCOMP = $(MATPATH)/bin/mex
  CCOMP = gcc
# Include "-DGPU" if compiling with cuda libraries
  CFLAGS= -c -fPIC -ansi -D_GNU_SOURCE -pthread 
#-I/var/MATLAB/R2011a/extern/include
#-D__CUDA_ARCH__
# Use this for release code
  COPTIMFLAGS= -O2 -DNDEBUG
# Use this for debugging code
#  COPTIMFLAGS = -g

 LMAIN  = libLucretiaMain.a
 LMATL  = libLucretiaMatlab.a

all : GetRmats.$(MEXEXT) GetTwiss.$(MEXEXT) RmatAtoB.$(MEXEXT) TrackThru.$(MEXEXT)  VerifyLattice.$(MEXEXT)

release::
	-mv GetRmats.$(MEXEXT) ../RMatrix
	-mv GetTwiss.$(MEXEXT) ../Twiss
	-mv RmatAtoB.$(MEXEXT) ../RMatrix
	-mv TrackThru.$(MEXEXT) ../Tracking
	-mv VerifyLattice.$(MEXEXT) ../LatticeVerification

# Mex files:
 GetRmats.$(MEXEXT) : $(LMAIN) $(LMATL) GetRmats.c
	$(MCOMP) GetRmats.c -L. -lLucretiaMatlab -lLucretiaMain
 GetTwiss.$(MEXEXT) : $(LMAIN) $(LMATL) GetTwiss.c
	$(MCOMP) GetTwiss.c -L. -lLucretiaMatlab -lLucretiaMain
 RmatAtoB.$(MEXEXT) : $(LMAIN) $(LMATL) RmatAtoB.c
	$(MCOMP) RmatAtoB.c -L. -lLucretiaMatlab -lLucretiaMain
 TrackThru.$(MEXEXT) : $(LMAIN) $(LMATL) TrackThru.c
	$(MCOMP) -v TrackThru.c -L. -lLucretiaMatlab -lLucretiaMain
 VerifyLattice.$(MEXEXT) : $(LMAIN) $(LMATL) VerifyLattice.c
	$(MCOMP) VerifyLattice.c -L. -lLucretiaMatlab -lLucretiaMain

# LucretiaMain library:

 $(LMAIN) : $(LMATL)\
	$(LMAIN)(LucretiaCommon.o) \
	$(LMAIN)(LucretiaPhysics.o) 

# LucretiaMatlab library

 $(LMATL) : \
	$(LMATL)(LucretiaMatlab.o) \
	$(LMATL)(LucretiaMatlabErrMsg.o)

clean::
	-rm -f *.o
	-rm -f *.lst
	-rm -f *.a
	-rm -f *.$(MEXEXT)

# now for the next level of dependencies:

$(LMAIN)(LucretiaCommon.o) : \
   LucretiaCommon.c LucretiaCommon.h LucretiaPhysics.h \
   LucretiaGlobalAccess.h LucretiaVersionProto.h \
   LucretiaDictionary.h
	$(CCOMP) $(CFLAGS) $(COPTIMFLAGS) LucretiaCommon.c
	ar -ru $(LMAIN) LucretiaCommon.o
	rm LucretiaCommon.o
	sleep 1

$(LMAIN)(LucretiaPhysics.o) : \
   LucretiaPhysics.c LucretiaCommon.h LucretiaPhysics.h \
   LucretiaGlobalAccess.h LucretiaVersionProto.h 
	$(CCOMP) $(CFLAGS) $(COPTIMFLAGS) LucretiaPhysics.c
	ar -ru $(LMAIN) LucretiaPhysics.o
	rm LucretiaPhysics.o
	sleep 1

$(LMATL)(LucretiaMatlab.o) : \
   LucretiaMatlab.c LucretiaGlobalAccess.h LucretiaMatlab.h \
   LucretiaVersionProto.h
	$(CCOMP) $(CFLAGS) $(COPTIMFLAGS) -I$(MATPATH)/extern/include LucretiaMatlab.c
	ar -ru $(LMATL) LucretiaMatlab.o
	rm LucretiaMatlab.o
	sleep 1

$(LMATL)(LucretiaMatlabErrMsg.o) : \
   LucretiaMatlabErrMsg.c LucretiaGlobalAccess.h
	$(CCOMP) $(CFLAGS) $(COPTIMFLAGS) -I$(MATPATH)/extern/include LucretiaMatlabErrMsg.c
	ar -ru $(LMATL) LucretiaMatlabErrMsg.o
	rm LucretiaMatlabErrMsg.o
	sleep 1

