# g4track code - for interfacing GEANT4 tracking with Lucretia
# "make all" makes shared obj library for use when building Lucretia
# mex files
# "make clean" to clean up all object and library files
#
# =============================================================================
# ASSUMES GEANT4 (VERSION 10+) INSTALLED IN STANDARD WAY WITH PATH DEFINED
# IN SYSTEM ("geant4-config" script must run) AND ALL MATERIAL FILES INCLUDED
# =============================================================================
#
# GW: May 29, 2014

# Link to matlab installation directory
 MATPATH = /var/MATLAB/R2014a
 ARCHEXT = glnxa64

# CFLAGS := $(shell geant4-config --cflags) -O2 -fexceptions -fPIC -fno-omit-frame-pointer -pthread -DG4MULTITHREADED
 CFLAGS := $(shell geant4-config --cflags) -O2 -fexceptions -fPIC -fno-omit-frame-pointer -pthread
 CLIBS := $(shell geant4-config --libs) -I$(MATPATH)/extern/include -L$(MATPATH)/bin/$(ARCHEXT)/ -lmx -lmex -lmat
 LFLAGS := -O2 -pthread -shared -Wl,--version-script,$(MATPATH)/extern/lib/$(ARCHEXT)/mexFunction.map -Wl,--no-undefined

# g4track Library
all : g4track.o lucretiaManager.o geomConstruction.o actionInitialization.o eventAction.o primaryGeneratorAction.o runAction.o trackingAction.o FieldSetup.o GlobalField.o

# cleanup
clean::
	-rm -f *.o
	-rm -f *.so*
	-rm -f libg4track.a

g4track.o : ../LucretiaCommon.h g4track.h geomConstruction.hh actionInitialization.hh
	g++ $(CFLAGS) $(CLIBS) -c g4track.cpp
	ar -ru libg4track.a g4track.o

geomConstruction.o: geomConstruction.hh
	g++ $(CFLAGS) $(CLIBS) -c geomConstruction.cpp
	ar -ru libg4track.a geomConstruction.o

actionInitialization.o: actionInitialization.hh primaryGeneratorAction.hh runAction.hh eventAction.hh trackingAction.hh
	g++ $(CFLAGS) $(CLIBS) -c actionInitialization.cpp
	ar -ru libg4track.a actionInitialization.o

eventAction.o: eventAction.hh
	g++ $(CFLAGS) $(CLIBS) -c eventAction.cpp
	ar -ru libg4track.a eventAction.o

primaryGeneratorAction.o: primaryGeneratorAction.hh
	g++ $(CFLAGS) $(CLIBS) -c primaryGeneratorAction.cpp
	ar -ru libg4track.a primaryGeneratorAction.o

runAction.o: runAction.hh
	g++ $(CFLAGS) $(CLIBS) -c runAction.cpp
	ar -ru libg4track.a runAction.o

trackingAction.o: trackingAction.hh
	g++ $(CFLAGS) $(CLIBS) -c trackingAction.cpp
	ar -ru libg4track.a trackingAction.o

lucretiaManager.o: lucretiaManager.hh
	g++ $(CFLAGS) $(CLIBS) -c lucretiaManager.cpp
	ar -ru libg4track.a lucretiaManager.o

FieldSetup.o: FieldSetup.hh
	g++ $(CFLAGS) $(CLIBS) -c FieldSetup.cpp
	ar -ru libg4track.a FieldSetup.o

GlobalField.o: GlobalField.hh
	g++ $(CFLAGS) $(CLIBS) -c GlobalField.cpp
	ar -ru libg4track.a GlobalField.o
