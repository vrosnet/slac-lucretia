# g4track code - for interfacing GEANT4 tracking with Lucretia
# "make libg4track" makes shared obj library for use when building Lucretia
# mex files
# ASSUMES GEANT4 (VERSION 10+) INSTALLED IN STANDARD WAY WITH PATH DEFINED
# IN SYSTEM ("geant4-config" script must run) AND ALL MATERIAL FILES INCLUDED
#=============================================================================

# Link to matlab installation directory
 MATPATH = /var/MATLAB/R2013b
 ARCHEXT = glnxa64

# CFLAGS := $(shell geant4-config --cflags) -O2 -fexceptions -fPIC -fno-omit-frame-pointer -pthread -DG4MULTITHREADED
 CFLAGS := $(shell geant4-config --cflags) -O2 -fexceptions -fPIC -fno-omit-frame-pointer -pthread
 CLIBS := $(shell geant4-config --libs) -I$(MATPATH)/extern/include -L$(MATPATH)/bin/$(ARCHEXT)/ -lmx -lmex -lmat
 LFLAGS := -O2 -pthread -shared -Wl,--version-script,$(MATPATH)/extern/lib/$(ARCHEXT)/mexFunction.map -Wl,--no-undefined

test : g4track_test

clean::
	-rm -f *.o
	-rm -f *.so*
	-rm -f libg4track.a

#mex : g4track_mex.$(MEXEXT)

# Standalone executable test
g4track_test : geomConstruction.o actionInitialization.o eventAction.o primaryGeneratorAction.o runAction.o trackingAction.o
	g++ $(CFLAGS) $(CLIBS) -o g4track_test g4track_test.o geomConstruction.o actionInitialization.o eventAction.o primaryGeneratorAction.o runAction.o trackingAction.o

# g4track Library
libg4track : g4track.o geomConstruction.o actionInitialization.o eventAction.o primaryGeneratorAction.o runAction.o trackingAction.o lucretiaManager.o
#libg4track : g4track.o geomConstruction.o actionInitialization.o eventAction.o primaryGeneratorAction.o runAction.o trackingAction.o lucretiaManager.o
#	g++ $(LFLAGS) -o libg4track.so g4track.o geomConstruction.o actionInitialization.o eventAction.o primaryGeneratorAction.o runAction.o trackingAction.o lucretiaManager.o

# the code
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
