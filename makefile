#  Makefile pour des programmes utilitaires BAO/PhotoZ
#      R. Ansari , (C) LAL+UPS , Sep 2014    
include $(SOPHYABASE)/include/sophyamake.inc

#  Define our target list 
all : Objs/galcatext Objs/grid2pk

clean :
	rm -f Objs/*


######
Objs/galcatext : Objs/galcatext.o 
	$(CXXLINK) -o Objs/galcatext Objs/galcatext.o  $(SOPHYAEXTSLBLIST) 

Objs/galcatext.o : galcatext.cc 
	$(CXXCOMPILE)  -o Objs/galcatext.o  galcatext.cc

######
Objs/grid2pk : Objs/grid2pk.o Objs/gpkspec.o
	$(CXXLINK) -o Objs/grid2pk Objs/grid2pk.o Objs/gpkspec.o $(SOPHYAEXTSLBLIST) 

Objs/grid2pk.o : grid2pk.cc gpkspec.h 
	$(CXXCOMPILE)  -o Objs/grid2pk.o  grid2pk.cc

Objs/gpkspec.o : gpkspec.cc gpkspec.h 
	$(CXXCOMPILE)  -o Objs/gpkspec.o  gpkspec.cc
