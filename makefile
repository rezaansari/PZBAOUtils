#  Makefile pour des programmes utilitaires BAO/PhotoZ
#      R. Ansari , (C) LAL+UPS , Sep 2014    
include $(SOPHYABASE)/include/sophyamake.inc

#  Define our target list 
all : Objs/galcatext Objs/grid2pk Objs/tpk2d

clean :
	rm -f Objs/*


##############
Objs/galcatext : Objs/galcatext.o 
	$(CXXLINK) -o Objs/galcatext Objs/galcatext.o  $(SOPHYAEXTSLBLIST) 

Objs/galcatext.o : galcatext.cc 
	$(CXXCOMPILE)  -o Objs/galcatext.o  galcatext.cc

##############
Objs/grid2pk : Objs/grid2pk.o Objs/gpkspec.o
	$(CXXLINK) -o Objs/grid2pk Objs/grid2pk.o Objs/gpkspec.o $(SOPHYAEXTSLBLIST) 

Objs/grid2pk.o : grid2pk.cc gpkspec.h 
	$(CXXCOMPILE)  -o Objs/grid2pk.o  grid2pk.cc
##############
Objs/tpk2d : Objs/tpk2d.o Objs/gpkspec.o Objs/myinteg2d.o 
	$(CXXLINK) -o Objs/tpk2d Objs/tpk2d.o Objs/gpkspec.o Objs/myinteg2d.o $(SOPHYAEXTSLBLIST) 

Objs/tpk2d.o : tpk2d.cc gpkspec.h corfunc.h myinteg2d.h hsplfit.h 
	$(CXXCOMPILE)  -o Objs/tpk2d.o  tpk2d.cc
###
Objs/myinteg2d.o : myinteg2d.cc myinteg2d.h 
	$(CXXCOMPILE)  -o Objs/myinteg2d.o  myinteg2d.cc

###
Objs/gpkspec.o : gpkspec.cc gpkspec.h 
	$(CXXCOMPILE)  -o Objs/gpkspec.o  gpkspec.cc
