#  Makefile pour des programmes utilitaires BAO/PhotoZ
#      R. Ansari , (C) LAL+UPS , Sep 2014    
include $(SOPHYABASE)/include/sophyamake.inc

#  Define our target list 
all : Objs/galcatext

clean :
	rm -f Objs/*


######

Objs/galcatext : Objs/galcatext.o 
	$(CXXLINK) -o Objs/galcatext Objs/galcatext.o  $(SOPHYAEXTSLBLIST) 

Objs/galcatext.o : galcatext.cc 
	$(CXXCOMPILE)  -o Objs/galcatext.o  galcatext.cc

