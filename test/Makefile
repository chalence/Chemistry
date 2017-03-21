

CHEM_DIR = /Users/ctss/Science/Chemistry
INCLUDE_DIRS = -I. -I$(CHEM_DIR)

CC = g++
FC = gfortran
CFLAGS = -O3 
FFLAGS = -lifcore -O3 -ffree-line-length-none -extend_source -fdefault-real-8

HEADERS = 
OBJ = chemtest.o

#include $(CHEM_DIR)/makefile
VPATH = $(CHEM_DIR)

CHEM_OBJ     += chemF.o ChemDriver.o ChemBurner.o ChemIntegrate.o Jacobian.o ChemNetwork.o
OBJ          += $(CHEM_OBJ)
CHEM_HEADERS += ChemDriver.H ChemBurner.H ChemGlobals.H ChemIntegrate.H ChemExternal.H ChemNetwork.H

HEADERS      += $(CHEM_HEADERS)

CFLAGS += -I$(CHEM_DIR)
FFLAGS += -I$(CHEM_DIR)


#print:
#	@echo $(HEADERS)


chemtest: $(OBJ)
	$(CC)  $(OBJ)  -o $@ 

.SUFFIXES: .F90

.F90.o:
	$(FC) -c $(FFLAGS)  $<

.cpp.o: 
	$(CC) -c $(CFLAGS)  $<

$(OBJ):  $(HEADERS)

clean:
	rm -f $(RESULT) *.o *.mod




