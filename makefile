#VPATH        += $(SRC)/Chemistry
INCLUDE_DIRS += -I$(SRC)/Chemistry

CHEM_OBJ     += chemF.o ChemDriver.o ChemBurner.o ChemIntegrate.o Jacobian.o
OBJ          += $(CHEM_OBJ)
CHEM_HEADERS += ChemDriver.H ChemBurner.H ChemGlobals.H ChemIntegrate.H ChemExternal.H

HEADERS      += $(CHEM_HEADERS)

$(CHEM_OBJ):  $(HEADERS)

