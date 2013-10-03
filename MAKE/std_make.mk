srcDIR = $(DIR)

srcOBJGQS = \
Computers/cijkl.o \
Computers/UmfpackInterface.o \
Resources/Resources.o \
Resources/Exception.o \
Resources/index_map.o \
Resources/RNG_taus.o

#Compiler
CPP=icpc
CC=icc
FF=ifort

#Compiler flags
I_ACCURACY_FLAGS = -fp-model precise -prec-div -prec-sqrt
I_OPTIMIZATION_FLAGS = -O3 -xhost -ipo
FLAGS = $(I_OPTIMIZATION_FLAGS) $(I_ACCURACY_FLAGS) #-opt-report-file opt_report.txt
FFLAGS = $(FLAGS)
CFLAGS = $(FLAGS) -wr1125 -wr2196 -wr2536 #-wr21 -wr279 -wr1125 #-wr418 
LinkFLAGS = $(FLAGS)


ifeq ($(COMPUTER_NAME),$(Walnut_NAME))
include $(DIR)/MAKE/walnut_make.mk
else
ifeq ($(COMPUTER_NAME),$(Fiji_NAME))
include $(DIR)/MAKE/fiji_make.mk
endif
endif




FRULE = $(FF)  $(FFLAGS) $(INCLUDE) -c -o $@ $<
CRULE = $(CPP) $(CFLAGS) $(INCLUDE) -c -o $@ $<
ORULE = $(CPP) $(CFLAGS) -o $@ $(OBJGQS) $(LIBRARY) $(LINK)


StandardDependencies = \
	$(srcDIR)/Boundaries/*.h \
	$(srcDIR)/Computers/*.h \
	$(srcDIR)/Database/*.h \
	$(srcDIR)/Minimization/*.h \
	$(srcDIR)/Potentials/*.h \
	$(srcDIR)/Resources/*.h \
	$(srcDIR)/State/*.h



#    $@  means "the target"
#    $<  means "whatever the dependencies are"



