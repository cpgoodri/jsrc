#DIR = /data0/home/cpgoodri/jcode/ssdirectional/jsrc
DIR = /home/cpgoodri/jsrc

name=test4
obj=$(name).o

#define some computer names
Walnut_NAME=walnut
Fiji_NAME=fiji

#set this to the computer name you are using
COMPUTER_NAME=walnut


prjDIR = .
prjOBJGQS = \
$(obj)

#include $(DIR)/MAKE/simple_make.mk
include $(DIR)/MAKE/std_make.mk

OBJGQS = \
$(patsubst %,$(srcDIR)/%,$(srcOBJGQS)) \
$(patsubst %,$(prjDIR)/%,$(prjOBJGQS))





.f.o: 
	$(FRULE)

.cpp.o: 
	$(CRULE)

$(name).out: $(OBJGQS)
	$(ORULE)

#if any header file is changed, the project file gets recompiled.
$(obj): $(StandardDependencies)


clean:
	\rm $(OBJGQS)



