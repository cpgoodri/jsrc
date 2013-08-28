DIR = /data0/home/cpgoodri/jcode/test_project/jsrc
#DIR = /home/cpgoodri/projects/test3/jsrc

name=test
obj=$(name).o


Walnut_NAME=walnut
Fiji_NAME=fiji
#COMPUTER_NAME=walnut
COMPUTER_NAME=fiji


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



