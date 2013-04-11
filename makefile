DIR = /data0/home/cpgoodri/jcode/test_project/jsrc

name=Test_Computer
obj=$(name).o



prjDIR = .
prjOBJGQS = \
$(obj)

include $(DIR)/std_make.mk

OBJGQS = \
$(patsubst %,$(srcDIR)/%,$(srcOBJGQS)) \
$(patsubst %,$(prjDIR)/%,$(prjOBJGQS))




.f.o:
	$(FRULE)

.cpp.o:
	$(CRULE)

$(name).out: $(OBJGQS) $(StandardDependencies)
	$(ORULE)





