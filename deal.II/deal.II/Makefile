# $Id$
# Copyright W. Bangerth, University of Heidelberg, 1998


all: deal.II examples

clean: lib-clean source-clean examples-clean


deal.II:
	cd source ; make

examples:
	cd examples ; make

lib-clean:
	-rm -f *~ include/*~ include/*/*~ include/basic/forward-declarations.h*
	cd lib ; make clean

examples-clean:
	cd examples ; make clean

source-clean:
	cd source ; make clean
TAGS:
	etags include/*/*.h source/*/*.cc


.PHONY: all deal.II examples TAGS
.PHONY: clean lib-clean examples-clean source-clean
