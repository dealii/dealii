# $Id$
# Copyright W. Bangerth, University of Heidelberg, 1998


all: deal.II examples

clean: source-clean examples-clean lib-clean 


deal.II:
	cd source ; $(MAKE)

examples:
	cd examples ; $(MAKE)

lib-clean:
	-rm -f *~ include/*~ include/*/*~ include/basic/forward-declarations.h*
	cd lib ; $(MAKE) clean

examples-clean:
	cd examples ; $(MAKE) clean

source-clean:
	cd source ; $(MAKE) clean

TAGS:
	etags include/*/*.h source/*/*.cc


.PHONY: all deal.II examples TAGS
.PHONY: clean lib-clean examples-clean source-clean
