# $Id$
# Copyright W. Bangerth, University of Heidelberg, 1998


all: deal.II doc examples

clean: lib-clean doc-clean source-clean examples-clean


deal.II:
	cd source ; make

doc:
	cd doc ; make

examples:
	cd examples ; make

lib-clean:
	-rm -f *~ include/*~ include/*/*~
	cd lib ; make clean


doc-clean:
	cd doc ; make clean

examples-clean:
	cd examples ; make clean

source-clean:
	cd source ; make clean
TAGS:
	etags include/*/* source/*/*


.PHONY: all deal.II doc examples TAGS
.PHONY: clean lib-clean doc-clean examples-clean source-clean
