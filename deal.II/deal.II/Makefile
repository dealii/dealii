# $Id$


all: deal.II doc examples

clean: lib-clean examples-clean doc-clean


deal.II:
	cd source ; make

doc:
	cd doc ; make

examples:
	cd examples ; make


lib-clean:
	cd lib ; rm lib*

doc-clean:
	cd doc ; make clean

examples-clean:
	cd examples ; make clean



.PHONY: all deal.II doc examples
.PHONY: clean lib-clean doc-clean examples-clean
