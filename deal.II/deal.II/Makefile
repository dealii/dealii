# $Id$


all: deal.II deal doc examples

clean: lib-clean doc-clean examples-clean


deal.II:
	cd source ; make

deal:
	cd deal ; make

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



.PHONY: all deal.II deal doc examples
.PHONY: clean lib-clean doc-clean examples-clean
