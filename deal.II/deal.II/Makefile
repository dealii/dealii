# $Id$


all: deal.II doc examples

clean: examples-clean doc-clean source-clean
	rm -f *~ include/*~ include/*/*~ lib/lib*


deal.II:
	cd source ; make

doc:
	cd doc ; make

examples:
	cd examples ; make


doc-clean:
	cd doc ; make clean

examples-clean:
	cd examples ; make clean

source-clean:
	cd source ; make clean



.PHONY: all deal.II doc examples
.PHONY: clean doc-clean examples-clean source-clean
