# $Id$


all: lac doc 

clean: doc-clean source-clean
	rm -f *~ include/*~ include/*/*~ lib/lib*

lac:
	cd source ; make

doc:
	cd doc ; make


doc-clean:
	cd doc ; make clean

source-clean:
	cd source ; make clean


.PHONY: all deal.II doc
.PHONY: clean doc-clean source-clean
