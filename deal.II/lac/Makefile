# $Id$


all: lac doc 

clean: lib-clean doc-clean


lac:
	cd source ; make

doc:
	cd doc ; make


lib-clean:
	cd lib ; rm lib*

doc-clean:
	cd doc ; make clean


.PHONY: all deal.II doc
.PHONY: clean lib-clean doc-clean
