# $Id$


all: lac

clean: source-clean
	rm -f *~ include/*~ include/*/*~ lib/lib* lib/Makefile.dep

lac:
	cd source ; make

source-clean:
	cd source ; make clean


.PHONY: all deal.II
.PHONY: clean source-clean
