# $Id$


all: lac

clean: source-clean
	rm -f *~ include/*~ include/*/*~ lib/lib* lib/o/*.o lib/go/*.go lib/Makefile.dep include/lac/forward-declarations.h*

lac:
	cd source ; make

source-clean:
	cd source ; make clean


.PHONY: all deal.II
.PHONY: clean source-clean
