# $Id$


all: lac

clean: source-clean
	rm -f *~ include/*~ include/*/*~ lib/lib* lib/o/*.o lib/go/*.go include/base/forward-declarations.h*

lac:
	cd source ; $(MAKE)

source-clean:
	cd source ; $(MAKE) clean


.PHONY: all
.PHONY: clean source-clean
