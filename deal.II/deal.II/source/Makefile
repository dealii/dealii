# $Id$

#this file uses many small rules to allow
#   make -j4
#on multiprocessor machines. Feel free to use this feature,
#it gives you quicker compile time and annoys those FEATFLOW
#people ;-)


all: lib.g.a lib.a

lib.g.a: libbasic.g.a libgrid.g.a libfe.g.a libnumerics.g.a
lib.a:   libbasic.a libgrid.a libfe.a libnumerics.a

libbasic.g.a:
	@cd basic   ; $(MAKE) lib.g.a

libgrid.g.a:
	@cd grid    ; $(MAKE) lib.g.a

libfe.g.a:
	@cd fe      ; $(MAKE) lib.g.a

libnumerics.g.a:
	@cd numerics; $(MAKE) lib.g.a

libbasic.a:
	@cd basic   ; $(MAKE) lib.a

libgrid.a:
	@cd grid    ; $(MAKE) lib.a

libfe.a:
	@cd fe      ; $(MAKE) lib.a

libnumerics.a:
	@cd numerics; $(MAKE) lib.a


clean:
	@cd basic   ; $(MAKE) clean
	@cd grid    ; $(MAKE) clean
	@cd fe      ; $(MAKE) clean
	@cd numerics; $(MAKE) clean



.PHONY: all lib.g.a lib.a clean
.PHONY: libbasic.g.a libgrid.g.a libfe.g.a libnumerics.g.a
.PHONY: libbasic.a libgrid.a libfe.a libnumerics.a



