# $Id$
#
# Makefile for the whole library


help:
	@echo Possible Targets
	@echo " (1) baseg lacg: debug libraries base and lac"
	@echo " (2) base lac: debug and optimized libraries base and lac"
	@echo " (3) 1dg 2dg 3dg: (1) and deal.II debug version for specified dimension"
	@echo " (4) 1d 2d 3d: (2) and deal.II for specified dimension"
	@echo " (5) all: 1d, 2d and 3d"
	@echo
#	@echo The following additional parameters may be specified to make:
#	@echo 
#	@echo

baseg:
	cd ../base/source ; $(MAKE) lib.g.a

base:
	cd ../base/source ; $(MAKE) lib.g.a lib.a

lacg: baseg
	cd ../lac/source ; $(MAKE) lib.g.a

lac: base
	cd ../lac/source ; $(MAKE) lib.g.a lib.a

1d: baseg lacg base lac
	cd ../deal.II/source ; $(MAKE) 1d

2d: base lac
	cd ../deal.II/source ; $(MAKE) 2d

3d: base lac
	cd ../deal.II/source ; $(MAKE) 3d

1dg: baseg lacg
	cd ../deal.II/source ; $(MAKE) 1dg

2dg: baseg lacg
	cd ../deal.II/source ; $(MAKE) 2dg

3dg: baseg lacg
	cd ../deal.II/source ; $(MAKE) 3dg

all: base lac 2d 1d 3d
