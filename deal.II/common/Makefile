# $Id$
#
# Makefile for the whole library


help:
	@echo "========================================================================="
	@echo "=Global Makefile for the deal.II libraries                              ="
	@echo "========================================================================="
	@echo "=Possible Targets                                                       ="
	@echo "= (1) baseg lacg: debug libraries base and lac                          ="
	@echo "= (2) base lac: debug and optimized libraries base and lac              ="
	@echo "= (3) 1dg 2dg 3dg: (1) and deal.II debug version for specified dimension="
	@echo "= (4) 1d 2d 3d: (2) and deal.II for specified dimension                 ="
	@echo "= (5) all: 1d, 2d and 3d                                                ="
	@echo "= (6) online-doc: generate the documentation in HTML format             ="
	@echo "= (7) printable-doc: generate the documentation in Postscript format.   ="
	@echo "=     You need to have Doc++ and LaTeX installed for that.              ="
	@echo "= (8) tex-doc: same, but only generate the LaTeX file, don't compile it ="
	@echo "=                                                                       ="
	@echo "=The following additional parameters may be specified to make:          ="
	@echo "=  PAR=-jn where n is the number of parallel processes for parallel     ="
	@echo "=          builds.                                                      ="
	@echo "========================================================================="


MAKEOPTIONS = $(PAR)


baseg:
	cd ../base/source ; $(MAKE) $(MAKEOPTIONS) libg

base:
	cd ../base/source ; $(MAKE) $(MAKEOPTIONS) libg lib

lacg: baseg
	cd ../lac/source ; $(MAKE) $(MAKEOPTIONS) libg

lac: base
	cd ../lac/source ; $(MAKE) $(MAKEOPTIONS) libg lib

1d: baseg lacg base lac
	cd ../deal.II/source ; $(MAKE) $(MAKEOPTIONS) 1d

2d: base lac
	cd ../deal.II/source ; $(MAKE) $(MAKEOPTIONS) 2d

3d: base lac
	cd ../deal.II/source ; $(MAKE) $(MAKEOPTIONS) 3d

1dg: baseg lacg
	cd ../deal.II/source ; $(MAKE) $(MAKEOPTIONS) 1dg

2dg: baseg lacg
	cd ../deal.II/source ; $(MAKE) $(MAKEOPTIONS) 2dg

3dg: baseg lacg
	cd ../deal.II/source ; $(MAKE) $(MAKEOPTIONS) 3dg

all: base lac 2d 1d 3d

online-doc:
	cd ../doc/auto ; make

printable-doc:
	cd ../doc/auto ; make latex/deal_II_technical_reference.ps

tex-doc:
	cd ../doc/auto ; make latex/deal_II_technical_reference.tex


TAGS:
	etags ../base/source/*.cc ../base/include/*/*.h ../lac/source/*.cc ../lac/include/*/*.h ../deal.II/source/*/*.cc ../deal.II/include/*/*.h

.PHONY: base baseg lac lacg 1d 2d 3d 1dg 2dg 3dg all online-doc printable-doc tex-doc TAGS
