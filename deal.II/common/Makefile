# $Id$
#
# Makefile for the whole library


help:
	@echo "======================================================================"
	@echo "=Global Makefile for the deal.II libraries                           ="
	@echo "======================================================================"
	@echo "=Possible Targets                                                    ="
	@echo "=  baseg lacg: debug libraries base and lac                          ="
	@echo "=  base lac: debug and optimized libraries base and lac              ="
	@echo "=  1dg 2dg 3dg: (1) and deal.II debug version for specified dimension="
	@echo "=  1d 2d 3d: (2) and deal.II for specified dimension                 ="
	@echo "=  all: 1d, 2d and 3d                                                ="
	@echo "=  online-doc: generate the documentation in HTML format             ="
	@echo "=  printable-doc: generate the documentation in Postscript format.   ="
	@echo "=      You need to have Doc++ and LaTeX installed for that.          ="
	@echo "=  tex-doc: same, but only generate the LaTeX file, don't compile it ="
	@echo "=  TODO: create a TODO file from TODO lines in the source files      ="
	@echo "=  clean: removes all object files, libraries, etc in all subdirs    ="
	@echo "=                                                                    ="
	@echo "=The following additional parameters may be specified to make:       ="
	@echo "=  PAR=-jn where n is the number of parallel processes for parallel  ="
	@echo "=          builds.                                                   ="
	@echo "======================================================================"


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
	cd ../doc/auto ; $(MAKE)
	@echo
	@echo
	@echo =======================================================
	@echo "The online documentation can now be found in"
	@echo "  ../doc/auto/kdoc/base/index.html"
	@echo "  ../doc/auto/kdoc/lac/index.html"
	@echo "  ../doc/auto/kdoc/deal.II/index.html"
	@echo =======================================================
	@echo
	@echo

printable-doc:
	cd ../doc/auto ; $(MAKE) latex/deal_II_technical_reference.ps
	@echo
	@echo
	@echo =======================================================
	@echo "The generated postscript file can be found here:"
	@echo "  ../doc/auto/latex/deal_II_technical_reference.ps"
	@echo =======================================================
	@echo
	@echo

tex-doc:
	cd ../doc/auto ; $(MAKE) latex/deal_II_technical_reference.tex
	@echo
	@echo
	@echo =======================================================
	@echo "The generated LaTeX file can be found here:"
	@echo "  ../doc/auto/latex/deal_II_technical_reference.tex"
	@echo =======================================================
	@echo
	@echo

TODO:
	@cat ../base/Todo ../lac/Todo ../deal.II/Todo > TODO
	@grep //TODO ../base/include/*/*.h ../base/source/*.cc  ../lac/include/*/*.h ../lac/source/*.cc  ../deal.II/include/*/*.h ../deal.II/source/*/*.cc >> TODO
#	@grep \#TODO Makefile Make.g* Makefile.template ../*/source/Makefile ../*/Makefile >> TODO

clean:
	cd ../base ; $(MAKE) clean
	cd ../lac ; $(MAKE) clean
	cd ../deal.II ; $(MAKE) clean
	cd ../doc/auto ; $(MAKE) clean
	-rm TODO TAGS


TAGS:
	etags ../base/source/*.cc ../base/include/*/*.h ../lac/source/*.cc ../lac/include/*/*.h ../deal.II/source/*/*.cc ../deal.II/include/*/*.h

.PHONY: base baseg lac lacg 1d 2d 3d 1dg 2dg 3dg all online-doc printable-doc tex-doc clean TAGS TODO
