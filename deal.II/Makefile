# $Id$
#
# Makefile for the whole library

# include global options and pathes
include Make.global_options


help:
	@echo "======================================================================"
	@echo "=Global Makefile for the deal.II libraries                           ="
	@echo "======================================================================"
	@echo "=Possible Targets                                                    ="
	@echo "=  baseg lacg: debug libraries base and lac                          ="
	@echo "=  base lac: debug and optimized libraries base and lac              ="
	@echo "=  1dg 2dg 3dg: (1) and deal.II debug version for specified dimension="
	@echo "=  1d 2d 3d: (2) and deal.II for specified dimension                 ="
	@echo "=  all: base, lac 1d, 2d and 3d                                      ="
	@echo "=  online-doc: generate the documentation in HTML format             ="
	@echo "=  printable-doc: generate the documentation in Postscript format.   ="
	@echo "=      You need to have Doc++ and LaTeX installed for that.          ="
	@echo "=  tex-doc: same, but only generate the LaTeX file, don't compile it ="
	@echo "=  TODO: create a TODO file from TODO lines in the source files      ="
	@echo "=  TAGS: create a TAGS file from include and source files            ="
	@echo "=  clean: removes all object files, libraries, etc in all subdirs    ="
	@echo "=                                                                    ="
	@echo "=The following additional parameters may be specified to make:       ="
	@echo "=  PAR=-jn where n is the number of parallel processes for parallel  ="
	@echo "=          builds. If you use GNU make 3.78 or higher, you can       ="
	@echo "=          call 'make -jn' directly. This parameter is no more       ="
	@echo "=          necessary then.                                           ="
	@echo "======================================================================"


MAKEOPTIONS = $(PAR)


baseg:
	cd $D/base ; $(MAKE) $(MAKEOPTIONS) libg

base:
	cd $D/base ; $(MAKE) $(MAKEOPTIONS) libg
	cd $D/base ; $(MAKE) $(MAKEOPTIONS) lib

lacg: 
	cd $D/lac ; $(MAKE) $(MAKEOPTIONS) libg

lac:
	cd $D/lac ; $(MAKE) $(MAKEOPTIONS) libg
	cd $D/lac ; $(MAKE) $(MAKEOPTIONS) lib

1dg:
	cd $D/deal.II ; $(MAKE) $(MAKEOPTIONS) 1dg

2dg:
	cd $D/deal.II ; $(MAKE) $(MAKEOPTIONS) 2dg

3dg:
	cd $D/deal.II ; $(MAKE) $(MAKEOPTIONS) 3dg

1d:
	cd $D/deal.II ; $(MAKE) $(MAKEOPTIONS) 1dg
	cd $D/deal.II ; $(MAKE) $(MAKEOPTIONS) 1d

2d:
	cd $D/deal.II ; $(MAKE) $(MAKEOPTIONS) 2dg
	cd $D/deal.II ; $(MAKE) $(MAKEOPTIONS) 2d

3d:
	cd $D/deal.II ; $(MAKE) $(MAKEOPTIONS) 3dg
	cd $D/deal.II ; $(MAKE) $(MAKEOPTIONS) 3d

all:
	cd $D/base ; $(MAKE) $(MAKEOPTIONS) libg
	cd $D/base ; $(MAKE) $(MAKEOPTIONS) lib
	cd $D/lac ; $(MAKE) $(MAKEOPTIONS) libg
	cd $D/lac ; $(MAKE) $(MAKEOPTIONS) lib
	cd $D/deal.II ; $(MAKE) $(MAKEOPTIONS) 2dg
	cd $D/deal.II ; $(MAKE) $(MAKEOPTIONS) 2d
	cd $D/deal.II ; $(MAKE) $(MAKEOPTIONS) 1dg
	cd $D/deal.II ; $(MAKE) $(MAKEOPTIONS) 1d
	cd $D/deal.II ; $(MAKE) $(MAKEOPTIONS) 3dg
	cd $D/deal.II ; $(MAKE) $(MAKEOPTIONS) 3d


online-doc:
	cd $D/doc ; $(MAKE)
	@echo
	@echo
	@echo =======================================================
	@echo "The online documentation can now be accessed through"
	@echo "  $D/doc/index.html"
	@echo =======================================================
	@echo
	@echo

printable-doc:
	cd $D/doc/auto ; $(MAKE) latex/deal_II_technical_reference.ps
	@echo
	@echo
	@echo =======================================================
	@echo "The generated postscript file can be found here:"
	@echo "  $D/doc/auto/latex/deal_II_technical_reference.ps"
	@echo =======================================================
	@echo
	@echo

tex-doc:
	cd $D/doc/auto ; $(MAKE) latex/deal_II_technical_reference.tex
	@echo
	@echo
	@echo =======================================================
	@echo "The generated LaTeX file can be found here:"
	@echo "  $D/doc/auto/latex/deal_II_technical_reference.tex"
	@echo =======================================================
	@echo
	@echo

TODO:
	@cat $D/base/Todo $D/lac/Todo $D/deal.II/Todo > $@
	@grep //TODO $D/base/include/*/*.h $D/base/source/*.cc \
                     $D/lac/include/*/*.h $D/lac/source/*.cc   \
                     $D/deal.II/include/*/*.h $D/deal.II/source/*/*.cc \
             >> $@


TAGS: $(filter-out %/forward_declarations.h,                                \
                  $(shell echo $D/base/source/*.cc $D/base/include/*/*.h    \
                               $D/lac/source/*.cc $D/lac/include/*/*.h      \
                               $D/deal.II/source/*/*.cc $D/deal.II/include/*/*.h \
                               $D/tests/*/*.cc $D/tests/*/*.h ))
	@etags $^


clean: clean-base clean-lac clean-dealII clean-doc
	-rm -f TODO TAGS

clean-base:
	cd $D/base ; $(MAKE) clean

clean-lac:
	cd $D/lac ; $(MAKE) clean

clean-dealII:
	cd $D/deal.II ; $(MAKE) clean

clean-doc:
	cd $D/doc ; $(MAKE) clean



.PHONY: base baseg lac lacg 1d 2d 3d 1dg 2dg 3dg all \
        online-doc printable-doc tex-doc \
        clean clean-base clean-lac clean-dealII clean-doc TAGS TODO
