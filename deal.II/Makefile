# $Id$
#
# Makefile for the whole library

# include global options and pathes
include common/Make.global_options


help:
	@echo "========================================================================"
	@echo "=              Global Makefile for the deal.II libraries               ="
	@echo "========================================================================"
	@echo "=                                                                      ="
	@echo "= The following targets exist:                                         ="
	@echo "=    baseg lacg : debug libraries \"base\" and \"lac\"                     ="
	@echo "=    baseo laco : optimized libraries \"base\" and \"lac\"                 ="
	@echo "=    base lac   : debug and optimized libraries \"base\" and \"lac\"       ="
	@echo "=    1dg 2dg 3dg: deal.II debug version for specified dimension        ="
	@echo "=    1do 2do 3do: deal.II optimized version for specified dimension    ="
	@echo "=    1d 2d 3d   : deal.II debug and optimized for specified dimension  ="
	@echo "=                                                                      ="
	@echo "=    contrib    : additional libraries in contrib, if there are any    ="
	@echo "=                                                                      ="
	@echo "=    all        : base, lac, 1d, 2d and 3d, contrib                    ="
	@echo "=    debug      : baseg, lacg, 1dg, 2dg and 3dg, contrib               ="
	@echo "=    optimized  : base, lac, 1d, 2d, and 3d optimized libraries only   ="
	@echo "=                                                                      ="
	@echo "=    online-doc : generate the documentation in HTML format            ="
	@echo "=    TODO       : create a \"TODO\" file from the source files           ="
	@echo "=    TAGS       : create a TAGS file from include and source files     ="
	@echo "=                                                                      ="
	@echo "=    clean      : removes all object files in subdirs                  ="
	@echo "=    distclean  : removes all object files, libraries, etc in subdirs  ="
	@echo "========================================================================"



baseg:
	cd $D/base && $(MAKE) $(MAKEOPTIONS) libg

baseo:
	cd $D/base && $(MAKE) $(MAKEOPTIONS) libo

lacg: baseg contrib
	cd $D/lac && $(MAKE) $(MAKEOPTIONS) libg

laco: baseo contrib
	cd $D/lac && $(MAKE) $(MAKEOPTIONS) libo

1dg: baseg lacg
	cd $D/deal.II && $(MAKE) $(MAKEOPTIONS) 1dg

2dg: baseg lacg
	cd $D/deal.II && $(MAKE) $(MAKEOPTIONS) 2dg

3dg: baseg lacg
	cd $D/deal.II && $(MAKE) $(MAKEOPTIONS) 3dg

1do:baseo laco
	cd $D/deal.II && $(MAKE) $(MAKEOPTIONS) 1d

2do:baseo laco
	cd $D/deal.II && $(MAKE) $(MAKEOPTIONS) 2d

3do: baseo laco
	cd $D/deal.II && $(MAKE) $(MAKEOPTIONS) 3d

deps: 
	cd $D/base && $(MAKE) $(MAKEOPTIONS) Makefile.dep
	cd $D/lac && $(MAKE) $(MAKEOPTIONS) Makefile.dep
	cd $D/deal.II && $(MAKE) $(MAKEOPTIONS) Makefile.dep

base: baseg baseo
lac: lacg laco
1d: 1dg 1do
2d: 2dg 2do
3d: 3dg 3do

$(LIBDIR):
	mkdir $@

baseg base lacg lac 1dg 1d 2dg 2d 3dg 3d all debug contrib : $(LIBDIR)

all: deps
	$(MAKE) $(MAKEOPTIONS) contrib base lac 2d 1d 3d

debug: deps
	$(MAKE) $(MAKEOPTIONS) contrib baseg lacg 2dg 1dg 3dg

optimized: deps
	$(MAKE) $(MAKEOPTIONS) contrib baseo laco 2do 1do 3do

# when HSL functions are used, we also compile the detached_ma27 program
# which itself uses the optimized base library
contrib: $(subst no,,$(subst yes,baseo,$(USE_CONTRIB_HSL)))
	cd $D/contrib && $(MAKE)


online-doc:
	cd $D/doc && $(MAKE)
	@echo
	@echo
	@echo =======================================================
	@echo "The online documentation can now be accessed through"
	@echo "  $D/doc/index.html"
	@echo =======================================================
	@echo
	@echo

TODO:
	@egrep -i '// *todo' $D/base/include/*/*.h $D/base/source/*.cc \
                     $D/lac/include/*/*.h $D/lac/source/*.cc   \
                     $D/deal.II/include/*/*.h $D/deal.II/source/*/*.cc \
	| perl -pi -e 's#$D/?#\n#;' \
	| perl -pi -e 's#(.*)(TODO:?)(\[.+\])#\3 \1\2\3#i;' \
	| perl -pi -e 's#\s*//\s*TODO:?\s*(\[.+\]):?\s*#\n\1     #i;' \
             > $@


TAGS: $(shell echo $D/base/source/*.cc $D/base/include/*/*.h    \
                   $D/lac/source/*.cc $D/lac/include/*/*.h      \
                   $D/deal.II/source/*/*.cc $D/deal.II/include/*/*.h \
                   $D/tests/*/*.cc $D/tests/*/*.h )
	@cd $D/common ; etags $^


clean: clean-contrib clean-base clean-lac clean-dealII clean-doc clean-examples clean-lib clean-tests
	-rm -f TODO common/TAGS

clean-contrib:
	cd $D/contrib && $(MAKE) clean

clean-base:
	cd $D/base && $(MAKE) clean

clean-lac:
	cd $D/lac && $(MAKE) clean

clean-dealII:
	cd $D/deal.II && $(MAKE) clean

clean-doc:
	cd $D/doc && $(MAKE) clean

clean-examples:
	cd $D/examples && $(MAKE) clean

clean-lib:
	cd $D/lib && $(MAKE) clean

# if that directory exists, go into 'tests' and clean up there as well
clean-tests:
	-cd tests && $(MAKE) clean

distclean: clean
	cd lib && $(MAKE) distclean

.PHONY: base baseg lac lacg 1d 2d 3d 1dg 2dg 3dg all \
        online-doc printable-doc tex-doc contrib \
        clean clean-contrib clean-base clean-lac clean-dealII \
	clean-doc clean-examples clean-lib clean-tests distclean \
	TAGS TODO
