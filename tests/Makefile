# Generated automatically from Makefile.in by configure.
# Makefile,v 1.14 2002/09/19 13:53:47 wolf Exp
# Common Makefile for all test directories

all: compare base lac fe deal.II bits

base: compare abort.o
	cd base ; $(MAKE)

lac: compare abort.o
	cd lac ; $(MAKE)

fe: compare abort.o
	cd fe ; $(MAKE)

bits: compare abort.o
	cd bits ; $(MAKE)

deal.II: compare abort.o
	cd deal.II ; $(MAKE)

report: compare abort.o
	@for i in base lac fe deal.II bits ; do \
		echo =======Report: $$i ======= ; \
		cd $$i ; $(MAKE) report ; cd .. ; \
	 done
	@echo =======Report: all-headers =======
	@today=`date -u +"%Y-%m-%d %H:%M"` ; \
	 if (($(MAKE) all-headers 2>&1) > /dev/null) ; then \
		echo $$today + all-headers ; \
	 else \
		echo $$today - all-headers ; \
	 fi

report+mail:
	@$(MAKE) report | tee testresults
	@cat testresults | mail regression-tests@dealii.org
	@rm testresults


# collect all header files, and for each try to do the following:
# generate a file tmp.cc that contains nothing but "#include <xxx>"
# with this particular header file. then try to compile it. this 
# sometimes fails since this header does not include all other headers
# it actually needs. we don't usually see this since we may be 
# including these missing headers into our files before, so that nobody
# notices the problem. expose it this way.
#
# this all happens in a subdir tmp-all-headers, which is removed
# afterwards again.
all-headers:
	@rm -rf tmp-all-headers
	@mkdir tmp-all-headers 
	@for i in base lac deal.II ; do \
		pushd . > /dev/null ; \
		headers=`cd $D/$$i/include ; \
			 for j in */*.h ; do echo $$j ; done` ; \
		popd > /dev/null ; \
		cd tmp-all-headers ; \
		cp ../Makefile.headers Makefile ; \
		for j in $$headers ; do \
			echo "=====Header======== $$j " ; \
			echo "#include <$$j>" > tmp.cc ; \
			if $(MAKE) -s tmp.g.$(OBJEXT) ; then : ; else exit 1 ; fi ; \
			rm -f tmp.* ; \
		done ; \
		cd .. ; \
	 done
	@rm -rf tmp-all-headers

# Generate a link to system specific test program output

compare:
	cd results; $(MAKE)
	ln -s results/$(TARGET)+$(GXX-VERSION) compare

resultname:
	@echo $(TARGET)+$(GXX-VERSION)

clean: clean-base clean-lac clean-deal.II clean-fe clean-bits
distclean: distclean-base distclean-lac distclean-deal.II distclean-fe distclean-bits

clean-base:
	cd base ; $(MAKE) clean

distclean-base:
	cd base ; $(MAKE) distclean

clean-lac:
	cd lac ; $(MAKE) clean

distclean-lac:
	cd lac ; $(MAKE) distclean

clean-fe:
	cd fe ; $(MAKE) clean

distclean-fe:
	cd fe ; $(MAKE) distclean

clean-deal.II:
	cd deal.II ; $(MAKE) clean

distclean-deal.II:
	cd deal.II ; $(MAKE) distclean

clean-bits:
	cd bits ; $(MAKE) clean

distclean-bits:
	cd bits ; $(MAKE) distclean

abort.o: abort.cc
	 @$(CXX) -ggdb -o $@ -c $<

.PHONY : all base lac fe deal.II bits \
	clean clean-base clean-lac clean-fe clean-deal.II clean-bits \
	distclean distclean-base distclean-lac distclean-fe distclean-deal.II distclean-bits

include Makefile.paths
include $D/common/Make.global_options

