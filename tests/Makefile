# Generated automatically from Makefile.in by configure.
# Makefile,v 1.14 2002/09/19 13:53:47 wolf Exp
# Common Makefile for all test directories

all: base lac fe deal.II bits

base:
	cd base ; $(MAKE)

lac:
	cd lac ; $(MAKE)

fe:
	cd fe ; $(MAKE)

bits:
	cd bits ; $(MAKE)

deal.II:
	cd deal.II ; $(MAKE)

report:
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


clean: clean-base clean-lac clean-deal.II clean-fe clean-bits
veryclean: veryclean-base veryclean-lac veryclean-deal.II veryclean-fe veryclean-bits

clean-base:
	cd base ; $(MAKE) clean

veryclean-base:
	cd base ; $(MAKE) veryclean

clean-lac:
	cd lac ; $(MAKE) clean

veryclean-lac:
	cd lac ; $(MAKE) veryclean

clean-fe:
	cd fe ; $(MAKE) clean

veryclean-fe:
	cd fe ; $(MAKE) veryclean

clean-deal.II:
	cd deal.II ; $(MAKE) clean

veryclean-deal.II:
	cd deal.II ; $(MAKE) veryclean

clean-bits:
	cd bits ; $(MAKE) clean

veryclean-bits:
	cd bits ; $(MAKE) veryclean


.PHONY : all base lac fe deal.II bits \
	clean clean-base clean-lac clean-fe clean-deal.II clean-bits \
	veryclean veryclean-base veryclean-lac veryclean-fe veryclean-deal.II veryclean-bits

include Makefile.paths
include $D/common/Make.global_options

