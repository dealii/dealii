# Generated automatically from Makefile.in by configure.
# $Id$
# Common Makefile for all test directories

all: base lac fe deal.II big-tests

base:
	cd base ; $(MAKE)

lac:
	cd lac ; $(MAKE)

fe:
	cd fe ; $(MAKE)

deal.II:
	cd deal.II ; $(MAKE)

report:
	@for i in base lac fe deal.II ; do \
		echo =======Report: $$i ======= ; \
		cd $$i ; $(MAKE) report ; cd .. ; \
	done


clean: clean-base clean-lac clean-deal.II clean-fe
veryclean: veryclean-base veryclean-lac veryclean-deal.II veryclean-fe

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

clean-big:
	cd big-tests ; $(MAKE) clean

veryclean-big:
	cd big-tests ; $(MAKE) veryclean

.PHONY : all base lac fe deal.II big-tests \
	clean clean-base clean-lac clean-fe clean-deal.II clean-big \
	veryclean veryclean-base veryclean-lac veryclean-fe veryclean-deal.II veryclean-big
