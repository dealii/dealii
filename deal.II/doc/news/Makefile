# $Id$

# include defaults file
D = ../..
include $D/common/Make.global_options


validate-xrefs:
	$(PERL) $D/doc/auto/scripts/validate-xrefs.pl *.html
	for dir in 200* ; do cd $$dir ; $(MAKE) validate-xrefs ; cd .. ; done


.PHONY: validate-xrefs


