# $Id$

# include defaults file
D = ../..
include $D/common/Make.global_options


default: index.html

index.html: extract.pl glossary.data
	$(PERL) $^ > $@


validate-xrefs:
	$(PERL) $D/doc/auto/scripts/validate-xrefs.pl index.html


clean:
	-rm -f index.html



.PHONY: clean default validate-xrefs
