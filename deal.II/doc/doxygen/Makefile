# $Id$
# $Name$
#
# Copyright (c) the deal.II authors 2003

html: deal.tag

private:
	perl -pi -e 's/(EXTRACT_PRIVATE\s*=\s*)NO/$$1YES/;s/(INTERNAL_DOCS\s*=\s*)NO/$$1YES/;' options.dox
	$(MAKE) html

%.tag: %.dox
	doxygen $<

base.tag: options.dox base.dox ../../base/include/*/*.h headers/base/*.h
lac.tag:  options.dox lac.dox base.tag ../../lac/include/*/*.h headers/lac/*.h
deal.tag: options.dox deal.dox lac.tag  ../../deal.II/include/*/*.h headers/deal.II/*.h

base:    base.tag
lac:     lac.tag
deal.II: deal.tag

.PHONY: html private base lac deal.II
