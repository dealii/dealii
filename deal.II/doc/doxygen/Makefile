# $Id$
# $Name$
#
# Copyright (c) the deal.II authors 2003

html: deal.tag

private:
	perl -pi -e 's/(EXTRACT_PRIVATE\s*=\s*)NO/$$1YES/' options.dox
	$(MAKE) html

%.tag: %.dox
	doxygen $<

base.tag: ../../base/include/*/*.h headers/base/*.h
lac.tag:  base.tag ../../lac/include/*/*.h headers/lac/*.h
deal.tag: lac.tag  ../../deal.II/include/*/*.h headers/deal.II/*.h

base:    base.tag
lac:     lac.tag
deal.II: deal.tag

.PHONY: html private base lac deal.II
