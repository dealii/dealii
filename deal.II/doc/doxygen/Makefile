# $Id$
# $Name$
#
# Copyright (c) the deal.II authors 2003

html: deal.tag

%.tag: %.dox
	doxygen $<

base.tag: ../../base/include/*/*.h headers/base/*.h
lac.tag:  base.tag ../../lac/include/*/*.h headers/lac/*.h
deal.tag: lac.tag  ../../deal.II/include/*/*.h headers/deal.II/*.h

base:    base.tag
lac:     lac.tag
deal.II: deal.tag

.PHONY: html base lac deal.II
