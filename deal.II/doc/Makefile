# $Id$

default: autogen-doc tutorial


autogen-doc:
	cd auto ; $(MAKE)


autogen-doc-all: 
	cd auto ; $(MAKE) all


tutorial:
	cd tutorial ; $(MAKE)


clean:
	cd auto ; $(MAKE) clean
	cd tutorial ; $(MAKE) clean

.PHONY: default autogen-doc autogen-doc-all tutorial clean
