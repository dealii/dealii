# $Id$


# generic targets
default: autogen-doc tutorial

all: default autogen-doc-all


# specific targets
autogen-doc:
	cd auto ; $(MAKE)


autogen-doc-all: 
	cd auto ; $(MAKE) all


tutorial:
	cd tutorial ; $(MAKE)


clean:
	cd auto ; $(MAKE) clean
	cd tutorial ; $(MAKE) clean

.PHONY: default all autogen-doc autogen-doc-all tutorial clean
