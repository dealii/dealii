# $Id$


# generic targets
default: autogen-doc tutorial development

all: default autogen-doc-all


# specific targets
autogen-doc:
	cd auto ; $(MAKE)


autogen-doc-all: 
	cd auto ; $(MAKE) all


tutorial:
	cd tutorial ; $(MAKE)

development:
	cd development ; $(MAKE)

clean:
	cd auto ; $(MAKE) clean
	cd tutorial ; $(MAKE) clean
	cd development ; $(MAKE) clean

.PHONY: default all autogen-doc autogen-doc-all tutorial clean development
