# $Id$


# generic targets
default: autogen-doc tutorial development validate-xrefs glossary

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

glossary: glossary/index.html

validate-xrefs:
	@echo "Validating cross-references in HTML files."
	cd tutorial ; $(MAKE) validate-xrefs
	cd development ; $(MAKE) validate-xrefs

clean:
	cd auto ; $(MAKE) clean
	cd tutorial ; $(MAKE) clean
	cd development ; $(MAKE) clean

.PHONY: default all autogen-doc autogen-doc-all tutorial clean development validate-xrefs glossary

glossary/index.html: glossary/extract.pl glossary/glossary.data
	perl $^ > $@
