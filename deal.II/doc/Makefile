# $Id$


D = ..
include $D/common/Make.global_options



# generic targets
default: autogen-doc tutorial development glossary validate-xrefs

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

glossary:
	cd glossary ; $(MAKE)


validate-xrefs:
	@echo "Validating cross-references in HTML files."
	$(PERL) $D/doc/auto/scripts/validate-xrefs.pl *html
	cd tutorial ; $(MAKE) validate-xrefs
	cd development ; $(MAKE) validate-xrefs
	cd glossary ; $(MAKE) validate-xrefs
	cd publications ; $(MAKE) validate-xrefs


clean:
	cd auto ; $(MAKE) clean
	cd tutorial ; $(MAKE) clean
	cd development ; $(MAKE) clean
	cd glossary ; $(MAKE) clean


.PHONY: default all autogen-doc autogen-doc-all tutorial clean development validate-xrefs glossary

