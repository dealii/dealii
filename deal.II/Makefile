default:

style-h-files:= $(wildcard include/deal.II/*/*.h)
style-cc-files:= $(wildcard source/*/*.cc)
style-examples:= $(wildcard examples/*/*.cc)

indent:
	@echo "============ Indenting all files"
	@for f in  $(style-h-files) $(style-cc-files) $(style-examples) ;\
	do\
	  echo $$f ;\
	  astyle --options=common/astyle.rc $$f ;\
	done

.PHONY: indent

