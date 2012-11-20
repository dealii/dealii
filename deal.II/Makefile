default:

style-h-files:= $(wildcard include/deal.II/*/*.h)
style-cc-files:= $(wildcard source/*/*.cc)
style-examples:= $(wildcard examples/*/*.cc)

indent:
        @echo "============ Indenting all files"
        astyle --options=contrib/utilities/astyle.rc $(style-h-files) $(style-cc-files) $(style-examples)

