D = ../..

dealparse: lex.yy.c dealparse.tab.c
	g++ -o $@ $^ -lfl

dealparse.tab.h: dealparse.y
	bison -d -t -v $<

lex.yy.c: dealparse.l dealparse.tab.h
	flex $<

database: dealparse
	cat $D/base/include/*/*.h $D/lac/include/*/*.h $D/deal.II/include/*/*.h | dealparse > $@

postgres: database
	perl pass2 database
