dealparse: lex.yy.c dealparse.tab.c
	g++ -o $@ $^ -lfl

dealparse.tab.h: dealparse.y
	bison++ -d -t -v $<

lex.yy.c: dealparse.l dealparse.tab.h
	flex $<
