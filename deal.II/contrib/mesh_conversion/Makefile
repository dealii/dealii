SHELL = /bin/sh

#MACROS
CPP = g++ -g -ansi -pedantic -Wall
PROG = convert_mesh

#SUFFIX RULES
.SUFFIXES: .o .cc
.cc.o:
	$(CPP) -c $<

all: $(PROG)

$(PROG): Main.o MeshConversion.o
	$(CPP)	-o	$(PROG)	Main.o MeshConversion.o
Main.o:	Main.cc 	MeshConversion.o
MeshConversion.o:	MeshConversion.cc

clean:
	@echo "-------------------------------------------------"
	@echo "REMOVING EXECUTABLES, *.o, a.out AND FILE BACKUPS"
	-rm -f  ${PROG} *.o *.*~ makefile~ *~
	@echo "-------------------------------------------------"

clean-temp:
	@echo "-------------------------------------------------"
	@echo "REMOVING EXECUTABLES, *.o, a.out AND FILE BACKUPS"
	-rm -f  *.o *.*~ makefile~ *~
	@echo "-------------------------------------------------"
