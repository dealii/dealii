# $Id$

all: lib.g.a lib.a

lib.g.a:
	@cd basic   ; make lib.g.a
	@cd grid    ; make lib.g.a
	@cd fe      ; make lib.g.a
	@cd numeric ; make lib.g.a


lib.a:
	@cd basic   ; make lib.a
	@cd grid    ; make lib.a
	@cd fe      ; make lib.a
	@cd numeric ; make lib.a

clean:
	@cd basic   ; make clean
	@cd grid    ; make clean
	@cd fe      ; make clean
	@cd numeric ; make clean



.PHONY: all lib.g.a lib.a clean
