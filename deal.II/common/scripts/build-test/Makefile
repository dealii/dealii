sources=$(wildcard configurations/*.cfg)
#targets=$(patsubst %.cfg,%.bash,$(notdir $(sources)))
targets=$(patsubst %.cfg,build-%,$(notdir $(sources)))

vpath %.cfg configurations

all: $(targets)

build-%: %.cfg template
	perl -p $< template > $@
