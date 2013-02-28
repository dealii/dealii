NPROCS:=$(shell grep -c ^processor /proc/cpuinfo)
ifeq ($(NPROCS), )
  NPROCS:=1
endif

bin/cmake:
	wget -c "http://www.cmake.org/files/v2.8/cmake-2.8.8.tar.gz"
	tar -xf cmake-2.8.8.tar.gz
	cd cmake-2.8.8 && ./configure --parallel=$(NPROCS) --prefix="$(PWD)"
	cd cmake-2.8.8 && make -j $(NPROCS) install

.PHONY: clean
clean:
	rm -rf bin cmake-2.8.8 cmake-2.8.8.tar.gz doc man share

all: bin/cmake
