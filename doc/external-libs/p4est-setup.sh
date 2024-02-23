#! /bin/bash
#
# SPDX-License-Identifier: GPL-2.0-or-later
#
# This file is part of p4est [1].
# p4est is a C library to manage a collection (a forest) of multiple
# connected adaptive quadtrees or octrees in parallel.
#
# Copyright (C) 2010 The University of Texas System
# Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac
# Modified 2010 by Wolfgang Bangerth
# Modified 2010 by Timo Heister
# Modified 2013 by Matthias Maier
#
# p4est is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# [1] https://www.p4est.org
#

# This program comes with ABSOLUTELY NO WARRANTY.

# error message when zlib is not found
MISSING_ZLIB_MESSAGE="deal.II requires that p4est be built with zlib support. Please \
either ensure that zlib is installed in a standard location or add appropriate \
flags to LDFLAGS and CPPFLAGS to both calls to configure to describe where zlib's \
shared object files and headers are (e.g., LDFLAGS=\"-L/path/to/shared-objects/\" \
and CPPFLAGS=\"-DSC_LOG_PRIORITY=SC_LP_ESSENTIAL -I/path/to/headers/\")."

# unpack under current directory
UNPACK=`pwd`

# choose names for fast and debug compilation directories
BUILD_DIR="$UNPACK/p4est-build"
BUILD_FAST="$BUILD_DIR/FAST"
BUILD_DEBUG="$BUILD_DIR/DEBUG"

function busage() {
        echo "Usage: `basename $0` <p4est_tar.gz_file> [<install location>]"
        echo "   or: `basename $0` /path/to/p4est-src/ [<install location>]"
}
function bdie () {
        echo "Error: $@"
        exit 1
}

if test -z "$CFLAGS" -a -z "$P4EST_CFLAGS_FAST" ; then
        export CFLAGS_FAST="-O2"
else
        export CFLAGS_FAST="$CFLAGS $P4EST_CFLAGS_FAST"
fi
echo "CFLAGS_FAST: $CFLAGS_FAST"
if test -z "$CFLAGS" -a -z "$P4EST_CFLAGS_DEBUG" ; then
        export CFLAGS_DEBUG="-O0 -g"
else
        export CFLAGS_DEBUG="$CFLAGS $P4EST_CFLAGS_DEBUG"
fi
echo "CFLAGS_DEBUG: $CFLAGS_DEBUG"

TGZ="$1"; shift
if test -d "$TGZ" ; then
  SRCDIR="$TGZ"  
  echo "using existing source dir '$SRCDIR'"
else
    if test ! -f "$TGZ" ; then
        busage
        bdie "File not found"
    fi
    if ! (echo "$TGZ" | grep -q 'p4est.*.tar.gz') ; then
        busage
        bdie "File name mismatch"
    fi
fi

# choose names for fast and debug installation directories
INSTALL_DIR="$1"; shift
if test -z "$INSTALL_DIR" ; then
        INSTALL_DIR="$UNPACK/p4est-install"
fi
INSTALL_FAST="$INSTALL_DIR/FAST"
INSTALL_DEBUG="$INSTALL_DIR/DEBUG"

echo
echo "This script tries to unpack, configure and build the p4est library."
echo "Build FAST: $BUILD_FAST"
echo "Build DEBUG: $BUILD_DEBUG"
echo "Install FAST: $INSTALL_FAST"
echo "Install DEBUG: $INSTALL_DEBUG"
echo "Checking environment: CFLAGS P4EST_CFLAGS_FAST P4EST_CFLAGS_DEBUG"

# remove old versions
if test -d "$BUILD_DIR" ; then
        rm -rf "$BUILD_DIR"
fi

if test -f "$TGZ" ; then
    DIR=`echo "$TGZ" | sed 's/\(p4est.*\).tar.gz/\1/'`
    DIR=`basename $DIR`
    echo "Unpack directory: $UNPACK/$DIR"
    if test -d "$UNPACK/$DIR" ; then
        echo "Source directory found (remove it to unpack anew)"
    else
        echo -n "Unpacking... "
        tar -xvz -f "$TGZ" -C "$UNPACK" >/dev/null
        echo "done"
    fi
    SRCDIR=$UNPACK/$DIR
fi
test -f "$SRCDIR/src/p4est.h" || bdie "Main header file missing"
test -f "$SRCDIR/configure" || bdie "Configure script missing"

echo "See output in files .../config.output and .../make.output"
echo
echo "Build FAST version in $BUILD_FAST"
mkdir -p "$BUILD_FAST"
cd "$BUILD_FAST"
"$SRCDIR/configure" --enable-mpi --enable-shared \
        --disable-vtk-binary --without-blas \
        --prefix="$INSTALL_FAST" CFLAGS="$CFLAGS_FAST" \
        CPPFLAGS="-DSC_LOG_PRIORITY=SC_LP_ESSENTIAL" \
        "$@" > config.output || bdie "Error in configure"
make -C sc -j 8 > make.output || bdie "Error in make sc"
make -j 8 >> make.output || bdie "Error in make p4est"
# ensure that we built p4est with zlib
find "$BUILD_FAST" -name "p4est_config.h" -type f -exec \
    grep -q "P4EST_HAVE_ZLIB *1" {} \; \
    || bdie "$MISSING_ZLIB_MESSAGE"
make install >> make.output || bdie "Error in make install"
echo "FAST version installed in $INSTALL_FAST"

echo
echo "Build DEBUG version in $BUILD_DEBUG"
mkdir -p "$BUILD_DEBUG"
cd "$BUILD_DEBUG"
"$SRCDIR/configure" --enable-debug --enable-mpi --enable-shared \
        --disable-vtk-binary --without-blas \
        --prefix="$INSTALL_DEBUG" CFLAGS="$CFLAGS_DEBUG" \
        CPPFLAGS="-DSC_LOG_PRIORITY=SC_LP_ESSENTIAL" \
        "$@" > config.output || bdie "Error in configure"
make -C sc -j 8 > make.output || bdie "Error in make sc"
make -j 8 >> make.output || bdie "Error in make p4est"
# ensure that we built p4est with zlib
find "$BUILD_DEBUG" -name "p4est_config.h" -type f -exec \
    grep -q "P4EST_HAVE_ZLIB *1" {} \; \
    || bdie "$MISSING_ZLIB_MESSAGE"
make install >> make.output || bdie "Error in make install"
echo "DEBUG version installed in $INSTALL_DEBUG"
echo
