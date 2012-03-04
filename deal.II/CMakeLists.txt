# -----------------------------------------------------------------
#      Configuration and build system description for deal.II
#
# This file is subject to QPL and may not be  distributed
# without copyright and license information. Please refer
# to the file deal.II/doc/license.html for the  text  and
# further information on this license.
#
# Copyright (C) 2012 by the deal.II authors
#
# $Id: aclocal.m4 24924 2012-01-25 12:35:17Z kormann $


# ###########################################################
#                Configuration of deal.II                   #
# ###########################################################

message ("")
message ("=========================================")
message ("==         Configuring deal.II         ==")
message ("=========================================")
message ("")

message ("------ Configuring the compilers")

cmake_minimum_required(VERSION 2.8)
project(deal.II)


# -----------------------------------------------------------
# Find a few programs that we need throughout this script and
# deal.II.
# -----------------------------------------------------------

message ("")
message ("------ Searching for external programs")
find_program (perl perl)
if (${perl} STREQUAL perl-NOTFOUND)
  message (FATAL_ERROR "** Searching for perl: not found")
else()
  message ("-- Searching for perl: " ${perl})
endif()

find_program (doxygen doxygen)
if (${doxygen} STREQUAL doxygen-NOTFOUND)
  message (FATAL_ERROR "** Searching for doxygen: not found")
else()
  message ("-- Searching for doxygen: " ${doxygen})
endif()



# ###########################################################
#                 Build system of deal.II                   #
# ###########################################################

# common defines
include_directories(include)
include_directories(contrib/boost-1.46.1/include)
include_directories(contrib/tbb/tbb30_104oss/include)
include_directories(contrib)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_CONFIG_H -DHAVE_ISNAN -ggdb  -DBOOST_NO_HASH -DBOOST_NO_SLIST -DDEBUG -Wall -W -Wpointer-arith -Wwrite-strings -Wsynth -Wsign-compare -Wswitch -ftemplate-depth-128 -Wno-long-long -std=c++0x -Wa,--compress-debug-sections -pthread -D_REENTRANT -fPIC")

set (BUILD_SHARED_LIBS on)


# define what's needed for the program that builds the .inst files from
# the .inst.in files
set (EXPAND_INSTANTIATIONS_SOURCES common/scripts/expand_instantiations.cc)
add_executable (expand_instantiations
                ${EXPAND_INSTANTIATIONS_SOURCES})
get_target_property (expand_instantiations_exe
                     expand_instantiations
                     LOCATION)


# we then need to describe how we can build the .inst files from the
# .inst.in files
file (GLOB_RECURSE inst_in_files source/*.inst.in)
foreach (inst_file ${inst_in_files})
  add_custom_command (OUTPUT ${PROJECT_BINARY_DIR}/${inst_file}
                      DEPENDS expand_instantiations ${inst_file}.in
                      COMMAND ${expand_instantiations_exe}
                      ARGS common/template-arguments < ${inst_file}.in
                                                     > ${PROJECT_BINARY_DIR}/${inst_file})
  set (inst_files ${inst_files} ${PROJECT_BINARY_DIR}/${inst_file})
endforeach()

# then do the same for the .inst.pl files
file (GLOB_RECURSE inst_pl_files source/*.inst.pl)
foreach (inst_file ${inst_pl_files})
  add_custom_command (OUTPUT ${PROJECT_BINARY_DIR}/${inst_file}
                      DEPENDS ${inst_file}.pl
                      COMMAND ${perl}
                      ARGS -I common ${inst_file}.pl
                                     > ${PROJECT_BINARY_DIR}/${inst_file})
  set (inst_files ${inst_files} ${PROJECT_BINARY_DIR}/${inst_file})
endforeach ()


# then define what's needed for the deal.II library
file (GLOB_RECURSE cc_files source/*.cc)
add_library(deal_II ${inst_files} ${cc_files})
