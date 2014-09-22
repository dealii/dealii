## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2014 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------

#
# General setup for the Intel C++ Compiler
#
# Please read the fat note in setup_compiler_flags.cmake prior to
# editing this file.
#

IF(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "10.0" )
  MESSAGE(WARNING "\n"
    "You're using an old version of the Intel C++ Compiler (icc/icpc)!\n"
    "It is strongly recommended to use at least version 10.\n"
    )
ENDIF()


########################
#                      #
#    General setup:    #
#                      #
########################

#
# Set the pic flag.
#
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-fpic")

#
# Check whether the -as-needed flag is available. If so set it to link
# the deal.II library with it.
#
ENABLE_IF_LINKS(DEAL_II_LINKER_FLAGS "-Wl,--as-needed")

#
# Set ansi mode:
#
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-ansi")

#
# Enable verbose warnings:
#
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-w2")

#
# Disable some warnings that lead to a lot of false positives:
#
#   -w68   integer conversion resulted in a change of sign
#          (triggers a lot in functionparser)
#   -w175  subscript out of range
#   -w135  class template "dealii::FE_Q_Base<POLY, dim, spacedim>"
#          has no member "Implementation"
#          (the compiler is objectively wrong since the warning
#           triggers also on code of the form
#           class FE_Q_Base {
#             struct Implementation; // forward declaration
#             friend struct Implementation;
#           };)
#   -w177  declared but not referenced
#   -w191  type qualifier is meaningless on cast type
#          Warnings from this warn about code like this:
#             static_cast<T const * const>(p)
#          There are many places in boost that do this kind of stuff
#   -w193  zero used for undefined preprocessing identifier "..."
#          This happens when using undefined preprocessor names in
#          conditions such as
#            #if (abc && def)
#          instead of
#            #if (defined(abc) && defined(def))
#          The standard says that in such cases, the undefined symbol
#          is assumed to be zero. The warning is in principle
#          useful, but the pattern appears exceedingly often in the TBB
#   -w279  controlling expression is constant
#   -w327  NULL reference is not allowed
#          (the compiler is correct here in that statements like
#           *static_cast<int*>(0) are not allowed to initialize
#           references; however, it's the only useful way to do
#           so if you need an invalid value for a reference)
#   -w383  value copied to temporary, reference to temporary used
#   -w981  operands are evaluated in unspecified order
#   -w1418 external function definition with no prior declaration
#          (happens in boost)
#   -w1478 deprecation warning
#   -w1572 floating-point equality and inequality comparisons are unreliable
#   -w2259 non-pointer conversion from "double" to "float" may
#          lose significant bits
#   -w21   type qualifiers are meaningless in this declaration
#   -w2536 type qualifiers are meaningless here
#
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-wd68")
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-wd135")
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-wd175")
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-wd177")
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-wd191")
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-wd193")
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-wd279")
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-wd327")
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-wd383")
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-wd981")
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-wd1418")
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-wd1478")
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-wd1572")
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-wd2259")
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-wd21")
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-wd2536")


#
# Also disable the following warnings that we frequently
# trigger writing dimension independent code:
#   -w111 statement is unreachable
#         Happens in code that is guarded by a check on 'dim'
#   -w128 loop is not reachable from preceding
#         Same as above
#   -w185 dynamic initialization in unreachable code
#         When initializing a local variable in code
#         that is executed only for one specific dimension
#   -w280 selector expression is constant
#         When writing 'switch(dim)'
#
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-wd111")
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-wd128")
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-wd185")
ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS "-wd280")


IF(DEAL_II_STATIC_EXECUTABLE)
  #
  # To produce a static executable, we have to statically link intel's
  # support libraries:
  #
  ENABLE_IF_SUPPORTED(DEAL_II_LINKER_FLAGS "-static")
  ENABLE_IF_SUPPORTED(DEAL_II_LINKER_FLAGS "-static-intel")
  ENABLE_IF_SUPPORTED(DEAL_II_LINKER_FLAGS "-static-gcc")
  ENABLE_IF_SUPPORTED(DEAL_II_LINKER_FLAGS "-pthread")
ELSE()
  #
  # Explicitly link intel support libraries dynamically:
  #
  ENABLE_IF_SUPPORTED(DEAL_II_LINKER_FLAGS "-shared-intel")
ENDIF()


#############################
#                           #
#    For Release target:    #
#                           #
#############################

IF (CMAKE_BUILD_TYPE MATCHES "Release")
  #
  # General optimization flags:
  #

  ADD_FLAGS(DEAL_II_CXX_FLAGS_RELEASE "-O2")

  # equivalent to -fno-strict-aliasing:
  ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS_RELEASE "-no-ansi-alias")

  ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS_RELEASE "-ip")

  ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS_RELEASE "-funroll-loops")
ENDIF()


###########################
#                         #
#    For Debug target:    #
#                         #
###########################

IF (CMAKE_BUILD_TYPE MATCHES "Debug")
  LIST(APPEND DEAL_II_DEFINITIONS_DEBUG "DEBUG")
  LIST(APPEND DEAL_II_USER_DEFINITIONS_DEBUG "DEBUG")

  ADD_FLAGS(DEAL_II_CXX_FLAGS_DEBUG "-O0")

  ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS_DEBUG "-g")
  ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS_DEBUG "-gdwarf-2")
  ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS_DEBUG "-grecord-gcc-switches")
ENDIF()

