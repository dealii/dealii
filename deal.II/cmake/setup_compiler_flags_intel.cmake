#####
##
## Copyright (C) 2012, 2013 by the deal.II authors
##
## This file is part of the deal.II library.
##
## <TODO: Full License information>
## This file is dual licensed under QPL 1.0 and LGPL 2.1 or any later
## version of the LGPL license.
##
## Author: Matthias Maier <matthias.maier@iwr.uni-heidelberg.de>
##
#####

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
ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-fpic")

#
# Check whether the -as-needed flag is available. If so set it to link
# the deal.II library with it.
#
ENABLE_IF_LINKS(DEAL_II_LINKER_FLAGS "-Wl,--as-needed")

#
# Set ansi mode:
#
ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-ansi")

#
# Enable verbose warnings:
#
ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-w2")

#
# Disable some warnings that lead to a lot of false positives:
#
#   -w68   integer conversion resulted in a change of sign
#          (triggers a lot in functionparser)
#   -w175  subscript out of range
#   -w177  declared but not referenced
#   -w279  controlling expression is constant
#   -w383  value copied to temporary, reference to temporary used
#   -w981  operands are evaluated in unspecified order
#   -w1418 external function definition with no prior declaration
#          (happens in boost)
#   -w1478 deprecation warning
#   -w1572 floating-point equality and inequality comparisons are unreliable
#
ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-wd68")
ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-wd175")
ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-wd177")
ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-wd279")
ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-wd383")
ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-wd981")
ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-wd1418")
ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-wd1478")
ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-wd1572")


IF(DEAL_II_STATIC_EXECUTABLE)
  #
  # To produce a static executable, we have to statically link intel's
  # support libraries:
  #
  ENABLE_IF_SUPPORTED(DEAL_II_LINKER_FLAGS "-static")
  ENABLE_IF_SUPPORTED(DEAL_II_LINKER_FLAGS "-static-intel")
  ENABLE_IF_SUPPORTED(DEAL_II_LINKER_FLAGS "-static-gcc")
  ENABLE_IF_SUPPORTED(DEAL_II_LINKER_FLAGS "-pthread")
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
  ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS_RELEASE "-ip")

  ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS_RELEASE "-funroll-loops")
  ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS_RELEASE "-ansi-alias") # equiv. -fstrict-aliasing
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

