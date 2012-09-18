#
# Setup default compiler flags:
#
# FAT NOTE:
#
# TODO
#

#
# TODO: For the moment assume that CC and CXX are the same compilers...
#




#######################################################
#                                                     #
#    General setup for GCC and GCC like compilers:    #
#                                                     #
#######################################################

IF( CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR
    CMAKE_CXX_COMPILER_ID MATCHES "Clang" )

  ########################
  #                      #
  #    General setup:    #
  #                      #
  ########################

  #
  # Set the pic flag.
  # On some systems, -fpic/PIC is implied, so don't set anything to avoid
  # a warning. (TODO)
  #
  ENABLE_IF_AVAILABLE(CMAKE_C_FLAGS "-fpic")
  ENABLE_IF_AVAILABLE(CMAKE_CXX_FLAGS "-fpic")

  #
  # Check whether the -as-needed flag is available. If so set it to link
  # the deal.II library with it.
  #
  ENABLE_IF_AVAILABLE(CMAKE_SHARED_LINKER_FLAGS "-Wl,-as-needed")

  #
  # Set -pedantic if the compiler supports it.
  # Do not enable -pedantic for gcc-4.4, though...
  #
  IF(NOT ( CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND CMAKE_CXX_COMPILER_VERSION MATCHES "4.4." ))
    ENABLE_IF_AVAILABLE(CMAKE_C_FLAGS "-pedantic")
    ENABLE_IF_AVAILABLE(CMAKE_CXX_FLAGS "-pedantic")
  ENDIF()

  ENABLE_IF_AVAILABLE(CMAKE_C_FLAGS "-Wall")
  ENABLE_IF_AVAILABLE(CMAKE_CXX_FLAGS "-Wall")
  ENABLE_IF_AVAILABLE(CMAKE_C_FLAGS "-Wpointer-arith")
  ENABLE_IF_AVAILABLE(CMAKE_CXX_FLAGS "-Wpointer-arith")
  ENABLE_IF_AVAILABLE(CMAKE_C_FLAGS "-Wwrite-strings")
  ENABLE_IF_AVAILABLE(CMAKE_CXX_FLAGS "-Wwrite-strings")
  ENABLE_IF_AVAILABLE(CMAKE_C_FLAGS "-Wsynth")
  ENABLE_IF_AVAILABLE(CMAKE_CXX_FLAGS "-Wsynth")
  ENABLE_IF_AVAILABLE(CMAKE_C_FLAGS "-Wsign-compare")
  ENABLE_IF_AVAILABLE(CMAKE_CXX_FLAGS "-Wsign-compare")
  ENABLE_IF_AVAILABLE(CMAKE_C_FLAGS "-Wswitch")
  ENABLE_IF_AVAILABLE(CMAKE_CXX_FLAGS "-Wswitch")

  #################################
  #                               #
  #    For the Release target:    #
  #                               #
  #################################

  #
  # General optimization flags:
  #
  ADD_FLAGS(CMAKE_C_FLAGS_RELEASE "-O2")
  ADD_FLAGS(CMAKE_CXX_FLAGS_RELEASE "-O2")

  ENABLE_IF_AVAILABLE(CMAKE_C_FLAGS_RELEASE "-funroll-loops")
  ENABLE_IF_AVAILABLE(CMAKE_CXX_FLAGS_RELEASE "-funroll-loops")
  ENABLE_IF_AVAILABLE(CMAKE_C_FLAGS_RELEASE "-fstrict-aliasing")
  ENABLE_IF_AVAILABLE(CMAKE_CXX_FLAGS_RELEASE "-fstrict-aliasing")
  ENABLE_IF_AVAILABLE(CMAKE_C_FLAGS_RELEASE "-felide-constructors")
  ENABLE_IF_AVAILABLE(CMAKE_CXX_FLAGS_RELEASE "-felide-constructors")

  ###############################
  #                             #
  #    For the Debug target:    #
  #                             #
  ###############################

  IF (CMAKE_BUILD_TYPE MATCHES "Debug")
    ADD_DEFINITIONS("-DDEBUG")
  ENDIF()

  ADD_FLAGS(CMAKE_C_FLAGS_DEBUG "-O0")
  ADD_FLAGS(CMAKE_CXX_FLAGS_DEBUG "-O0")

  ENABLE_IF_AVAILABLE(CMAKE_C_FLAGS_DEBUG "-ggdb")
  ENABLE_IF_AVAILABLE(CMAKE_CXX_FLAGS_DEBUG "-ggdb")
  ENABLE_IF_AVAILABLE(CMAKE_SHARED_LINKER_FLAGS "-ggdb")
  IF(NOT DEAL_II_HAVE_FLAG_-ggdb)
    # If -ggdb is not available, fall back to -g:
    ENABLE_IF_AVAILABLE(CMAKE_C_FLAGS_DEBUG "-g")
    ENABLE_IF_AVAILABLE(CMAKE_CXX_FLAGS_DEBUG "-g")
    ENABLE_IF_AVAILABLE(CMAKE_SHARED_LINKER_FLAGS "-g")
  ENDIF()

ELSE()

  MESSAGE(WARNING "Unrecognized compiler. Please set the relevant compiler options by hand.")
ENDIF()

