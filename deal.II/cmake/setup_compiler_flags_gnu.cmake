#
# General setup for GCC and compilers sufficiently close to GCC
#
# Please read the fat note in setup_compiler_flags.cmake prior to
# editing this file.
#


########################
#                      #
#    General setup:    #
#                      #
########################

#
# Set the pic flag.
# TODO: On some systems, -fpic/PIC is implied, so don't set anything to
# avoid a warning.
#
ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-fpic")

#
# Check whether the -as-needed flag is available. If so set it to link
# the deal.II library with it.
#
ENABLE_IF_LINKS(CMAKE_SHARED_LINKER_FLAGS "-Wl,--as-needed")


#
# Set -pedantic if the compiler supports it.
#
ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-pedantic")

#
# Setup various warnings:
#
ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-Wall")
ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-Wpointer-arith")
ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-Wwrite-strings")
ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-Wsynth")
ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-Wsign-compare")
ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-Wswitch")


#
# Newer versions of gcc have a flag -Wunused-local-typedefs that, though in
# principle a good idea, triggers a lot in BOOST in various places.
# Unfortunately, this warning is included in -W/-Wall, so disable it if the
# compiler supports it.
#
ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-Wno-unused-local-typedefs")


IF(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  #
  # Like many other compilers, clang produces warnings for array
  # accesses out of bounds, even if they are in code that's dead
  # for this dimension. Suppress this.
  #
  # There are a number of other warnings we get that can't easily
  # be worked around and that are definitely not useful. Suppress
  # those too.
  #
  ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-Wno-array-bounds")
  ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-Wno-parentheses")
  ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-Wno-delete-non-virtual-dtor")
  ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-Wno-unneeded-internal-declaration")
  ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-Wno-unused-function")
  ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-Wno-unused-variable")
  ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-Wno-long-long")
ENDIF()


IF(CMAKE_SYSTEM_NAME MATCHES "Darwin")
  #
  # Use -Wno-long-long on Apple Darwin to avoid some unnecessary
  # warnings. However, newer gccs on that platform do not have
  # this flag any more, so check whether we can indeed do this
  #
  ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-Wno-long-double")
ENDIF()


IF(CMAKE_SYSTEM_NAME MATCHES "CYGWIN") # TODO: Check for correct name
  #
  # TODO: Is this flag still necessary? Sound pretty invasive...
  #
  # On Cygwin, when using shared libraries, there might occur
  # difficulties when linking libraries for several dimensions,
  # as some symbols are defined in all of them. This leads to a
  # linker error. We force the linker to ignore multiple symbols,
  # but of course this might lead to strange program behaviour if
  # you accidentally defined one symbol multiple times...
  # (added 2005/07/13, Ralf B. Schulz)
  # (modified 2005/12/20, Ralf B. Schulz)
  ENABLE_IF_LINKS(CMAKE_SHARED_LINKER_FLAGS
    "-Xlinker --allow-multiple-definition"
    )
ENDIF()


##############################
#                            #
#    For Release targets:    #
#                            #
##############################

IF (CMAKE_BUILD_TYPE MATCHES "Release")
  #
  # General optimization flags:
  #
  ADD_FLAGS(DEAL_II_CXX_FLAGS_RELEASE "-O2")

  ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS_RELEASE "-funroll-loops")
  ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS_RELEASE "-fstrict-aliasing")
  ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS_RELEASE "-felide-constructors")

  ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS_RELEASE "-Wno-unused")
ENDIF()


###########################
#                         #
#    For Debug target:    #
#                         #
###########################

IF (CMAKE_BUILD_TYPE MATCHES "Debug")

  LIST(APPEND DEAL_II_DEFINITIONS_DEBUG "-DDEBUG")
  LIST(APPEND DEAL_II_USER_DEFINTIONS_DEBUG "-DDEBUG")

  ADD_FLAGS(DEAL_II_CXX_FLAGS_DEBUG "-O0")

  ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS_DEBUG "-ggdb")
  ENABLE_IF_SUPPORTED(DEAL_II_SHARED_LINKER_FLAGS_DEBUG "-ggdb")
  #
  # If -ggdb is not available, fall back to -g:
  #
  IF(NOT DEAL_II_HAVE_FLAG_-ggdb)
    ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS_DEBUG "-g")
    ENABLE_IF_SUPPORTED(DEAL_II_SHARED_LINKER_FLAGS_DEBUG "-g")
  ENDIF()
ENDIF()


#########################
#                       #
#    Set up C FLAGS:    #
#                       #
#########################

#
# For the moment we assume that CC and CXX are the same compiler and that
# we can set (almost) the same default flags for both:
#
SET(CMAKE_C_FLAGS ${CMAKE_CXX_FLAGS})
SET(DEAL_II_C_FLAGS_RELEASE ${DEAL_II_CXX_FLAGS_RELEASE})
SET(DEAL_II_C_FLAGS_DEBUG ${DEAL_II_CXX_FLAGS_DEBUG})

#
# OK, touché, touché. We have to strip flags not supported by a C target:
#
STRIP_FLAG(CMAKE_C_FLAGS "-Wsynth")
STRIP_FLAG(DEAL_II_C_FLAGS_RELEASE "-felide-constructors")

#
# and disable some warnings:
#
STRIP_FLAG(CMAKE_C_FLAGS "-Wall") # There is no other way to disable -Wunknown-pragma atm...
STRIP_FLAG(CMAKE_C_FLAGS "-Wsign-compare")
STRIP_FLAG(CMAKE_C_FLAGS "-Wwrite-strings")

