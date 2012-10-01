#####
##
## Copyright (C) 2012 by the deal.II authors
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

IF(DEAL_II_USE_MT_POSIX)
  SET(src_boost_thread
    pthread/once.cpp
    pthread/thread.cpp
    )
ELSE()
  SET(src_boost_thread
    win32/once.cpp
    win32/thread.cpp
    )
ENDIF()

DEAL_II_ADD_LIBRARY(obj_boost_thread OBJECT ${src_boost_thread})

IF(DEAL_II_USE_MT_POSIX)

  # TODO: Necessary?
    # include <boost/cstdint.hpp>

    # ifndef UINTMAX_C
    #   define UINTMAX_C(x) x ## ULL
    # endif

  DEAL_II_ADD_DEFINITIONS(obj_boost_thread
    "BOOST_THREAD_BUILD_LIB" "BOOST_ALL_NO_LIB"
    )
ELSE()
  DEAL_II_ADD_DEFINITIONS(obj_boost_thread
    "BOOST_THREAD_BUILD_LIB" "BOOST_ALL_NO_LIB" "BOOST_THREAD_POSIX"
    )
ENDIF()

