## ---------------------------------------------------------------------
## $Id$
##
## Copyright (C) 2012 - 2019 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE.md at
## the top level directory of deal.II.
##
## ---------------------------------------------------------------------

if(DEAL_II_USE_MT_POSIX)
  set(src_boost_thread
    pthread/once.cpp
    pthread/thread.cpp
    )
else()
  set(src_boost_thread
    win32/thread.cpp
    win32/tss_dll.cpp
    win32/tss_pe.cpp
    )
endif()

deal_ii_add_library(obj_boost_thread OBJECT ${src_boost_thread})
deal_ii_add_definitions(obj_boost_thread
  "BOOST_THREAD_BUILD_LIB" "BOOST_ALL_NO_LIB"
  )

if(CMAKE_SYSTEM_NAME MATCHES "Windows")
  deal_ii_add_definitions(obj_boost_thread "BOOST_HAS_WINTHREADS")
else()
  deal_ii_add_definitions(obj_boost_thread "BOOST_THREAD_POSIX")
endif()

