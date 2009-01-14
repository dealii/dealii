//---------------------------------------------------------------------------
//    $Id: thread_management.cc 18205 2009-01-13 13:53:08Z bangerth $
//    Version: $Name$
//
//    Copyright (C) 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


// include all the files that form BOOST's thread implementation so that we
// don't have to build BOOST itself only to get at this small part of it. it
// also ensures that we use the correct compiler and flags

#define BOOST_THREAD_POSIX
#define BOOST_THREAD_BUILD_LIB 1

#include <../libs/thread/src/pthread/once.cpp>
#include <../libs/thread/src/pthread/exceptions.cpp>
#include <../libs/thread/src/pthread/thread.cpp>
