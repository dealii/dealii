// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// common include file for all serialization tests

#include <deal.II/base/logstream.h>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <fstream>
#include <iomanip>
#include <sstream>

#include "../tests.h"


// compare objects for equality and pointers for equality of the object
// pointed to
template <typename T>
bool
compare(const T &t1, const T &t2)
{
  return t1 == t2;
}

template <typename T>
bool
compare(T *t1, T *t2)
{
  return *t1 == *t2;
}


template <typename T>
void
verify(const T &t1, T &t2)
{
  // save data to archive
  std::ostringstream oss;
  {
    boost::archive::text_oarchive oa(oss, boost::archive::no_header);
    oa << t1;
    // archive and stream closed when
    // destructors are called
  }
  deallog << oss.str() << std::endl;

  // verify correctness of the
  // serialization
  {
    std::istringstream            iss(oss.str());
    boost::archive::text_iarchive ia(iss, boost::archive::no_header);


    ia >> t2;

    AssertThrow(compare(t1, t2), ExcInternalError());
  }
}
