// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2010 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// common include file for all serialization tests

#include <deal.II/base/logstream.h>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
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
  // First test serialization to a text archive. This one we can actually
  // output into the log file and use as a reference for future tests.
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

  // Now do the same thing again, but with a binary archive. This one we
  // can't output into the log, but we can at least verify that it works
  // correctly. (In practice, one might then also want to compress
  // the binary archive, but that is immaterial to us here. The point
  // is just to verify that serialization to a binary archive works
  // correctly.)
  {
    std::ostringstream oss;
    {
      boost::archive::binary_oarchive oa(oss, boost::archive::no_header);
      oa << t1;
    }
    {
      std::istringstream              iss(oss.str());
      boost::archive::binary_iarchive ia(iss, boost::archive::no_header);

      ia >> t2;

      AssertThrow(compare(t1, t2), ExcInternalError());
    }
  }
}
