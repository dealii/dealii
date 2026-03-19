// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// Check serialization for AlignedVector. There was a bug for the case
// where the loading object had nonzero size, and the loaded object
// has zero size. We forgot to call resize() in that case.

#include <deal.II/base/aligned_vector.h>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <sstream>
#include <string>

#include "../tests.h"


template <typename T>
void
test_serialization()
{
  // Create an empty vector and serialize it.
  AlignedVector<T> vec;

  std::ostringstream              out;
  boost::archive::binary_oarchive oarchive(out);
  oarchive << vec;
  const auto serialization = out.str();

  // Then resize it to something nonzero
  vec.resize(10);

  // Deserialize into the same vector. It should have size zero again
  // (but that was broken at the time of writing the test).
  std::istringstream              in(serialization);
  boost::archive::binary_iarchive iarchive(in);
  iarchive >> vec;

  AssertThrow(vec.size() == 0, ExcInternalError());
  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();

  // Test with trivial and non-trivial types:
  test_serialization<int>();
  test_serialization<std::string>();
}
