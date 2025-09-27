// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Make sure that we don't get into trouble for a vector that is empty
// -- the AlignedVector class uses a custom deleter that aborts if the
// object is not storing any data, but this deleter should not be
// called in that case.

#include <deal.II/base/aligned_vector.h>

#include "../tests.h"


void
test()
{
  {
    AlignedVector<int> i;

    // Do nothing -- the object just goes out of scope again.
  }

  {
    AlignedVector<int> i;
    i.resize(0); // Do not allocate any elements here either

    // Do nothing more -- the object just goes out of scope again.
  }

  {
    AlignedVector<int> i(0);

    // Do nothing more -- the object just goes out of scope again.
  }

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  test();
}
