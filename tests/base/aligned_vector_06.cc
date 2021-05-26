// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


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
