// ---------------------------------------------------------------------
//
// Copyright (C) 2015-2016 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// test that a default-constructed ArrayView can actually be copied

#include "../tests.h"
#include <iomanip>

#include <deal.II/base/array_view.h>


void test ()
{
  ArrayView<int> a;
  ArrayView<int> b = a;

  deallog << "OK" << std::endl;
}




int main()
{
  initlog();

  test ();
}
