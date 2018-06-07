// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2017 by the deal.II authors
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


// ensure that we can refer to symbols using the old namespace
// std_cxx1x even if we only include the newer base/std_cxx11/*h
// include files (e.g., if we use the old namespace name but don't
// #include anything at all directly -- in which case we get things
// from indirect #includes made in the library which only include
// base/std_cxx11/*h)


#include <deal.II/base/std_cxx11/shared_ptr.h>

#include "../tests.h"


int
main()
{
  initlog();

  std_cxx1x::shared_ptr<int> p;
  deallog << "OK" << std::endl;
}
