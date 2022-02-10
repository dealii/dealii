// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2022 by the deal.II authors
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


// test the testsuite framework. this test is supposed to run successfully
// and the diff stage must be skipped

#include "../tests.h"


int
main()
{
  initlog();

  deallog << "OUTPUT!" << std::endl;
}
