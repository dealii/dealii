// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Check that it is possible to put SmartPointer objects into a
// std::any object.


#include <deal.II/base/smartpointer.h>
#include <deal.II/base/subscriptor.h>

#include <any>
#include <iostream>

#include "../tests.h"


class Test : public Subscriptor
{};


int
main()
{
  initlog();

  Test               t;
  SmartPointer<Test> r(&t);
  std::any           a = r;
}
