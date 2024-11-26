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



// Check that it is possible to put ObserverPointer objects into a
// std::any object.


#include <deal.II/base/enable_observer_pointer.h>
#include <deal.II/base/observer_pointer.h>

#include <any>
#include <iostream>

#include "../tests.h"


class Test : public EnableObserverPointer
{};


int
main()
{
  initlog();

  Test                  t;
  ObserverPointer<Test> r(&t);
  std::any              a = r;
}
