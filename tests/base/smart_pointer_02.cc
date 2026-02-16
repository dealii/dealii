// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// Ensure that the old SmartPointer class name continues to be
// available after the renaming to ObserverPointer.


#include <deal.II/base/enable_observer_pointer.h>
#include <deal.II/base/observer_pointer.h>
#include <deal.II/base/smartpointer.h>

#include <any>
#include <iostream>

#include "../tests.h"


class Test : public EnableObserverPointer
{};


int
main()
{
  static_assert(std::is_same_v<ObserverPointer<Test>, SmartPointer<Test>>);
  static_assert(
    std::is_same_v<ObserverPointer<Test, int>, SmartPointer<Test, int>>);
}
