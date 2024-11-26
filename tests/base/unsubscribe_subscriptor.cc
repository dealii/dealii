// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check that unsubscribing with a wrong id is handled correctly


#include <deal.II/base/enable_observer_pointer.h>
#include <deal.II/base/observer_pointer.h>

#include <iostream>
#include <vector>

#include "../tests.h"

class Test : public EnableObserverPointer
{};

int
main()
{
  deal_II_exceptions::disable_abort_on_exception();

  initlog();

  EnableObserverPointer subscriptor;
  std::atomic<bool>     dummy_a;
  const char           *foo_a = "a";
  const char           *foo_b = "b";
  subscriptor.subscribe(&dummy_a, foo_a);
  subscriptor.unsubscribe(&dummy_a, foo_b);
  std::atomic<bool> dummy_b;
  subscriptor.unsubscribe(&dummy_b, foo_a);
  subscriptor.unsubscribe(&dummy_a, foo_a);

  return 0;
}
