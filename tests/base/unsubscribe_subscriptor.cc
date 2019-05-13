// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2019 by the deal.II authors
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



// check that unsubscribing with a wrong id is handled correctly


#include <deal.II/base/smartpointer.h>
#include <deal.II/base/subscriptor.h>

#include <iostream>
#include <vector>

#include "../tests.h"

class Test : public Subscriptor
{};

int
main()
{
  deal_II_exceptions::disable_abort_on_exception();

  initlog();

  Subscriptor       subscriptor;
  std::atomic<bool> dummy_a;
  const char *      foo_a = "a";
  const char *      foo_b = "b";
  subscriptor.subscribe(&dummy_a, foo_a);
  subscriptor.unsubscribe(&dummy_a, foo_b);
  std::atomic<bool> dummy_b;
  subscriptor.unsubscribe(&dummy_b, foo_a);
  subscriptor.unsubscribe(&dummy_a, foo_a);

  return 0;
}
