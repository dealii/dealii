// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check MGLevelObject::apply()

#include <deal.II/base/mg_level_object.h>

#include <algorithm>

#include "../tests.h"



template <class T>
void
check(MGLevelObject<T> &obj)
{
  obj.apply([&](const unsigned int lvl, T &value) { value = (T)lvl; });

  obj.apply([&](const unsigned int lvl, T &value) {
    deallog << "lvl: " << lvl << " value: " << value << std::endl;
  });
}

int
main()
{
  initlog();

  MGLevelObject<double> o(2, 4);
  check(o);
  o.resize(0, 1);
  check(o);
}
