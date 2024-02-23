// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test the functions of ConditionalOStream


#include <deal.II/base/conditional_ostream.h>

#include <limits>

#include "../tests.h"


int
main()
{
  initlog();

  ConditionalOStream o(deallog.get_file_stream(), true);
  o << "Yes" << std::endl;
  deallog << o.is_active() << std::endl;

  o.set_condition(false);
  o << "No" << std::endl;
  deallog << o.is_active() << std::endl;
}
