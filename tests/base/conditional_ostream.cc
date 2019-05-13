// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2018 by the deal.II authors
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
