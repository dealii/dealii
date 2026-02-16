// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// Test constructor of hp::MappingCollection that accepts multiple mappings.


#include <deal.II/fe/mapping_q.h>

#include <deal.II/hp/mapping_collection.h>

#include "../tests.h"



template <int dim>
void
test()
{
  hp::MappingCollection<dim> mappings(MappingQ<dim>(1), MappingQ<dim>(2));

  deallog << mappings.size() << std::endl;
}



int
main()
{
  initlog();
  deallog.get_file_stream().precision(2);

  test<1>();
  test<2>();
  test<3>();

  deallog << "OK" << std::endl;
}
