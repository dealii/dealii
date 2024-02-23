// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// a test that shows that mapping_collection_0[1-3] really is due to
// the fact that MappingQ has a dysfunctional copy constructor...


#include <deal.II/fe/mapping_q.h>

#include <deal.II/hp/mapping_collection.h>

#include "../tests.h"



template <int dim>
void
test()
{
  MappingQ<dim> mapping(2);
  {
    deallog << "Copying..." << std::endl;
    MappingQ<dim> copy(mapping);
    deallog << "Deleting clone..." << std::endl;
  }
  deallog << "Destroying original..." << std::endl;
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
