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



// test that MappingCollection objects are copyable without running into
// troubles when the copy is destroyed earlier than the original
// object


#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/hp/mapping_collection.h>

#include "../tests.h"



template <int dim>
void
test()
{
  hp::MappingCollection<dim> mapping_collection;
  mapping_collection.push_back(MappingQ<dim>(2));
  mapping_collection.push_back(MappingQ<dim>(1));

  // now create a copy and make sure
  // it goes out of scope before the
  // original
  {
    hp::MappingCollection<dim> copy(mapping_collection);
  }
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
