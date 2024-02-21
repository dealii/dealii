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



// a test that triggers really hard to track down failures in
// mapping_collection_01 in a really simple way


#include <deal.II/fe/mapping_q.h>

#include <deal.II/hp/mapping_collection.h>

#include "../tests.h"



template <int dim>
void
test()
{
  hp::MappingCollection<dim> mapping_collection(MappingQ<dim>(2));
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
