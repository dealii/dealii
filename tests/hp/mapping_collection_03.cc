// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// a test that triggers really hard to track down failures in
// mapping_collection_01 in a really simple way; an even further
// simplified form of mapping_collection_02


#include <deal.II/fe/mapping_q.h>

#include <deal.II/hp/mapping_collection.h>

#include "../tests.h"



template <int dim>
void
test()
{
  MappingQ<dim> mapping(2);
  deallog << "Cloning..." << std::endl;
  std::unique_ptr<Mapping<dim>> copy = mapping.clone();
  deallog << "Deleting clone..." << std::endl;
  copy.reset();
  deallog << "Destroying original..." << std::endl;
}



int
main()
{
  std::ofstream logfile("output");
  logfile.precision(2);

  deallog.attach(logfile);

  test<1>();
  test<2>();
  test<3>();

  deallog << "OK" << std::endl;
}
