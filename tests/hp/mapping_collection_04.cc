// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2017 by the deal.II authors
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
  std::ofstream logfile("output");
  logfile.precision(2);

  deallog.attach(logfile);

  test<1>();
  test<2>();
  test<3>();

  deallog << "OK" << std::endl;
}
