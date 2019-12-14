// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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

// Test GridGenerator::generate_from_name_and_arguments.

#include <deal.II/base/tensor.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

template <int dim, int spacedim>
void
test(const std::string &name, const std::string &arguments)
{
  Triangulation<dim, spacedim> tria;

  deallog << "Generating Triangulation<" << dim << ", " << spacedim
          << "> : " << name << "(" << arguments << ")" << std::endl;

  GridGenerator::generate_from_name_and_arguments(tria, name, arguments);

  GridOut go;
  go.write_msh(tria, deallog.get_file_stream());

  std::ofstream ofile(name + "_" + std::to_string(dim) + "_" +
                      std::to_string(spacedim) + ".msh");
  go.write_msh(tria, ofile);
}


int
main()
{
  initlog();

  test<2, 2>("cheese", "3, 2");
  test<2, 3>("cheese", "3, 2");
  test<3, 3>("cheese", "4, 3, 2");

  test<1, 1>("general_cell", "0;     1    : true");
  test<1, 2>("general_cell", "0,0  ; 1,0  : true");
  test<1, 3>("general_cell", "0,0,0; 1,0,0: true");

  test<2, 2>("general_cell", "0,0; 1,0; 0,1; 1,1: true");
  test<2, 3>("general_cell", "0,0,0; 1,0,0; 0,1,0; 1,1,0: true");
  test<3, 3>("general_cell",
             "0,0,0; 1,0,0; 0,1,0; 1,1,0; "
             "0,0,1; 1,0,1; 0,1,1; 1,1,1: true");

  test<1, 1>("hyper_cube", "0: 1: true");
  test<1, 2>("hyper_cube", "0: 1: true");
  test<1, 3>("hyper_cube", "0: 1: true");
  test<2, 2>("hyper_cube", "0: 1: true");
  test<2, 3>("hyper_cube", "0: 1: true");
  test<3, 3>("hyper_cube", "0: 1: true");
}
