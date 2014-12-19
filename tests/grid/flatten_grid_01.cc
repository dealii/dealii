// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2013 by the deal.II authors
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

// Generate a grid, refine it once, flatten it and output the result.

#include "../tests.h"

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_generator.h>

template <int dim, int spacedim1, int spacedim2> 
void test() {
  deallog << "Testing <" << dim << "," << spacedim1
	  << "> VS <" << dim << "," << spacedim2 
	  << ">" << std::endl;
  
  Triangulation<dim, spacedim1> tria1;
  GridGenerator::hyper_cube(tria1);
  tria1.refine_global(1);
  
  Triangulation<dim, spacedim2> tria2;
  GridGenerator::flatten_triangulation(tria1, tria2);
  GridOut go;
  go.write_msh(tria2, deallog.get_file_stream());
}

int main()
{
  initlog();
  test<1,1,1>();
  test<1,1,2>();
  test<1,1,3>();
  // 
  test<1,2,1>();
  test<1,2,2>();
  test<1,2,3>();
  // 
  test<1,3,1>();
  test<1,3,2>();
  test<1,3,3>();
  //
  test<2,2,2>();
  test<2,2,3>();
  //
  test<2,3,2>();
  test<2,3,3>();
  //
  test<3,3,3>();
}
