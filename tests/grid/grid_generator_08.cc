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
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// Test grid generation functions  in GridGenerator for spacedim>dim

#include "../tests.h"
#include <deal.II/base/tensor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>



template <int dim, int spacedim>
void
test(std::ostream &out)
{
  Point<dim> p1;
  p1[0] = 2.;
  if (dim>1) p1[1] = -1.;
  Point<dim> p2;
  p2[0] = 3.;
  if (dim>1) p2[1] = 2.;
  Point<dim> p3;
  p3[0] = 2.;
  if (dim>1) p3[1] = 1.;

  GridOut go;
  GridOut::OutputFormat format = GridOut::msh;

  if (true)
    {
      deallog << "hyper_cube" << std::endl;
      Triangulation<dim,spacedim> tr;
      GridGenerator::hyper_cube(tr, 3., 7.);
      if (tr.n_cells() > 0)
        go.write(tr, out, format);
    }
  if (true)
    {
      deallog << "subdivided_hyper_cube" << std::endl;
      Triangulation<dim,spacedim> tr;
      GridGenerator::subdivided_hyper_cube(tr, 3, 1., 7.);
      if (tr.n_cells() > 0)
        go.write(tr, out, format);
    }
  if (true)
    {
      deallog << "hyper_rectangle" << std::endl;
      Triangulation<dim,spacedim> tr;
      GridGenerator::hyper_rectangle(tr, p1, p2, true);
      if (tr.n_cells() > 0)
        go.write(tr, out, format);
    }
  if (true)
    {
      deallog << "subdivided_hyper_rectangle" << std::endl;
      Triangulation<dim,spacedim> tr;
      std::vector<unsigned int> sub(dim);
      sub[0] = 2;
      if (dim>1) sub[1] = 3;
      if (dim>2) sub[2] = 4;
      GridGenerator::subdivided_hyper_rectangle(tr, sub, p1, p2, true);
      if (tr.n_cells() > 0)
        go.write(tr, out, format);
    }
}

int
main()
{
  initlog();

  deallog.push("1-2");
  test<1,2>(deallog.get_file_stream());
  deallog.pop();
  deallog.push("1-3");
  test<1,3>(deallog.get_file_stream());
  deallog.pop();
  deallog.push("2-3");
  test<2,3>(deallog.get_file_stream());
  deallog.pop();
}
