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



// Test GridGenerator::subdivided_hyper_rectangle with vector of step
// sizes.

#include <deal.II/base/tensor.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"



template <int dim>
void
test(std::ostream &out)
{
  Point<dim> p1;
  p1[0] = 2.;
  if (dim > 1)
    p1[1] = -1.;
  if (dim > 2)
    p1[2] = 0.;
  Point<dim> p2;
  p2[0] = 3.;
  if (dim > 1)
    p2[1] = 2.;
  if (dim > 2)
    p2[2] = 4.;
  Point<dim> p3;
  p3[0] = 2.;
  if (dim > 1)
    p3[1] = 1.;
  if (dim > 2)
    p3[2] = 4.;

  GridOut go;

  // uniformly subdivided mesh
  if (true)
    {
      deallog << "subdivided_hyper_rectangle" << std::endl;
      Triangulation<dim>               tr;
      std::vector<std::vector<double>> sub(dim);
      for (unsigned int i = 0; i < dim; ++i)
        sub[i] = std::vector<double>(i + 2, (p2[i] - p1[i]) / (i + 2));

      GridGenerator::subdivided_hyper_rectangle(tr, sub, p1, p2, true);
      if (tr.n_cells() > 0)
        go.write_gnuplot(tr, out);
    }


  // non-uniformly subdivided mesh
  if (true)
    {
      deallog << "subdivided_hyper_rectangle" << std::endl;
      Triangulation<dim>               tr;
      std::vector<std::vector<double>> sub(dim);
      for (unsigned int i = 0; i < dim; ++i)
        {
          sub[i] = std::vector<double>(i + 2, (p2[i] - p1[i]) / (i + 2));
          sub[i][0] /= 2;
          sub[i].back() *= 1.5;
        }

      GridGenerator::subdivided_hyper_rectangle(tr, sub, p1, p2, true);
      if (tr.n_cells() > 0)
        go.write_gnuplot(tr, out);
    }
}


int
main()
{
  initlog();

  deallog.push("1d");
  test<1>(deallog.get_file_stream());
  deallog.pop();
  deallog.push("2d");
  test<2>(deallog.get_file_stream());
  deallog.pop();
  deallog.push("3d");
  test<3>(deallog.get_file_stream());
  deallog.pop();
}
