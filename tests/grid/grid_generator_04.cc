// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2007 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test GridGenerator::hyper_cube_with_hole

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"



template <int dim, int spacdim = dim>
void
test(std::ostream &out)
{
  GridOut go;
  go.set_flags(GridOutFlags::Ucd(false, true));
  Triangulation<dim, spacdim> tr;

  std::vector<double> radii;
  radii.push_back(.2);
  radii.push_back(.3);
  radii.push_back(.4);


  std::vector<double> radiiext;
  radiiext.push_back(.3);
  radiiext.push_back(.4);
  radiiext.push_back(.5);

  std::vector<double> Ls;
  Ls.push_back(.1);
  Ls.push_back(1);
  Ls.push_back(10);

  std::vector<unsigned int> Nzs;
  Nzs.push_back(2);
  Nzs.push_back(3);
  Nzs.push_back(4);

  for (unsigned int i = 0; i < radii.size(); ++i)
    for (unsigned int k = 0; k < (dim == 2 ? 1 : Ls.size()); ++k)
      for (unsigned int l = 0; l < (dim == 2 ? 1 : Ls.size()); ++l)
        {
          out << "               ====================" << std::endl
              << "Inner radius = " << radii[i] << std::endl
              << "Outer radius = " << radiiext[i] << std::endl
              << "Depth        = " << Ls[k] << std::endl
              << "Nzs          = " << Nzs[l] << std::endl
              << "No colorize    ====================" << std::endl;

          // No colorize
          try
            {
              GridGenerator::hyper_cube_with_cylindrical_hole(
                tr, radii[i], radiiext[i], Ls[k], Nzs[l], false);
            }
          catch (...)
            {
              out << "Generation failed." << std::endl;
            }
          if (tr.n_cells() > 0)
            go.write_ucd(tr, out);

          tr.clear();

          out << "Colorize       ====================" << std::endl;
          try
            {
              GridGenerator::hyper_cube_with_cylindrical_hole(
                tr, radii[i], radiiext[i], Ls[k], Nzs[l], true);
            }
          catch (...)
            {
              out << "Generation failed." << std::endl;
            }
          if (tr.n_cells() > 0)
            go.write_ucd(tr, out);
          tr.clear();

          out << "               ====================" << std::endl;
        }
}


int
main()
{
  initlog();

  deallog.push("2d");
  test<2>(deallog.get_file_stream());
  deallog.pop();
  deallog.push("2d-3d");
  test<2, 3>(deallog.get_file_stream());
  deallog.pop();
  deallog.push("3d");
  test<3>(deallog.get_file_stream());
  deallog.pop();
}
