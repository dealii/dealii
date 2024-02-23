// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check that we can read in and edge orient meshes in 3d with a
// significant number of cells. this tests Mike's reorienting
// algorithm. while the grid_in_3d test checks artificial meshes, this
// one uses a mesh of the 747 nacelle

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"

// #include <deal.II/fe/fe_q.h>
// #include <deal.II/numerics/data_out.h>

void
test(const char *filename)
{
  const unsigned int dim = 3;

  Triangulation<3> tria;
  GridIn<3>        gi;
  gi.attach_triangulation(tria);
  std::ifstream in(filename);

  try
    {
      gi.read_ucd(in);
    }
  catch (const std::exception &exc)
    {
      deallog << "  caught exception:" << std::endl << exc.what() << std::endl;
      return;
    }

  deallog << "  " << tria.n_active_cells() << " active cells" << std::endl;

  int hash  = 0;
  int index = 0;
  for (Triangulation<dim>::active_cell_iterator c = tria.begin_active();
       c != tria.end();
       ++c, ++index)
    for (const unsigned int i : c->vertex_indices())
      hash += (index * i * c->vertex_index(i)) % (tria.n_active_cells() + 1);
  deallog << "  hash=" << hash << std::endl;

  //   DoFHandler<dim> dh (tria);
  //   FE_Q<dim> fe(1);
  //   dh.distribute_dofs (fe);

  //   DataOut<dim> d_o;
  //   d_o.attach_dof_handler (dh);
  //   Vector<double> zero (dh.n_dofs());
  //   d_o.add_data_vector (zero, "Zero");
  //   d_o.build_patches ();
  //   std::ofstream gnufile("xx.vtk");
  //   d_o.write_vtk (gnufile);
}


int
main()
{
  initlog();
  deallog << std::setprecision(2);

  test(SOURCE_DIR "/grid_in_3d_02/747.ucd");
}
