// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2004 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check GridTools::partition_triangulation. generate output in VTK format

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>

#include "../tests.h"


template <int dim>
void
test()
{
  deallog << "Dimension = " << dim << std::endl;

  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(4 - dim);
  for (unsigned int i = 0; i < 11 - 2 * dim; ++i)
    {
      typename Triangulation<dim>::active_cell_iterator cell =
        triangulation.begin_active();
      for (unsigned int index = 0; cell != triangulation.end(); ++cell, ++index)
        if (index % (3 * dim) == 0)
          cell->set_refine_flag();
      triangulation.execute_coarsening_and_refinement();
    }


  // subdivide into 5 subdomains
  GridTools::partition_triangulation(5, triangulation);

  // generate a field where the value equals
  // the subdomain number, and output it
  FE_DGQ<dim>     fe(0);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  Vector<double> partitions(triangulation.n_active_cells());
  {
    typename Triangulation<dim>::active_cell_iterator cell =
      triangulation.begin_active();
    for (unsigned int index = 0; cell != triangulation.end(); ++cell, ++index)
      partitions(index) = cell->subdomain_id();
  }

  DataOut<dim> data_out;
  data_out.attach_triangulation(triangulation);
  data_out.add_data_vector(partitions, "partitions");
  data_out.build_patches();

  data_out.write_vtk(deallog.get_file_stream());
}



int
main()
{
  initlog();

  try
    {
      test<1>();
      test<2>();
      test<3>();
    }
  catch (const std::exception &exc)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;

      return 1;
    }
  catch (...)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    };
}
