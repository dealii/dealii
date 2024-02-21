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



// check GridTools::partition_triangulation with cell weighting.
// generate output in VTK format.
// This test is based off of metis_03.cc.

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
test(const bool with_weighting, const bool write_to_vtk = false)
{
  deallog << "Dimension = " << dim << std::endl;
  deallog << "With weighting = " << std::boolalpha << with_weighting
          << std::endl;

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

  // cell weighting
  std::vector<unsigned int> cell_weighting;
  if (with_weighting)
    {
      cell_weighting =
        std::vector<unsigned int>(triangulation.n_active_cells(), 1u);

      const Point<dim>                                  origin;
      typename Triangulation<dim>::active_cell_iterator cell =
        triangulation.begin_active();
      for (unsigned int index = 0; cell != triangulation.end(); ++cell, ++index)
        {
          if (origin.distance(cell->center()) < 0.5)
            cell_weighting[index] = 10u;
        }
    }


  // subdivide into 5 subdomains
  GridTools::partition_triangulation(5, cell_weighting, triangulation);

  // generate a field where the value equals
  // the subdomain number, and output it
  FE_DGQ<dim>     fe(0);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  Vector<double> partitions(triangulation.n_active_cells());
  Vector<double> weights(triangulation.n_active_cells());
  {
    typename Triangulation<dim>::active_cell_iterator cell =
      triangulation.begin_active();
    for (unsigned int index = 0; cell != triangulation.end(); ++cell, ++index)
      {
        partitions(index) = cell->subdomain_id();
        weights(index)    = (with_weighting ? cell_weighting[index] : 1.0);
      }
  }

  DataOut<dim> data_out;
  data_out.attach_triangulation(triangulation);
  data_out.add_data_vector(partitions, "partitions");
  data_out.add_data_vector(weights, "weights");
  data_out.build_patches();

  data_out.write_vtk(deallog.get_file_stream());

  if (write_to_vtk)
    {
      std::stringstream filename;
      filename << "grid_" << dim << "d.with_weighting_" << with_weighting
               << ".vtk";
      std::ofstream outfile(filename.str());
      data_out.write_vtk(outfile);
    }
}



int
main()
{
  initlog();

  try
    {
      // Without weighting
      test<1>(false);
      test<2>(false);
      test<3>(false);

      // With weighting
      test<1>(true);
      test<2>(true);
      test<3>(true);
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
