// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// extracted from count_dofs_per_block that crashed at one point


#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/intergrid_map.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/numerics/data_out.h>

#include <numeric>

#include "../tests.h"


template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> triangulation(
    MPI_COMM_WORLD, Triangulation<dim>::limit_level_difference_at_vertices);

  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(1);

  {
    typename Triangulation<dim>::cell_iterator it =
      triangulation.begin_active();
    it->set_refine_flag();
    triangulation.execute_coarsening_and_refinement();

    it = triangulation.begin(1);
    ++it;
    ++it;
    ++it;
    it = triangulation.begin(1);
    for (unsigned int a = 0; a < 4; ++a)
      it->child(a)->set_coarsen_flag();

    triangulation.prepare_coarsening_and_refinement();

    {
      // data out
      FE_Q<dim>       fe(1);
      DoFHandler<dim> dof_handler(triangulation);

      dof_handler.distribute_dofs(fe);

      unsigned int  n_coarse = 0;
      Vector<float> subdomain(triangulation.n_active_cells());
      {
        unsigned int index = 0;

        for (typename Triangulation<dim>::active_cell_iterator cell =
               triangulation.begin_active();
             cell != triangulation.end();
             ++cell, ++index)
          {
            subdomain(index) = 0;

            if (cell->is_ghost() || cell->is_artificial())
              subdomain(index) = -4;

            if (cell->refine_flag_set())
              subdomain(index) += 1;
            if (cell->coarsen_flag_set())
              {
                subdomain(index) += 2;
                ++n_coarse;
              }
          }
      }
      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
        {
          deallog << "id=" << triangulation.locally_owned_subdomain()
                  << " n_coarsen=" << n_coarse << std::endl;
          for (unsigned int i = 0; i < subdomain.size(); ++i)
            deallog << subdomain(i) << std::endl;
        }

      triangulation.execute_coarsening_and_refinement();
    }
  }
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  if (myid == 0)
    {
      initlog();

      deallog.push("2d");
      test<2>();
      deallog.pop();

      deallog << "OK" << std::endl;
    }
  else
    {
      test<2>();
    }
}
