// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2013 by the deal.II authors
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



// extracted from count_dofs_per_block that crashed at one point


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/intergrid_map.h>
#include <deal.II/base/utilities.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <numeric>
#include <cstdlib>


template<int dim>
void test()
{
  parallel::distributed::Triangulation<dim>
  triangulation (MPI_COMM_WORLD,
                 Triangulation<dim>::limit_level_difference_at_vertices);

  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global (1);

  {
    typename Triangulation<dim>::cell_iterator
    it=triangulation.begin_active();
    it->set_refine_flag();
    triangulation.execute_coarsening_and_refinement ();

    it=triangulation.begin(1);
    ++it;
    ++it;
    ++it;
    it=triangulation.begin(1);
    for (unsigned int a=0; a<4; a++)
      it->child(a)->set_coarsen_flag();

    triangulation.prepare_coarsening_and_refinement ();

    {
      //data out
      FE_Q<dim>            fe(1);
      DoFHandler<dim>      dof_handler(triangulation);

      dof_handler.distribute_dofs(fe);

      unsigned int n_coarse=0;
      Vector<float> subdomain (triangulation.n_active_cells());
      {
        unsigned int index = 0;

        for (typename Triangulation<dim>::active_cell_iterator
             cell = triangulation.begin_active();
             cell != triangulation.end(); ++cell, ++index)
          {
            subdomain(index)=0;

            if (cell->is_ghost() || cell->is_artificial())
              subdomain(index)=-4;

            if (cell->refine_flag_set())
              subdomain(index)+=1;
            if (cell->coarsen_flag_set())
              {
                subdomain(index)+=2;
                ++n_coarse;
              }
          }
      }
      if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
        {
          deallog << "id=" << triangulation.locally_owned_subdomain()
                  << " n_coarsen=" << n_coarse << std::endl;
          for (unsigned int i=0; i<subdomain.size(); ++i)
            deallog << subdomain(i) << std::endl;
        }

      triangulation.execute_coarsening_and_refinement ();
    }
  }
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  if (myid == 0)
    {
      std::ofstream logfile("output");
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

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
