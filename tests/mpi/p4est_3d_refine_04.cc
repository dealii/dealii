// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2017 by the deal.II authors
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



// create a relatively simple 3d mesh and refine a few cells. this
// should produce the same output regardless of the number of
// processors
//
// this test is extracted from step-32 while searching for some
// apparently rather subtle and indeed pretty frustrating bug...

#include "../tests.h"
#include <deal.II/base/tensor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/base/utilities.h>

#include <ostream>


template <int dim>
void test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);

  if (true)
    {
      parallel::distributed::Triangulation<dim>
      triangulation(MPI_COMM_WORLD);

      // create a mesh for the
      // hypershell and refine it
      // twice globally
      const double R0 = 6371000.-2890000.;
      const double R1 = 6371000.-  35000.;
      GridGenerator::hyper_shell (triangulation,
                                  Point<dim>(),
                                  R0,
                                  R1,
                                  (dim==3) ? 96 : 12,
                                  true);
      triangulation.reset_manifold(0);
      triangulation.refine_global (2);

      // then flag all cells for
      // refinement that are close to
      // the north pole
      unsigned int x_flagged_cells = 0;
      for (typename Triangulation<dim>::active_cell_iterator
           cell = triangulation.begin_active();
           cell != triangulation.end(); ++cell)
        if (!cell->is_ghost() && !cell->is_artificial())
          if (cell->center()[2] > R1*0.75)
            {
              ++x_flagged_cells;
              cell->set_refine_flag();
            }

      // count how many cells
      // actually have the refinement
      // flag after preparing for
      // refinement, and ensure that
      // this number is still the
      // same as before
      triangulation.prepare_coarsening_and_refinement ();
      {
        unsigned int n_flagged_cells = 0;
        for (typename Triangulation<dim>::active_cell_iterator
             cell = triangulation.begin_active();
             cell != triangulation.end(); ++cell)
          if (!cell->is_ghost() && !cell->is_artificial())
            if (cell->refine_flag_set())
              ++n_flagged_cells;

        Assert (n_flagged_cells == x_flagged_cells,
                ExcInternalError());

        unsigned int global_f_c = 0;
        MPI_Allreduce (&n_flagged_cells, &global_f_c, 1, MPI_UNSIGNED,
                       MPI_SUM, MPI_COMM_WORLD);

        if (myid == 0)
          deallog << "# flagged cells = " << global_f_c << std::endl;
      }

      triangulation.execute_coarsening_and_refinement ();

      if (myid == 0)
        deallog << "#cells = " << triangulation.n_global_active_cells()
                << std::endl;
    }

  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << "OK" << std::endl;
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);

  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      initlog();

      deallog.push("3d");
      test<3>();
      deallog.pop();
    }
  else
    test<3>();

}
