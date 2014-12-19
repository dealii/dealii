// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
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



// like _01, but say that we don't want to coarsen any cells at
// all. make sure this does indeed happen
//
// the test is a bit fickle: the threshold for refinement should be
// right at 17, based on counting the error distribution. If we
// iterate 10 times in
// parallel::distributed::GridRefinement::refine_and_coarsen_fixed_fraction,
// we end up with a threshold just above that, yielding the original
// number of refined cells; however, in r21925 we changed the number
// of iterations to 25, which yields a refinement threshold just
// *under* 17, and a different number of cells to refine. I guess
// that's just how things go with these sort of operations

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/vector.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/base/utilities.h>


#include <fstream>


void test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);

  parallel::distributed::Triangulation<2> tr(MPI_COMM_WORLD);

  std::vector<unsigned int> sub(2);
  sub[0] = 5*Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);
  sub[1] = 1;
  GridGenerator::subdivided_hyper_rectangle(static_cast<Triangulation<2>&>(tr),
                                            sub, Point<2>(0,0), Point<2>(1,1));
  tr.refine_global (1);

  Vector<float> indicators (tr.dealii::Triangulation<2>::n_active_cells());
  {
    unsigned int cell_index = 0;
    unsigned int my_cell_index = 0;
    for (Triangulation<2>::active_cell_iterator
         cell = tr.begin_active(); cell != tr.end(); ++cell, ++cell_index)
      if (cell->subdomain_id() == myid)
        {
          ++my_cell_index;
          indicators(cell_index) = my_cell_index+1;
        }
    Assert (my_cell_index == 20, ExcInternalError());
  }

  parallel::distributed::GridRefinement
  ::refine_and_coarsen_fixed_fraction (tr, indicators,
                                       (74.+0.5)/210, 0);

  // now count number of cells
  // flagged for refinement and
  // coarsening. we have to
  // accumulate over all processors
  unsigned int my_refined   = 0,
               my_coarsened = 0;
  for (Triangulation<2>::active_cell_iterator
       cell = tr.begin_active(); cell != tr.end(); ++cell)
    if (cell->refine_flag_set())
      ++my_refined;
    else if (cell->coarsen_flag_set())
      ++my_coarsened;

  unsigned int n_refined   = 0,
               n_coarsened = 0;
  MPI_Reduce (&my_refined, &n_refined, 1, MPI_UNSIGNED, MPI_SUM, 0,
              MPI_COMM_WORLD);
  MPI_Reduce (&my_coarsened, &n_coarsened, 1, MPI_UNSIGNED, MPI_SUM, 0,
              MPI_COMM_WORLD);

  // make sure we have indeed flagged
  // exactly 20% of cells
  if (myid == 0)
    {
      deallog << "total active cells = "
              << tr.n_global_active_cells() << std::endl;
      deallog << "n_refined = " << n_refined << std::endl;
      // the following should be 0
      deallog << "n_coarsened = " << n_coarsened << std::endl;
    }

  tr.execute_coarsening_and_refinement ();
  if (myid == 0)
    deallog << "total active cells = "
            << tr.n_global_active_cells() << std::endl;
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile("output");
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      test();
    }
  else
    test();

}
