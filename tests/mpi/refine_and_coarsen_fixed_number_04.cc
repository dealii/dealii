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



// create a triangulation on a single processor and test that
// parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number
// does roughly the same as
// ::GridRefinement::refine_and_coarsen_fixed_number

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/vector.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/base/utilities.h>


#include <fstream>


void test()
{
  const unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);

  parallel::distributed::Triangulation<2> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(static_cast<Triangulation<2>&>(tr));
  tr.refine_global (4);

  Vector<float> indicators (tr.dealii::Triangulation<2>::n_active_cells());
  float min_indicator = tr.dealii::Triangulation<2>::n_active_cells(),
        max_indicator = 0;
  {
    unsigned int cell_index = 0;
    unsigned int my_cell_index = 0;
    for (Triangulation<2>::active_cell_iterator
         cell = tr.begin_active(); cell != tr.end(); ++cell, ++cell_index)
      if (cell->subdomain_id() == myid)
        {
          ++my_cell_index;
          indicators(cell_index) = my_cell_index;
          min_indicator = std::min (min_indicator, indicators(cell_index));
          max_indicator = std::max (max_indicator, indicators(cell_index));
        }
  }

  // use one strategy to compute
  // thresholds and obtain those
  // thresholds
  parallel::distributed::GridRefinement
  ::refine_and_coarsen_fixed_number (tr, indicators, 2./3, 1./6);
  {
    float coarsen_indicator = min_indicator-1,
          refine_indicator  = max_indicator+1;
    unsigned int cell_index=0;
    for (Triangulation<2>::active_cell_iterator
         cell = tr.begin_active(); cell != tr.end(); ++cell, ++cell_index)
      if (cell->refine_flag_set())
        refine_indicator = std::min (refine_indicator,
                                     indicators(cell_index));
      else if (cell->coarsen_flag_set())
        coarsen_indicator = std::max (coarsen_indicator,
                                      indicators(cell_index));
    if (myid == 0)
      {
        deallog << "thresholds = " << refine_indicator << ' '
                << coarsen_indicator << std::endl;
      }

    for (Triangulation<2>::active_cell_iterator
         cell = tr.begin_active(); cell != tr.end(); ++cell)
      {
        cell->clear_coarsen_flag();
        cell->clear_refine_flag();
      }
  }

  // now use the second strategy to
  // compute thresholds and obtain
  // those thresholds. note that this
  // only works because we are
  // working on only a single
  // processor
  dealii::GridRefinement
  ::refine_and_coarsen_fixed_number (tr, indicators, 2./3, 1./6);
  {
    float coarsen_indicator = min_indicator-1,
          refine_indicator  = max_indicator+1;
    unsigned int cell_index=0;
    for (Triangulation<2>::active_cell_iterator
         cell = tr.begin_active(); cell != tr.end(); ++cell, ++cell_index)
      if (cell->refine_flag_set())
        refine_indicator = std::min (refine_indicator,
                                     indicators(cell_index));
      else if (cell->coarsen_flag_set())
        coarsen_indicator = std::max (coarsen_indicator,
                                      indicators(cell_index));
    if (myid == 0)
      {
        deallog << "thresholds = " << refine_indicator << ' '
                << coarsen_indicator << std::endl;
      }

    for (Triangulation<2>::active_cell_iterator
         cell = tr.begin_active(); cell != tr.end(); ++cell)
      {
        cell->clear_coarsen_flag();
        cell->clear_refine_flag();
      }
  }

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
