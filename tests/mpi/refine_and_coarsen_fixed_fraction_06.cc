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



// like _04 but with indicators and a mesh that matches step-32

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/vector.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/base/utilities.h>


#include <fstream>


namespace EquationData
{
  const double R0      = 6371000.-2890000.;
  const double R1      = 6371000.-  35000.;
}

void test()
{
  const unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);

  parallel::distributed::Triangulation<2> tr(MPI_COMM_WORLD);

  const unsigned int initial_refinement = 3;
  GridGenerator::hyper_shell (tr,
                              Point<2>(),
                              EquationData::R0,
                              EquationData::R1,
                              12,
                              true);
  static HyperShellBoundary<2> boundary;
  tr.set_boundary (0, boundary);
  tr.set_boundary (1, boundary);

  tr.refine_global (initial_refinement);

  Assert (tr.dealii::Triangulation<2>::n_active_cells()==768,
          ExcInternalError());

  // now read indicators
  Vector<float> indicators (tr.dealii::Triangulation<2>::n_active_cells());
  {
    std::ifstream in(SOURCE_DIR "/refine_and_coarsen_fixed_fraction_06/indicators");
    Assert(in, ExcMessage("File missing"));
    for (unsigned int i=0; i<indicators.size(); ++i)
      in >> indicators(i);
  }

  double min_indicator = *std::min_element (indicators.begin(),
                                            indicators.end()),
                         max_indicator = *std::max_element (indicators.begin(),
                                                            indicators.end());



  // use one strategy to compute
  // thresholds and obtain those
  // thresholds
  parallel::distributed::GridRefinement
  ::refine_and_coarsen_fixed_fraction (tr, indicators, 0.6, 0.2);
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
  ::refine_and_coarsen_fixed_fraction (tr, indicators, 0.6, 0.2);
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
