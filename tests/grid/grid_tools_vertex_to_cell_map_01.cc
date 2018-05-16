//
// Copyright (C) 2015 - 2017 by the deal.II authors
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

// Test vertex_to_cell_map for 2D problem with mpi and hanging nodes.

#include "../tests.h"
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/cell_id.h>

#include <vector>

void
test()
{
  const MPI_Comm mpi_comm = MPI_COMM_WORLD;

  parallel::distributed::Triangulation<2> tria(mpi_comm);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);
  unsigned int rank = Utilities::MPI::this_mpi_process(mpi_comm);
  if (rank==0)
    {
      typename Triangulation<2>::active_cell_iterator cell = tria.begin_active(),
                                                      end_cell = tria.end();
      for (; cell!=end_cell; ++cell)
        if (cell->is_locally_owned())
          cell->set_refine_flag();
    }

  tria.execute_coarsening_and_refinement();

  std::vector<std::set<typename Triangulation<2>::active_cell_iterator> >
  vertex_to_cell = GridTools::vertex_to_cell_map(tria);

  AssertThrow(tria.n_vertices()==vertex_to_cell.size(),
              ExcMessage("Wrong number of vertices"));
  std::vector<unsigned int> n_cells;
  for (unsigned int i=0; i<vertex_to_cell.size(); ++i)
    n_cells.push_back(vertex_to_cell[i].size());

  if (rank==0)
    {
      std::vector<unsigned int> histogram(5,0);
      for (unsigned int i=0; i<n_cells.size(); ++i)
        histogram[n_cells[i]] += 1;

      AssertThrow(histogram[0]==0, ExcMessage("Wrong cell distribution"));
      AssertThrow(histogram[1]==4, ExcMessage("Wrong cell distribution"));
      AssertThrow(histogram[2]==20, ExcMessage("Wrong cell distribution"));
      AssertThrow(histogram[3]==4, ExcMessage("Wrong cell distribution"));
      AssertThrow(histogram[4]==27, ExcMessage("Wrong cell distribution"));
    }
  if (rank==1)
    {
      std::vector<unsigned int> histogram(5,0);
      for (unsigned int i=0; i<n_cells.size(); ++i)
        histogram[n_cells[i]] += 1;

      AssertThrow(histogram[0]==0, ExcMessage("Wrong cell distribution"));
      AssertThrow(histogram[1]==4, ExcMessage("Wrong cell distribution"));
      AssertThrow(histogram[2]==18, ExcMessage("Wrong cell distribution"));
      AssertThrow(histogram[3]==6, ExcMessage("Wrong cell distribution"));
      AssertThrow(histogram[4]==24, ExcMessage("Wrong cell distribution"));
    }
}


int
main (int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);
  MPILogInitAll log;

  test();

  deallog << "OK" << std::endl;

  return 0;
}
