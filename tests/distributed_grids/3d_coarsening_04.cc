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



// Like coarsening_03, but with a much more complex grid

#include "../tests.h"
#include "coarse_grid_common.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/intergrid_map.h>

#include <fstream>
#include <cstdlib>


template<int dim>
void test(std::ostream & /*out*/)
{
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);
  Triangulation<dim> tr2 (Triangulation<dim>::limit_level_difference_at_vertices);

  {
    GridIn<dim> gi;
    gi.attach_triangulation (tr);
    std::ifstream in (SOURCE_DIR "/../deal.II/grid_in_3d_02/747.ucd");
    gi.read (in);
    //tr.refine_global (1);
  }

  {
    GridIn<dim> gi;
    gi.attach_triangulation (tr2);
    std::ifstream in (SOURCE_DIR "/../deal.II/grid_in_3d_02/747.ucd");
    gi.read (in);
    //tr2.refine_global (1);
  }

  Assert (tr.n_active_cells() == tr2.n_active_cells(),
          ExcInternalError());
  deallog << " Number of cells: "
          << tr.n_active_cells() << ' '
          << tr2.n_active_cells()
          << std::endl;


  for (unsigned int i=0; i<2; ++i)
    {
      // refine one-tenth of cells randomly
      std::vector<bool> flags (tr.n_active_cells(), false);
      for (unsigned int k=0; k<flags.size()/30 + 1; ++k)
        flags[Testing::rand() % flags.size()] = true;
      // make sure there's at least one that
      // will be refined
      flags[0] = true;

      InterGridMap<Triangulation<dim> > intergrid_map;
      intergrid_map.make_mapping (tr, tr2);

      // refine tr and tr2
      unsigned int index=0;
      for (typename Triangulation<dim>::active_cell_iterator
           cell = tr.begin_active();
           cell != tr.end(); ++cell, ++index)
        if (flags[index])
          {
            cell->set_refine_flag();
            intergrid_map[cell]->set_refine_flag();
          }
      Assert (index == tr.n_active_cells(), ExcInternalError());

      // flag all other cells for coarsening
      // (this should ensure that at least
      // some of them will actually be
      // coarsened)
      index=0;
      for (typename Triangulation<dim>::active_cell_iterator
           cell = tr.begin_active();
           cell != tr.end(); ++cell, ++index)
        if (!flags[index])
          {
            cell->set_coarsen_flag();
            intergrid_map[cell]->set_coarsen_flag();
          }

      tr.execute_coarsening_and_refinement ();
      tr2.execute_coarsening_and_refinement ();

      deallog << i << " Number of cells: "
              << tr.n_active_cells() << ' '
              << tr2.n_active_cells()
              << std::endl;

      assert_tria_equal(tr, tr2);

    }
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  deallog.push("3d");
  test<3>(logfile);
  deallog.pop();


}
