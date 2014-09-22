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



// Like refinement_02, but with a complex grid

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
    std::ifstream in (SOURCE_DIR "/../deal.II/grid_in_02/2d.xda");
    try
      {
        gi.read_xda (in);
      }
    catch (const typename Triangulation<dim>::DistortedCellList &distorted_cells)
      {
        // ignore distorted cells
        deallog << distorted_cells.distorted_cells.size()
                << " distorted cells after creating mesh."
                << std::endl;
      }
  }

  {
    GridIn<dim> gi;
    gi.attach_triangulation (tr2);
    std::ifstream in (SOURCE_DIR "/../deal.II/grid_in_02/2d.xda");
    try
      {
        gi.read_xda (in);
      }
    catch (const typename Triangulation<dim>::DistortedCellList &distorted_cells)
      {
        // ignore distorted cells
        deallog << distorted_cells.distorted_cells.size()
                << " distorted cells after creating mesh."
                << std::endl;
      }
  }

  Assert (tr.n_active_cells() == tr2.n_active_cells(),
          ExcInternalError());


  for (unsigned int i=0; i<1; ++i)
    {
      std::vector<bool> flags (tr.n_active_cells());
      for (unsigned int j=0; j<flags.size(); ++j)
        flags[j] = (Testing::rand() < 0.2*RAND_MAX);

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

      try
        {
          tr.execute_coarsening_and_refinement ();
        }
      catch (const typename Triangulation<dim>::DistortedCellList &distorted_cells)
        {
          // ignore distorted cells
          deallog << distorted_cells.distorted_cells.size()
                  << " distorted cells after adaptively refining mesh."
                  << std::endl;
        }
      try
        {
          tr2.execute_coarsening_and_refinement ();
        }
      catch (const typename Triangulation<dim>::DistortedCellList &distorted_cells)
        {
          // ignore distorted cells
          deallog << distorted_cells.distorted_cells.size()
                  << " distorted cells after adaptively refining mesh."
                  << std::endl;
        }

//      write_vtk(tr, "1");
      deallog << std::endl;

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

  deallog.push("2d");
  test<2>(logfile);
  deallog.pop();


}
