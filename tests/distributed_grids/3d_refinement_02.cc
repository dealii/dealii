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



// Test interaction with p4est with a simple grid in 3d. here, we test
// that refining a mesh actually works, where we refine several times
// more or less randomly by choosing a single cell for refinement

// in addition to just refining, have a parallel triangulation that should
// look the same but isn't attached to a distribution manager, and compare
// with it in each step

// note that due to the way we set up triangulations in parallel, we
// can't expect that cells on the two meshes have the same numbers (or
// are in the same order, for that matter). we therefore use an
// IntergridMap. there should be matching cells, however

// note that p4est refines meshes in a way that is equivalent to
// specifying Triangulation<dim>::limit_level_difference_at_vertices,
// so this is what we also give as argument to the mesh with which we
// compare

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

  GridGenerator::hyper_cube(tr);
  tr.refine_global (1);

  GridGenerator::hyper_cube(tr2);
  tr2.refine_global (1);

  Assert (tr.n_active_cells() == tr2.n_active_cells(),
          ExcInternalError());


  for (unsigned int i=0; i<15-2*dim; ++i)
    {
      std::vector<bool> flags (tr.n_active_cells(), false);
      {
        const unsigned int x = Testing::rand() % flags.size();
        deallog << "Refining cell " << x << std::endl;
        flags[x] = true;
      }

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
      tr.execute_coarsening_and_refinement ();
      tr2.execute_coarsening_and_refinement ();

      write_vtk(tr, "1");
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

  deallog.push("3d");
  test<3>(logfile);
  deallog.pop();


}
