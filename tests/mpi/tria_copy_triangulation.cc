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



// Test interaction with p4est with a simple grid in 2d. here, we test
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

#include <fstream>
#include <cstdlib>


template<int dim>
void test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);

  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridGenerator::subdivided_hyper_cube(tr, 6-dim);

  parallel::distributed::Triangulation<dim> new_tr(MPI_COMM_WORLD);

  new_tr.copy_triangulation (tr);

  typename Triangulation<dim,dim>::active_cell_iterator cell1, cell2;

  for (cell1 = tr.begin_active(), cell2 = new_tr.begin_active();
       cell1 != tr.end(), cell2 != new_tr.end();
       ++cell1, ++cell2)
    {
      if (cell1->is_locally_owned ())
        {
          if (myid == 0)
            {
              deallog << "triangulation::cell     " << cell1 << ": locally owned, subdomain_id = " << cell1->subdomain_id () << std::endl;
              deallog << "new_triangulation::cell " << cell2 << ": locally owned, subdomain_id = " << cell2->subdomain_id () << std::endl;
            };

          Assert (cell2->is_locally_owned (), ExcInternalError());
          Assert (cell1->subdomain_id () == cell2->subdomain_id (), ExcInternalError());

          for (unsigned int vertex=0;
               vertex<GeometryInfo<dim>::vertices_per_cell;
               ++vertex)
            {
              Assert (cell1->vertex(vertex).distance (cell2->vertex(vertex)) < 1.e-14, ExcInternalError());

              if (myid == 0)
                deallog << "  vertices " << vertex << " coincide" << std::endl;
            };
        }
      else if (cell1->is_ghost())
        {
          if (myid == 0)
            {
              deallog << "triangulation::cell     " << cell1 << ": ghost, subdomain_id = " << cell1->subdomain_id () << std::endl;
              deallog << "new_triangulation::cell " << cell2 << ": ghost, subdomain_id = " << cell2->subdomain_id () << std::endl;
            };

          Assert (cell2->is_ghost (), ExcInternalError());
        }
      else if (cell1->is_artificial())
        {
          if (myid == 0)
            {
              deallog << "triangulation::cell     " << cell1 << ": artificial" << std::endl;
              deallog << "new_triangulation::cell " << cell2 << ": artificial" << std::endl;
            };

          Assert (cell2->is_artificial (), ExcInternalError());
        }
      else
        Assert (false, ExcInternalError());
    };

//   assert_tria_equal("tria_copy_triangulation", tr, new_tr);

  Assert (tr.n_active_cells() == new_tr.n_active_cells(), ExcInternalError());
  Assert (tr.n_levels() == new_tr.n_levels(), ExcInternalError());
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

      deallog.push("2d");
      test<2>();
      deallog.pop();
    }
  else
    test<2>();

}
