// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Same as tria_copy_triangulation_01, but for locally refined meshes.


#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>

#include "../tests.h"


template <int dim>
void
test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  parallel::distributed::Triangulation<dim> tria1(MPI_COMM_WORLD);
  {
    GridGenerator::hyper_shell(tria1, Point<dim>(), 1., 2., 0, true);
    tria1.begin_active()->set_refine_flag();
    tria1.execute_coarsening_and_refinement();
  }

  parallel::distributed::Triangulation<dim> tria2(MPI_COMM_WORLD);
  tria2.copy_triangulation(tria1);

  Assert(tria1.n_active_cells() == tria2.n_active_cells(), ExcInternalError());
  Assert(tria1.n_levels() == tria2.n_levels(), ExcInternalError());

  typename Triangulation<dim, dim>::active_cell_iterator cell1, cell2;
  for (cell1 = tria1.begin_active(), cell2 = tria2.begin_active();
       cell1 != tria1.end();
       ++cell1, ++cell2)
    {
      for (unsigned int f = 0; f < cell1->n_faces(); ++f)
        if (cell1->face(f)->at_boundary())
          Assert(cell1->face(f)->boundary_id() == cell2->face(f)->boundary_id(),
                 ExcInternalError());
      Assert(cell1->manifold_id() == cell2->manifold_id(), ExcInternalError());

      if (cell1->is_locally_owned())
        {
          if (myid == 0)
            {
              deallog << "triangulation::cell     " << cell1
                      << ": locally owned, subdomain_id = "
                      << cell1->subdomain_id() << std::endl;
              deallog << "new_triangulation::cell " << cell2
                      << ": locally owned, subdomain_id = "
                      << cell2->subdomain_id() << std::endl;
            };

          Assert(cell2->is_locally_owned(), ExcInternalError());
          Assert(cell1->subdomain_id() == cell2->subdomain_id(),
                 ExcInternalError());

          for (const unsigned int vertex : GeometryInfo<dim>::vertex_indices())
            {
              Assert(cell1->vertex(vertex).distance(cell2->vertex(vertex)) <
                       1.e-14,
                     ExcInternalError());

              if (myid == 0)
                deallog << "  vertices " << vertex << " coincide" << std::endl;
            };
        }
      else if (cell1->is_ghost())
        {
          if (myid == 0)
            {
              deallog << "triangulation::cell     " << cell1
                      << ": ghost, subdomain_id = " << cell1->subdomain_id()
                      << std::endl;
              deallog << "new_triangulation::cell " << cell2
                      << ": ghost, subdomain_id = " << cell2->subdomain_id()
                      << std::endl;
            };

          Assert(cell2->is_ghost(), ExcInternalError());
        }
      else if (cell1->is_artificial())
        {
          if (myid == 0)
            {
              deallog << "triangulation::cell     " << cell1 << ": artificial"
                      << std::endl;
              deallog << "new_triangulation::cell " << cell2 << ": artificial"
                      << std::endl;
            };

          Assert(cell2->is_artificial(), ExcInternalError());
        }
      else
        Assert(false, ExcInternalError());
    };
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      initlog();

      deallog.push("2d");
      test<2>();
      deallog.pop();
    }
  else
    test<2>();
}
