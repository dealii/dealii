// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



// check CellId

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include <sstream>

#include "../tests.h"


template <int dim>
void
check(Triangulation<dim> &tr)
{
  typename Triangulation<dim>::cell_iterator cell = tr.begin(), endc = tr.end();


  for (; cell != endc; ++cell)
    {
      deallog << cell->level() << " " << cell->index() << std::endl;

      // Store the CellId, convert it to a binary representation,
      // create a new CellId from that, and create a cell iterator
      // pointing to the same cell

      const CellId cid = cell->id();

      const CellId::binary_type cids = cid.to_binary<dim>();

      const CellId cid2(cids);

      typename Triangulation<dim>::cell_iterator cell2 = cid2.to_cell(tr);

      Assert(cid2 == cid, ExcInternalError());
      Assert(cell2 == cell, ExcInternalError());

      deallog << cell2->level() << " " << cell2->index() << std::endl;
    }

  deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  // Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv,
  // testing_max_num_threads());

  initlog();
  deal_II_exceptions::disable_abort_on_exception();

  Triangulation<2> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(3);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  check(tria);

  Triangulation<3> tria2;
  GridGenerator::hyper_cube(tria2);
  tria2.refine_global(1);
  tria2.begin_active()->set_refine_flag();
  tria2.execute_coarsening_and_refinement();
  check(tria2);
}
