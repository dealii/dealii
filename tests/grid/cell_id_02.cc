// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2013 by the deal.II authors
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



// check CellId

#include "../tests.h"

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>

#include <fstream>
#include <sstream>

template <class TRIA>
void check (TRIA &tr)
{
  typename TRIA::cell_iterator cell = tr.begin(),
                               endc = tr.end();


  for (; cell!=endc; ++cell)
    {
      std::ostringstream outb;
      outb << cell->id();
      CellId tmp;
      std::istringstream in(outb.str());
      in >> tmp;
      deallog << cell->level() << " " << cell->index() << " " << cell->id() << " " << tmp << std::endl;
    }


  CellId empty;

  Assert(tr.begin()->id() != tr.begin_active()->id(), ExcInternalError());
  Assert(tr.begin()->id() != empty, ExcInternalError());
  Assert(tr.begin()->id() == tr.begin()->id(), ExcInternalError());

  deallog << "OK" << std::endl;
}


int main (int argc, char *argv[])
{
  // Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);

  initlog();
  deal_II_exceptions::disable_abort_on_exception();

  Triangulation<2> tria;
  GridGenerator::hyper_cube (tria);
  tria.refine_global (2);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  check(tria);
}



