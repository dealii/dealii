//----------------------------  refine_and_coarsen_3d.cc  ---------------------------
//    $Id: refine_and_coarsen_for_parents_02.cc 23710 2011-05-17 04:50:10Z bangerth $
//    Version: $Name$ 
//
//    Copyright (C) 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------  refine_and_coarsen_for_parents_02.cc  ----------------------


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
void check (TRIA & tr)
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
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
  
  initlog(__FILE__);
  deal_II_exceptions::disable_abort_on_exception();

  Triangulation<2> tria;
  GridGenerator::hyper_cube (tria);
  tria.refine_global (2);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  check(tria);
}

  
  
