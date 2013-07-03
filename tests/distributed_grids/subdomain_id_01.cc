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


// check (level)subdomain_id(). Should be =0.

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



template <class TRIA>
void check (TRIA & tr)
{
  typename TRIA::cell_iterator cell = tr.begin(),
        endc = tr.end();
  
  for (; cell!=endc; ++cell)
    {
      deallog << cell->level_subdomain_id() << " ";
      try
	{
	  deallog << cell->subdomain_id();
	}
      catch (...)
	{
	  deallog << ".";
	  
	}
      deallog << std::endl;
    }
  
  deallog << "OK" << std::endl;
}


int main (int argc, char *argv[]) 
{
  deal_II_exceptions::disable_abort_on_exception();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  std::ofstream logfile("subdomain_id_01/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  Triangulation<2> tria;
  parallel::distributed::Triangulation<2> tria2(MPI_COMM_WORLD);
  GridGenerator::hyper_cube (tria);
  tria.refine_global (2);
  GridGenerator::hyper_cube (tria2);
  tria2.refine_global (2);
  check(tria);
  check(tria2);
}

  
  
