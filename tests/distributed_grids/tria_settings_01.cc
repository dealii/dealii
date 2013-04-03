//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


// Test the setting
// parallel::distributed::Triangulation<dim>::mesh_reconstruction_after_repartitioning
// by looking which order the cells get enumerated. With the flag, the cells
// stay in z-order. Without it, the order in which the cells are created
// matters.

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

#include <fstream>

template<int dim>
void testit(parallel::distributed::Triangulation<dim> & tr)
{
  GridGenerator::hyper_cube(tr);
  tr.refine_global(1);

  (++tr.begin_active())->set_refine_flag();
  tr.execute_coarsening_and_refinement ();

  tr.begin_active()->set_refine_flag();
  tr.execute_coarsening_and_refinement ();

  typename  parallel::distributed::Triangulation<dim>::active_cell_iterator it=tr.begin_active();
  for (;it!=tr.end();++it)
    {
      deallog << it->center() << ", ";
    }
    deallog << std::endl;
}


template<int dim>
void test(std::ostream& /*out*/)
{
  {
    
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD,
					       dealii::Triangulation<dim>::none,
					       parallel::distributed::Triangulation<dim>::default_setting);
  testit(tr);
  
  }
  
  {
    
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD,
					       dealii::Triangulation<dim>::none,
					       parallel::distributed::Triangulation<dim>::mesh_reconstruction_after_repartitioning);
  testit(tr);
  
  }


  
  

}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  std::ofstream logfile("tria_settings_01/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  deallog.push("2d");
  test<2>(logfile);
  deallog.pop();
}
