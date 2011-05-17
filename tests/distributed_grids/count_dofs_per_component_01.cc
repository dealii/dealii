//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2008, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


// Test DoFTools::count_dofs_per_component


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/intergrid_map.h>
#include <deal.II/base/utilities.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>

#include <fstream>
#include <numeric>
#include <cstdlib>


template<int dim>
void test()
{
  parallel::distributed::Triangulation<dim>
    triangulation (MPI_COMM_WORLD,
		   Triangulation<dim>::limit_level_difference_at_vertices);

  FESystem<dim> fe (FE_Q<dim>(3),2,
		    FE_DGQ<dim>(1),1);

  DoFHandler<dim> dof_handler (triangulation);

  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global (2);
  dof_handler.distribute_dofs (fe);

  std::vector<unsigned int> dofs_per_component (fe.n_components());
  DoFTools::count_dofs_per_component (dof_handler, dofs_per_component);

  Assert (std::accumulate (dofs_per_component.begin(), dofs_per_component.end(), 0U)
	  == dof_handler.n_dofs(),
	  ExcInternalError());
      
  unsigned int myid = Utilities::System::get_this_mpi_process (MPI_COMM_WORLD);
  if (myid == 0)
    {
      deallog << "Total number of dofs: " << dof_handler.n_dofs() << std::endl;
      for (unsigned int i=0; i<dofs_per_component.size(); ++i)
	deallog << "Block " << i << " has " << dofs_per_component[i] << " global dofs"
		<< std::endl;
    }
}


int main(int argc, char *argv[])
{
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
  MPI_Init (&argc,&argv);
#else
  (void)argc;
  (void)argv;
#endif

  unsigned int myid = Utilities::System::get_this_mpi_process (MPI_COMM_WORLD);
  if (myid == 0)
    {
      std::ofstream logfile("count_dofs_per_component_01/output");
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      deallog.push("2d");
      test<2>();
      deallog.pop();

      deallog.push("3d");
      test<3>();
      deallog.pop();
    }
  else
    {
      test<2>();
      test<3>();
    }

#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
  MPI_Finalize();
#endif
}
