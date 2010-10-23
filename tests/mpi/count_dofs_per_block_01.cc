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


// test DoFTools::count_dofs_per_block in parallel


#include "../tests.h"
#include <base/logstream.h>
#include <base/tensor.h>
#include <distributed/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_generator.h>
#include <grid/intergrid_map.h>
#include <base/utilities.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_tools.h>
#include <fe/fe_system.h>
#include <fe/fe_q.h>
#include <fe/fe_dgq.h>

#include <fstream>
#include <cstdlib>
#include <numeric>


template<int dim>
void test()
{
  parallel::distributed::Triangulation<dim>
    triangulation (MPI_COMM_WORLD,
		   Triangulation<dim>::limit_level_difference_at_vertices);

  FESystem<dim> fe (FE_Q<dim>(3),2,
		    FE_DGQ<dim>(1),1);

  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global (4);

  DoFHandler<dim> dof_handler (triangulation);
  dof_handler.distribute_dofs (fe);

  std::vector<unsigned int> dofs_per_block (fe.n_blocks());
  DoFTools::count_dofs_per_block (dof_handler, dofs_per_block);

  if (Utilities::System::get_this_mpi_process (MPI_COMM_WORLD) == 0)
    for (unsigned int i=0; i<fe.n_blocks(); ++i)
      deallog << "Block " << i << " has " << dofs_per_block[i] << " dofs" << std::endl;
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
      std::ofstream logfile(output_file_for_mpi("count_dofs_per_block_01").c_str());
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
