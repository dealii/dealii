//---------------------------------------------------------------------------
//    $Id: 2d_coarse_grid_01.cc 17444 2008-10-31 19:35:14Z bangerth $
//    Version: $Name$
//
//    Copyright (C) 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


// create a parallel DoFHandler on a 2d mesh and check componentwise
// renumbering


#include "../tests.h"
#include "coarse_grid_common.h"
#include <base/logstream.h>
#include <base/tensor.h>
#include <grid/tria.h>
#include <distributed/tria.h>
#include <grid/grid_generator.h>
#include <grid/grid_out.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_tools.h>
#include <dofs/dof_renumbering.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>

#include <fe/fe_q.h>
#include <fe/fe_system.h>

#include <fstream>


template<int dim>
void test()
{
  unsigned int myid = Utilities::System::get_this_mpi_process (MPI_COMM_WORLD);
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tr);

  const FE_Q<dim> fe_q(1);
  FESystem<dim> fe(fe_q,2);
  {
    tr.refine_global (2);
    tr.begin_active()->set_refine_flag();
    tr.execute_coarsening_and_refinement();

    DoFHandler<dim> dofh(tr);
    dofh.distribute_dofs (fe);

    {
      IndexSet dof_set;
      DoFTools::extract_locally_active_dofs (dofh, dof_set);
				       //if (myid==0)
      if (Utilities::System::get_this_mpi_process (MPI_COMM_WORLD) == 0)
	dof_set.print(deallog);
    }

    if (Utilities::System::get_this_mpi_process (MPI_COMM_WORLD) == 0)
      deallog << "****" << std::endl;
    DoFRenumbering::component_wise(dofh);
    {
      IndexSet dof_set;
      DoFTools::extract_locally_active_dofs (dofh, dof_set);
				       //if (myid==0)
      if (Utilities::System::get_this_mpi_process (MPI_COMM_WORLD) == 0)
	dof_set.print(deallog);

      if (myid==0)
	{

	  std::vector<unsigned int> local_dof_indices;
	typename DoFHandler<dim>::active_cell_iterator
	  cell, endc = dofh.end();

	if (Utilities::System::get_this_mpi_process (MPI_COMM_WORLD) == 0)
	  for (cell = dofh.begin_active(); cell != endc; ++cell)
	    if (!cell->is_artificial() && !cell->is_ghost())
	      {
		local_dof_indices.resize (cell->get_fe().dofs_per_cell);
		cell->get_dof_indices (local_dof_indices);
		for (unsigned int i=0;i< cell->get_fe().dofs_per_cell;++i)
		  deallog << local_dof_indices[i] << " ";
		deallog << std::endl;
	      }
	}

    }
  }

  if (0){
    tr.refine_global (1);
    DoFHandler<dim> dofh(tr);
    dofh.distribute_dofs (fe);

  }

  if (Utilities::System::get_this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << "OK" << std::endl;
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


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile(output_file_for_mpi("p4est_2d_renumber_01").c_str());
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
      deallog.push("2d");
      test<2>();
      deallog.pop();

      deallog.push("3d");
      test<3>();
      deallog.pop();
    }

#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
  MPI_Finalize();
#endif
}
