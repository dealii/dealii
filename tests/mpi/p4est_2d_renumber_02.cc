//---------------------------------------------------------------------------
//    $Id$
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
//
// this is like the _01 testcase but with a simpler mesh for which I
// can actually count DoFs and decide who owns what


#include "../tests.h"
#include "coarse_grid_common.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <fstream>


template<int dim>
void test()
{
  unsigned int myid = Utilities::System::get_this_mpi_process (MPI_COMM_WORLD);
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  std::vector<unsigned int> sub(2);
  sub[0] = Utilities::System::get_n_mpi_processes (MPI_COMM_WORLD);
  sub[1] = 1;
  GridGenerator::subdivided_hyper_rectangle(static_cast<Triangulation<dim>&>(tr),
					    sub, Point<2>(0,0), Point<2>(1,1));

  const FE_Q<dim> fe_q(1);
  FESystem<dim> fe(fe_q,2);
  DoFHandler<dim> dofh(tr);
  dofh.distribute_dofs (fe);

  if (Utilities::System::get_this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << "Total dofs=" << dofh.n_dofs() << std::endl;

  {
    IndexSet dof_set;
    DoFTools::extract_locally_active_dofs (dofh, dof_set);
    if (Utilities::System::get_this_mpi_process (MPI_COMM_WORLD) == 0)
      dof_set.print(deallog);
  }

  if (Utilities::System::get_this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << "****" << std::endl;

  DoFRenumbering::component_wise(dofh);
  {
    IndexSet dof_set;
    DoFTools::extract_locally_active_dofs (dofh, dof_set);

    if (Utilities::System::get_this_mpi_process (MPI_COMM_WORLD) == 0)
      {
	dof_set.print(deallog);
	for (unsigned int i=0; i<Utilities::System::get_n_mpi_processes (MPI_COMM_WORLD); ++i)
	  {
	    deallog << "Dofs owned by processor " << i << ": ";
	    dofh.locally_owned_dofs_per_processor()[i].print(deallog);
	    deallog << std::endl;
	  }
      }

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
      std::ofstream logfile(output_file_for_mpi("p4est_2d_renumber_02").c_str());
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      deallog.push("2d");
      test<2>();
      deallog.pop();
    }
  else
    {
      deallog.push("2d");
      test<2>();
      deallog.pop();
    }

#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
  MPI_Finalize();
#endif
}
