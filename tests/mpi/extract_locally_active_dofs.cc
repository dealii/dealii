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


// Test DoFTools::extract_locally_active_dofs and ensure that it returns the
// same result as DoFTools::extract_dofs_with_subdomain_association()

#include "../tests.h"
#include "coarse_grid_common.h"
#include <base/logstream.h>
#include <base/tensor.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <distributed/tria.h>
#include <grid/grid_generator.h>
#include <grid/grid_out.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_tools.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>

#include <fe/fe_q.h>

#include <fstream>


template<int dim>
void test()
{
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_ball(tr);
  tr.refine_global (2);

  const FE_Q<dim> fe(2);
  DoFHandler<dim> dofh(tr);
  dofh.distribute_dofs (fe);

  IndexSet locally_active;
  DoFTools::extract_locally_active_dofs (dofh, locally_active);

  Assert (locally_active ==
	  DoFTools::dof_indices_with_subdomain_association (dofh,
							    tr.locally_owned_subdomain()),
	  ExcInternalError());
  // Assert (locally_active.n_elements() ==
  // 	  DoFTools::count_dofs_with_subdomain_association (dofh,
  // 							   tr.locally_owned_subdomain()),
  // 	  ExcInternalError());

  if (Utilities::System::get_this_mpi_process (MPI_COMM_WORLD) == 0)
    {
      deallog << locally_active.size() << ' ' << locally_active.n_elements()
	      << std::endl;

      for (unsigned int i=0; i<locally_active.size(); ++i)
	if (locally_active.is_element(i))
	  deallog << i << ' ';
      deallog << "OK" << std::endl;
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
      std::ofstream logfile(output_file_for_mpi("extract_locally_active_dofs").c_str());
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
