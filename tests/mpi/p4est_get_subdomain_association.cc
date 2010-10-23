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


// we couldn't call DoFTools::dof_indices_with_subdomain_association
// without causing an abort
//
// in this test, let processor 1 output stuff since its index range
// is non-contiguous

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
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>

#include <fe/fe_q.h>

#include <fstream>


template<int dim>
void test()
{
  unsigned int myid = Utilities::System::get_this_mpi_process (MPI_COMM_WORLD);

  if (Utilities::System::get_this_mpi_process (MPI_COMM_WORLD) == 1)
    deallog << "hyper_cube" << std::endl;

  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tr);
  tr.refine_global (3);
  DoFHandler<dim> dofh(tr);

  static const FE_Q<dim> fe(2);
  dofh.distribute_dofs (fe);

   if (myid==1)
     {
       deallog << "dofh.n_dofs() " << dofh.n_locally_owned_dofs_per_processor() << std::endl;
       deallog << "dofh.n_locally_owned_dofs() " << dofh.n_locally_owned_dofs() << std::endl;

       const IndexSet set
	 = DoFTools::dof_indices_with_subdomain_association
	 (static_cast<const DoFHandler<dim>&>(dofh),
	  tr.locally_owned_subdomain());

       deallog << set.n_elements() << std::endl;
       for (unsigned int i=0; i<set.n_elements(); ++i)
	 deallog << "   " << set.nth_index_in_set(i) << std::endl;
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

  if (myid == 1)
    {
      std::ofstream logfile(output_file_for_mpi("p4est_get_subdomain_association").c_str());
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      deallog.push("2d");
      test<2>();
      deallog.pop();
    }
  else
    test<2>();

#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
  MPI_Finalize();
#endif
}
