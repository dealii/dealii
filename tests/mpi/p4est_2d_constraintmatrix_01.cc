//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2009, 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


// check ConstraintMatrix for a distributed mesh,
// also compare with/without sparse line_cache via IndexSet.
// Simple mesh.

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
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_q.h>

#include <fstream>
#include <sstream>

template<int dim>
void test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << "hyper_cube" << std::endl;

  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tr);
  tr.refine_global (1);
  tr.begin_active()->set_refine_flag();
  tr.execute_coarsening_and_refinement();

  DoFHandler<dim> dofh(tr);

  static const FE_Q<dim> fe(1);
  dofh.distribute_dofs (fe);

  IndexSet dof_set;
  DoFTools::extract_locally_active_dofs (dofh, dof_set);

  ConstraintMatrix cm;
  DoFTools::make_hanging_node_constraints (dofh, cm);
  ConstraintMatrix cm2(dof_set);
  DoFTools::make_hanging_node_constraints (dofh, cm2);

  {
    std::ofstream file((std::string("p4est_2d_constraintmatrix_01/ncpu_") + Utilities::int_to_string(Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD)) + "/dat." + Utilities::int_to_string(myid)).c_str());
    file << "**** proc " << myid << std::endl;
    cm.print(file);
    file << "****" << std::endl;
    cm2.print(file);
  }

  MPI_Barrier(MPI_COMM_WORLD);




  if (myid==0)
    {
      for (unsigned int i=0;i<numproc;++i)
	{
	  cat_file((std::string("p4est_2d_constraintmatrix_01/ncpu_") + Utilities::int_to_string(Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD)) + "/dat." + Utilities::int_to_string(i)).c_str());
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

  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile(output_file_for_mpi("p4est_2d_constraintmatrix_01").c_str());
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
