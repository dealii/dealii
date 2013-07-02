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


// test FETools::interpolation_difference

#include "../tests.h"
#include "coarse_grid_common.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
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
#include <deal.II/fe/fe_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/lac/petsc_parallel_vector.h>

#include <fstream>




template<int dim>
void test()
{
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tr);

//tr.refine_global (2);

  const FE_Q<dim> fe(2);
  DoFHandler<dim> dofh1(tr);
  DoFHandler<dim> dofh2(tr);
  dofh1.distribute_dofs (fe);
  dofh2.distribute_dofs (fe);

  ConstraintMatrix cm1;
  cm1.close();
  ConstraintMatrix cm2;
  cm2.close();
  
  IndexSet  dof1_locally_owned_dofs = dofh1.locally_owned_dofs();
  IndexSet  dof2_locally_owned_dofs = dofh2.locally_owned_dofs();
  IndexSet  dof1_locally_relevant_dofs;
  IndexSet  dof2_locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs (dofh2,
					   dof2_locally_relevant_dofs);
  DoFTools::extract_locally_relevant_dofs (dofh1,
					   dof1_locally_relevant_dofs);
  
  PETScWrappers::MPI::Vector u1(MPI_COMM_WORLD, dof1_locally_owned_dofs, dof1_locally_relevant_dofs);

  PETScWrappers::MPI::Vector out(MPI_COMM_WORLD, dof1_locally_owned_dofs);

  FETools::interpolation_difference
    (dofh1, cm1, u1, dofh2, cm2, out);

  double norm = out.l2_norm();

  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << "norm = " << norm
	    << std::endl;
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);

  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile(output_file_for_mpi("interpolate_04").c_str());
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      deallog.push("2d");
      test<2>();
      deallog.pop();

      deallog.push("3d");
				       //    test<3>();
      deallog.pop();
    }
  else
    {
      deallog.push("2d");
      test<2>();
      deallog.pop();

      deallog.push("3d");
//      test<3>();
      deallog.pop();
    }
}
