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


// check whether all constraints are set on locally active set
// Mesh: hypercube where 7 of 8 octants are refined.

#include "../tests.h"
#include "coarse_grid_common.h"
#include <deal.II/distributed/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/numerics/data_out.h>


template <int dim>
void test ()
{
  parallel::distributed::Triangulation<dim> tria (MPI_COMM_WORLD);
  FE_Q<dim> fe(2);
  DoFHandler<dim> dof (tria);

  GridGenerator::hyper_cube (tria, 0, 1);
  tria.refine_global(1);

				// refine all cells except the one located at
				// [0.25, 0.25, 0.25]
  for (typename Triangulation<dim>::active_cell_iterator
	 cell = tria.begin_active();
       cell != tria.end(); ++cell)
    {
      Point<dim> diff(0.25,0.25,0.25);
      diff -= cell->center();
      if (diff.norm() > 0.25 && !(cell->is_ghost() || cell->is_artificial()))
	cell->set_refine_flag();
    }
  tria.execute_coarsening_and_refinement();
  dof.distribute_dofs(fe);

				// build constraint matrix
  IndexSet locally_relevant (dof.n_dofs());
  DoFTools::extract_locally_relevant_dofs (dof, locally_relevant);
  ConstraintMatrix constraints (locally_relevant);
  DoFTools::make_hanging_node_constraints (dof, constraints);
  constraints.close();

				// print out constraints for each
				// processor.
  const unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  const unsigned int numproc = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

  IndexSet locally_active (dof.n_dofs());
  DoFTools::extract_locally_active_dofs (dof, locally_active);
  std::ofstream file((std::string("p4est_3d_constraintmatrix_04/ncpu_") + Utilities::int_to_string(Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD)) + "/dat." + Utilities::int_to_string(myid)).c_str());
  file << "**** proc " << myid << ": \n\n";
  file << "Constraints:\n";
  constraints.print(file);
  file << std::endl;

  file << "My active dofs: ";
  locally_active.print(file);
  file << std::endl;
  MPI_Barrier(MPI_COMM_WORLD);

  if (myid==0)
    {
      for (unsigned int i=0;i<numproc;++i)
	{
	  cat_file((std::string("p4est_3d_constraintmatrix_04/ncpu_") + Utilities::int_to_string(Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD)) + "/dat." + Utilities::int_to_string(i)).c_str());
	}
    }
}

int main(int argc, char** argv)
{
  Utilities::System::MPI_InitFinalize mpi_initialization(argc, argv);
  const unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);

  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile(output_file_for_mpi("p4est_3d_constraintmatrix_04").c_str());
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      deallog.push("3d");
      test<3>();
      deallog.pop();
    }
  else
    test<3>();
}
