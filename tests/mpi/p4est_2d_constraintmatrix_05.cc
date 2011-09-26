//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


// check that ConstraintMatrix.distribute() is not doing anything in a
// distributed computation for a vector that already has the entries set
// correctly

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/numerics/vectors.h>

#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/fe/fe_q.h>


template<int dim>
void test()
{
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tr);

  tr.refine_global (1);
  tr.begin_active()->set_refine_flag ();
  tr.execute_coarsening_and_refinement ();

  DoFHandler<dim> dofh(tr);

  static FE_Q<dim> fe(1);

  dofh.distribute_dofs (fe);

  IndexSet owned_set = dofh.locally_owned_dofs();

  IndexSet relevant_set;
  DoFTools::extract_locally_relevant_dofs (dofh, relevant_set);

  TrilinosWrappers::MPI::Vector x_ref;
  x_ref.reinit(owned_set, MPI_COMM_WORLD);
  VectorTools::interpolate(dofh,
			   ConstantFunction<dim> (1.),
			   x_ref);
  x_ref.compress();

  TrilinosWrappers::MPI::Vector x1 (x_ref);

				// we have interpolated values, so
				// ConstraintMatrix::distribute should not do
				// anything
  x1 = x_ref;
  ConstraintMatrix cm(relevant_set);
  DoFTools::make_hanging_node_constraints (dofh, cm);
  cm.close ();
  cm.distribute(x1);

  x1 -= x_ref;
  double err = x1.linfty_norm();
  if (err>1.0e-12)
    if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
      deallog << "err:" << err << std::endl;

				// now test the same thing with a fresh vector
				// that we manually fill with ones, not by a
				// function in interpolate
  TrilinosWrappers::MPI::Vector x2 (owned_set, MPI_COMM_WORLD);
  x2 = 1;
  cm.distribute(x2);
  x2 -= x_ref;
  err = x2.linfty_norm();
  if (err>1.0e-12)
    if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
      deallog << "err:" << err << std::endl;

  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
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

  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile(output_file_for_mpi("p4est_2d_constraintmatrix_05").c_str());
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
