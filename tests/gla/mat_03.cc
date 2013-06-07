//---------------------------------------------------------------------------
//    $Id: ghost_01.cc 28440 2013-02-17 14:35:08Z heister $
//    Version: $Name$ 
//
//    Copyright (C) 2004, 2005, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


// document bug in assembling a LA::MPI::SparseMatrix

#include "../tests.h"
#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/base/index_set.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/sparsity_tools.h>

#include "gla.h"

template <class LA, int dim>
void test ()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);
  
  if (myid==0)
    deallog << "numproc=" << numproc << std::endl;

  ConstraintMatrix cm;
  cm.close();

  parallel::distributed::Triangulation<dim>
    triangulation (MPI_COMM_WORLD,
                     typename Triangulation<dim>::MeshSmoothing
                     (Triangulation<dim>::smoothing_on_refinement |
                      Triangulation<dim>::smoothing_on_coarsening));
  const double R0      = 6371000.-2890000;
  const double R1      = 6371000.-  35000.;
  GridGenerator::hyper_shell (triangulation,
                                  Point<dim>(),
                                  R0,
                                  R1,
                                  (dim==3) ? 96 : 12,
                                  true);

  FE_Q<dim>                                 temperature_fe(2);

  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs (temperature_fe);

  IndexSet owned = dof_handler.locally_owned_dofs();
  IndexSet relevant;
  DoFTools::extract_locally_relevant_dofs (dof_handler, relevant);

  CompressedSimpleSparsityPattern sp (owned);
  typename LA::MPI::SparseMatrix matrix;
  DoFTools::make_sparsity_pattern (dof_handler, sp,
                                       cm, false,
                                       Utilities::MPI::
                                       this_mpi_process(MPI_COMM_WORLD));
  SparsityTools::distribute_sparsity_pattern(sp,
      dof_handler.n_locally_owned_dofs_per_processor(),
      MPI_COMM_WORLD, relevant);
  sp.compress();
  matrix.reinit (owned, owned, sp, MPI_COMM_WORLD);

  FullMatrix<double> local_mat(9,9);
  for (unsigned int i=0;i<9;++i)
    for (unsigned int j=0;j<9;++j)
      local_mat(i,j)=1 + i + 100*j;
  std::vector<unsigned int> local_dofs;

  if (myid==0)
    {
      unsigned int v[]={0, 1, 2, 3, 4, 5, 6, 7, 8};
      local_dofs.assign(v,v+9);
    }
  else
    {
      unsigned int v[]={21, 39, 22, 40, 23, 41, 42, 43, 44};
      local_dofs.assign(v,v+9);
    }
  cm.distribute_local_to_global (local_mat, local_dofs, matrix);
  matrix.compress(VectorOperation::add);

				   // done
  if (myid==0)
    deallog << "OK" << std::endl;
}



int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll log(__FILE__);

  {	
    deallog.push("PETSc");
    test<LA_PETSc,2>();
    deallog.pop();	
    deallog.push("Trilinos");
    test<LA_Trilinos,2>();
    deallog.pop();	
  }
  
  // compile, don't run
  //if (myid==9999)
    //  test<LA_Dummy>();
  

}
