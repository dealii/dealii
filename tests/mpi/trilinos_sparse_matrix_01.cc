//----------------------  trilinos_matvec_01.cc  ---------------------------
//    $Id: trilinos_matvec_01.cc 28697 2013-03-01 16:48:48Z heister $
//    Version: $Name$
//
//    Copyright (C) 2004, 2005, 2008, 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------  trilinos_matvec_01.cc  ---------------------------


// Test whether TrilinosWrappers::SparseMatrix::vmult gives same result with
// Trilinos vector and parallel distributed vector

#include "../tests.h"
#include <deal.II/base/utilities.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/parallel_vector.h>
#include <fstream>
#include <iostream>
#include <vector>


void test ()
{
  const unsigned int n_procs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  const unsigned int my_id = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  const unsigned int n_rows = 3;
  const unsigned int n_cols = 4;

  IndexSet row_partitioning (n_rows);
  IndexSet col_partitioning (n_cols);

  if (n_procs == 1)
    {
      row_partitioning.add_range(0, n_rows);
      col_partitioning.add_range(0, n_cols);
    }
  else if (n_procs == 2)
    {
				       // row_partitioning should be { [0, 2), [2, n_rows) }
				       // col_partitioning should be { [0, 2), [2, n_cols) }
				       // col_relevant_set should be { [0, 3), [1, n_cols) }
      if (my_id == 0)
	{
	  row_partitioning.add_range(0, 2);
	  col_partitioning.add_range(0, 2);
	}
      else if (my_id == 1)
	{
	  row_partitioning.add_range(2, n_rows);
	  col_partitioning.add_range(2, n_cols);
	}
    }
  else
    Assert (false, ExcNotImplemented());

  TrilinosWrappers::SparsityPattern sp (row_partitioning,
					col_partitioning, MPI_COMM_WORLD);
  if (my_id == 0)
    {
      sp.add (0, 0);
      sp.add (0, 2);
    }
  if ((n_procs == 1) || (my_id == 1))
    sp.add(2,3);
  sp.compress();

  TrilinosWrappers::SparseMatrix A;
  A.clear ();
  A.reinit (sp);
  if (my_id==0)
    {
      A.set(0, 0, 0.1);
      A.set(0, 2, 0.2);
    }
  if ((n_procs == 1) || (my_id == 1))
    A.set(2,3, 0.3);

  A.compress(VectorOperation::insert);
  
  if (my_id == 0)
    {
      //A.add (0, 0, 1.33);
      //A.add (0, 2, 1.0);
    }
  if ((n_procs == 1) || (my_id == 1))
    A.add(0,0, 1.67);
    //A.begin()->value() += 1.67; // crashes
  
  A.compress(VectorOperation::add);

  if (my_id==0)
    {
      deallog << "A(0,0)=" << A(0,0) << std::endl;
    }
  
 
    if (my_id == 0) deallog << "OK" << std::endl;
}



int main (int argc, char **argv)
{  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);

  const unsigned int n_procs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile(output_file_for_mpi("trilinos_sparse_matrix_01").c_str());
      deallog.attach(logfile);
      deallog << std::setprecision(4);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      test();
    }
  else
    test();

}
