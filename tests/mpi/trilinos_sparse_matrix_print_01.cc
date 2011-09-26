//----------------------------  trilinos_vector_equality_4.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2004, 2005, 2008, 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  trilinos_vector_equality_4.cc  ---------------------------


// TrilinosWrappers::SparseMatrix::print got column indices wrong

#include "../tests.h"
#include <deal.II/base/utilities.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>
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
  if ((n_procs == 1) || (my_id == 1))
    sp.add(2,3);
  sp.compress();

  TrilinosWrappers::SparseMatrix A;
  A.clear ();
  A.reinit (sp);
  if ((n_procs == 1) || (my_id == 1))
    A.add(2,3, 2.0);
  A.compress();

  if ((n_procs == 1) || (my_id == 1))
    A.print(deallog.get_file_stream());
}



int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);

  const unsigned int n_procs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  deallog.push(Utilities::int_to_string(myid));

				   // let processor 1 speak if we run
				   // in parallel
  if ((n_procs == 1) || (myid == 1))
    {
      std::ofstream logfile(output_file_for_mpi("trilinos_sparse_matrix_print_01").c_str());
      deallog.attach(logfile);
      deallog << std::setprecision(4);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      test();
    }
  else
    test();

}
