//----------------------------  trilinos_sparse_matrix_iterator_10.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2004, 2005, 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  trilinos_sparse_matrix_iterator_10.cc  ---------------------------


// this test, extracted from dof_constraints_09 and trilinos_sparse_matrix_iterator_10,
// used to fail with aborts. this test is equivalent to
// trilinos_sparse_matrix_iterator_09, except that we use the postfix operator++
// instead of the prefix form

#include "../tests.h"
#include <deal.II/lac/trilinos_sparsity_pattern.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <fstream>
#include <iomanip>


void test ()
{
                                   // create a sparsity pattern with totally
                                   // empty lines (not even diagonals, since
                                   // not quadratic)
  TrilinosWrappers::SparsityPattern sparsity(4,5,1);
  sparsity.add (1,1);
  sparsity.add (3,1);
  sparsity.compress ();

                                   // attach a sparse matrix to it
  TrilinosWrappers::SparseMatrix A(sparsity);

                                   // and loop over the elements of it
  for (TrilinosWrappers::SparseMatrix::const_iterator k=A.begin();
       k!=A.end(); k++)
    deallog << k->row() << ' ' << k->column() << ' ' << k->value()
            << std::endl;
}



int main (int argc, char **argv)
{
  std::ofstream logfile("trilinos_sparse_matrix_iterator_10/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv);

  try
    {
      test ();
    }
  catch (std::exception &exc)
    {
      deallog << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      deallog << "Exception on processing: " << std::endl
		<< exc.what() << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;

      return 1;
    }
  catch (...)
    {
      deallog << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      deallog << "Unknown exception!" << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      return 1;
    };
}
