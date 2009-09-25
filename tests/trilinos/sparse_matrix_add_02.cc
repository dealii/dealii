//-----------------  trilinos_sparse_matrix_add_02.cc  -------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2004, 2005, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------  trilinos_sparse_matrix_add_02.cc  -------------------------


// compare collective adding of elements in a trilinos matrix using
// TrilinosWrappers::SparseMatrix::add() and a FullMatrix<double> with
// setting the same elements on an element-by-element level. Use the entries
// as they would result from a Laplace operator in 1D.

#include "../tests.h" 
#include <base/utilities.h>
#include <lac/trilinos_sparse_matrix.h>
#include <lac/full_matrix.h>
#include <fstream>
#include <iostream>


void test (TrilinosWrappers::SparseMatrix &m)
{
  TrilinosWrappers::SparseMatrix m2(m.m(), m.n(), 3);

                                   // first set a few entries one-by-one and
                                   // initialize the sparsity pattern for m2
  for (unsigned int i=0; i<m.m(); ++i)
    for (unsigned int j=0; j<m.n(); ++j)
      if (std::fabs((double)i-j) < 2)
	{
	  double value;
	  if (i == j) 
	    if (i>0 && i<m.m()-1)
	      value = 2.;
	    else 
	      value = 1.;
	  else
	    value = -1.;

	  m.set (i,j, value);
	  m2.set (i,j, 0.);
	}

  m.compress ();
  m2.compress();

				   // now add the same elements from a full
				   // matrix
  {
    FullMatrix<double> full_matrix(2,2);
    full_matrix(0,0) = full_matrix(1,1) = 1.;
    full_matrix(0,1) = full_matrix(1,0) = -1.;
    std::vector<unsigned int> local_indices (2);

    for (unsigned int i=0; i<m.m()-1; ++i)
      {
	local_indices[0] = i;
	local_indices[1] = i+1;

	m2.add (local_indices, local_indices, full_matrix);
      }
  }

  m2.compress();

				   // subtract the matrix m from this one,
				   // we should get a zero matrix
  m2.add(-1.0, m);
				   // calculate the Frobenius norm of the
				   // matrix in order to check whether all
				   // elements really are zero
  double norm = m2.frobenius_norm();
  Assert (norm == 0, ExcInternalError());

  deallog << "OK" << std::endl;
}



int main (int argc,char **argv) 
{
  std::ofstream logfile("sparse_matrix_add_02/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10); 

  Utilities::System::MPI_InitFinalize mpi_initialization (argc, argv);

  try
    {
      {
        TrilinosWrappers::SparseMatrix m (16,16,3);
        test (m);
      }
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      std::cerr << "Exception on processing: " << std::endl
		<< exc.what() << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      
      return 1;
    }
  catch (...) 
    {
      std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      std::cerr << "Unknown exception!" << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      return 1;
    };
}
