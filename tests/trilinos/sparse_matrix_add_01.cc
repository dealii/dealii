//-----------------  trilinos_sparse_matrix_add_01.cc  -------------------------
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
//-----------------  trilinos_sparse_matrix_add_01.cc  -------------------------


// compare collective adding of elements in a trilinos matrix using
// TrilinosWrappers::SparseMatrix::add() with element-wise setting

#include "../tests.h" 
#include <deal.II/base/utilities.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <fstream>
#include <iostream>


void test (TrilinosWrappers::SparseMatrix &m)
{
  TrilinosWrappers::SparseMatrix m2(m.m(), m.n(), m.m()/3+1);

                                   // first set a few entries one-by-one
  for (unsigned int i=0; i<m.m(); ++i)
    for (unsigned int j=0; j<m.n(); ++j)
      if ((i+2*j+1) % 3 == 0)
	{
	  m.set (i,j, i*j*.5+.5);
	  m2.set (i,j, 0.);
	}

  m.compress ();
  m2.compress();

				   // now add the same elements row-wise
  {
    std::vector<unsigned int> col_indices (m.n()/3+1);
    std::vector<double> col_values (m.n()/3+1);
    for (unsigned int i=0; i<m.m(); ++i)
      {
	unsigned int col_index = 0;
				   // count the number of elements in this
				   // row
	for (unsigned int j=0; j<m.n(); ++j)
	  if ((i+2*j+1) % 3 == 0)
	    ++col_index;

	col_indices.resize(col_index);
	col_values.resize(col_index);
	col_index = 0;

				   // extract column values
	for (unsigned int j=0; j<m.n(); ++j)
	  if ((i+2*j+1) % 3 == 0)
	    {
	      col_indices[col_index] = j;
	      col_values[col_index] =  i*j*.5+.5;
	      col_index++;
	    }

	m2.add (i, col_indices, col_values);
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
  std::ofstream logfile("sparse_matrix_add_01/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10); 

  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv);

  try
    {
      {
        TrilinosWrappers::SparseMatrix m (5U,5U,3U);
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
