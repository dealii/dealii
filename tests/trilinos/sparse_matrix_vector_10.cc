//----------------------------  trilinos_sparse_matrix_vector_10.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2013 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  trilinos_sparse_matrix_vector_10.cc  ---------------------------


// check SparseMatrix::vmult, vmult_add with distributed deal.II vector (but
// without using distributed things)

#include "../tests.h" 
#include <deal.II/base/utilities.h>
#include <deal.II/lac/parallel_vector.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <fstream>
#include <iostream>
#include <vector>


void test (parallel::distributed::Vector<double> &v,
           parallel::distributed::Vector<double> &w)
{
  TrilinosWrappers::SparseMatrix m(w.size(),v.size(),v.size());
  for (unsigned int i=0; i<m.m(); ++i)
    for (unsigned int j=0; j<m.n(); ++j)
        m.set (i,j, i+2*j);

  for (unsigned int i=0; i<v.size(); ++i)
    v(i) = i;
  
  m.compress (VectorOperation::insert);

                                   // w:=Mv
  m.vmult (w,v);

                                   // make sure we get the expected result
  for (unsigned int i=0; i<m.m(); ++i)
    {
      double result = 0;
      for (unsigned int j=0; j<m.n(); ++j)
        result += (i+2*j)*j;
      Assert (w(i) == result, ExcInternalError());
    }

  m.vmult_add (w, v);
                                   // make sure we get the expected result
  for (unsigned int i=0; i<m.m(); ++i)
    {
      double result = 0;
      for (unsigned int j=0; j<m.n(); ++j)
        result += (i+2*j)*j;
      Assert (w(i) == result+result, ExcInternalError());
    }

  deallog << "OK" << std::endl;
}



int main (int argc, char **argv)
{
  std::ofstream logfile("sparse_matrix_vector_10/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10); 

  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv);


  try
    {
      {
        parallel::distributed::Vector<double> v (100);
        parallel::distributed::Vector<double> w (95);
        test (v,w);
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
