//----------------------------  petsc_sparse_matrix_vector_03.cc  ---------------------------
//    vector_11.cc,v 1.4 2003/07/03 10:31:46 guido Exp
//    Version: 
//
//    Copyright (C) 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  petsc_sparse_matrix_vector_03.cc  ---------------------------


// check SparseMatrix::vmult_add

#include "../tests.h"
#include <lac/petsc_vector.h>
#include <lac/petsc_sparse_matrix.h>
#include <fstream>
#include <iostream>
#include <vector>


void test (petsc_wrappers::Vector &v,
           petsc_wrappers::Vector &w)
{
  petsc_wrappers::SparseMatrix m(v.size(),v.size(),v.size());
  for (unsigned int i=0; i<m.m(); ++i)
    for (unsigned int j=0; j<m.m(); ++j)
        m.set (i,j, i+2*j);

  for (unsigned int i=0; i<v.size(); ++i)
    {
      v(i) = i;
      w(i) = i;
    }
  
  m.compress ();
  v.compress ();
  w.compress ();

                                   // w:=Mv
  m.vmult_add (w,v);

                                   // make sure we get the expected result
  for (unsigned int i=0; i<v.size(); ++i)
    {
      Assert (v(i) == i, ExcInternalError());

      double result = 0;
      for (unsigned int j=0; j<m.m(); ++j)
        result += (i+2*j)*j;
      Assert (w(i) == i+result, ExcInternalError());
    }

  deallog << "OK" << std::endl;
}



int main (int argc, char **argv)
{
  std::ofstream logfile("petsc_sparse_matrix_vector_03.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  try
    {
      PetscInitialize(&argc,&argv,0,0);
      {
        petsc_wrappers::Vector v (100);
        petsc_wrappers::Vector w (100);
        test (v,w);
      }
      PetscFinalize ();
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
