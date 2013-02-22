//----------------------------  petsc_sparse_matrix_vector_07.cc  ---------------------------
//    $Id: sparse_matrix_vector_07.cc 24924 2012-01-25 12:35:17Z kormann $
//    Version: $Name$ 
//
//    Copyright (C) 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  petsc_sparse_matrix_vector_07.cc  ---------------------------


// check SparseMatrix::add(other, factor)

#include "../tests.h"
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <fstream>
#include <iostream>
#include <vector>


void test ()
{
  const unsigned int s = 10;
  PETScWrappers::SparseMatrix m(s,s,s);
  for (unsigned int i=0; i<m.m(); ++i)
    for (unsigned int j=0; j<=i; ++j)
        m.set (i,j, i+2*j);
      
  m.compress ();

  PETScWrappers::SparseMatrix m2(s,s,s);
  m2.set(0,1,5.0);
  for (unsigned int i=0; i<m2.m(); ++i)
    m2.set(i,i,1.0+ i);
  m2.compress();

  // we now print the matrix m,
  // print after adding m2, and then subtract m2 again
  // to get the original matrix back.
  
  deallog << "before: " << m(0,1) << std::endl;
  for (unsigned int i=0;i<s;++i)
    deallog << m(i,i) << " ";
  deallog << std::endl;

  m.add(m2,1.0);

  deallog << "after: " << m(0,1) << std::endl;
  for (unsigned int i=0;i<s;++i)
    deallog << m(i,i) << " ";
  deallog << std::endl;

  m.add(m2,-1.0);

  deallog << "back to original: " << m(0,1) << std::endl;
  for (unsigned int i=0;i<s;++i)
    deallog << m(i,i) << " ";
  deallog << std::endl;
  
  deallog << "OK" << std::endl;
}



int main (int argc, char **argv)
{
  std::ofstream logfile("sparse_matrix_01/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  PetscInitialize(&argc,&argv,0,0);
  test ();
  PetscFinalize ();
}
