//----------------------------  symmetric_tensor_19.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  symmetric_tensor_19.cc  ---------------------------

// check inversion of rank-4 tensor

#include "../tests.h"
#include <base/symmetric_tensor.h>
#include <base/logstream.h>
#include <fstream>
#include <iostream>




template <int dim>
void check (const SymmetricTensor<4,dim> &A)
{
  const SymmetricTensor<4,dim> B = invert (A);

                                   // check left inverse
  deallog << "    checking left inverse" << std::endl;
  const SymmetricTensor<4,dim> T_left = B*A;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
	for (unsigned int l=0; l<dim; ++l)
	  {
	    deallog << "      "
                    << A[i][j][k][l] << ' '
                    << B[i][j][k][l] << ' '
                    << T_left[i][j][k][l] << std::endl;

	    Assert (std::fabs(T_left[i][j][k][l] -
                              identity_tensor<dim>()[i][j][k][l])
                    < 1e-10,
		    ExcInternalError());
	  }

                                   // check left inverse
  deallog << "    checking right inverse" << std::endl;
  const SymmetricTensor<4,dim> T_right = A*B;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
	for (unsigned int l=0; l<dim; ++l)
	  {
	    deallog << "      "
                    << A[i][j][k][l] << ' '
                    << B[i][j][k][l] << ' '
                    << T_right[i][j][k][l] << std::endl;

	    Assert (std::fabs(T_right[i][j][k][l] -
                              identity_tensor<dim>()[i][j][k][l])
                    < 1e-10,
		    ExcInternalError());
	  }
}

  

template <int dim>
void test ()
{
  deallog << "dim=" << dim << std::endl;

  deallog << "  unit tensor" << std::endl;
  check (identity_tensor<dim>());

                                   // do something with a more complicated
                                   // tensor. make sure it is not
                                   // rank-deficient, so choose elements by
                                   // hand
  deallog << "  complicated tensor" << std::endl;
  SymmetricTensor<4,dim> A;
  switch (dim)
    {
      case 1:
            A[0][0][0][0] = 2;
            break;
            
      case 2:
            A[0][0][0][0] = 2;
            A[0][0][1][1] = 4;
            A[0][0][0][1] = 8;
            A[1][1][0][0] = 4;
            A[1][1][1][1] = 6;
            A[1][1][0][1] = 10;
            A[0][1][0][0] = 6;
            A[0][1][1][1] = 10;
            A[0][1][0][1] = 16;
            break;

      case 3:
                                             // I'm too lazy to code something
                                             // up by hand here
            for (unsigned int i=0; i<3; ++i)
              for (unsigned int j=0; j<3; ++j)
                for (unsigned int k=0; k<3; ++k)
                  for (unsigned int l=0; l<3; ++l)
                    A[i][j][k][l] = 1.*rand()/RAND_MAX;
            break;
            
      default:
            Assert (false, ExcNotImplemented());
    }
  check (A);
}



int main ()
{
  std::ofstream logfile("symmetric_tensor_19.output");
  logfile.precision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1> ();
  test<2> ();
  test<3> ();
  
  deallog << "OK" << std::endl;
}
