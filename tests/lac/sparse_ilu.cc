//----------------------------  sparse_ilu.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  sparse_ilu.cc  ---------------------------


// make sure that the SparseILU applied with infinite fill-in
// generates the exact inverse matrix

#include <cmath>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include "testmatrix.h"
#include <base/logstream.h>
#include <lac/sparse_matrix.h>
#include <lac/sparse_ilu.h>
#include <lac/vector.h>



int main()
{
  ofstream logfile("sparse_ilu.output");
  logfile.setf(ios::fixed);
  logfile.precision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  

  for (unsigned int size=4; size <= 16; size *= 2)
    {
      unsigned int dim = (size-1)*(size-1);

      deallog << "Size " << size << " Unknowns " << dim << std::endl;
      
				       // Make matrix
      FDMatrix testproblem(size, size);
      SparsityPattern structure(dim, dim, 5);
      testproblem.build_structure(structure);
      structure.compress();
      SparseMatrix<double>  A(structure);
      testproblem.laplacian(A);

      
      for (unsigned int test=0; test<2; ++test)
	{
	  deallog << "Test " << test << std::endl;
	  
					   // generate sparse ILU.
					   //
					   // for test 1, test with
					   // full pattern.  for test
					   // 2, test with same
					   // pattern as A
	  SparsityPattern ilu_pattern (dim, dim,
				       (test==0 ? dim : 5));
	  switch (test)
	    {
	      case 0:
		    for (unsigned int i=0; i<dim; ++i)
		      for (unsigned int j=0; j<dim; ++j)
			ilu_pattern.add(i,j);
		    break;

	      case 1:
		    for (unsigned int i=0; i<dim; ++i)
		      for (unsigned int j=0; j<dim; ++j)
			if (structure(i,j) != SparsityPattern::invalid_entry)
			  ilu_pattern.add(i,j);
		    break;

	      default:
		    Assert (false, ExcNotImplemented());
	    };
	  ilu_pattern.compress();
	  SparseILU<double> ilu (ilu_pattern);
	  ilu.decompose (A);
	  
					   // now for three test vectors v
					   // determine norm of
					   // (I-BA)v, where B is ILU of A.
					   // since matrix is symmetric,
					   // likewise test for right
					   // preconditioner
	  Vector<double> v(dim);
	  Vector<double> tmp1(dim), tmp2(dim);
	  for (unsigned int i=0; i<3; ++i)
	    {
	      for (unsigned int j=0; j<dim; ++j)
		v(j) = 1. * rand()/RAND_MAX;
	      
	      A.vmult (tmp1, v);
	      ilu.vmult (tmp2, tmp1);
	      tmp2 -= v;
	      const double left_residual = tmp2.l2_norm();
	      
	      ilu.vmult (tmp1, v);
	      A.vmult (tmp2, tmp1);
	      tmp2 -= v;
	      const double right_residual = tmp2.l2_norm();
	      
	      
	      deallog << "Residual with test vector " << i << ":  "
		      << " left=" << left_residual
		      << ", right=" << right_residual
		      << std::endl;
	    };
	};
      
    };
};

