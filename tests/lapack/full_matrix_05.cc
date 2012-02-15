//--------------------------------------------------------------------
//    $Id: full_matrix_05.cc 23710 2011-05-17 04:50:10Z bangerth $
//    Version: $Name$ 
//
//    Copyright (C) 2005, 2006, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//--------------------------------------------------------------------

// Tests generalized eigenvalues of FullMatrix

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <fstream>
#include <iostream>
#include <vector>

const double left[] =
{
      4., -1., -1., -1.,
      -1., 4., -1., -1.,
      -1., -1., 4., -1.,
      -1., -1., -1., 4.
};

const double right[] =
{
      4., -1., -1., -1.,
      -1., 5., -1., -1.,
      -1., -1., 6., -1.,
      -1., -1., -1., 7.
};



int main()
{
  const std::string logname = JobIdentifier::base_name(__FILE__) + std::string("/output");
  std::ofstream logfile(logname.c_str());
  logfile.precision(1);
  deallog.attach(logfile);
  deallog.depth_console(0);
  
  FullMatrix<double> A(4,4,left),
                     B(4,4,right);
  LAPACKFullMatrix<double> LA(4,4),
                           LB(4,4);
  for(unsigned int itype=1;itype<=3; ++itype)
    {                     
      deallog << std::endl 
              << "generalized eigenvalue problem of type " 
	      << itype << std::endl;
      LA = A;
      LB = B;
      std::vector<Vector<double> > eigenvectors(A.m());
      LA.compute_generalized_eigenvalues_symmetric (LB, eigenvectors, itype);
      
      for (unsigned int i=0;i<A.m();++i)
	{
	  std::complex<double> lambda = LA.eigenvalue(i);
	  deallog << "generalized eigenvalue "
		  << std::scientific << lambda.real() << '\t'
		  << std::scientific << lambda.imag() << std::endl
		  << "generalized eigenvector ";
	  for (unsigned int j=0;j<A.m();++j)
	    {
	      deallog << std::scientific 
	              << eigenvectors[i](j)/eigenvectors[i](0) 
		      << '\t';
	    }
	  deallog << std::endl;
	}
    }
}
