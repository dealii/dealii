//--------------------------------------------------------------------
//    $Id$
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

// Tests eigenvalues of FullMatrix

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <float.h>
#include <fstream>
#include <iostream>
#include <vector>

const double left[] =
{
  
   1.75, -0.433012701892219, 0.0, 0.0,
   -0.433012701892219, 1.25, 0.0, 0.0,
   0.0, 0.0, 3.5, -0.5,
   0.0, 0.0, -0.5, 3.5
};



int main()
{
  const std::string logname = JobIdentifier::base_name(__FILE__) + std::string("/output");
  std::ofstream logfile(logname.c_str());
  logfile.precision(1);
  deallog.attach(logfile);
  deallog.depth_console(0);
  
  FullMatrix<double> A(4,4,left);
  LAPACKFullMatrix<double> LA(4,4);                   
  LA = A;
  FullMatrix<double> eigenvectors;
  Vector<double> eigenvalues(0);
  
  LA.compute_eigenvalues_symmetric (0.5, 2.5,
				    2.0*DBL_MIN,
				    eigenvalues,
				    eigenvectors);
  
  for (unsigned int i=0;i<eigenvalues.size();++i)
    {
      deallog << "eigenvalue "
	      << std::scientific << eigenvalues(i) << std::endl
	      << "eigenvector ";
      for (unsigned int j=0;j<A.m();++j)
	{
	  deallog << std::scientific 
		  << eigenvectors(j,i)/eigenvectors(0,i)
		  << '\t';
	}
      deallog << std::endl;
    }
}
