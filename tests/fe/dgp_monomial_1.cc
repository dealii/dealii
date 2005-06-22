//----------------------------  dgp_monomial_1.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003, 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    fuqher information on this license.
//
//----------------------------  dgp_monomial_1.cc  ---------------------------


// Just output the embedding matrices of the FE_DGPMonomial
// element. Test introduced when we started to compute them on the
// fly, rather than precomputing them for a number of elements and
// storing them in a table

#include "../tests.h"
#include <base/logstream.h>
#include <fe/fe_dgp_monomial.h>

#include <fstream>
#include <string>

#define PRECISION 3



template<int dim>
void
test(const unsigned int degree)
{
  deallog << "FE_DGPMonomial<" << dim << "> (" << degree << ")"
	  << std::endl;
  
  FE_DGPMonomial<dim> fe_q(degree);

  for (unsigned int c=0; c<GeometryInfo<dim>::children_per_cell; ++c)
    {
      const FullMatrix<double> & m = fe_q.get_prolongation_matrix(c);

      for (unsigned int i=0; i<m.m(); ++i)
	{
	  for (unsigned int j=0; j<m.n(); ++j)
	    deallog << m(i,j) << ' ';
	  deallog << std::endl;
	}
      
      deallog << std::endl;
    }
  
  deallog << std::endl;
}


int
main()
{
  std::ofstream logfile ("dgp_monomial_1.output");
  logfile.precision (PRECISION);
  logfile.setf(std::ios::fixed);  
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

                                   // we had the matrices precomputed
                                   // up to Q4 for 1d, Q3 for 2d and
                                   // Q2 for 3d
  for (unsigned int degree=1; degree<=4; ++degree)
    test<1>(degree);

  for (unsigned int degree=1; degree<=4; ++degree)
    test<2>(degree);

  for (unsigned int degree=1; degree<=2; ++degree)
    test<3>(degree);
  
  return 0;
}



