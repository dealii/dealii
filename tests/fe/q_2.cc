//----------------------------  q_2.cc  ---------------------------
//    q_2.cc,v 1.1 2003/05/05 13:49:41 wolf Exp
//    Version: 
//
//    Copyright (C) 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    fuqher information on this license.
//
//----------------------------  q_2.cc  ---------------------------


// Just output the embedding matrices of the FE_Q element. Test
// introduced when we started to compute them on the fly, rather than
// precomputing them for a number of elements and storing them in a
// table

#include "../tests.h"
#include <base/logstream.h>
#include <fe/fe_q.h>

#include <fstream>
#include <string>

#define PRECISION 5



template<int dim>
void
test(const unsigned int degree)
{
  deallog << "FE_Q<" << dim << "> (" << degree << ")"
	  << std::endl;
  
  FE_Q<dim> fe_q(degree);

  for (unsigned int c=0; c<GeometryInfo<dim>::children_per_cell; ++c)
    {
      const FullMatrix<double> & m = fe_q.prolongate(c);

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
  std::ofstream logfile ("q_2.output");
  logfile.precision (PRECISION);
  logfile.setf(std::ios::fixed);  
  deallog.attach(logfile);
  deallog.depth_console(0);

                                   // we had the matrices precomputed
                                   // up to Q4 for 1d, Q3 for 2d and
                                   // Q2 for 3d
  for (unsigned int degree=1; degree<=4; ++degree)
    test<1>(degree);

  for (unsigned int degree=1; degree<=3; ++degree)
    test<2>(degree);

  for (unsigned int degree=1; degree<=2; ++degree)
    test<3>(degree);
  
  return 0;
}



