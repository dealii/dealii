//----------------------------  q_1.cc  ---------------------------
//    q_1.cc,v 1.1 2003/05/05 13:49:41 wolf Exp
//    Version: 
//
//    Copyright (C) 2003, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    fuqher information on this license.
//
//----------------------------  q_1.cc  ---------------------------


// Just output the constraint matrices of the FE_Q element. Test
// introduced when we started to compute them on the fly, rather than
// precomputing them for a number of elements and storing them in a
// table

#include "../tests.h"
#include <base/logstream.h>
#include <fe/fe_q.h>

#include <fstream>
#include <string>

#define PRECISION 2



template<int dim>
void
test(const unsigned int degree)
{
  deallog << "FE_Q<" << dim << "> (" << degree << ")"
	  << std::endl;
  
  FE_Q<dim> fe_q(degree);
  const FullMatrix<double> & constraints = fe_q.constraints();

  for (unsigned int i=0; i<constraints.m(); ++i)
    {
      for (unsigned int j=0; j<constraints.n(); ++j)
	deallog << 100*constraints(i,j) << ' ';
      deallog << std::endl;
    }
  
  deallog << std::endl;
}


int
main()
{
  std::ofstream logfile ("q_1.output");
  logfile.precision (PRECISION);
  logfile.setf(std::ios::fixed);  
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

                                   // no constraints in 1d, but we had
                                   // the matrices precomputed up to
                                   // Q4 for 2d and Q2 for 3d
  for (unsigned int degree=1; degree<=4; ++degree)
    test<2>(degree);

  for (unsigned int degree=1; degree<=2; ++degree)
    test<3>(degree);

  return 0;
}



