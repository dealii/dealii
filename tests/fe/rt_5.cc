//----------------------------  rt_5.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  rt_5.cc  ---------------------------


// Just output the restriction matrices of the RT element

#include <base/logstream.h>
#include <fe/fe_raviart_thomas.h>

#include <fstream>
#include <string>

#define PRECISION 5



template<int dim>
void
test(const unsigned int degree)
{
  deallog << "FE_RaviartThomas<" << dim << "> (" << degree << ")"
	  << std::endl;
  
  FE_RaviartThomas<dim> fe_rt(degree);

  for (unsigned int c=0; c<GeometryInfo<dim>::children_per_cell; ++c)
    {
      const FullMatrix<double> & m = fe_rt.restrict(c);

      for (unsigned int i=0; i<m.m(); ++i)
	{
	  for (unsigned int j=0; j<m.n(); ++j)
	    deallog << m(i,j) << ' ';
	  deallog << std::endl;
	}
      
      deallog << std::endl;
    }
}


int
main()
{
  std::ofstream logfile ("rt_5.output");
  logfile.precision (PRECISION);
  logfile.setf(std::ios::fixed);  
  deallog.attach(logfile);
  deallog.depth_console(0);

  for (unsigned int degree=0; degree<1; ++degree)
    test<2>(degree);
//  test<3>(degree);
  
  return 0;
}



