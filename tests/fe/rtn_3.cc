//----------------------------  rtn_3.cc  ---------------------------
//    rtn_3.cc,v 1.3 2003/06/09 16:00:38 wolf Exp
//    Version: 
//
//    Copyright (C) 2003, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  rtn_3.cc  ---------------------------

// output the number of dofs per cell for the RTN element; that seemed wrong
// in the past. compare it also with the number of dofs per cell of the RT
// element

#include "../tests.h"
#include <deal.II/fe/fe_raviart_thomas.h>

#include <vector>
#include <fstream>
#include <string>

#define PRECISION 2



template<int dim>
void test ()
{
  for (unsigned int degree=0; degree<9-2*dim; ++degree)
    {
      FE_RaviartThomasNodal<dim> fe_rtn(degree);
      deallog << fe_rtn.get_name() << ' ' << fe_rtn.dofs_per_cell << std::endl;

      if (dim != 3)
        {
          FE_RaviartThomas<dim> fe_rt(degree);
          deallog << fe_rt.get_name() << ' ' << fe_rt.dofs_per_cell << std::endl;
        }
    }
}


int
main()
{
  std::ofstream logfile ("rtn_3/output");
  deallog << std::setprecision(PRECISION);
  deallog << std::fixed;  
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<2>();
  test<3>();
  
  return 0;
}



