//--------------------------------------------------------------------------
//    fe_tools.cc,v 1.1 2003/11/28 15:03:26 guido Exp
//    Version: 
//
//    Copyright (C) 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//--------------------------------------------------------------------------

// Test FETools::get_fe_by_name

#include "../tests.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include <deal.II/base/logstream.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/base/quadrature.h>

template <int dim>
void test_fe(const char* name)
{
  FiniteElement<dim>* fe = FETools::get_fe_from_name<dim>(std::string(name));

  deallog << fe->get_name() << std::endl
	  << '\t' << fe->dofs_per_cell
	  << '\t' << fe->dofs_per_vertex
	  << '\t' << fe->dofs_per_line
	  << '\t' << fe->dofs_per_quad
	  << '\t' << fe->dofs_per_hex
	  << std::endl;
}


int
main()
{
  std::ofstream logfile("fe_tools_01/output");
  deallog.attach(logfile);
  deallog.depth_console(0);

				   // These names are all correct.
  test_fe<1>("FE_Q(1)");
  test_fe<1>("FE_Q(2)");
  test_fe<1>("FE_Q<1>(2)");
  test_fe<1>("FE_Q(QGaussLobatto(3))");
  test_fe<1>("FE_Q(QGaussLobatto(4))");
  test_fe<2>("FE_Q(1)");
  test_fe<2>("FE_Q<2>(2)");
  test_fe<2>("FE_Q<dim>(2)");
  test_fe<2>("FE_DGQ(1)");
  test_fe<2>("FE_RaviartThomas(1)");
  test_fe<2>("FE_Q(QGaussLobatto(3))");
  test_fe<2>("FE_Q(QGaussLobatto(4))");
  test_fe<3>("FE_Q(1)");
  test_fe<1>("FESystem<1>[FE_Q<dim>(2)^dim-FE_DGQ<d>(1)]");
  test_fe<2>("FESystem<2>[FE_Q<2>(2)^dim-FE_DGQ<2>(1)]");
  test_fe<2>("FESystem[FESystem<2>[FE_Q<2>(2)^2-FE_DGQ<2>(1)]^2-FE_Q(1)]");
  test_fe<2>("FESystem[FESystem[FESystem[FE_Q(1)^2-FE_Q(1)]^2]-FESystem[FE_Q(1)^2]-FESystem[FE_Q(1)-FE_DGP(0)]]");

				   // Now set up a list of malformed
				   // names
  std::vector<const char*> names;
//  names.push_back("FE_Q[2]");

  for(unsigned int i=0;i<names.size();++i)
    {
      try
	{
	  test_fe<2>(names[i]);
	}
      catch(ExceptionBase& e)
	{
	  logfile << e.what();
	}
    }
}
