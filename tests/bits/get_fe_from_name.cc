// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


#include "../tests.h"
#include <fstream>
#include <deal.II/base/logstream.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_tools.h>

class Test
{
public:

  void generate(const char *myname)
  {
    std::string name = myname;

    fe1 = FETools::get_fe_from_name<1>(name);
    fe2 = FETools::get_fe_from_name<2>(name);
    fe3 = FETools::get_fe_from_name<3>(name);

    deallog << "Read " << name << std::endl;
    deallog << "Generated :" << std::endl;
    deallog << fe1->get_name()  << std::endl;
    deallog << fe2->get_name()  << std::endl;
    deallog << fe3->get_name()  << std::endl;

    delete fe1;
    delete fe2;
    delete fe3;
  }

  FiniteElement<1> *fe1;
  FiniteElement<2> *fe2;
  FiniteElement<3> *fe3;

};

int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);

  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  Test gen;

  gen.generate("FE_Q_Hierarchical(1)");
  gen.generate("FE_DGPNonparametric(1)");
  gen.generate("FE_DGQ(1)");
  gen.generate("FE_Q(1)");

  gen.generate("FE_Q_Hierarchical(2)");
  gen.generate("FE_DGPNonparametric(2)");
  gen.generate("FE_DGQ(2)");
  gen.generate("FE_Q(2)");

  gen.generate("FESystem[FE_Q_Hierarchical(1)^2-FE_Q_Hierarchical(1)]");
  gen.generate("FESystem[FE_DGPNonparametric(1)^2-FE_Q_Hierarchical(1)]");
  gen.generate("FESystem[FE_DGQ(1)^2-FE_Q_Hierarchical(1)]");
  gen.generate("FESystem[FE_Q(1)^2-FE_Q_Hierarchical(1)]");

  gen.generate("FESystem[FE_Q_Hierarchical(1)^2-FE_DGPNonparametric(1)]");
  gen.generate("FESystem[FE_DGPNonparametric(1)^2-FE_DGPNonparametric(1)]");
  gen.generate("FESystem[FE_DGQ(1)^2-FE_DGPNonparametric(1)]");
  gen.generate("FESystem[FE_Q(1)^2-FE_DGPNonparametric(1)]");

  gen.generate("FESystem[FE_Q_Hierarchical(1)^2-FE_DGQ(1)]");
  gen.generate("FESystem[FE_DGPNonparametric(1)^2-FE_DGQ(1)]");
  gen.generate("FESystem[FE_DGQ(1)^2-FE_DGQ(1)]");
  gen.generate("FESystem[FE_Q(1)^2-FE_DGQ(1)]");


  gen.generate("FESystem[FE_Q_Hierarchical(1)^2-FE_Q(1)]");
  gen.generate("FESystem[FE_DGPNonparametric(1)^2-FE_Q(1)]");
  gen.generate("FESystem[FE_DGQ(1)^2-FE_Q(1)]");
  gen.generate("FESystem[FE_Q(1)^2-FE_Q(1)]");


  return 0;
}

