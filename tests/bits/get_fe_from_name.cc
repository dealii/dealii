//----------------------------  get_fe_from_name.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005 by the deal.II authors 
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  get_fe_from_name.cc  ---------------------------
// Test the get_fe_from_name function of the fe_tools. 

#include "../tests.h"
#include <fstream>
#include <base/logstream.h>
#include <fe/fe.h>
#include <fe/fe_tools.h>

class Test {
public:
  
  void generate(const char * myname) {
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
  
  FiniteElement<1> * fe1;
  FiniteElement<2> * fe2;
  FiniteElement<3> * fe3;
 
};

int main () 
{  
  std::ofstream logfile("get_fe_from_name.output");
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
    
