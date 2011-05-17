//----------------------------  sparsity_pattern.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  sparsity_pattern.cc  ---------------------------

// check serialization for SparsityPattern

#include "serialization.h"
#include "../lac/testmatrix.h"
#include <deal.II/lac/sparsity_pattern.h>


void test ()
{
  const unsigned int N1 = 5;
  SparsityPattern sp1((N1-1)*(N1-1), (N1-1)*(N1-1), 5);
  FDMatrix(N1,N1).five_point_structure (sp1);
  sp1.compress ();
 
  const unsigned int N2 = 3;
  SparsityPattern sp2((N2-1)*(N2-1), (N2-1)*(N2-1), 5);
  FDMatrix(N2,N2).five_point_structure (sp2);
  sp2.compress ();

  SparsityPattern sp3;
  
  verify (sp1, sp2);
  
  verify (sp1, sp3);
}


int main ()
{
  std::ofstream logfile("sparsity_pattern/output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();

  deallog << "OK" << std::endl;
}
