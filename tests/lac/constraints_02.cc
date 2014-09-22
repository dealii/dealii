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
#include  <deal.II/lac/constraint_matrix.h>

#include <fstream>
#include <vector>

void test()
{
  ConstraintMatrix constraints;

  constraints.add_line(1);
  constraints.add_entry(1,2,1.);
  constraints.add_entry(1,3,1.);
  
  constraints.add_line(3);
  constraints.add_entry(3,4,1.);
  constraints.add_entry(3,5,1.);

  constraints.add_line(5);
  constraints.add_entry(5,0,1.);

  constraints.close();

  std::vector<types::global_dof_index> indices(4);
  indices[0] = 1;
  indices[1] = 2;
  indices[2] = 3;
  indices[3] = 5;
  
  constraints.resolve_indices(indices);

  for (unsigned int i=0; i<indices.size(); ++i)
    deallog<<"Index: "<<indices[i]<<std::endl;
}

int main()
{
  std::ofstream logfile("output");
  logfile.precision(2);

  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test();

  deallog<<"OK"<<std::endl;
}
