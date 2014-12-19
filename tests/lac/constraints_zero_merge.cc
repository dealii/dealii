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



// generate the constraints for a case where there are nodes that have
// a constraint x[i]=0, i.e. where the right hand side is a trivial
// linear combination of other degrees of freedom. then merge two such
// constraint matrices


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/constraint_matrix.h>

#include <fstream>


void test ()
{

  ConstraintMatrix cm1, cm2;

  // a "regular" and a singular
  // constraint to each of the
  // constraint matrices
  cm1.add_line (1);
  cm1.add_entry (1, 2, 42.);
  cm1.add_line (4);

  cm2.add_line (11);
  cm2.add_entry (11, 22, 42.);
  cm2.add_line (14);

  cm2.merge (cm1);

  deallog << "CM1" << std::endl;
  cm1.print (deallog.get_file_stream());

  deallog << "CM2" << std::endl;
  cm2.print (deallog.get_file_stream());
}


int main ()
{
  std::ofstream logfile("output");
  logfile.precision(2);

  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();

  deallog << "OK" << std::endl;
}
