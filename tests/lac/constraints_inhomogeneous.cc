// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
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



// test: set inhomogeneous constraints and indirectly apply those to other
// constrained nodes.


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/constraint_matrix.h>

#include <fstream>


void test ()
{

  ConstraintMatrix cm;

  // an inhomogeneous constraint
  cm.add_line (4);
  cm.set_inhomogeneity (4, 3.14159);

  // a homogeneous constraint that is
  // constrained to the inhomogeneous one
  cm.add_line (1);
  cm.add_entry (1, 2, 42.);
  cm.add_entry (1, 4, 1.);

  // and a standard homogeneous constraint
  cm.add_line (17);
  cm.add_entry(17, 6, 2.);
  cm.add_entry(17, 15, 3.);

  // a "singular" constraint
  cm.add_line (3);

  // now close the constraint matrix
  cm.close();

  cm.print (deallog.get_file_stream());
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
