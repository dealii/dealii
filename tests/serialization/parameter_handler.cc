// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2013 by the deal.II authors
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



// check serialization for Parameter_handler

#include "serialization.h"
#include <deal.II/base/parameter_handler.h>
#include <boost/property_tree/ptree_serialization.hpp>
#include <boost/serialization/vector.hpp>


void test ()
{
  ParameterHandler prm1;
  prm1.declare_entry ("int1",
                      "1",
                      Patterns::Integer(),
                      "doc 1");
  prm1.declare_entry ("int2",
                      "2",
                      Patterns::Integer(),
                      "doc 2");
  prm1.enter_subsection ("ss1");
  {
    prm1.declare_entry ("double 1",
                        "1.234",
                        Patterns::Double(),
                        "doc 3");

    prm1.enter_subsection ("ss2");
    {
      prm1.declare_entry ("double 2",
                          "4.321",
                          Patterns::Double(),
                          "doc 4");
    }
    prm1.leave_subsection ();
  }
  prm1.leave_subsection ();

  // things with strange characters
  prm1.enter_subsection ("Testing%testing");
  {
    prm1.declare_entry ("string&list",
                        "< & > ; /",
                        Patterns::Anything(),
                        "docs 1");
    prm1.declare_entry ("int*int",
                        "2",
                        Patterns::Integer());
    prm1.declare_entry ("double+double",
                        "6.1415926",
                        Patterns::Double(),
                        "docs 3");
  }
  prm1.leave_subsection ();

  ParameterHandler prm2;
  prm2.enter_subsection ("s234s1");
  {
    prm2.declare_entry ("dummy 1",
                        "19.4",
                        Patterns::Double(),
                        "paper 4");

    prm2.enter_subsection ("s2skjds2");
    {
      prm2.declare_entry ("var 2",
                          "6.321",
                          Patterns::Double(),
                          "invalid");
    }
    prm2.leave_subsection ();
  }
  prm2.leave_subsection ();

  // things with strange characters
  prm2.enter_subsection ("Testing%testing");
  {
    prm2.declare_entry ("int*int",
                        "2",
                        Patterns::Integer());
    prm2.declare_entry ("string&list",
                        "< & > ; /",
                        Patterns::Anything(),
                        "docs 1");
  }
  prm2.leave_subsection ();

  prm2.declare_entry ("int1",
                      "1.",
                      Patterns::Double(),
                      "doc 1");
  prm2.declare_entry ("int2",
                      "2",
                      Patterns::Anything(),
                      "doc 2");

  ParameterHandler prm3;

  verify (prm1, prm2);

  verify (prm1, prm3);
}


int main ()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();

  deallog << "OK" << std::endl;
}


