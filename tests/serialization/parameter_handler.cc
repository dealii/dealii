//----------------------------  parameter_handler.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2002, 2003, 2004, 2005, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  parameter_handler.cc  ---------------------------


// check serialization for Parameter_handler

#include "serialization.h"
#include <base/parameter_handler.h>
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
  std::ofstream logfile("parameter_handler/output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();

  deallog << "OK" << std::endl;
}


