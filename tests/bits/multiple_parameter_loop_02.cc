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



// check the MultipleParameterLoop class

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/parameter_handler.h>
#include <fstream>


// Take the example from the
// documentation but instead of
// running anything simply print
// parameters
class HelperClass : public MultipleParameterLoop::UserClass
{
public:
  virtual void create_new (unsigned int run_no)
  {
    this->run_no = run_no;
  }
  virtual void declare_parameters (ParameterHandler &prm);
  virtual void run (ParameterHandler &prm);
  unsigned int run_no;
};


void HelperClass::declare_parameters (ParameterHandler &prm)
{
  prm.enter_subsection ("Testing");
  prm.declare_entry ("string list",
                     "a",
                     Patterns::List(Patterns::Selection("a|b|c|d|e|f|g|h")),
                     "docs 1");
  prm.declare_entry ("int",
                     "1",
                     Patterns::Integer());
  prm.declare_entry ("double",
                     "3.1415926",
                     Patterns::Double(),
                     "docs 3");
  prm.leave_subsection ();
}


void HelperClass::run (ParameterHandler &prm)
{
  deallog << "Number of run: " << run_no << std::endl;

  prm.print_parameters (deallog.get_file_stream(), ParameterHandler::Text);
}



void check (const char *p)
{
  class MultipleParameterLoop prm;
  HelperClass h;

  h.declare_parameters (prm);
  prm.read_input (p);
  prm.loop (h);
}


int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check (SOURCE_DIR "/prm/multiple_parameter_loop_02.prm");

  return 0;
}
