// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check the MultipleParameterLoop class

#include <deal.II/base/parameter_handler.h>

#include "../tests.h"


// Take the example from the
// documentation but instead of
// running anything simply print
// parameters
class HelperClass : public MultipleParameterLoop::UserClass
{
public:
  virtual void
  create_new(unsigned int run_no)
  {
    this->run_no = run_no;
  }
  virtual void
  declare_parameters(ParameterHandler &prm);
  virtual void
               run(ParameterHandler &prm);
  unsigned int run_no;
};


void
HelperClass::declare_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("Testing");
  prm.declare_entry("string list",
                    "a",
                    Patterns::List(Patterns::Selection("a|b|c|d|e|f|g|h")),
                    "docs 1");
  prm.declare_entry("int", "1", Patterns::Integer());
  prm.declare_entry("double", "3.1415926", Patterns::Double(), "docs 3");
  prm.leave_subsection();
}


void
HelperClass::run(ParameterHandler &prm)
{
  deallog << "Number of run: " << run_no << std::endl;

  prm.print_parameters(deallog.get_file_stream(), ParameterHandler::PRM);
}



void
check(const char *p)
{
  class MultipleParameterLoop prm;
  HelperClass                 h;

  h.declare_parameters(prm);
  prm.parse_input(p);
  prm.loop(h);
}


int
main()
{
  initlog();

  check(SOURCE_DIR "/prm/multiple_parameter_loop_02.prm");

  return 0;
}
