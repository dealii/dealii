// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test escaping of underscores in ParameterHandler::print_parameters(LaTeX)

#include <deal.II/base/parameter_handler.h>

#include "../tests.h"


int
main()
{
  initlog();

  // latex header to make it easy to check the .tex output:
  deallog.get_file_stream()
    << "\\documentclass{article}\n"
       "\\usepackage{amsmath}\n"
       "\\usepackage{imakeidx}\n"
       "\\makeindex[name=prmindex, title=Index of run-time parameter entries]\n"
       "\\makeindex[name=prmindexfull, title=Index of run-time parameters with section names]\n"
       "\\usepackage[colorlinks,linkcolor=blue,urlcolor=blue,citecolor=blue,baseurl=../]{hyperref}\n"
       "\\begin{document}\n"
    << std::endl;

  ParameterHandler prm;
  prm.enter_subsection("section_1");
  {
    prm.enter_subsection("section_2&3");
    {
      prm.declare_entry("list",
                        "default_1",
                        Patterns::List(
                          Patterns::Selection("default_1|b|c|d|e|f|g|h")),
                        "docs 1");
      prm.declare_entry("int_hello", "1", Patterns::Integer());
      prm.declare_entry("text",
                        "default text with _ ~ \\ # & { }",
                        Patterns::Anything(),
                        "documentation with formulas: $x_3$");
    }
    prm.leave_subsection();
    prm.declare_entry("some other entry with {%}",
                      "3.1415926",
                      Patterns::Double(),
                      "documentation");
  }
  prm.leave_subsection();

  prm.print_parameters(deallog.get_file_stream(), ParameterHandler::LaTeX);

  deallog.get_file_stream() << "\n\n"
                               "\\printindex[prmindex]\n"
                               "\\printindex[prmindexfull]\n"
                               "\\end{document}"
                            << std::endl;
}
