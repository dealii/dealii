// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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
      prm.declare_entry(
        "list",
        "default_1",
        Patterns::List(Patterns::Selection("default_1|b|c|d|e|f|g|h")),
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
