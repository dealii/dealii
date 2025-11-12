// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2025 - by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// This program tests that we ignore a trailing comma in function constants
// and a trailing semicolon in function expressions.

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/base/point.h>
#include <deal.II/base/utilities.h>

#include "../tests.h"

template <int dim>
void
Test()
{
  // A parameter handler
  ParameterHandler prm;

  // Test vector declaration
  for (unsigned int i = 0; i < dim; ++i)
    {
      const std::string id = "Function " + Utilities::int_to_string(dim) +
                             " - " + Utilities::int_to_string(i);
      prm.enter_subsection(id);

      Functions::ParsedFunction<dim>::declare_parameters(prm, i + 1);
      prm.set("Function constants",
              "f=" + Utilities::int_to_string(i + 1) + ",");

      std::string expr;

      for (unsigned int j = 0; j < i + 1; ++j)
        {
          expr += "f + 0.5;";
        }

      prm.set("Function expression", expr);

      Functions::ParsedFunction<dim> function(i + 1);
      function.parse_parameters(prm);

      prm.leave_subsection();

      deallog << "Dim: " << dim << ". Components: " << i
              << ". Value: " << function.value(Point<dim>()) << std::endl;
    }

  deallog << "Tested on: " << std::endl;
  prm.log_parameters(deallog);
}

int
main()
{
  initlog();

  Test<1>();
  Test<2>();
  Test<3>();
}
