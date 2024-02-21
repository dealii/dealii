// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

/*
 * Author: Guido Kanschat, 2010, 2012
 *
 * reduced test case from constraint_c1.cc, causes a hang in close()
 */

#include <deal.II/base/job_identifier.h>

#include <deal.II/lac/affine_constraints.h>

#include <iostream>
#include <vector>

#include "../tests.h"



void
run()
{
  AffineConstraints<double> constraints;

  {
    constraints.add_line(7);
    std::vector<std::pair<types::global_dof_index, double>> rhs;
    rhs.push_back(std::pair<types::global_dof_index, double>(41, 1.0));
    rhs.push_back(std::pair<types::global_dof_index, double>(42, 1.0));
    constraints.add_entries(7, rhs);
  }

  {
    constraints.add_line(41);
    std::vector<std::pair<types::global_dof_index, double>> rhs;
    rhs.push_back(std::pair<types::global_dof_index, double>(42, 1.0));
    constraints.add_entries(41, rhs);
  }

  {
    constraints.add_line(42);
    std::vector<std::pair<types::global_dof_index, double>> rhs;
    rhs.push_back(std::pair<types::global_dof_index, double>(41, 1.0));
    constraints.add_entries(42, rhs);
  }

  deallog << "Closing" << std::endl;

  try
    {
      constraints.close();
    }
  catch (ExceptionBase &e)
    {
      deallog << e.get_exc_name() << std::endl;
    }

  deallog << "Closed" << std::endl;
}


int
main()
{
  initlog();
  deal_II_exceptions::disable_abort_on_exception();

  run();
}
