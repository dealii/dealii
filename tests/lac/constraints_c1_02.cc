// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

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
