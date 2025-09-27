// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check that adding time information to a .visit file works as intended


#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>

#include "../tests.h"



void
test()
{
  const unsigned int number_of_time_steps = 3;
  std::vector<std::pair<double, std::vector<std::string>>>
    times_and_piece_names(number_of_time_steps);

  times_and_piece_names[0].first = 0.0;
  times_and_piece_names[0].second.push_back("subdomain-01.time_step_0.vtk");
  times_and_piece_names[0].second.push_back("subdomain-02.time_step_0.vtk");

  times_and_piece_names[1].first = 0.5;
  times_and_piece_names[1].second.push_back("subdomain-01.time_step_1.vtk");
  times_and_piece_names[1].second.push_back("subdomain-02.time_step_1.vtk");

  times_and_piece_names[2].first = 1.0;
  times_and_piece_names[2].second.push_back("subdomain-01.time_step_2.vtk");
  times_and_piece_names[2].second.push_back("subdomain-02.time_step_2.vtk");
  DataOutBase::write_visit_record(deallog.get_file_stream(),
                                  times_and_piece_names);

  deallog << "OK" << std::endl;
}


int
main()
{
  deal_II_exceptions::disable_abort_on_exception();
  initlog();
  test();
}
