// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


/* clang-format off */
// find the most dominating FE from a set of FEs on faces (codim=1).
// for this task we call FECollection::find_dominated_fe_extended(), which
// concatenates the two functions FECollection::find_enclosing_fes() and
// FECollection::find_dominating_fe().
// we test the results for the following collections:
//   {Q1, Q2, Q3, Q4}             with {2,3} => Q4          3
//   {Q5xQ5, Q4xQ4, Q3xQ4, Q4xQ3} with {2,3} => Q4xQ4       1
//   {Q5xQ5, Q3xQ4, Q4xQ3}        with {2,3} => Q5xQ5       0
//   {Q1x0, 0xQ1, 0x0, 0x0}       with {2,3} => none        invalid_fe_index
//   {Q1x0, 0xQ1, 0x0, 0x0}       with {2,3} => 0x0         2   (with dominating FE_Nothing)
//   {Q2xQ2, Q2xQ2, Q2xQ1, Q1xQ2} with {2,3} => Q2xQ2       0
//   {Q2xQ2, Q3xQ3, Q3xQ4, Q4xQ3} with {2,3} => none        invalid_fe_index
//   {Q1, Q2, Q4, Q3}             with {3}   => Q3          3
//   {Q3, Q4, Q1, Q1}             with {2,3} => Q1          2   (self-domination)
/* clang-format on */


#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/hp/fe_collection.h>

#include "../tests.h"


template <int dim>
void
test()
{
  std::set<unsigned int> fes;
  fes.insert(2);
  fes.insert(3);

  // {Q1, Q2, Q3, Q4}
  {
    hp::FECollection<dim> fe_collection;
    fe_collection.push_back(FE_Q<dim>(1));
    fe_collection.push_back(FE_Q<dim>(2));
    fe_collection.push_back(FE_Q<dim>(3));
    fe_collection.push_back(FE_Q<dim>(4));
    deallog << fe_collection.find_dominated_fe_extended(fes, /*codim=*/1)
            << std::endl;
  }

  // {Q5xQ5, Q4xQ4, Q3xQ4, Q4xQ3}
  {
    hp::FECollection<dim> fe_collection;
    fe_collection.push_back(FESystem<dim>(FE_Q<dim>(5), 2));
    fe_collection.push_back(FESystem<dim>(FE_Q<dim>(4), 2));
    fe_collection.push_back(FESystem<dim>(FE_Q<dim>(3), 1, FE_Q<dim>(4), 1));
    fe_collection.push_back(FESystem<dim>(FE_Q<dim>(4), 1, FE_Q<dim>(3), 1));
    deallog << fe_collection.find_dominated_fe_extended(fes, /*codim=*/1)
            << std::endl;
  }

  // {Q5xQ5, Q3xQ4, Q4xQ3}
  {
    hp::FECollection<dim> fe_collection;
    fe_collection.push_back(FESystem<dim>(FE_Q<dim>(5), 2));
    fe_collection.push_back(FESystem<dim>(FE_Q<dim>(3), 1, FE_Q<dim>(4), 1));
    fe_collection.push_back(FESystem<dim>(FE_Q<dim>(4), 1, FE_Q<dim>(3), 1));
    std::set<unsigned int> fes;
    fes.insert(1);
    fes.insert(2);
    deallog << fe_collection.find_dominated_fe_extended(fes, /*codim=*/1)
            << std::endl;
  }

  // {Q1x0, 0xQ1, 0x0, 0x0}
  {
    hp::FECollection<dim> fe_collection;
    fe_collection.push_back(
      FESystem<dim>(FE_Q<dim>(1), 1, FE_Nothing<dim>(), 1));
    fe_collection.push_back(
      FESystem<dim>(FE_Nothing<dim>(), 1, FE_Q<dim>(1), 1));
    fe_collection.push_back(FESystem<dim>(FE_Nothing<dim>(), 2));
    fe_collection.push_back(FESystem<dim>(FE_Nothing<dim>(), 2));
    const unsigned int ind =
      fe_collection.find_dominated_fe_extended(fes, /*codim=*/1);
    if (ind == numbers::invalid_fe_index)
      deallog << "numbers::invalid_fe_index" << std::endl;
    else
      deallog << ind << std::endl;
  }

  // dominating FE_Nothing
  // {Q1x0, 0xQ1, 0x0, 0x0}
  {
    hp::FECollection<dim> fe_collection;
    fe_collection.push_back(
      FESystem<dim>(FE_Q<dim>(1), 1, FE_Nothing<dim>(1, true), 1));
    fe_collection.push_back(
      FESystem<dim>(FE_Nothing<dim>(1, true), 1, FE_Q<dim>(1), 1));
    fe_collection.push_back(FESystem<dim>(FE_Nothing<dim>(1, true), 2));
    fe_collection.push_back(FESystem<dim>(FE_Nothing<dim>(1, true), 2));
    deallog << fe_collection.find_dominated_fe_extended(fes, /*codim=*/1)
            << std::endl;
  }


  // {Q2xQ2, Q2xQ2, Q2xQ1, Q1xQ2}
  {
    hp::FECollection<dim> fe_collection;
    fe_collection.push_back(FESystem<dim>(FE_Q<dim>(2), 2));
    fe_collection.push_back(FESystem<dim>(FE_Q<dim>(2), 2));
    fe_collection.push_back(FESystem<dim>(FE_Q<dim>(2), 1, FE_Q<dim>(1), 1));
    fe_collection.push_back(FESystem<dim>(FE_Q<dim>(1), 1, FE_Q<dim>(2), 1));
    deallog << fe_collection.find_dominated_fe_extended(fes, /*codim=*/1)
            << std::endl;
  }

  // {Q2xQ2, Q3xQ3, Q3xQ4, Q4xQ3}
  {
    hp::FECollection<dim> fe_collection;
    fe_collection.push_back(FESystem<dim>(FE_Q<dim>(2), 2));
    fe_collection.push_back(FESystem<dim>(FE_Q<dim>(3), 2));
    fe_collection.push_back(FESystem<dim>(FE_Q<dim>(3), 1, FE_Q<dim>(4), 1));
    fe_collection.push_back(FESystem<dim>(FE_Q<dim>(4), 1, FE_Q<dim>(3), 1));
    const unsigned int ind =
      fe_collection.find_dominated_fe_extended(fes, /*codim=*/1);
    if (ind == numbers::invalid_fe_index)
      deallog << "numbers::invalid_fe_index" << std::endl;
    else
      deallog << ind << std::endl;
  }

  // {Q1, Q2, Q4, Q3}
  {
    hp::FECollection<dim> fe_collection;
    fe_collection.push_back(FE_Q<dim>(1));
    fe_collection.push_back(FE_Q<dim>(2));
    fe_collection.push_back(FE_Q<dim>(4));
    fe_collection.push_back(FE_Q<dim>(3));
    std::set<unsigned int> fes;
    fes.insert(3);
    deallog << fe_collection.find_dominated_fe_extended(fes, /*codim=*/1)
            << std::endl;
  }

  // {Q3, Q4, Q1, Q1}
  {
    hp::FECollection<dim> fe_collection;
    fe_collection.push_back(FE_Q<dim>(3));
    fe_collection.push_back(FE_Q<dim>(4));
    fe_collection.push_back(FE_Q<dim>(1));
    fe_collection.push_back(FE_Q<dim>(1));
    deallog << fe_collection.find_dominated_fe_extended(fes, /*codim=*/1)
            << std::endl;
  }
}

int
main()
{
  initlog();
  deallog.get_file_stream().precision(2);

  deallog.push("2D");
  test<2>();
  deallog.pop();
  deallog.push("3D");
  test<3>();
  deallog.pop();

  deallog << "OK" << std::endl;
}
