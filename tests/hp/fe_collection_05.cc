// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2017 by the deal.II authors
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



// test the results of FECollection::find_least_face_dominating_fe(), namely for:
// {Q1,Q2,Q3,Q4}                with {2,3} => Q3          2
// {Q1xQ1, Q2xQ2, Q3xQ4, Q4xQ3} with {2,3} => Q2xQ2       1
// {Q1xQ1, Q3xQ4, Q4xQ3}        with {1,2} => Q1xQ1       0
// {0x0, 0x0, Q1x0, 0xQ1}       with {2,3} => none        invalid_unsigned_int
// {0x0, 0x0, Q1x0, 0xQ1}       with {2,3} => 0x0         0   (with dominating FE_Nothing)
// {Q1xQ1,Q1xQ1,Q2xQ1,Q1,Q2}    with {2,3} => Q1          0
// {Q4xQ4, Q5xQ5, Q3xQ4, Q4xQ3} with {2,3} => none        invalid_unsigned_int
// {Q1,Q2,Q4,Q3}                with {3}   => Q3          3
// {Q3,Q4,Q1,Q1}                with {2,3} => Q1          2    // self-domination


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

  // {Q1,Q2,Q3,Q4}
  {
    hp::FECollection<dim> fe_collection;
    fe_collection.push_back(FE_Q<dim>(1));
    fe_collection.push_back(FE_Q<dim>(2));
    fe_collection.push_back(FE_Q<dim>(3));
    fe_collection.push_back(FE_Q<dim>(4));
    deallog << fe_collection.find_least_face_dominating_fe(fes) << std::endl;
  }

  // {Q1xQ1, Q2xQ2, Q3xQ4, Q4xQ3}
  {
    hp::FECollection<dim> fe_collection;
    fe_collection.push_back(FESystem<dim>(FE_Q<dim>(1), 1, FE_Q<dim>(1), 1));
    fe_collection.push_back(FESystem<dim>(FE_Q<dim>(2), 1, FE_Q<dim>(2), 1));
    fe_collection.push_back(FESystem<dim>(FE_Q<dim>(3), 1, FE_Q<dim>(4), 1));
    fe_collection.push_back(FESystem<dim>(FE_Q<dim>(4), 1, FE_Q<dim>(3), 1));
    deallog << fe_collection.find_least_face_dominating_fe(fes) << std::endl;
  }

  // {Q1xQ1, Q3xQ4, Q4xQ3}
  {
    hp::FECollection<dim> fe_collection;
    fe_collection.push_back(FESystem<dim>(FE_Q<dim>(1), 1, FE_Q<dim>(1), 1));
    fe_collection.push_back(FESystem<dim>(FE_Q<dim>(3), 1, FE_Q<dim>(4), 1));
    fe_collection.push_back(FESystem<dim>(FE_Q<dim>(4), 1, FE_Q<dim>(3), 1));
    std::set<unsigned int> fes;
    fes.insert(1);
    fes.insert(2);
    deallog << fe_collection.find_least_face_dominating_fe(fes) << std::endl;
  }

  // {0x0, 0x0, Q1x0, 0xQ1}
  {
    hp::FECollection<dim> fe_collection;
    fe_collection.push_back(
      FESystem<dim>(FE_Nothing<dim>(), 1, FE_Nothing<dim>(), 1));
    fe_collection.push_back(
      FESystem<dim>(FE_Nothing<dim>(), 1, FE_Nothing<dim>(), 1));
    fe_collection.push_back(
      FESystem<dim>(FE_Q<dim>(1), 1, FE_Nothing<dim>(), 1));
    fe_collection.push_back(
      FESystem<dim>(FE_Nothing<dim>(), 1, FE_Q<dim>(1), 1));
    const unsigned int ind = fe_collection.find_least_face_dominating_fe(fes);
    if (ind == numbers::invalid_unsigned_int)
      deallog << "numbers::invalid_unsigned_int" << std::endl;
    else
      deallog << ind << std::endl;
  }

  // dominating FE_Nothing
  // {0x0, 0x0, Q1x0, 0xQ1}
  {
    hp::FECollection<dim> fe_collection;
    fe_collection.push_back(
      FESystem<dim>(FE_Nothing<dim>(1, true), 1, FE_Nothing<dim>(1, true), 1));
    fe_collection.push_back(
      FESystem<dim>(FE_Nothing<dim>(1, true), 1, FE_Nothing<dim>(1, true), 1));
    fe_collection.push_back(
      FESystem<dim>(FE_Q<dim>(1), 1, FE_Nothing<dim>(1, true), 1));
    fe_collection.push_back(
      FESystem<dim>(FE_Nothing<dim>(1, true), 1, FE_Q<dim>(1), 1));
    deallog << fe_collection.find_least_face_dominating_fe(fes) << std::endl;
  }


  // {Q1xQ1,Q1xQ1,Q2xQ1,Q1,Q2}
  {
    hp::FECollection<dim> fe_collection;
    fe_collection.push_back(FESystem<dim>(FE_Q<dim>(1), 1, FE_Q<dim>(1), 1));
    fe_collection.push_back(FESystem<dim>(FE_Q<dim>(1), 1, FE_Q<dim>(1), 1));
    fe_collection.push_back(FESystem<dim>(FE_Q<dim>(2), 1, FE_Q<dim>(1), 1));
    fe_collection.push_back(FESystem<dim>(FE_Q<dim>(1), 1, FE_Q<dim>(2), 1));
    deallog << fe_collection.find_least_face_dominating_fe(fes) << std::endl;
  }

  // {Q4xQ4, Q5xQ5, Q3xQ4, Q4xQ3
  {
    hp::FECollection<dim> fe_collection;
    fe_collection.push_back(FESystem<dim>(FE_Q<dim>(4), 1, FE_Q<dim>(4), 1));
    fe_collection.push_back(FESystem<dim>(FE_Q<dim>(5), 1, FE_Q<dim>(5), 1));
    fe_collection.push_back(FESystem<dim>(FE_Q<dim>(3), 1, FE_Q<dim>(4), 1));
    fe_collection.push_back(FESystem<dim>(FE_Q<dim>(4), 1, FE_Q<dim>(3), 1));
    const unsigned int ind = fe_collection.find_least_face_dominating_fe(fes);
    if (ind == numbers::invalid_unsigned_int)
      deallog << "numbers::invalid_unsigned_int" << std::endl;
    else
      deallog << ind << std::endl;
  }

  // {Q1,Q2,Q4,Q3}
  {
    hp::FECollection<dim> fe_collection;
    fe_collection.push_back(FE_Q<dim>(1));
    fe_collection.push_back(FE_Q<dim>(2));
    fe_collection.push_back(FE_Q<dim>(4));
    fe_collection.push_back(FE_Q<dim>(3));
    std::set<unsigned int> fes;
    fes.insert(3);
    deallog << fe_collection.find_least_face_dominating_fe(fes) << std::endl;
  }

  {
    hp::FECollection<dim> fe_collection;
    fe_collection.push_back(FE_Q<dim>(3));
    fe_collection.push_back(FE_Q<dim>(4));
    fe_collection.push_back(FE_Q<dim>(1));
    fe_collection.push_back(FE_Q<dim>(1));
    deallog << fe_collection.find_least_face_dominating_fe(fes) << std::endl;
  }
}

int
main()
{
  std::ofstream logfile("output");
  logfile.precision(2);

  deallog.attach(logfile);

  test<2>();
  test<3>();

  deallog << "OK" << std::endl;
}
