// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test internal typetraits used in FEEvaluation

#include <deal.II/matrix_free/fe_evaluation.h>

#include "../tests.h"

// dummy class we use to check typetraits and internal function.
// this one mimics LA::d::Vec
template <typename Number>
class Dummy
{
public:
  using value_type = Number;

  Number
  local_element(const unsigned int) const
  {
    deallog << "Dummy::local_element() const" << std::endl;
    return Number();
  }

  void
  set_local_element(const unsigned int, const Number)
  {
    deallog << "Dummy::set_local_element()" << std::endl;
  }

  void
  add_local_element(const unsigned int, const Number)
  {
    deallog << "Dummy::add_local_element()" << std::endl;
  }

  Number
  operator()(const unsigned int) const
  {
    deallog << "Dummy::operator() const" << std::endl;
    return Number();
  }
};


template <typename Number>
class Dummy2
{
public:
  using value_type = Number;

  Number
  local_element(const unsigned int) const
  {
    deallog << "Dummy2::local_element() const" << std::endl;
    return Number();
  }

  Number &
  local_element(const unsigned int)
  {
    deallog << "Dummy2::local_element()" << std::endl;
    return dummy;
  }

  Number
  operator()(const unsigned int) const
  {
    deallog << "Dummy2::operator() const" << std::endl;
    return Number();
  }

  Number &
  operator()(const unsigned int)
  {
    deallog << "Dummy2::operator()" << std::endl;
    return dummy;
  }

private:
  Number dummy;
};


int
main()
{
  initlog();

  Dummy<double>  dummy;
  Dummy2<double> dummy2;

  deallog << "has_add_local_element:" << std::endl
          << "Dummy = "
          << internal::has_add_local_element<Dummy<double>> << std::endl
          << "Dummy2 = "
          << internal::has_add_local_element<Dummy2<double>> << std::endl
          << "has_set_local_element:" << std::endl
          << "Dummy = "
          << internal::has_set_local_element<Dummy<double>> << std::endl
          << "Dummy2 = "
          << internal::has_set_local_element<Dummy2<double>> << std::endl;

  // now check internal::vector_access_[set|add] wrappers
  deallog << "internal::vector_access_set:" << std::endl;
  internal::vector_access_set(dummy, 0, 0);
  internal::vector_access_set(dummy2, 0, 0);

  deallog << "internal::vector_access_add:" << std::endl;
  internal::vector_access_add(dummy, 0, 0);
  internal::vector_access_add(dummy2, 0, 0);

  deallog << "OK" << std::endl;
}
