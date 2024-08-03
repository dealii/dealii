// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2024 by the deal.II authors
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

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/trilinos_vector.h>

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

  Number &
  local_element(const unsigned int)
  {
    deallog << "Dummy::local_element()" << std::endl;
    return dummy;
  }

  Number
  operator()(const unsigned int) const
  {
    deallog << "Dummy::operator() const" << std::endl;
    return Number();
  }

  Number &
  operator()(const unsigned int)
  {
    deallog << "Dummy::operator()" << std::endl;
    return dummy;
  }

private:
  Number dummy;
};


template <typename Number>
class Dummy2
{
public:
  using value_type = Number;

  Number
  operator[](const unsigned int) const
  {
    deallog << "Dummy2::operator() const" << std::endl;
    return Number();
  }

  Number &
  operator[](const unsigned int)
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

  deallog
    << "has_local_element:" << std::endl
    << "LinearAlgebra::distributed::Vector = "
    << internal::has_local_element<
         LinearAlgebra::distributed::Vector<double>> << std::endl
    << "TrilinosWrappers::MPI::Vector = "
    << internal::has_local_element<TrilinosWrappers::MPI::Vector> << std::endl
    << "Dummy = " << internal::has_local_element<Dummy<double>> << std::endl
    << "Dummy2 = " << internal::has_local_element<Dummy2<double>> << std::endl
    << "Vector = " << internal::has_local_element<Vector<double>> << std::endl;

  // now check internal::vector_access wrapper
  // we expect local_element() to be called
  deallog << "internal::vector_access:" << std::endl;
  internal::vector_access(dummy, 0);
  // internal::vector_access(dummy2, 0);


  // now check has_partitioners_are_compatible:
  deallog
    << "has_partitioners_are_compatible:" << std::endl
    << "LinearAlgebra::distributed::Vector = "
    << internal::has_partitioners_are_compatible<
         LinearAlgebra::distributed::Vector<double>> << std::endl
    << "TrilinosWrappers::MPI::Vector = "
    << internal::has_partitioners_are_compatible<
         TrilinosWrappers::MPI::Vector> << std::endl
    << "Vector = "
    << internal::has_partitioners_are_compatible<Vector<double>> << std::endl;

  // check has_begin:
  deallog << "has_begin:" << std::endl
          << "LinearAlgebra::distributed::Vector = "
          << internal::has_begin<
               LinearAlgebra::distributed::Vector<double>> << std::endl
          << "TrilinosWrappers::MPI::Vector = "
          << internal::has_begin<TrilinosWrappers::MPI::Vector> << std::endl
          << "Vector = " << internal::has_begin<Vector<double>> << std::endl
          << "Dummy = " << internal::has_begin<Dummy<double>> << std::endl;
  // check is_vectorizable:
  deallog
    << "is_vectorizable:" << std::endl
    << "LinearAlgebra::distributed::Vector<double> && double = "
    << internal::is_vectorizable<LinearAlgebra::distributed::Vector<double>,
                                 double>::value
    << std::endl
    << "LinearAlgebra::distributed::Vector<double> && float = "
    << internal::is_vectorizable<LinearAlgebra::distributed::Vector<double>,
                                 float>::value
    << std::endl
    << "TrilinosWrappers::MPI::Vector && double = "
    << internal::is_vectorizable<TrilinosWrappers::MPI::Vector, double>::value
    << std::endl
    << "Vector = " << internal::is_vectorizable<Vector<double>, double>::value
    << std::endl;

  deallog << "OK" << std::endl;
}
