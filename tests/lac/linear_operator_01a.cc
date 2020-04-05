// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2018 by the deal.II authors
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

// Variant of linear_operator_01 that uses std::complex<double> as
// value_type: Test the LinearOperator template on a trivial vector
// implementation :: RightVector -> LeftVector with complex numbers.

#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/vector_memory.templates.h>

#include <complex>

#include "../tests.h"



// Dummy vectors with different, non convertible types:

struct LeftVector
{
  typedef std::complex<double> value_type;
  value_type                   value;

  LeftVector &
  operator=(value_type new_value)
  {
    value = new_value;
    return *this;
  }
  LeftVector &
  operator*=(value_type scale)
  {
    value *= scale;
    return *this;
  }
  LeftVector &
  operator/=(value_type scale)
  {
    value /= scale;
    return *this;
  }
  LeftVector &
  operator+=(const LeftVector &u)
  {
    value += u.value;
    return *this;
  }
  int
  size() const
  {
    return 1;
  }
  std::size_t
  memory_consumption() const
  {
    return 1;
  }
};

struct RightVector
{
  typedef std::complex<double> value_type;
  value_type                   value;

  RightVector &
  operator=(value_type new_value)
  {
    value = new_value;
    return *this;
  }
  RightVector &
  operator*=(value_type scale)
  {
    value *= scale;
    return *this;
  }
  RightVector &
  operator/=(value_type scale)
  {
    value /= scale;
    return *this;
  }
  RightVector &
  operator+=(const RightVector &u)
  {
    value += u.value;
    return *this;
  }
  int
  size() const
  {
    return 1;
  }
  std::size_t
  memory_consumption() const
  {
    return 1;
  }
};

int
main()
{
  initlog();

  // Create to distinct linear operators:

  typedef dealii::internal::LinearOperatorImplementation::EmptyPayload Payload;
  LinearOperator<LeftVector, RightVector, Payload> multiply2;
  multiply2.vmult = [](LeftVector &v, const RightVector &u) {
    v.value = 2. * u.value;
  };
  multiply2.vmult_add = [](LeftVector &v, const RightVector &u) {
    v.value += 2. * u.value;
  };
  multiply2.Tvmult = [](RightVector &v, const LeftVector &u) {
    v.value = 2. * u.value;
  };
  multiply2.Tvmult_add = [](RightVector &v, const LeftVector &u) {
    v.value += 2. * u.value;
  };
  multiply2.reinit_range_vector  = [](LeftVector &, bool) {};
  multiply2.reinit_domain_vector = [](RightVector &, bool) {};

  auto multiply4  = multiply2;
  multiply4.vmult = [](LeftVector &v, const RightVector &u) {
    v.value = 4. * u.value;
  };
  multiply4.vmult_add = [](LeftVector &v, const RightVector &u) {
    v.value += 4. * u.value;
  };
  multiply4.Tvmult = [](RightVector &v, const LeftVector &u) {
    v.value = 4. * u.value;
  };
  multiply4.Tvmult_add = [](RightVector &v, const LeftVector &u) {
    v.value += 4. * u.value;
  };


  // Small unit tests for all functions:

  RightVector u = {4.};
  LeftVector  v = {0.};

  // vmult, vmult_add

  multiply2.vmult(v, u);
  deallog << "2 * " << u.value << " = " << v.value << std::endl;

  multiply4.vmult_add(v, u);
  deallog << "... + 4 * " << u.value << " = " << v.value << std::endl;

  multiply4.vmult(v, u);
  deallog << "4 * " << u.value << " = " << v.value << std::endl;

  multiply2.vmult_add(v, u);
  deallog << "... + 2 * " << u.value << " = " << v.value << std::endl;

  // Tvmult, Tvmult_add

  v.value = 4.;

  multiply2.Tvmult(u, v);
  deallog << "2 * " << v.value << " = " << u.value << std::endl;

  multiply4.Tvmult_add(u, v);
  deallog << "... + 4 * " << v.value << " = " << u.value << std::endl;

  multiply4.Tvmult(u, v);
  deallog << "4 * " << v.value << " = " << u.value << std::endl;

  multiply2.Tvmult_add(u, v);
  deallog << "... + 2 * " << v.value << " = " << u.value << std::endl;

  // operator+, operator-, operator+=, operator-=

  auto test = multiply2 + multiply4;
  test.vmult(v, u);
  deallog << "(2 + 4) * " << u.value << " = " << v.value << std::endl;

  test = multiply2 - multiply4;
  test.vmult(v, u);
  deallog << "(2 - 4) * " << u.value << " = " << v.value << std::endl;

  test += multiply2;
  test.vmult(v, u);
  deallog << "(2 - 4 + 2) * " << u.value << " = " << v.value << std::endl;

  test -= multiply4;
  test.vmult(v, u);
  deallog << "(2 - 4 + 2 - 4) * " << u.value << " = " << v.value << std::endl;

  // operator* with scalar

  test = 4. * multiply4;
  test.vmult(v, u);
  deallog << "(4 * 4) * " << u.value << " = " << v.value << std::endl;

  test = multiply4;
  test *= std::complex<double>(4.);
  test.vmult(v, u);
  deallog << "(4 * 4) * " << u.value << " = " << v.value << std::endl;

  test = multiply4 * 4.;
  test.vmult(v, u);
  deallog << "(4 * 4) * " << u.value << " = " << v.value << std::endl;

  // operator* and transpose

  auto        test2 = transpose_operator(multiply2) * multiply4;
  RightVector w     = {0.};
  test2.vmult(w, u);
  deallog << "(2 * 4) * " << u.value << " = " << w.value << std::endl;

  test2 *= identity_operator(test2.reinit_range_vector);
  test2.vmult(w, u);
  deallog << "(2 * 4) * 1 * " << u.value << " = " << w.value << std::endl;

  // identity

  auto test3 = identity_operator(test2.reinit_range_vector) + test2;
  test3.vmult(w, u);
  deallog << "(1 + 2 * 4) * " << u.value << " = " << w.value << std::endl;

  // null operator

  auto test4 = null_operator(test2);
  test4.vmult(w, u);
  deallog << " 0 * " << u.value << " = " << w.value << std::endl;
}
