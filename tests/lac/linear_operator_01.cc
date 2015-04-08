// ---------------------------------------------------------------------
//
// Copyright (C) 2015 by the deal.II authors
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

// Test the LinearOperator template on a trivial vector implementation
// :: RightVector -> LeftVector

#include "../tests.h"

#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/vector_memory.templates.h>

using namespace dealii;

// Dummy vectors with different, non convertible types:

struct LeftVector
{
  typedef double value_type;
  double value;

  LeftVector & operator *= (double scale)
  {
    value *= scale;
    return *this;
  }
  LeftVector & operator /= (double scale)
  {
    value /= scale;
    return *this;
  }
  LeftVector & operator += (const LeftVector &u)
  {
    value += u.value;
    return *this;
  }
  int size() const
  {
    return 1;
  }
  std::size_t memory_consumption () const
  {
    return 1;
  }
};

struct RightVector
{
  typedef double value_type;
  double value;

  RightVector & operator *= (double scale)
  {
    value *= scale;
    return *this;
  }
  RightVector & operator /= (double scale)
  {
    value /= scale;
    return *this;
  }
  RightVector & operator += (const RightVector &u)
  {
    value += u.value;
    return *this;
  }
  int size() const
  {
    return 1;
  }
  std::size_t memory_consumption () const
  {
    return 1;
  }
};

int main()
{
  initlog();

  // Create to distinct linear operators:

  LinearOperator<LeftVector, RightVector> multiply2;
  multiply2.vmult = [](LeftVector &v, const RightVector &u)
  {
    v.value = 2 * u.value;
  };
  multiply2.vmult_add = [](LeftVector &v, const RightVector &u)
  {
    v.value += 2 * u.value;
  };
  multiply2.Tvmult = [](RightVector &v, const LeftVector &u)
  {
    v.value = 2 * u.value;
  };
  multiply2.Tvmult_add = [](RightVector &v, const LeftVector &u)
  {
    v.value += 2 * u.value;
  };
  multiply2.reinit_range_vector = [](LeftVector &, bool)
  {
  };
  multiply2.reinit_domain_vector = [](RightVector &, bool)
  {
  };

  auto multiply4 = multiply2;
  multiply4.vmult = [](LeftVector &v, const RightVector &u)
  {
    v.value = 4 * u.value;
  };
  multiply4.vmult_add = [](LeftVector &v, const RightVector &u)
  {
    v.value += 4 * u.value;
  };
  multiply4.Tvmult = [](RightVector &v, const LeftVector &u)
  {
    v.value = 4 * u.value;
  };
  multiply4.Tvmult_add = [](RightVector &v, const LeftVector &u)
  {
    v.value += 4 * u.value;
  };


  // Small unit tests for all functions:

  RightVector u = { 4. };
  LeftVector v = { 0. };

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

  test = multiply4 * 4.;
  test.vmult(v, u);
  deallog << "(4 * 4) * " << u.value << " = " << v.value << std::endl;

  // operator* and transpose

  auto test2 = transpose_linop(multiply2) * multiply4;
  RightVector w = { 0. };
  test2.vmult(w, u);
  deallog << "(2 * 4) * " << u.value << " = " << w.value << std::endl;

  // identity

  auto test3 = identity_linop(test2.reinit_range_vector) + test2;
  test3.vmult(w, u);
  deallog << "(1 + 2 * 4) * " << u.value << " = " << w.value << std::endl;
}
