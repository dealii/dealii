// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check ROLAdaptor with an IndexSet mask.

#include <deal.II/base/index_set.h>
#include <deal.II/base/signaling_nan.h>

#include <deal.II/lac/vector.h>

#include <deal.II/trilinos/rol_adaptor.h>

#include <ROL_Vector.hpp>

#include "../tests.h"


template <typename VectorType>
void
test()
{
  const unsigned int size = 5;

  VectorType u(size);
  VectorType v(size);
  VectorType w(size);

  // mask excludes first and last indices
  IndexSet mask(size);
  mask.add_range(1, size - 1);

  // fill masked parts of vectors with some numbers, and the other parts with
  // signaling nans
  unsigned int tmp = 0;
  for (unsigned int i = 0; i < size; ++i)
    {
      if (mask.is_element(i))
        {
          u[i] = ++tmp;
          v[i] = ++tmp;
          w[i] = ++tmp;
        }
      else
        {
          u[i] = numbers::signaling_nan<typename VectorType::value_type>();
          v[i] = numbers::signaling_nan<typename VectorType::value_type>();
          w[i] = numbers::signaling_nan<typename VectorType::value_type>();
        }
    }

  // wrap for ROL
  TrilinosWrappers::ROLAdaptor<VectorType> u_rol(ROL::makePtrFromRef(u), mask);
  TrilinosWrappers::ROLAdaptor<VectorType> v_rol(ROL::makePtrFromRef(v), mask);
  TrilinosWrappers::ROLAdaptor<VectorType> w_rol(ROL::makePtrFromRef(w), mask);

  // let ROL do some checks
  u_rol.checkVector(v_rol, w_rol, true, deallog.get_file_stream());
}


int
main(int argc, char **argv)
{
  initlog();

  test<Vector<double>>();
}
