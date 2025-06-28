// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Verify that Tensor is trivially copyable.

// TODO not all compilers that support enough of a subset of C++11 to compile
// the library (notably GCC 4.8) implement std::is_trivially_copyable. At some
// point in the future we should use that instead of the boost equivalent.

#include <deal.II/base/tensor.h>

#include <boost/type_traits.hpp>

#include <complex>

#include "../tests.h"

template <typename Number>
void
test()
{
  deallog << "Tensor<0, 1> is trivially copyable: "
          << boost::has_trivial_copy<Tensor<0, 1, Number>>::value << std::endl;
  deallog << "Tensor<0, 2> is trivially copyable: "
          << boost::has_trivial_copy<Tensor<0, 2, Number>>::value << std::endl;

  deallog << "Tensor<1, 1> is trivially copyable: "
          << boost::has_trivial_copy<Tensor<1, 1, Number>>::value << std::endl;
  deallog << "Tensor<1, 2> is trivially copyable: "
          << boost::has_trivial_copy<Tensor<1, 2, Number>>::value << std::endl;
  deallog << "Tensor<1, 3> is trivially copyable: "
          << boost::has_trivial_copy<Tensor<1, 3, Number>>::value << std::endl;

  deallog << "Tensor<2, 1> is trivially copyable: "
          << boost::has_trivial_copy<Tensor<2, 1, Number>>::value << std::endl;
  deallog << "Tensor<2, 2> is trivially copyable: "
          << boost::has_trivial_copy<Tensor<2, 2, Number>>::value << std::endl;
  deallog << "Tensor<2, 3> is trivially copyable: "
          << boost::has_trivial_copy<Tensor<2, 3, Number>>::value << std::endl;

  deallog << "Tensor<3, 1> is trivially copyable: "
          << boost::has_trivial_copy<Tensor<3, 1, Number>>::value << std::endl;
  deallog << "Tensor<3, 2> is trivially copyable: "
          << boost::has_trivial_copy<Tensor<3, 2, Number>>::value << std::endl;
  deallog << "Tensor<3, 3> is trivially copyable: "
          << boost::has_trivial_copy<Tensor<3, 3, Number>>::value << std::endl;
}

int
main()
{
  initlog();

  deallog << std::boolalpha;
  deallog << "testing float" << std::endl;
  test<float>();

  deallog << "testing double" << std::endl;
  test<double>();

  deallog << "testing std::complex<float>" << std::endl;
  test<std::complex<float>>();

  deallog << "testing std::complex<double>" << std::endl;
  test<std::complex<double>>();
}
