// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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


// Verify that SymmetricTensor is trivially copyable.

// TODO not all compilers that support enough of a subset of C++11 to compile
// the library (notably GCC 4.8) implement std::is_trivally_copyable. At some
// point in the future we should use that instead of the boost equivalent.

#include <deal.II/base/symmetric_tensor.h>

#include <boost/type_traits.hpp>

#include <complex>

#include "../tests.h"

template <typename Number>
void
test()
{
  deallog << "SymmetricTensor<2, 1> is trivially copyable: "
          << boost::has_trivial_copy<SymmetricTensor<2, 1, Number>>::value
          << std::endl;
  deallog << "SymmetricTensor<2, 2> is trivially copyable: "
          << boost::has_trivial_copy<SymmetricTensor<2, 2, Number>>::value
          << std::endl;
  deallog << "SymmetricTensor<2, 3> is trivially copyable: "
          << boost::has_trivial_copy<SymmetricTensor<2, 3, Number>>::value
          << std::endl;

  deallog << "SymmetricTensor<4, 1> is trivially copyable: "
          << boost::has_trivial_copy<SymmetricTensor<4, 1, Number>>::value
          << std::endl;
  deallog << "SymmetricTensor<4, 2> is trivially copyable: "
          << boost::has_trivial_copy<SymmetricTensor<4, 2, Number>>::value
          << std::endl;
  deallog << "SymmetricTensor<4, 3> is trivially copyable: "
          << boost::has_trivial_copy<SymmetricTensor<4, 3, Number>>::value
          << std::endl;
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
