// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
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


// Check that some of the common tensor operations can be performed with
// tensors of expressions

#include <deal.II/base/numbers.h>

#include <deal.II/differentiation/sd.h>

#include "../tests.h"

namespace SD = Differentiation::SD;

template <int dim, typename NumberType>
Tensor<2, dim, NumberType>
make_tensor(const NumberType &val)
{
  Tensor<2, dim, NumberType> out;
  for (unsigned int i = 0; i < out.n_independent_components; ++i)
    out[out.unrolled_to_component_indices(i)] = NumberType(i) + val;
  return out;
}

template <int dim, typename NumberType>
SymmetricTensor<2, dim, NumberType>
make_symm_tensor(const NumberType &val)
{
  SymmetricTensor<2, dim, NumberType> out;
  for (unsigned int i = 0; i < out.n_independent_components; ++i)
    out[out.unrolled_to_component_indices(i)] = NumberType(i) + val;
  return out;
}

template <int dim, typename NumberType>
void
test_tensor()
{
  typedef SD::Expression                       SD_number_t;
  typedef Tensor<2, dim, SD_number_t>          SD_tensor_t;
  typedef SymmetricTensor<2, dim, SD_number_t> SD_symm_tensor_t;

  // --- Values ---
  deallog.push("Initialize");
  SD_tensor_t      a(make_tensor<dim>(NumberType(10)));
  SD_tensor_t      b(make_tensor<dim>(NumberType(10)));
  SD_symm_tensor_t as(make_symm_tensor<dim>(NumberType(10)));
  SD_symm_tensor_t bs(make_symm_tensor<dim>(NumberType(10)));
  deallog << "a: " << a << std::endl;
  deallog << "b: " << b << std::endl;
  deallog << "as: " << as << std::endl;
  deallog << "bs: " << bs << std::endl;
  deallog.pop();

  deallog.push("Relational operators 1");
  deallog << "a == b: " << (a == b) << std::endl;
  deallog << "a != b: " << (a != b) << std::endl;
  if (numbers::NumberTraits<NumberType>::is_complex == false)
    deallog << "a.norm() < b.norm(): " << (a.norm() < b.norm()) << std::endl;
  deallog.pop();

  deallog.push("Set new values");
  a  = make_tensor<dim>(NumberType(5));
  b  = make_tensor<dim>(NumberType(9));
  as = make_symm_tensor<dim>(NumberType(5));
  bs = make_symm_tensor<dim>(NumberType(7));
  deallog << "a: " << a << std::endl;
  deallog << "b: " << b << std::endl;
  deallog << "as: " << as << std::endl;
  deallog << "bs: " << bs << std::endl;
  deallog.pop();
  //
  deallog.push("Relational operators 2");
  deallog << "a == b: " << (a == b) << std::endl;
  deallog << "a != b: " << (a != b) << std::endl;
  if (numbers::NumberTraits<NumberType>::is_complex == false)
    deallog << "a.norm() < b.norm(): " << (a.norm() < b.norm()) << std::endl;
  deallog.pop();

  deallog.push("Math operators");
  deallog.push("Tensor");
  deallog << "a+b: " << (a + b) << std::endl;
  deallog << "a-b: " << (a - b) << std::endl;
  deallog << "a*b: " << (a * b) << std::endl;
  deallog.pop();
  deallog.push("SymmetricTensor");
  deallog << "as+bs: " << (as + bs) << std::endl;
  deallog << "as-bs: " << (as - bs) << std::endl;
  deallog << "as*bs: " << (as * bs) << std::endl;
  const SD_tensor_t as_ns = as;
  deallog << "as: " << as << std::endl;
  deallog << "as_ns: " << as_ns << std::endl;
  const SD_tensor_t as_ns2(as);
  deallog << "as+b: " << (as + b) << std::endl;
  deallog << "as-b: " << (as - b) << std::endl;
  deallog << "as*b: " << (static_cast<SD_tensor_t>(as) * b) << std::endl;
  deallog.pop();
  deallog.pop();

  // --- Normal numbers ---
  deallog.push("Copy constructor");
  SD_tensor_t      c(a);   // Copy constructor
  SD_symm_tensor_t cs(as); // Copy constructor
  deallog << "c: " << c << std::endl;
  deallog << "cs: " << cs << std::endl;
  deallog.pop();

  deallog.push("Math operators"); // including those with arithmetic numbers
  const NumberType s(2);
  deallog.push("Tensor");
  c += b;
  deallog << "c+=b: " << c << std::endl;
  c *= s;
  deallog << "c*=s: " << c << std::endl;
  c -= b;
  deallog << "c-=b: " << c << std::endl;
  c /= s;
  deallog << "c/=s: " << c << std::endl;
  deallog.pop();
  deallog.push("SymmetricTensor");
  cs += bs;
  deallog << "cs+=bs: " << cs << std::endl;
  cs *= s;
  deallog << "cs*=s: " << cs << std::endl;
  cs -= bs;
  deallog << "cs-=bs: " << cs << std::endl;
  cs /= s;
  deallog << "cs/=s: " << cs << std::endl;
  deallog.pop();
  deallog.pop();

  // --- Symbols ---
  deallog.push("Symbols");
  const SD_number_t      y("y");
  const SD_tensor_t      x(SD::make_tensor_of_symbols<2, dim>("x"));
  const SD_symm_tensor_t z(SD::make_symmetric_tensor_of_symbols<2, dim>("z"));
  deallog << "x: " << x << std::endl;
  deallog << "y: " << y << std::endl;
  deallog << "z: " << z << std::endl;
  deallog.pop();

  deallog.push("Math operators"); // including those with arithmetic numbers
  deallog.push("Tensor");
  c += x;
  deallog << "c+=x: " << c << std::endl;
  c *= y;
  deallog << "c*=y: " << c << std::endl;
  c *= s;
  deallog << "c*=s: " << c << std::endl;
  c -= x;
  deallog << "c-=x: " << c << std::endl;
  c /= y;
  deallog << "c/=y: " << c << std::endl;
  c /= s;
  deallog << "c/=s: " << c << std::endl;
  deallog.pop();
  deallog.push("SymmetricTensor");
  cs += z;
  deallog << "cs+=z: " << cs << std::endl;
  cs *= y;
  deallog << "cs*=y: " << cs << std::endl;
  cs *= s;
  deallog << "cs*=s: " << cs << std::endl;
  cs -= z;
  deallog << "cs-=z: " << cs << std::endl;
  cs /= y;
  deallog << "cs/=y: " << cs << std::endl;
  cs /= s;
  deallog << "cs/=s: " << cs << std::endl;
  deallog.pop();
  deallog.pop();

  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();

  const int dim = 2;

  deallog.push("Integer");
  test_tensor<dim, int>();
  deallog.pop();

  deallog.push("Float");
  test_tensor<dim, float>();
  deallog.pop();

  deallog.push("Double");
  test_tensor<dim, double>();
  deallog.pop();

  deallog.push("Complex float");
  test_tensor<dim, std::complex<float>>();
  deallog.pop();

  deallog.push("Complex double");
  test_tensor<dim, std::complex<double>>();
  deallog.pop();

  deallog << "OK" << std::endl;
}
