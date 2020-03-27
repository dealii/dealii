// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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


// Check that symbol substitution works for tensors

#include <deal.II/differentiation/sd.h>

#include "../tests.h"

namespace SD = Differentiation::SD;

using SD_number_t = SD::Expression;

template <int rank, int dim, template <int, int, typename> class TensorType>
void
test(const TensorType<rank, dim, double>      a,
     const TensorType<rank, dim, double>      b,
     const TensorType<rank, dim, SD_number_t> symb_a,
     const TensorType<rank, dim, SD_number_t> symb_b)
{
  using Tensor_SD_number_t = TensorType<rank, dim, SD_number_t>;
  using Tensor_t           = TensorType<rank, dim, double>;

  deallog.push("Substitution: Individual");
  {
    const Tensor_SD_number_t symb_a_plus_b = symb_a + symb_b;
    deallog << "Symbolic a+b: " << symb_a_plus_b << std::endl;

    const Tensor_SD_number_t symb_a_plus_b_1 =
      SD::substitute(symb_a_plus_b, symb_a, a);
    deallog << "Symbolic a+b (a=1): " << symb_a_plus_b_1 << std::endl;

    const Tensor_SD_number_t symb_a_plus_b_subs =
      SD::substitute(symb_a_plus_b_1, symb_b, b);
    deallog << "Symbolic a+b (a=1, b=2): " << symb_a_plus_b_subs << std::endl;
  }
  deallog.pop();

  deallog.push("Substitution: Vector");
  {
    const Tensor_SD_number_t symb_a_plus_b = symb_a + symb_b;

    const std::vector<std::pair<Tensor_SD_number_t, Tensor_t>> symbol_values{
      std::make_pair(symb_a, a), std::make_pair(symb_b, b)};

    const Tensor_SD_number_t symb_a_plus_b_subs =
      SD::substitute(symb_a_plus_b, symbol_values);
    deallog << "Symbolic a+b (a=1, b=2): " << symb_a_plus_b_subs << std::endl;
  }
  deallog.pop();

  deallog.push("Substitution: Map");
  {
    const Tensor_SD_number_t symb_a_plus_b = symb_a + symb_b;

    const SD::types::substitution_map substitution_map =
      SD::make_substitution_map(std::make_pair(symb_a, a),
                                std::make_pair(symb_b, b));

    const Tensor_SD_number_t symb_a_plus_b_subs =
      SD::substitute(symb_a_plus_b, substitution_map);
    deallog << "Symbolic a+b (a=1, b=2): " << symb_a_plus_b_subs << std::endl;
  }
  deallog.pop();

  deallog.push("Substitution with evaluation");
  {
    const Tensor_SD_number_t symb_a_plus_b = symb_a + symb_b;

    const SD::types::substitution_map substitution_map =
      SD::make_substitution_map(std::make_pair(symb_a, a),
                                std::make_pair(symb_b, b));

    const Tensor_t symb_a_plus_b_subs =
      SD::substitute_and_evaluate<double>(symb_a_plus_b, substitution_map);
    Assert(symb_a_plus_b_subs == (a + b), ExcInternalError());
  }
  deallog.pop();

  deallog << "OK" << std::endl << std::endl;
}


template <int rank, int dim>
void
test_tensor()
{
  deallog << "Tensor: Rank " << rank << ", dim " << dim << std::endl;

  using Tensor_SD_number_t = Tensor<rank, dim, SD_number_t>;
  using Tensor_t           = Tensor<rank, dim, double>;

  Tensor_t t_a, t_b;
  for (auto it = t_a.begin_raw(); it != t_a.end_raw(); ++it)
    *it = 1.0;
  for (auto it = t_b.begin_raw(); it != t_b.end_raw(); ++it)
    *it = 2.0;

  const Tensor_SD_number_t symb_t_a =
    SD::make_tensor_of_symbols<rank, dim>("a");
  const Tensor_SD_number_t symb_t_b =
    SD::make_tensor_of_symbols<rank, dim>("b");

  test(t_a, t_b, symb_t_a, symb_t_b);
}


template <int rank, int dim>
void
test_symmetric_tensor()
{
  deallog << "SymmetricTensor: Rank " << rank << ", dim " << dim << std::endl;

  using Tensor_SD_number_t = SymmetricTensor<rank, dim, SD_number_t>;
  using Tensor_t           = SymmetricTensor<rank, dim, double>;

  Tensor_t t_a, t_b;
  for (auto it = t_a.begin_raw(); it != t_a.end_raw(); ++it)
    *it = 1.0;
  for (auto it = t_b.begin_raw(); it != t_b.end_raw(); ++it)
    *it = 2.0;

  const Tensor_SD_number_t symb_t_a =
    SD::make_symmetric_tensor_of_symbols<rank, dim>("a");
  const Tensor_SD_number_t symb_t_b =
    SD::make_symmetric_tensor_of_symbols<rank, dim>("b");

  test(t_a, t_b, symb_t_a, symb_t_b);
}


int
main()
{
  initlog();

  const int dim = 2;

  test_tensor<0, dim>();
  test_tensor<1, dim>();
  test_tensor<2, dim>();
  test_tensor<3, dim>();
  test_tensor<4, dim>();

  test_symmetric_tensor<2, dim>();
  test_symmetric_tensor<4, dim>();

  deallog << "OK" << std::endl;
}
