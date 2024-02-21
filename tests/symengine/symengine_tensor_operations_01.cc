// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check that symbol tensors and symbolic function tensors can be created
// using the utility functions

#include <deal.II/differentiation/sd.h>

#include "../tests.h"

namespace SD = Differentiation::SD;


template <int dim>
void
test()
{
  deallog.push("dim " + Utilities::to_string(dim));

  using SD_number_t = SD::Expression;

  const Tensor<1, dim, SD_number_t> symb_v =
    SD::make_vector_of_symbols<dim>("v");
  const Tensor<0, dim, SD_number_t> symb_T0 =
    SD::make_tensor_of_symbols<0, dim>("T0");
  const Tensor<1, dim, SD_number_t> symb_T1 =
    SD::make_tensor_of_symbols<1, dim>("T1");
  const Tensor<2, dim, SD_number_t> symb_T2 =
    SD::make_tensor_of_symbols<2, dim>("T2");
  const Tensor<3, dim, SD_number_t> symb_T3 =
    SD::make_tensor_of_symbols<3, dim>("T3");
  const Tensor<4, dim, SD_number_t> symb_T4 =
    SD::make_tensor_of_symbols<4, dim>("T4");
  const SymmetricTensor<2, dim, SD_number_t> symb_ST2 =
    SD::make_symmetric_tensor_of_symbols<2, dim>("ST2");
  const SymmetricTensor<4, dim, SD_number_t> symb_ST4 =
    SD::make_symmetric_tensor_of_symbols<4, dim>("ST4");

  deallog << "Symbol vector: " << symb_v << std::endl;
  deallog << "Symbol tensor (rank-0): " << symb_T0 << std::endl;
  deallog << "Symbol tensor (rank-1): " << symb_T1 << std::endl;
  deallog << "Symbol tensor (rank-2): " << symb_T2 << std::endl;
  deallog << "Symbol tensor (rank-3): " << symb_T3 << std::endl;
  deallog << "Symbol tensor (rank-4): " << symb_T4 << std::endl;
  deallog << "Symbol symmetric tensor (rank-2): " << symb_ST2 << std::endl;
  deallog << "Symbol symmetric tensor (rank-4): " << symb_ST4 << std::endl;

  SD::types::substitution_map sub_map;
  sub_map[symb_T0] = 0.0;
  for (unsigned int i = 0; i < dim; ++i)
    {
      sub_map[symb_v[i]] = 0.0;
    }

  const Tensor<1, dim, SD_number_t> symb_func_v =
    SD::make_vector_of_symbolic_functions<dim>("f_v", sub_map);
  const Tensor<0, dim, SD_number_t> symb_func_T0 =
    SD::make_tensor_of_symbolic_functions<0, dim>("f_T0", sub_map);
  const Tensor<1, dim, SD_number_t> symb_func_T1 =
    SD::make_tensor_of_symbolic_functions<1, dim>("f_T1", sub_map);
  const Tensor<2, dim, SD_number_t> symb_func_T2 =
    SD::make_tensor_of_symbolic_functions<2, dim>("f_T2", sub_map);
  const Tensor<3, dim, SD_number_t> symb_func_T3 =
    SD::make_tensor_of_symbolic_functions<3, dim>("f_T3", sub_map);
  const Tensor<4, dim, SD_number_t> symb_func_T4 =
    SD::make_tensor_of_symbolic_functions<4, dim>("f_T4", sub_map);
  const SymmetricTensor<2, dim, SD_number_t> symb_func_ST2 =
    SD::make_symmetric_tensor_of_symbolic_functions<2, dim>("f_ST2", sub_map);
  const SymmetricTensor<4, dim, SD_number_t> symb_func_ST4 =
    SD::make_symmetric_tensor_of_symbolic_functions<4, dim>("f_ST4", sub_map);

  deallog << "Symbol function vector: " << symb_func_v << std::endl;
  deallog << "Symbol function tensor (rank-0): " << symb_func_T0 << std::endl;
  deallog << "Symbol function tensor (rank-1): " << symb_func_T1 << std::endl;
  deallog << "Symbol function tensor (rank-2): " << symb_func_T2 << std::endl;
  deallog << "Symbol function tensor (rank-3): " << symb_func_T3 << std::endl;
  deallog << "Symbol function tensor (rank-4): " << symb_func_T4 << std::endl;
  deallog << "Symbol function symmetric tensor (rank-2): " << symb_func_ST2
          << std::endl;
  deallog << "Symbol function symmetric tensor (rank-4): " << symb_func_ST4
          << std::endl;

  deallog << "OK" << std::endl;
  deallog.pop();
}


int
main()
{
  initlog();

  test<1>();
  test<2>();
  test<3>();

  deallog << "OK" << std::endl;
}
