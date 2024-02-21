// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test internal::GenericDoFsPerObject::generate().

#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_wedge_p.h>

#include "../tests.h"


template <int dim, int spacedim>
void
test(const FiniteElement<dim, spacedim> &fe)
{
  const auto dpo = internal::GenericDoFsPerObject::generate(fe);

  deallog << fe.get_name() << std::endl;

  for (const auto &row : dpo.dofs_per_object_exclusive)
    for (const auto &i : row)
      deallog << i << " ";
  deallog << std::endl;

  for (const auto &row : dpo.dofs_per_object_inclusive)
    for (const auto &i : row)
      deallog << i << " ";
  deallog << std::endl;

  for (const auto &row : dpo.object_index)
    for (const auto &i : row)
      deallog << i << " ";
  deallog << std::endl;

  for (const auto &row : dpo.first_object_index_on_face)
    for (const auto &i : row)
      deallog << i << " ";
  deallog << std::endl;
  deallog << std::endl;
}

int
main()
{
  initlog();

  test(FE_Q<1>(1));
  test(FE_Q<1>(2));
  test(FE_Q<2>(1));
  test(FE_Q<2>(2));
  test(FE_Q<3>(1));
  test(FE_Q<3>(2));

  test(FE_SimplexP<2>(1));
  test(FE_SimplexP<2>(2));
  test(FE_SimplexP<3>(1));
  test(FE_SimplexP<3>(2));

  test(FE_PyramidP<3>(1));

  test(FE_WedgeP<3>(1));
  test(FE_WedgeP<3>(2));
}
