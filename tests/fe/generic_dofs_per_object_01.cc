// ---------------------------------------------------------------------
//
// Copyright (C) 2022 by the deal.II authors
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

// Test internal::GenericDoFsPerObject::generate().

#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_wedge_p.h>

#include "../tests.h"

using namespace dealii;

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
