// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/lac/sparsity_pattern.h>

#include "../tests.h"

#include "fe_tools_common.h"

// check
//   FETools::back_interpolate(6)
// for hp::DoFHandler without hanging nodes



template <int dim>
void
check_this(const FiniteElement<dim> &fe1, const FiniteElement<dim> &fe2)
{
  // only check if both elements have
  // support points. otherwise,
  // interpolation doesn't really
  // work
  if ((fe1.get_unit_support_points().size() == 0) ||
      (fe2.get_unit_support_points().size() == 0))
    return;
  //  likewise for non-primitive elements
  if (!fe1.is_primitive() || !fe2.is_primitive())
    return;
  // we need to have dof_constraints
  // for this test
  if (!fe1.constraints_are_implemented() || !fe2.constraints_are_implemented())
    return;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, 0., 1.);
  tria.refine_global(2);

  std::unique_ptr<DoFHandler<dim>> dof1(make_dof_handler(tria, fe1));

  Vector<double> in(dof1->n_dofs());
  for (unsigned int i = 0; i < in.size(); ++i)
    in(i) = i;
  Vector<double> out(dof1->n_dofs());

  FETools::back_interpolate(*dof1, in, fe2, out);
  output_vector(out);
}
