// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/lac/vector.h>

#include "../tests.h"

#include "dof_tools_common.h"

// check
//   FiniteElement::n_dofs_per_object



template <int dim>
void
check_this(const DoFHandler<dim> &dof_handler)
{
  const FiniteElement<dim> &fe = dof_handler.get_fe();
  deallog << fe.dofs_per_vertex << ' ' << fe.dofs_per_line << ' '
          << fe.dofs_per_quad << ' ' << fe.dofs_per_hex << std::endl;
  deallog << fe.template n_dofs_per_object<0>() << ' '
          << fe.template n_dofs_per_object<1>() << ' '
          << fe.template n_dofs_per_object<2>() << ' '
          << fe.template n_dofs_per_object<3>() << std::endl;
}
