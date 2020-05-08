// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2020 by the deal.II authors
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


#include <deal.II/lac/vector.h>

#include "../tests.h"

#include "dof_tools_common.h"
#include "dof_tools_common_fake_hp.h"

// check
//   FiniteElement::n_dofs_per_object



template <typename DoFHandlerType>
void
check_this(const DoFHandlerType &dof_handler)
{
  const FiniteElement<DoFHandlerType::dimension> &fe = dof_handler.get_fe();
  deallog << fe.dofs_per_vertex << ' ' << fe.dofs_per_line << ' '
          << fe.dofs_per_quad << ' ' << fe.dofs_per_hex << std::endl;
  deallog << fe.template n_dofs_per_object<0>() << ' '
          << fe.template n_dofs_per_object<1>() << ' '
          << fe.template n_dofs_per_object<2>() << ' '
          << fe.template n_dofs_per_object<3>() << std::endl;
}
