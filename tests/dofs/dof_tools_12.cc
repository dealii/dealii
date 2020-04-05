// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2018 by the deal.II authors
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


#include "../tests.h"

#include "dof_tools_common.h"
//#include "dof_tools_common_fake_hp.h"

// check
//   DoFTools::extract_dofs



template <typename DoFHandlerType>
void
check_this(const DoFHandlerType &dof_handler)
{
  std::vector<bool> mask(dof_handler.get_fe(0).n_components(), false);

  // only select first component
  mask[0] = true;
  output_bool_vector(DoFTools::extract_dofs(dof_handler, ComponentMask(mask)));

  // also select last component
  mask.back() = true;
  output_bool_vector(DoFTools::extract_dofs(dof_handler, ComponentMask(mask)));
}
