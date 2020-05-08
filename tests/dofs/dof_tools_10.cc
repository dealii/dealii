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


#include <deal.II/fe/mapping_q.h>

#include "../tests.h"

#include "dof_tools_common.h"
#include "dof_tools_common_fake_hp.h"

// check
//   DoFTools::map_dofs_to_support_points (const Mapping<dim> &,
//                   const DoFHandler<dim> &,
//                   std::vector<Point<dim> > &)



template <typename DoFHandlerType>
void
check_this(const DoFHandlerType &dof_handler)
{
  // don't check if fe has no support
  // points
  if (dof_handler.get_fe().get_unit_support_points().size() == 0)
    return;

  std::vector<Point<DoFHandlerType::dimension>> map(dof_handler.n_dofs());
  MappingQ<DoFHandlerType::dimension>           mapping(2);

  DoFTools::map_dofs_to_support_points(mapping, dof_handler, map);

  // output every third element
  for (unsigned int i = 0; i < map.size(); i += 3)
    deallog << map[i] << " ";
  deallog << std::endl;
}
