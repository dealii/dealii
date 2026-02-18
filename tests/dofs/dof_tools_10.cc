// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2003 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


#include <deal.II/fe/mapping_q.h>

#include "../tests.h"

#include "dof_tools_common.h"

// check
//   DoFTools::map_dofs_to_support_points (const Mapping<dim> &,
//                   const DoFHandler<dim> &,
//                   std::vector<Point<dim> > &)



template <int dim>
void
check_this(const DoFHandler<dim> &dof_handler)
{
  // don't check if FE has no support
  // points
  if (dof_handler.get_fe().get_unit_support_points().size() == 0)
    return;

  std::vector<Point<dim>> map(dof_handler.n_dofs());
  MappingQ<dim>           mapping(2);

  DoFTools::map_dofs_to_support_points(mapping, dof_handler, map);

  // output every third element
  for (unsigned int i = 0; i < map.size(); i += 3)
    deallog << map[i] << ' ';
  deallog << std::endl;
}
