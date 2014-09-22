// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
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


#include "../tests.h"
#include "dof_tools_common.h"
#include <deal.II/fe/mapping_q.h>

// check
//   DoFTools::map_dofs_to_support_points (const Mapping<dim> &,
//                   const DoFHandler<dim> &,
//                   std::vector<Point<dim> > &)


std::string output_file_name = "output";


template <int dim>
void
check_this (const DoFHandler<dim> &dof_handler)
{
  // don't check if fe has no support
  // points
  if (dof_handler.get_fe().get_unit_support_points().size() == 0)
    return;

  std::vector<Point<dim> > map(dof_handler.n_dofs());
  MappingQ<dim> mapping(2);

  DoFTools::map_dofs_to_support_points (mapping, dof_handler, map);

  // output every third element
  for (unsigned int i=0; i<map.size(); i+=3)
    deallog << map[i] << " ";
  deallog << std::endl;
}
