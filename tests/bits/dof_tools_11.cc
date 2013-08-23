// ---------------------------------------------------------------------
// $Id$
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
//   DoFTools::map_support_points_to_dofs (const Mapping<dim> &,
//                   const DoFHandler<dim> &,
//                   std::map<Point<dim>, unsigned int, Comp> &)


std::string output_file_name = "dof_tools_11/output";

struct Comp
{
  template <int dim>
  bool operator() (const Point<dim> p, const Point<dim> q) const
  {
    // avoid distinguishing points by roundoff
    for (unsigned int d=0; d<dim; ++d)
      if (p(d) < q(d) - 1e-14)
        return true;
      else if (p(d) > q(d) + 1e-14)
        return false;

    return false;
  }
};


template <int dim>
void
check_this (const DoFHandler<dim> &dof_handler)
{
  // don't check if fe has no support
  // points
  if (dof_handler.get_fe().get_unit_support_points().size() == 0)
    return;

  std::map<Point<dim>, types::global_dof_index, Comp> map;
  MappingQ<dim> mapping(2);

  DoFTools::map_support_points_to_dofs (mapping, dof_handler, map);

  // output every second element
  unsigned int j=0;
  for (typename std::map<Point<dim>, types::global_dof_index, Comp>::const_iterator
       i=map.begin(); i!=map.end(); ++i,++j)
    if (j%2 == 0)
      deallog << i->first << " " << i->second << std::endl;
  deallog << std::endl;
}
