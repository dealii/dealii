//----------------------------  dof_tools_10.cc  ---------------------------
//    dof_tools_10.cc,v 1.1 2003/02/16 23:55:56 wolf Exp
//    Version: 
//
//    Copyright (C) 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  dof_tools_10.cc  ---------------------------

#include "../tests.h"
#include "dof_tools_common.cc"
#include <fe/mapping_q.h>

// check
//   DoFTools::map_dofs_to_support_points (const Mapping<dim> &,
//				           const DoFHandler<dim> &,
//				           std::vector<Point<dim> > &)


std::string output_file_name = "dof_tools_10.output";


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
