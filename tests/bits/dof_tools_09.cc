//----------------------------  dof_tools_9.cc  ---------------------------
//    $Id$
//    Version: 
//
//    Copyright (C) 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  dof_tools_9.cc  ---------------------------

#include "../tests.h"
#include "dof_tools_common.cc"

// check
//   DoFTools::map_dof_to_boundary_indices(const DoFHandler<int>     &,
//                                         const std::set<unsigned char> &,
//                                         std::vector<unsigned int> &)


std::string output_file_name = "dof_tools_09.output";


template <int dim>
void
check_this (const DoFHandler<dim> &dof_handler)
{
  std::vector<unsigned int> map(dof_handler.n_dofs());
  std::set<unsigned char> boundary_ids;

                                   // check for boundary id 0 alone
  boundary_ids.insert (0);
  DoFTools::map_dof_to_boundary_indices (dof_handler, map);
  for (unsigned int i=0; i<map.size(); ++i)
    deallog << (map[i] == DoFHandler<dim>::invalid_dof_index ?
                -1 : static_cast<signed int>(map[i]))
            << " ";
  deallog << std::endl;

                                   // check for boundary id 0 and 1
  boundary_ids.insert (1);
  DoFTools::map_dof_to_boundary_indices (dof_handler, map);
  for (unsigned int i=0; i<map.size(); ++i)
    deallog << (map[i] == DoFHandler<dim>::invalid_dof_index ?
                -1 : static_cast<signed int>(map[i]))
            << " ";
  deallog << std::endl;
}
