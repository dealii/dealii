//----------------------------  dof_tools_11.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  dof_tools_11.cc  ---------------------------

#include "../tests.h"
#include "dof_tools_common.cc"
#include <fe/mapping_q.h>

// check
//   DoFTools::map_support_points_to_dofs (const Mapping<dim> &,
//				           const DoFHandler<dim> &,
//				           std::map<Point<dim>, unsigned int, Comp> &)


std::string output_file_name = "dof_tools_11.output";

struct Comp
{
    template <int dim>
    bool operator() (const Point<dim> p, const Point<dim> q) const
      {
        for (unsigned int d=0; d<dim; ++d)
          if (p(d) < q(d))
            return true;
          else if (p(d) > q(d))
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

  std::map<Point<dim>, unsigned int, Comp> map;
  MappingQ<dim> mapping(2);

  DoFTools::map_support_points_to_dofs (mapping, dof_handler, map);

                                   // output every second element
  unsigned int j=0;
  for (typename std::map<Point<dim>, unsigned int, Comp>::const_iterator
         i=map.begin(); i!=map.end(); ++i,++j)
    if (j%2 == 0)
      deallog << i->first << " " << i->second << std::endl;
  deallog << std::endl;
}
