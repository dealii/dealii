//----------------------------  dof_tools_13.cc  ---------------------------
//    dof_tools_13.cc,v 1.1 2003/02/16 23:55:57 wolf Exp
//    Version: 
//
//    Copyright (C) 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  dof_tools_13.cc  ---------------------------

#include "../tests.h"
#include "dof_tools_common.cc"
#include <lac/vector.h>

// check
//   DoFTools::distribute_cell_to_dof_vector



std::string output_file_name = "dof_tools_13.output";


template <int dim>
void
check_this (const DoFHandler<dim> &dof_handler)
{
                                   // this doesn't make much sense if
                                   // the element is not primitive
  if (dof_handler.get_fe().is_primitive() == false)
    return;
  
  Vector<double> cell_data (dof_handler.get_tria().n_active_cells());
  for (unsigned int i=0; i<cell_data.size(); ++i)
    cell_data(i) = i;

                                   // distribute to first component
  Vector<double> dof_data (dof_handler.n_dofs());
  DoFTools::distribute_cell_to_dof_vector (dof_handler,
                                           cell_data,
                                           dof_data);
                                   // output every third element
  for (unsigned int i=0; i<dof_data.size(); i+=3)
    deallog << dof_data(i) << " ";
  deallog << std::endl;  


                                   // distribute to last component
  DoFTools::distribute_cell_to_dof_vector (dof_handler,
                                           cell_data,
                                           dof_data,
                                           dof_handler.get_fe().n_components()-1);
  for (unsigned int i=0; i<dof_data.size(); i+=3)
    deallog << dof_data(i) << " ";
  deallog << std::endl;  
}
