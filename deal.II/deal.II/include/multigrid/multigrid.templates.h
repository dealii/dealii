//----------------------------  multigrid.templates.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  multigrid.templates.h  ---------------------------
#ifndef __deal2__multigrid_templates_h
#define __deal2__multigrid_templates_h


#include <dofs/dof_constraints.h>
//#include <numerics/data_out.h>
#include <multigrid/multigrid.h>
#include <algorithm>
#include <fstream>

#include <lac/sparse_matrix.h>


template <int dim>
void
Multigrid<dim>::print_vector (const unsigned int /*level*/,
			      const Vector<double>& /*v*/,
			      const char* /*name*/) const
{
//    if (level!=maxlevel)
//      return;
  
//    const DoFHandler<dim>* dof = mg_dof_handler;
  
//    Vector<double> out_vector(dof->n_dofs());

//    out_vector = -10000;
  
//    const unsigned int dofs_per_cell = mg_dof_handler->get_fe().dofs_per_cell;

//    std::vector<unsigned int> global_dof_indices (dofs_per_cell);
//    std::vector<unsigned int> level_dof_indices (dofs_per_cell);

//    typename DoFHandler<dim>::active_cell_iterator
//      global_cell = dof->begin_active(level);
//    typename MGDoFHandler<dim>::active_cell_iterator
//      level_cell = mg_dof_handler->begin_active(level);
//    const typename MGDoFHandler<dim>::cell_iterator
//      endc = mg_dof_handler->end(level);

//  				   // traverse all cells and copy the
//  				   // data appropriately to the output
//  				   // vector
//    for (; level_cell != endc; ++level_cell, ++global_cell)
//      {
//        global_cell->get_dof_indices (global_dof_indices);
//        level_cell->get_mg_dof_indices(level_dof_indices);

//  				       // copy level-wise data to
//  				       // global vector
//        for (unsigned int i=0; i<dofs_per_cell; ++i)
//  	out_vector(global_dof_indices[i])
//  	  = v(level_dof_indices[i]);
//      }

//    std::ofstream out_file(name);
//    DataOut<dim> out;
//    out.attach_dof_handler(*dof);
//    out.add_data_vector(out_vector, "v");
//    out.build_patches(5);
//    out.write_gnuplot(out_file);
}

#endif
