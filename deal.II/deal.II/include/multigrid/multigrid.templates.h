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


#include <multigrid/multigrid.h>

/* --------------------------- inline functions --------------------- */


template <class VECTOR>
template <int dim>
Multigrid<VECTOR>::Multigrid (const MGDoFHandler<dim>& mg_dof_handler,
			      const MGMatrixBase<VECTOR>& matrix,
			      const MGCoarseGrid<VECTOR>& coarse,
			      const MGTransfer<VECTOR>& transfer,
			      const MGSmoother<VECTOR>& pre_smooth,
			      const MGSmoother<VECTOR>& post_smooth,
			      const unsigned int min_level,
			      const unsigned int max_level)
		:
		minlevel(min_level),
		maxlevel(std::min(mg_dof_handler.get_tria().n_levels()-1,
				  max_level)),
		defect(minlevel,maxlevel),
		solution(minlevel,maxlevel),
		t(minlevel,maxlevel),
		matrix(&matrix),
		coarse(&coarse),
		transfer(&transfer),
		pre_smooth(&pre_smooth),
		post_smooth(&post_smooth)
{
};


template <class VECTOR>
void
Multigrid<VECTOR>::level_mgstep(const unsigned int level)
{
#ifdef MG_DEBUG
  char *name = new char[100];
  sprintf(name, "MG%d-defect",level);
  print_vector(level, defect[level], name);
#endif

  solution[level] = 0.;
  
  if (level == minlevel)
    {
      (*coarse)(level, solution[level], defect[level]);
#ifdef MG_DEBUG
      sprintf(name, "MG%d-solution",level);
      print_vector(level, solution[level], name);
#endif
      return;
    }

			   // smoothing of the residual by modifying s
  pre_smooth->smooth(level, solution[level], defect[level]);

#ifdef MG_DEBUG
  sprintf(name, "MG%d-pre",level);
  print_vector(level, solution[level], name);
#endif
  
				   // t = -A*solution[level]
  matrix->vmult(level, t[level], solution[level]);
  t[level].scale(-1);
  
				   // make t rhs of lower level
				   // The non-refined parts of the
				   // coarse-level defect already contain
				   // the global defect, the refined parts
				   // its restriction.
  for (unsigned int l = level;l>minlevel;--l)
    {
      t[l-1] = 0.;
      transfer->restrict_and_add (l, t[l-1], t[l]);
      defect[l-1] += t[l-1];
    }

				   // add additional DG contribution
//  edge_vmult(level, defect[level-1], defect[level]);
  
				   // do recursion
  level_mgstep(level-1);

				   // reset size of the auxiliary
				   // vector, since it has been
				   // resized in the recursive call to
				   // level_mgstep directly above
  t[level] = 0.;

				   // do coarse grid correction

  transfer->prolongate(level, t[level], solution[level-1]);

#ifdef MG_DEBUG
  sprintf(name, "MG%d-cgc",level);
  print_vector(level, t, name);
#endif

  solution[level] += t[level];
  
				   // post-smoothing
  post_smooth->smooth(level, solution[level], defect[level]);

#ifdef MG_DEBUG
  sprintf(name, "MG%d-post",level);
  print_vector(level, solution[level], name);

  delete[] name;
#endif
}


template <class VECTOR>
void
Multigrid<VECTOR>::vcycle()
{
				   // The defect vector has been
				   // initialized by copy_to_mg. Now
				   // adjust the other vectors.
  solution.resize(minlevel, maxlevel);
  t.resize(minlevel, maxlevel);
  
  for (unsigned int level=minlevel; level<=maxlevel;++level)
    {
      solution[level].reinit(defect[level]);
      t[level].reinit(defect[level]);
    }

  level_mgstep (maxlevel);
//  abort ();
}



// template <int dim>
// void
// Multigrid<dim>::print_vector (const unsigned int /*level*/,
// 			      const Vector<double>& /*v*/,
// 			      const char* /*name*/) const
// {
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
//}

#endif
