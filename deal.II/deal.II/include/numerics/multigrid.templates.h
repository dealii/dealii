/*----------------------------   multigrid.templates.h     ---------------------------*/
/*      $Id$                 */
#ifndef __multigrid.templates_H
#define __multigrid.templates_H
/*----------------------------   multigrid.templates.h     ---------------------------*/


#include <grid/dof_constraints.h>
#include <numerics/multigrid.h>
#include <algorithm>
#include <fstream>




template <int dim>
template <typename number>
void
MG<dim>::copy_to_mg(const Vector<number>& src)
{
  const unsigned int fe_dofs = mg_dof_handler->get_fe().total_dofs;
  const unsigned int face_dofs = mg_dof_handler->get_fe().dofs_per_face;

				   // set the elements of the vectors
				   // on all levels to zero
  defect.clear();
//  constraints->condense(src);
  
  vector<int> global_dof_indices (fe_dofs);
  vector<int> level_dof_indices (fe_dofs);
  vector<int> level_face_indices (face_dofs);

				   // traverse the grid top-down
				   // (i.e. starting with the most
				   // refined grid). this way, we can
				   // always get that part of one
				   // level of the output vector which
				   // corresponds to a region which is
				   // more refined, by restriction of
				   // the respective vector on the
				   // next finer level, which we then
				   // already have built.
  for (int level=maxlevel; level>=static_cast<int>(minlevel); --level)
    {
      DoFHandler<dim>::active_cell_iterator
	global_cell = mg_dof_handler->DoFHandler<dim>::begin_active(level);
      MGDoFHandler<dim>::active_cell_iterator
	level_cell = mg_dof_handler->begin_active(level);
      const MGDoFHandler<dim>::active_cell_iterator
	level_end  = mg_dof_handler->end_active(level);

      for (; level_cell!=level_end; ++level_cell, ++global_cell)
	{
					   // get the dof numbers of
					   // this cell for the global
					   // and the level-wise
					   // numbering
	  global_cell->get_dof_indices(global_dof_indices);
	  level_cell->get_mg_dof_indices (level_dof_indices);

					   // transfer the global
					   // defect in the vector
					   // into the level-wise one
	  for (unsigned int i=0; i<fe_dofs; ++i)
	    defect[level](level_dof_indices[i]) = src(global_dof_indices[i]);

//TODO: what happens here?
					   // Delete values on refinement edge
	  for (unsigned int face_n=0; face_n<GeometryInfo<dim>::faces_per_cell; ++face_n)
	    {
	      const MGDoFHandler<dim>::face_iterator face = level_cell->face(face_n);
	      if (face->has_children())
		{
		  face->get_mg_dof_indices(level_face_indices);
		  for (unsigned int i=0; i<face_dofs; ++i)
		    defect[level](level_face_indices[i]) = 0.;
		};
	    };
	};

				       // for that part of the level
				       // which is further refined:
				       // get the defect by
				       // restriction of the defect on
				       // one level higher
      if (static_cast<unsigned int>(level) < maxlevel)
	transfer->restrict_and_add (level+1, defect[level], defect[level+1]);
    };
};



template <int dim>
template <typename number>
void
MG<dim>::copy_from_mg(Vector<number> &dst) const
{
  const unsigned int fe_dofs = mg_dof_handler->get_fe().total_dofs;

  vector<int> global_dof_indices (fe_dofs);
  vector<int> level_dof_indices (fe_dofs);

  DoFHandler<dim>::active_cell_iterator
    global_cell = mg_dof_handler->DoFHandler<dim>::begin_active();
  MGDoFHandler<dim>::active_cell_iterator
    level_cell = mg_dof_handler->begin_active();
  const MGDoFHandler<dim>::active_cell_iterator
    endc = mg_dof_handler->end();

				   // traverse all cells and copy the
				   // data appropriately to the output
				   // vector
  for (; level_cell != endc; ++level_cell, ++global_cell)
    {
      const unsigned int level = level_cell->level();
      
				       // get the dof numbers of
				       // this cell for the global
				       // and the level-wise
				       // numbering
      global_cell->get_dof_indices (global_dof_indices);
      level_cell->get_mg_dof_indices(level_dof_indices);

				       // copy level-wise data to
				       // global vector
      for (unsigned int i=0; i<fe_dofs; ++i)
	dst(global_dof_indices[i]) = solution[level](level_dof_indices[i]);
    };

				   // clear constrained nodes
  constraints->set_zero(dst);
}




/*----------------------------   multigrid.templates.h     ---------------------------*/
/* end of #ifndef __multigrid.templates_H */
#endif
/*----------------------------   multigrid.templates.h     ---------------------------*/
