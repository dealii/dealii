//----------------------------  multigrid.templates.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
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
#include <numerics/data_out.h>
#include <multigrid/multigrid.h>
#include <algorithm>
#include <fstream>

#include <lac/sparse_matrix.h>


//TODO: This function needs to be specially implemented, since in 2d mode we use faces
#if deal_II_dimension == 1

template <int dim>
template <typename number>
void
Multigrid<dim>::copy_to_mg (const Vector<number>&)
{
  Assert(false, ExcNotImplemented());
}

#else

template <int dim>
template <typename number>
void
Multigrid<dim>::copy_to_mg (const Vector<number>& osrc)
{
				   // Make src a real finite element function
  Vector<number> src = osrc;
//  constraints->distribute(src);

  const unsigned int dofs_per_cell = mg_dof_handler->get_fe().dofs_per_cell;
  const unsigned int dofs_per_face = mg_dof_handler->get_fe().dofs_per_face;

				   // set the elements of the vectors
				   // on all levels to zero
  defect.clear();
  
  vector<unsigned int> global_dof_indices (dofs_per_cell);
  vector<unsigned int> level_dof_indices  (dofs_per_cell);
  vector<unsigned int> level_face_indices (dofs_per_face);

				   // initialize the objects with
				   // their correct size
  for (unsigned int level=minlevel; level<=maxlevel ; ++level)
    {
      const unsigned int system_size = (*level_matrices)[level].m();
      
      solution[level].reinit(system_size);
      defect[level].reinit(system_size);
    };

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
      MGDoFHandler<dim>::active_cell_iterator
	level_cell = mg_dof_handler->begin_active(level);
      const MGDoFHandler<dim>::active_cell_iterator
	level_end  = mg_dof_handler->end_active(level);

//TODO: Treat hanging nodes properly
// The single-level vector is not an FE-function, because the values at
// hanging nodes are set to zero. This should be treated before the restriction.

				       // Compute coarse level right hand side
				       // by restricting from fine level.
      for (; level_cell!=level_end; ++level_cell)
	{
	  DoFObjectAccessor<dim, dim>& global_cell = *level_cell;
					   // get the dof numbers of
					   // this cell for the global
					   // and the level-wise
					   // numbering
	  global_cell.get_dof_indices(global_dof_indices);
	  level_cell->get_mg_dof_indices (level_dof_indices);

					   // transfer the global
					   // defect in the vector
					   // into the level-wise one
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    defect[level](level_dof_indices[i]) = src(global_dof_indices[i]);

	  for (unsigned int face_n=0; face_n<GeometryInfo<dim>::faces_per_cell; ++face_n)
	    {
	      const MGDoFHandler<dim>::face_iterator face = level_cell->face(face_n);
	      if (face->has_children())
		{
		  face->get_mg_dof_indices(level_face_indices);


						   // Delete values on refinement edge,
						   // since restriction will add them again.
		  for (unsigned int i=0; i<dofs_per_face; ++i)
		    defect[level](level_face_indices[i])
		      = 0.;
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

#endif


template <int dim>
template <typename number>
void
Multigrid<dim>::copy_from_mg(Vector<number> &dst) const
{
  const unsigned int dofs_per_cell = mg_dof_handler->get_fe().dofs_per_cell;

  vector<unsigned int> global_dof_indices (dofs_per_cell);
  vector<unsigned int> level_dof_indices (dofs_per_cell);

  MGDoFHandler<dim>::active_cell_iterator
    level_cell = mg_dof_handler->begin_active();
  const MGDoFHandler<dim>::active_cell_iterator
    endc = mg_dof_handler->end();

				   // traverse all cells and copy the
				   // data appropriately to the output
				   // vector

				   // Is the level monotonuosly increasing?

  for (; level_cell != endc; ++level_cell)
    {
      DoFObjectAccessor<dim, dim>& global_cell = *level_cell;
      const unsigned int level = level_cell->level();
      
				       // get the dof numbers of
				       // this cell for the global
				       // and the level-wise
				       // numbering
      global_cell.get_dof_indices (global_dof_indices);
      level_cell->get_mg_dof_indices(level_dof_indices);

				       // copy level-wise data to
				       // global vector
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	dst(global_dof_indices[i]) = solution[level](level_dof_indices[i]);
    };

				   // clear constrained nodes
  constraints->set_zero(dst);
}



template <int dim>
void
Multigrid<dim>::print_vector (const unsigned int level,
			      const Vector<double>& v,
			      const char* name) const
{
  if (level!=maxlevel)
    return;
  
  const DoFHandler<dim>* dof = mg_dof_handler;
  
  Vector<double> out_vector(dof->n_dofs());

  out_vector = -10000;
  
  const unsigned int dofs_per_cell = mg_dof_handler->get_fe().dofs_per_cell;

  vector<unsigned int> global_dof_indices (dofs_per_cell);
  vector<unsigned int> level_dof_indices (dofs_per_cell);

  DoFHandler<dim>::active_cell_iterator
    global_cell = dof->begin_active(level);
  MGDoFHandler<dim>::active_cell_iterator
    level_cell = mg_dof_handler->begin_active(level);
  const MGDoFHandler<dim>::cell_iterator
    endc = mg_dof_handler->end(level);

				   // traverse all cells and copy the
				   // data appropriately to the output
				   // vector
  for (; level_cell != endc; ++level_cell, ++global_cell)
    {
      global_cell->get_dof_indices (global_dof_indices);
      level_cell->get_mg_dof_indices(level_dof_indices);

				       // copy level-wise data to
				       // global vector
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	out_vector(global_dof_indices[i])
	  = v(level_dof_indices[i]);
    }

  ofstream out_file(name);
  DataOut<dim> out;
  out.attach_dof_handler(*dof);
  out.add_data_vector(out_vector, "v");
  out.build_patches(5);
  out.write_gnuplot(out_file);
}

#endif

