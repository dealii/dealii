// $Id$

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
  
  vector<int> index(fe_dofs);
  vector<int> mgindex(fe_dofs);
  vector<int> mgfaceindex(face_dofs);

  for (int level=maxlevel; level>=static_cast<int>(minlevel); --level)
    {
      DoFHandler<dim>::active_cell_iterator
	dc = mg_dof_handler->DoFHandler<dim>::begin_active(level);
      MGDoFHandler<dim>::active_cell_iterator
	level_cell = mg_dof_handler->begin_active(level);
      const MGDoFHandler<dim>::active_cell_iterator
	level_end  = mg_dof_handler->end_active(level);
	
//TODO: Document on difference between #end# and #end_level#

      for (; level_cell!=level_end; ++level_cell, ++dc)
	{
	  dc->get_dof_indices(index);
	  level_cell->get_mg_dof_indices (mgindex);
	  for (unsigned int i=0; i<fe_dofs; ++i)
	    defect[level](mgindex[i]) = src(index[i]);

//TODO: what happens here?
					   // Delete values on refinement edge
	  for (unsigned int face_n=0; face_n<GeometryInfo<dim>::faces_per_cell; ++face_n)
	    {
	      const MGDoFHandler<dim>::face_iterator face = level_cell->face(face_n);
	      if (face->has_children())
		{
		  face->get_mg_dof_indices(mgfaceindex);
		  for (unsigned int i=0; i<face_dofs; ++i)
		    defect[level](mgfaceindex[i]) = 0.;
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
MG<dim>::copy_from_mg(Vector<number>& dst) const
{
  const unsigned int fe_dofs = mg_dof_handler->get_fe().total_dofs;

  vector<int> index(fe_dofs);
  vector<int> mgindex(fe_dofs);

//   for (unsigned int i=minlevel; i<= maxlevel; ++i)
//     {
//       char name[10];
//       sprintf(name,"CF%d",i);
//       ofstream of(name);
//       write_gnuplot(*dofs, s[i], i, of);
//     }

  DoFHandler<dim>::active_cell_iterator dc = mg_dof_handler->DoFHandler<dim>::begin_active();
  for (MGDoFHandler<dim>::active_cell_iterator c = mg_dof_handler->begin_active();
       c != mg_dof_handler->end(); ++c, ++dc)
    {
      dc->get_dof_indices(index);
      c->get_mg_dof_indices(mgindex);
      for (unsigned int i=0;i<fe_dofs;++i)
	dst(index[i]) = solution[c->level()](mgindex[i]);
    } 
  constraints->set_zero(dst);
//   {
//     ofstream of("CF");
//     DataOut<2> out;
//     out.attach_dof_handler(*dofs);
//     out.add_data_vector(dst,"solution");
//     out.write_gnuplot(of,1);
//   }
}

