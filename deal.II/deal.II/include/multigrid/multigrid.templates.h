// $Id$

#include <numerics/multigrid.h>
#include <algorithm>
#include <fstream>




template <int dim>
MG<dim>::MG(const MGDoFHandler<dim>& dofs,
	    ConstraintMatrix& hanging_nodes,
	    const MGMatrix<SparseMatrix<double> >& matrices,
	    const MGTransferBase& transfer,
	    unsigned int minl, unsigned int maxl)
		:
		MGBase(transfer, minl,
		       min(dofs.get_tria().n_levels()-1, maxl)),
		dofs(&dofs),
		matrices(&matrices),
		hanging_nodes(hanging_nodes)
{
  for (unsigned int level = minlevel; level<=maxlevel ; ++level)
    {
      solution[level].reinit(matrices[level].m());
      defect[level].reinit(matrices[level].m());
    };
};



template <int dim>
template <typename number>
void
MG<dim>::copy_to_mg(const Vector<number>& src)
{
  const unsigned int fe_dofs = dofs->get_fe().total_dofs;
  const unsigned int face_dofs = dofs->get_fe().dofs_per_face;

				   // set the elements of the vectors
				   // on all levels to zero
  defect.clear();
//  hanging_nodes.condense(src);
  
  vector<int> index(fe_dofs);
  vector<int> mgindex(fe_dofs);
  vector<int> mgfaceindex(face_dofs);
  
//   {
//     ofstream of("CT");
//     DataOut<2> out;
//     out.attach_dof_handler(*dofs);
//     out.add_data_vector(src,"solution");
//     out.write_gnuplot(of,1);
//   }

  for (int level = maxlevel; level >= static_cast<int>(minlevel); --level)
    {
      DoFHandler<dim>::active_cell_iterator dc = dofs->DoFHandler<dim>::begin_active(level);
      MGDoFHandler<dim>::active_cell_iterator c = dofs->begin_active(level);

//TODO: Document on difference between #end# and #end_level#

      for (; c != dofs->end_active(level) ; ++c, ++dc)
	{
	  dc->get_dof_indices(index);
	  c->get_mg_dof_indices(mgindex);
	  for (unsigned int i=0;i<fe_dofs;++i)
	    defect[c->level()](mgindex[i]) = src(index[i]);

					   // Delete values on refinement edge
	  for (unsigned int face_n = 0; face_n< GeometryInfo<dim>::faces_per_cell; ++face_n)
	    {
	      MGDoFHandler<dim>::face_iterator face = c->face(face_n);
	      if (face->has_children())
		{
		  face->get_mg_dof_indices(mgfaceindex);
		  for (unsigned int i=0;i<face_dofs;++i)
		    defect[c->level()](mgfaceindex[i]) = 0.;
		}
	    }
	}
//      if (level > (int) minlevel)
//        transfer->restrict(level, d[level-1], d[level]);
      if (level < (int) maxlevel)
	transfer->restrict(level+1, defect[level], defect[level+1]);
    }
  
//   for (unsigned int i=minlevel; i<= maxlevel; ++i)
//     {
//       char name[10];
//       sprintf(name,"CT%d",i);
//       ofstream of(name);
//       write_gnuplot(*dofs, d[i], i, of);
//     }
}

template <int dim>
template <typename number>
void
MG<dim>::copy_from_mg(Vector<number>& dst) const
{
  const unsigned int fe_dofs = dofs->get_fe().total_dofs;

  vector<int> index(fe_dofs);
  vector<int> mgindex(fe_dofs);

//   for (unsigned int i=minlevel; i<= maxlevel; ++i)
//     {
//       char name[10];
//       sprintf(name,"CF%d",i);
//       ofstream of(name);
//       write_gnuplot(*dofs, s[i], i, of);
//     }

  DoFHandler<dim>::active_cell_iterator dc = dofs->DoFHandler<dim>::begin_active();
  for (MGDoFHandler<dim>::active_cell_iterator c = dofs->begin_active();
       c != dofs->end(); ++c, ++dc)
    {
      dc->get_dof_indices(index);
      c->get_mg_dof_indices(mgindex);
      for (unsigned int i=0;i<fe_dofs;++i)
	dst(index[i]) = solution[c->level()](mgindex[i]);
    } 
  hanging_nodes.set_zero(dst);
//   {
//     ofstream of("CF");
//     DataOut<2> out;
//     out.attach_dof_handler(*dofs);
//     out.add_data_vector(dst,"solution");
//     out.write_gnuplot(of,1);
//   }
}



template <int dim>
void
MG<dim>::level_vmult(unsigned int l,
			Vector<double> &r,
			const Vector<double> &u,
			const Vector<double> &)
{
  (*matrices)[l].vmult(r,u);
  r.scale(-1.);
}

/*
template <int dim>
MG::
template <int dim>
MG::
template <int dim>
MG::

*/
