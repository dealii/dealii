// $Id$

#include <numerics/multigrid.h>

template <int dim>
MG<dim>::MG(const MGDoFHandler<dim>& dofs,
	    const MGMatrix<SparseMatrix<float> >& matrices,
	    const MGTransferBase& transfer)
		:
		MGBase(transfer, 0, dofs.get_tria().n_levels()-1),
		dofs(&dofs),
		matrices(&matrices)
{
  for (unsigned int level = minlevel; level<=maxlevel ; ++level)
    {
      s[level].reinit(matrices[level].m());
      d[level].reinit(matrices[level].m());
    }
}

template <int dim>
template <typename number>
void
MG<dim>::copy_to_mg(const Vector<number>& src)
{
  const unsigned int fe_dofs = dofs->get_fe().total_dofs;
  
  vector<int> index(fe_dofs);
  vector<int> mgindex(fe_dofs);

  DoFHandler<dim>::active_cell_iterator dc = dofs->DoFHandler<dim>::begin_active();
  for (MGDoFHandler<dim>::active_cell_iterator c = dofs->begin_active()
				     ; c != dofs->end() ; ++c, ++dc)
    {
      dc->get_dof_indices(index);
      c->get_mg_dof_indices(mgindex);
      for (unsigned int i=0;i<fe_dofs;++i)
	d[maxlevel](mgindex[i]) = src(index[i]);
    } 
}

template <int dim> template <typename number>
void
MG<dim>::copy_from_mg(Vector<number>& dst) const
{
  const unsigned int fe_dofs = dofs->get_fe().total_dofs;

  vector<int> index(fe_dofs);
  vector<int> mgindex(fe_dofs);

  DoFHandler<dim>::active_cell_iterator dc = dofs->DoFHandler<dim>::begin_active();
  for (MGDoFHandler<dim>::active_cell_iterator c = dofs->begin_active()
				     ; c != dofs->end() ; ++c, ++dc)
    {
      dc->get_dof_indices(index);
      c->get_mg_dof_indices(mgindex);
      for (unsigned int i=0;i<fe_dofs;++i)
	dst(index[i]) = s[maxlevel](mgindex[i]);
    }
}

template <int dim>
void
MG<dim>::level_residual(unsigned int l,
			Vector<float> &r,
			const Vector<float> &u,
			const Vector<float> &b)
{
  (*matrices)[l].residual(r,u,b);
}

/*
template <int dim>
MG::
template <int dim>
MG::
template <int dim>
MG::

*/
