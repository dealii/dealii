// $Id$

#include <numerics/multigrid.h>

/*
template <int dim>
MultiGrid::
*/

template <int dim>
MultiGrid<dim>::MultiGrid(const MGDoFHandler<dim>& dofs, const MGTransferBase& transfer)
		:
		MultiGridBase(transfer, 0, dofs.get_tria().n_levels()-1),
		dofs(&dofs),
		matrix_structures(maxlevel-minlevel),
		matrices(maxlevel-minlevel)
{
  for(unsigned l = minlevel; l < maxlevel; ++l)
    {
      dofs.make_sparsity_pattern(l, matrix_structures[l-minlevel]);
      matrices[l-minlevel].reinit(matrix_structures[l-minlevel]);
    }
}

template <int dim> template <typename number>
void
MultiGrid<dim>::copy_to_mg(vector<Vector<float> >& dst,
			   const Vector<number>& src) const
{

}

template <int dim> template <typename number>
void
MultiGrid<dim>::copy_from_mg(Vector<number>& dst,
			     const vector<Vector<float> >& src) const
{
// active iterator
}


/*
template <int dim>
MultiGrid::
template <int dim>
MultiGrid::
template <int dim>
MultiGrid::
template <int dim>
MultiGrid::

*/
