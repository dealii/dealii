/* $Id$ */


#include <lac/sparsity_pattern.h>
#include <multigrid/mg_dof_handler.h>
#include <multigrid/mg_dof_accessor.h>
#include <grid/tria_iterator.h>
#include <multigrid/mg_dof_tools.h>
#include <fe/fe.h>


template <int dim>
void MGDoFTools::make_sparsity_pattern (const MGDoFHandler<dim> &mg_dof_handler,
					const unsigned int       level,
					SparsityPattern         &sparsity)
{
  Assert (sparsity.n_rows() == mg_dof_handler.n_dofs(level),
	  ExcDimensionMismatch (sparsity.n_rows(), mg_dof_handler.n_dofs(level)));
  Assert (sparsity.n_cols() == mg_dof_handler.n_dofs(level),
	  ExcDimensionMismatch (sparsity.n_cols(), mg_dof_handler.n_dofs(level)));

  const unsigned int dofs_per_cell = mg_dof_handler.get_fe().dofs_per_cell;
  vector<unsigned int> dofs_on_this_cell(dofs_per_cell);
  MGDoFHandler<dim>::cell_iterator cell = mg_dof_handler.begin(level),
				   endc = mg_dof_handler.end(level);
  for (; cell!=endc; ++cell) 
    {
      cell->get_mg_dof_indices (dofs_on_this_cell);
				       // make sparsity pattern for this cell
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=0; j<dofs_per_cell; ++j)
	  sparsity.add (dofs_on_this_cell[i],
			dofs_on_this_cell[j]);
    };
};



// explicit instantiations
template void MGDoFTools::make_sparsity_pattern (const MGDoFHandler<deal_II_dimension> &,
						 const unsigned int,
						 SparsityPattern &);

