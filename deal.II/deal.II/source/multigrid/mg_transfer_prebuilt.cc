//-----------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------

#include <base/logstream.h>

#include <lac/vector.h>
#include <lac/block_vector.h>
#include <lac/sparse_matrix.h>
#include <lac/block_sparse_matrix.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_tools.h>
#include <fe/fe.h>
#include <multigrid/mg_dof_handler.h>
#include <multigrid/mg_dof_accessor.h>
#include <multigrid/mg_transfer.h>
#include <multigrid/mg_transfer.templates.h>


template <typename number>
template <int dim>
void MGTransferPrebuilt<number>::build_matrices (
  const MGDoFHandler<dim> &mg_dof)
{
  const unsigned int n_levels      = mg_dof.get_tria().n_levels();
  const unsigned int dofs_per_cell = mg_dof.get_fe().dofs_per_cell;
  
  sizes.resize(n_levels);
  for (unsigned int l=0;l<n_levels;++l)
    sizes[l] = mg_dof.n_dofs(l);
  
				   // reset the size of the array of
				   // matrices. call resize(0) first,
				   // in order to delete all elements
				   // and clear their memory. then
				   // repopulate these arrays
				   //
				   // note that on resize(0), the
				   // shared_ptr class takes care of
				   // deleting the object it points to
				   // by itself
  prolongation_matrices.resize (0);
#ifndef DEAL_PREFER_MATRIX_EZ
  prolongation_sparsities.resize (0);
#endif

  for (unsigned int i=0; i<n_levels-1; ++i)
    {
#ifndef DEAL_PREFER_MATRIX_EZ
      prolongation_sparsities
	.push_back (boost::shared_ptr<SparsityPattern> (new SparsityPattern));
      prolongation_matrices
	.push_back (boost::shared_ptr<SparseMatrix<double> > (new SparseMatrix<double>));
#else
      prolongation_matrices
	.push_back (boost::shared_ptr<SparseMatrix<double> > (new SparseMatrixEZ<double>));
#endif
    }
  
				   // two fields which will store the
				   // indices of the multigrid dofs
				   // for a cell and one of its children
  std::vector<unsigned int> dof_indices_parent (dofs_per_cell);
  std::vector<unsigned int> dof_indices_child (dofs_per_cell);
  
				   // for each level: first build the sparsity
				   // pattern of the matrices and then build the
				   // matrices themselves. note that we only
				   // need to take care of cells on the coarser
				   // level which have children
  for (unsigned int level=0; level<n_levels-1; ++level)
    {
				       // reset the dimension of the structure.
				       // note that for the number of entries
				       // per row, the number of parent dofs
				       // coupling to a child dof is
				       // necessary. this, of course, is the
				       // number of degrees of freedom per
				       // cell
#ifdef DEAL_PREFER_MATRIX_EZ      
      prolongation_matrices[level]->reinit (sizes[level+1],
					    sizes[level],
					    dofs_per_cell);
#else
				       // increment dofs_per_cell
				       // since a useless diagonal
				       // element will be stored
      prolongation_sparsities[level]->reinit (sizes[level+1],
					      sizes[level],
					      dofs_per_cell+1);
      
      for (typename MGDoFHandler<dim>::cell_iterator cell=mg_dof.begin(level);
	   cell != mg_dof.end(level); ++cell)
	if (cell->has_children())
	  {
	    cell->get_mg_dof_indices (dof_indices_parent);

	    for (unsigned int child=0;
		 child<GeometryInfo<dim>::children_per_cell; ++child)
	      {
						 // set an alias to the
						 // prolongation matrix for
						 // this child
		const FullMatrix<double> &prolongation
		  = mg_dof.get_fe().get_prolongation_matrix (child);

		Assert (prolongation.n() != 0, ExcNoProlongation());
		
		cell->child(child)->get_mg_dof_indices (dof_indices_child);

						 // now tag the entries in the
						 // matrix which will be used
						 // for this pair of parent/child
		for (unsigned int i=0; i<dofs_per_cell; ++i)
		  for (unsigned int j=0; j<dofs_per_cell; ++j)
		    if (prolongation(i,j) != 0)
		      {
			prolongation_sparsities[level]->add (dof_indices_child[i],
							     dof_indices_parent[j]);
		      };
	      };
	  };
      prolongation_sparsities[level]->compress ();

      prolongation_matrices[level]->reinit (*prolongation_sparsities[level]);
#endif

				       // now actually build the matrices
      for (typename MGDoFHandler<dim>::cell_iterator cell=mg_dof.begin(level);
	   cell != mg_dof.end(level); ++cell)
	if (cell->has_children())
	  {
	    cell->get_mg_dof_indices (dof_indices_parent);

	    for (unsigned int child=0;
		 child<GeometryInfo<dim>::children_per_cell; ++child)
	      {
						 // set an alias to the
						 // prolongation matrix for
						 // this child
		const FullMatrix<double> &prolongation
		  = mg_dof.get_fe().get_prolongation_matrix (child);
	    
		cell->child(child)->get_mg_dof_indices (dof_indices_child);

						 // now set the entries in the
						 // matrix
		for (unsigned int i=0; i<dofs_per_cell; ++i)
		  for (unsigned int j=0; j<dofs_per_cell; ++j)
		    if (prolongation(i,j) != 0)
		      prolongation_matrices[level]->set (dof_indices_child[i],
							 dof_indices_parent[j],
							 prolongation(i,j));
	      };
	  };
    };
}



template
void MGTransferPrebuilt<float>::build_matrices<deal_II_dimension>
(const MGDoFHandler<deal_II_dimension> &mg_dof);

template
void MGTransferPrebuilt<double>::build_matrices<deal_II_dimension>
(const MGDoFHandler<deal_II_dimension> &mg_dof);

template void
MGTransferPrebuilt<float>::copy_to_mg (
  const MGDoFHandler<deal_II_dimension>&,
  MGLevelObject<Vector<float> >&,
  const Vector<double>&) const;
template void
MGTransferPrebuilt<float>::copy_to_mg (
  const MGDoFHandler<deal_II_dimension>&,
  MGLevelObject<Vector<float> >&,
  const BlockVector<double>&) const;
template void
MGTransferPrebuilt<float>::copy_from_mg (
  const MGDoFHandler<deal_II_dimension>&,
  Vector<double>&,
  const MGLevelObject<Vector<float> >&) const;
template void
MGTransferPrebuilt<float>::copy_from_mg (
  const MGDoFHandler<deal_II_dimension>&,
  BlockVector<double>&,
  const MGLevelObject<Vector<float> >&) const;

template void
MGTransferPrebuilt<double>::copy_to_mg (
  const MGDoFHandler<deal_II_dimension>&,
  MGLevelObject<Vector<double> >&,
  const Vector<double>&) const;
template void
MGTransferPrebuilt<double>::copy_to_mg (
  const MGDoFHandler<deal_II_dimension>&,
  MGLevelObject<Vector<double> >&,
  const BlockVector<double>&) const;
template void
MGTransferPrebuilt<double>::copy_from_mg (
  const MGDoFHandler<deal_II_dimension>&,
  Vector<double>&,
  const MGLevelObject<Vector<double> >&) const;
template void
MGTransferPrebuilt<double>::copy_from_mg (
  const MGDoFHandler<deal_II_dimension>&,
  BlockVector<double>&,
  const MGLevelObject<Vector<double> >&) const;

template void
MGTransferPrebuilt<double>::copy_from_mg_add (
  const MGDoFHandler<deal_II_dimension>&,
  Vector<double>&,
  const MGLevelObject<Vector<double> >&) const;
template void
MGTransferPrebuilt<double>::copy_from_mg_add (
  const MGDoFHandler<deal_II_dimension>&,
  BlockVector<double>&,
  const MGLevelObject<Vector<double> >&) const;
