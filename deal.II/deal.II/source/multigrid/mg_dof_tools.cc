//----------------------------  mg_dof_tools.cc  ---------------------------
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
//----------------------------  mg_dof_tools.cc  ---------------------------


#include <base/multithread_info.h>
#include <base/thread_management.h>
#include <lac/sparsity_pattern.h>
#include <lac/block_sparsity_pattern.h>
#include <lac/compressed_sparsity_pattern.h>
#include <lac/sparsity_pattern.h>
#include <lac/block_vector.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <multigrid/mg_dof_handler.h>
#include <multigrid/mg_dof_accessor.h>
#include <multigrid/mg_dof_tools.h>
#include <multigrid/mg_base.h>
#include <dofs/dof_tools.h>
#include <fe/fe.h>

//#include <multigrid/mg_tools.templates.h>

#include <algorithm>
#include <numeric>

template <int dim, class SparsityPattern>
void MGTools::make_sparsity_pattern (const MGDoFHandler<dim> &dof,
					SparsityPattern         &sparsity,
					const unsigned int       level)
{
  const unsigned int n_dofs = dof.n_dofs(level);

  Assert (sparsity.n_rows() == n_dofs,
	  ExcDimensionMismatch (sparsity.n_rows(), n_dofs));
  Assert (sparsity.n_cols() == n_dofs,
	  ExcDimensionMismatch (sparsity.n_cols(), n_dofs));

  const unsigned int dofs_per_cell = dof.get_fe().dofs_per_cell;
  std::vector<unsigned int> dofs_on_this_cell(dofs_per_cell);
  typename MGDoFHandler<dim>::cell_iterator cell = dof.begin(level),
					    endc = dof.end(level);
  for (; cell!=endc; ++cell) 
    {
      cell->get_mg_dof_indices (dofs_on_this_cell);
				       // make sparsity pattern for this cell
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=0; j<dofs_per_cell; ++j)
	  sparsity.add (dofs_on_this_cell[i],
			dofs_on_this_cell[j]);
    }
}



template<int dim, class SparsityPattern>
void
MGTools::make_flux_sparsity_pattern (const MGDoFHandler<dim> &dof,
					SparsityPattern       &sparsity,
					const unsigned int level)
{
  const unsigned int n_dofs = dof.n_dofs(level);
  
  Assert (sparsity.n_rows() == n_dofs,
	  ExcDimensionMismatch (sparsity.n_rows(), n_dofs));
  Assert (sparsity.n_cols() == n_dofs,
	  ExcDimensionMismatch (sparsity.n_cols(), n_dofs));

  const unsigned int dofs_per_cell = dof.get_fe().dofs_per_cell;
  std::vector<unsigned int> dofs_on_this_cell(dofs_per_cell);
  std::vector<unsigned int> dofs_on_other_cell(dofs_per_cell);
  typename MGDoFHandler<dim>::cell_iterator cell = dof.begin(level),
					    endc = dof.end(level);
  for (; cell!=endc; ++cell)
    {
      cell->get_mg_dof_indices (dofs_on_this_cell);
				       // make sparsity pattern for this cell
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=0; j<dofs_per_cell; ++j)
	  sparsity.add (dofs_on_this_cell[i],
			dofs_on_this_cell[j]);

				       // Loop over all interior neighbors
      for (unsigned int face = 0;
	   face < GeometryInfo<dim>::faces_per_cell;
	   ++face)
	{
	  if ( (! cell->at_boundary(face)) &&
	       (static_cast<unsigned int>(cell->neighbor_level(face)) == level) )
	    {
	      typename MGDoFHandler<dim>::cell_iterator
		neighbor = cell->neighbor(face);
	      neighbor->get_mg_dof_indices (dofs_on_other_cell);
					       // only add one direction
					       // The other is taken care of
					       // by neighbor.
	      for (unsigned int i=0; i<dofs_per_cell; ++i)
		{
		  for (unsigned int j=0; j<dofs_per_cell; ++j)
		    {
		      sparsity.add (dofs_on_this_cell[i],
				    dofs_on_other_cell[j]);
		    }
		}
	    }
	} 
    }
}



template<int dim, class SparsityPattern>
void
MGTools::make_flux_sparsity_pattern_edge (const MGDoFHandler<dim> &dof,
					     SparsityPattern       &sparsity,
					     const unsigned int level)
{
  Assert ((level>=1) && (level<dof.get_tria().n_levels()),
	  ExcIndexRange(level, 1, dof.get_tria().n_levels()));

  const unsigned int fine_dofs = dof.n_dofs(level);
  const unsigned int coarse_dofs = dof.n_dofs(level-1);
  
  // Matrix maps from fine level to coarse level

  Assert (sparsity.n_rows() == coarse_dofs,
	  ExcDimensionMismatch (sparsity.n_rows(), coarse_dofs));
  Assert (sparsity.n_cols() == fine_dofs,
	  ExcDimensionMismatch (sparsity.n_cols(), fine_dofs));

  const unsigned int dofs_per_cell = dof.get_fe().dofs_per_cell;
  std::vector<unsigned int> dofs_on_this_cell(dofs_per_cell);
  std::vector<unsigned int> dofs_on_other_cell(dofs_per_cell);
  typename MGDoFHandler<dim>::cell_iterator cell = dof.begin(level),
					    endc = dof.end(level);
  for (; cell!=endc; ++cell)
    {
      cell->get_mg_dof_indices (dofs_on_this_cell);
				       // Loop over all interior neighbors
      for (unsigned int face = 0;
	   face < GeometryInfo<dim>::faces_per_cell;
	   ++face)
	{
	  // Neighbor is coarser

	  if ( (! cell->at_boundary(face)) &&
	       (static_cast<unsigned int>(cell->neighbor_level(face)) != level) )
	    {
	      typename MGDoFHandler<dim>::cell_iterator
		neighbor = cell->neighbor(face);
	      neighbor->get_mg_dof_indices (dofs_on_other_cell);

	      for (unsigned int i=0; i<dofs_per_cell; ++i)
		{
		  for (unsigned int j=0; j<dofs_per_cell; ++j)
		    {
		      sparsity.add (dofs_on_other_cell[i],
				    dofs_on_this_cell[j]);
		      sparsity.add (dofs_on_other_cell[j],
				    dofs_on_this_cell[i]);
		    }
		}
	    }
	} 
    }
}



//TODO: Replace FullMatrix by vector2d<bool>
template<int dim, class SparsityPattern>
void
MGTools::make_flux_sparsity_pattern (const MGDoFHandler<dim> &dof,
					SparsityPattern       &sparsity,
					const unsigned int level,
					const FullMatrix<double>& int_mask,
					const FullMatrix<double>& flux_mask)
{
  const unsigned int n_dofs = dof.n_dofs(level);
  const unsigned int n_comp = dof.get_fe().n_components();
  
  Assert (sparsity.n_rows() == n_dofs,
	  ExcDimensionMismatch (sparsity.n_rows(), n_dofs));
  Assert (sparsity.n_cols() == n_dofs,
	  ExcDimensionMismatch (sparsity.n_cols(), n_dofs));
  Assert (int_mask.m() == n_comp,
	  ExcDimensionMismatch (int_mask.m(), n_comp));
  Assert (int_mask.n() == n_comp,
	  ExcDimensionMismatch (int_mask.n(), n_comp));
  Assert (flux_mask.m() == n_comp,
	  ExcDimensionMismatch (flux_mask.m(), n_comp));
  Assert (flux_mask.n() == n_comp,
	  ExcDimensionMismatch (flux_mask.n(), n_comp));
  
  const unsigned int total_dofs = dof.get_fe().dofs_per_cell;
  std::vector<unsigned int> dofs_on_this_cell(total_dofs);
  std::vector<unsigned int> dofs_on_other_cell(total_dofs);
  typename MGDoFHandler<dim>::cell_iterator cell = dof.begin(level),
					    endc = dof.end(level);


  std::vector<std::vector<bool> > int_dof_mask(total_dofs,
				 std::vector<bool>(total_dofs, false));
  std::vector<std::vector<bool> > flux_dof_mask(total_dofs,
				 std::vector<bool>(total_dofs, false));
  for (unsigned int i=0; i<total_dofs; ++i)
    for (unsigned int j=0; j<total_dofs; ++j)
      {
	unsigned int ii = dof.get_fe().system_to_component_index(i).first;
	unsigned int jj = dof.get_fe().system_to_component_index(j).first;
	
	if (int_mask(ii,jj) != 0)
	  int_dof_mask[i][j] = true;
	if (flux_mask(ii,jj) != 0)
	  flux_dof_mask[i][j] = true;
      }
  
  (const_cast<Triangulation<dim>& > (dof.get_tria())).clear_user_flags();
  
  for (; cell!=endc; ++cell)
    {
      cell->get_mg_dof_indices (dofs_on_this_cell);
				       // make sparsity pattern for this cell
      for (unsigned int i=0; i<total_dofs; ++i)
	for (unsigned int j=0; j<total_dofs; ++j)
	  if (int_dof_mask[i][j])
	    sparsity.add (dofs_on_this_cell[i],
			  dofs_on_this_cell[j]);

				       // Loop over all interior neighbors
      for (unsigned int face = 0;
	   face < GeometryInfo<dim>::faces_per_cell;
	   ++face)
	{
	  typename MGDoFHandler<dim>::face_iterator cell_face = cell->face(face);
	  if (cell_face->user_flag_set ())
	    continue;

	  if (! cell->at_boundary (face) )
	    {
	      typename MGDoFHandler<dim>::cell_iterator
		neighbor = cell->neighbor(face);

	      if (neighbor->level() < cell->level())
		continue;

	      unsigned int neighbor_face = cell->neighbor_of_neighbor(face);

	      neighbor->get_mg_dof_indices (dofs_on_other_cell);
	      for (unsigned int i=0; i<total_dofs; ++i)
		{
		  for (unsigned int j=0; j<total_dofs; ++j)
		    if (flux_dof_mask[i][j])
		      {
			sparsity.add (dofs_on_this_cell[i],
				      dofs_on_other_cell[j]);
			sparsity.add (dofs_on_other_cell[i],
				      dofs_on_this_cell[j]);
		      }
		}
	      neighbor->face(neighbor_face)->set_user_flag (); 
	    }
	}
    }
}



template <int dim>
void
MGTools::count_dofs_per_component (const MGDoFHandler<dim>& dof_handler,
				      std::vector<std::vector<unsigned int> >& result)
{
  const unsigned int nlevels = dof_handler.get_tria().n_levels();
  
  Assert (result.size() == nlevels,
	  ExcDimensionMismatch(result.size(), nlevels));

  const unsigned int n_components = dof_handler.get_fe().n_components();
  for (unsigned int l=0;l<nlevels;++l)
    {
      result[l].resize (n_components);
      
				       // special case for only one
				       // component. treat this first
				       // since it does not require any
				       // computations
      if (n_components == 1)
	{
	  result[l][0] = dof_handler.n_dofs(l);
	} else {
					   // otherwise determine the number
					   // of dofs in each component
					   // separately. do so in parallel
	  std::vector<std::vector<bool> >
	    dofs_in_component (n_components,
			       std::vector<bool>(dof_handler.n_dofs(l),
						 false));
	  std::vector<std::vector<bool> >
	    component_select (n_components,
			      std::vector<bool>(n_components, false));
	  Threads::ThreadManager thread_manager;
	  for (unsigned int i=0; i<n_components; ++i)
	    {
	      void (*fun_ptr) (const unsigned int       level,
			       const MGDoFHandler<dim>    &,
			       const std::vector<bool>    &,
			       std::vector<bool>          &)
		= &DoFTools::template extract_level_dofs<dim>;
	      component_select[i][i] = true;
	      Threads::spawn (thread_manager,
			      Threads::encapsulate (fun_ptr)
			      .collect_args (l, dof_handler, component_select[i],
					     dofs_in_component[i]));
	    };
	  thread_manager.wait();
	  
					   // next count what we got
	  for (unsigned int i=0; i<n_components; ++i)
	      result[l][i] = std::count(dofs_in_component[i].begin(),
					dofs_in_component[i].end(),
					true);
	  
					   // finally sanity check
	  Assert (std::accumulate (result[l].begin(),
				   result[l].end(), 0U)
		  ==
		  dof_handler.n_dofs(l),
		  ExcInternalError());
	}
    }
}



template<int dim, typename number>
void
MGTools::reinit_vector (const MGDoFHandler<dim>& mg_dof,
			   MGLevelObject<Vector<number> >& v)
{
  for (unsigned int level=v.get_minlevel();
       level<=v.get_maxlevel();++level)
    {
      unsigned int n = mg_dof.n_dofs (level);
      v[level].reinit(n);
    }

}





template<int dim, typename number>
void
MGTools::reinit_vector (const MGDoFHandler<dim>& mg_dof,
			   MGLevelObject<BlockVector<number> >& v,
			   const std::vector<bool>& selected)
{
  const unsigned int ncomp = mg_dof.get_fe().n_components();
  
  if (selected.size() == 0)
    {
      selected.resize(ncomp);
      fill_n (selected.begin(), ncomp, true);
    }

  Assert (selected.size() == ncomp,
	  ExcDimensionMismatch(selected.size(), ncomp));

  unsigned int n_selected = std::accumulate (selected.begin(),
					     selected.end(),
					     0U);
			     
  std::vector<std::vector<unsigned int> >
    ndofs(mg_dof.n_levels(), ncomp);

  count_dofs_per_component (mg_dof, ndofs);
  
  for (unsigned int level=v.get_minlevel();
       level<=v.get_maxlevel();++level)
    {
      v[level].reinit(n_selected);
      unsigned int k=0;
      for (unsigned int i=0;i<ncomp;++i)
	if (selected[i])
	  v[level].block(k++).reinit(ndofs[level][i]);
    }
}





// explicit instantiations
template void
MGTools::make_sparsity_pattern (const MGDoFHandler<deal_II_dimension> &,
				SparsityPattern &,
				const unsigned int);

template void
MGTools::make_flux_sparsity_pattern (const MGDoFHandler<deal_II_dimension> &,
				     SparsityPattern &,
				     const unsigned int);

template void
MGTools::make_flux_sparsity_pattern_edge (const MGDoFHandler<deal_II_dimension> &,
					  SparsityPattern &,
					  const unsigned int);

template void
MGTools::make_sparsity_pattern (const MGDoFHandler<deal_II_dimension> &,
				CompressedSparsityPattern &,
				const unsigned int);

template void
MGTools::make_flux_sparsity_pattern (const MGDoFHandler<deal_II_dimension> &,
				     CompressedSparsityPattern &,
				     const unsigned int);

template void
MGTools::make_flux_sparsity_pattern_edge (const MGDoFHandler<deal_II_dimension> &,
					  CompressedSparsityPattern &,
					  const unsigned int);

template void
MGTools::make_sparsity_pattern (const MGDoFHandler<deal_II_dimension> &,
				BlockSparsityPattern &,
				const unsigned int);
template void
MGTools::make_flux_sparsity_pattern (const MGDoFHandler<deal_II_dimension> &,
				     BlockSparsityPattern &,
				     const unsigned int);

template void
MGTools::make_flux_sparsity_pattern_edge (const MGDoFHandler<deal_II_dimension> &,
					  BlockSparsityPattern &,
					  const unsigned int);

template void
MGTools::make_sparsity_pattern (const MGDoFHandler<deal_II_dimension> &,
				CompressedBlockSparsityPattern &,
				const unsigned int);

template void
MGTools::make_flux_sparsity_pattern (const MGDoFHandler<deal_II_dimension> &,
				     CompressedBlockSparsityPattern &,
				     const unsigned int);

template void
MGTools::make_flux_sparsity_pattern_edge (const MGDoFHandler<deal_II_dimension> &,
					  CompressedBlockSparsityPattern &,
					  const unsigned int);

#if deal_II_dimension > 1
template void
MGTools::make_flux_sparsity_pattern (const MGDoFHandler<deal_II_dimension> &,
				     SparsityPattern &,
				     const unsigned int,
				     const FullMatrix<double>&,
				     const FullMatrix<double>&);

template void
MGTools::make_flux_sparsity_pattern (const MGDoFHandler<deal_II_dimension> &,
				     CompressedSparsityPattern &,
				     const unsigned int,
				     const FullMatrix<double>&,
				     const FullMatrix<double>&);

template void
MGTools::make_flux_sparsity_pattern (const MGDoFHandler<deal_II_dimension> &,
				     BlockSparsityPattern &,
				     const unsigned int,
				     const FullMatrix<double>&,
				     const FullMatrix<double>&);

template void
MGTools::make_flux_sparsity_pattern (const MGDoFHandler<deal_II_dimension> &,
				     CompressedBlockSparsityPattern &,
				     const unsigned int,
				     const FullMatrix<double>&,
				     const FullMatrix<double>&);

#endif

template void MGTools::reinit_vector (const MGDoFHandler<deal_II_dimension>&,
				      MGLevelObject<Vector<double> >&);
template void MGTools::reinit_vector (const MGDoFHandler<deal_II_dimension>&,
				      MGLevelObject<Vector<float> >&);


//  template void MGTools::copy_to_mg (const MGDoFHandler<deal_II_dimension>&,
//  				   const MGTransferBase<Vector<double> >&,
//  				   MGLevelObject<Vector<double> >&,
//  				   const Vector<double>&);

//  template void MGTools::copy_from_mg(const MGDoFHandler<deal_II_dimension>&,
//  				    Vector<double>&,
//  				    const MGLevelObject<Vector<double> >&);

//  template void MGTools::copy_from_mg_add(const MGDoFHandler<deal_II_dimension>&,
//  					Vector<double>&,
//  					const MGLevelObject<Vector<double> >&);

//  template void MGTools::copy_to_mg (const MGDoFHandler<deal_II_dimension>&,
//  				   const MGTransferBase<Vector<double> >&,
//  				   MGLevelObject<Vector<double> >&,
//  				   const Vector<float>&);

//  template void MGTools::copy_from_mg(const MGDoFHandler<deal_II_dimension>&,
//  				    Vector<double>&,
//  				    const MGLevelObject<Vector<float> >&);

//  template void MGTools::copy_from_mg_add(const MGDoFHandler<deal_II_dimension>&,
//  					Vector<double>&,
//  					const MGLevelObject<Vector<float> >&);

//  template void MGTools::copy_to_mg (const MGDoFHandler<deal_II_dimension>&,
//  				   const MGTransferBase<Vector<double> >&,
//  				   MGLevelObject<Vector<double> >&,
//  				   const BlockVector<double>&);

//  template void MGTools::copy_from_mg(const MGDoFHandler<deal_II_dimension>&,
//  				    BlockVector<double>&,
//  				    const MGLevelObject<Vector<double> >&);

//  template void MGTools::copy_from_mg_add(const MGDoFHandler<deal_II_dimension>&,
//  					BlockVector<double>&,
//  					const MGLevelObject<Vector<double> >&);

//  template void MGTools::copy_to_mg (const MGDoFHandler<deal_II_dimension>&,
//  				   const MGTransferBase<Vector<double> >&,
//  				   MGLevelObject<Vector<double> >&,
//  				   const BlockVector<float>&);

//  template void MGTools::copy_from_mg(const MGDoFHandler<deal_II_dimension>&,
//  				    BlockVector<double>&,
//  				    const MGLevelObject<Vector<float> >&);

//  template void MGTools::copy_from_mg_add(const MGDoFHandler<deal_II_dimension>&,
//  					BlockVector<double>&,
//  					const MGLevelObject<Vector<float> >&);

template void MGTools::count_dofs_per_component (const MGDoFHandler<deal_II_dimension>&,
						 std::vector<std::vector<unsigned int> >&);
