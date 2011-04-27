//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <base/multithread_info.h>
#include <base/thread_management.h>
#include <base/quadrature_lib.h>
#include <base/table.h>
#include <base/template_constraints.h>
#include <base/utilities.h>
#include <lac/sparsity_pattern.h>
#include <lac/compressed_sparsity_pattern.h>
#include <lac/compressed_set_sparsity_pattern.h>
#include <lac/compressed_simple_sparsity_pattern.h>
#include <lac/trilinos_sparsity_pattern.h>
#include <lac/block_sparsity_pattern.h>
#include <lac/vector.h>
#include <lac/constraint_matrix.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <grid/intergrid_map.h>
#include <grid/grid_tools.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <fe/fe.h>
#include <fe/fe_values.h>
#include <fe/fe_tools.h>
#include <hp/fe_collection.h>
#include <dofs/dof_tools.h>
#include <numerics/vectors.h>

#include <multigrid/mg_dof_handler.h>
#include <multigrid/mg_dof_accessor.h>

#include <algorithm>
#include <numeric>

DEAL_II_NAMESPACE_OPEN



template <class DH, class SparsityPattern>
void
DoFTools::make_sparsity_pattern (const DH               &dof,
				 SparsityPattern        &sparsity,
				 const ConstraintMatrix &constraints,
				 const bool              keep_constrained_dofs,
				 const types::subdomain_id_t subdomain_id)
{
  const unsigned int n_dofs = dof.n_dofs();

  Assert (sparsity.n_rows() == n_dofs,
	  ExcDimensionMismatch (sparsity.n_rows(), n_dofs));
  Assert (sparsity.n_cols() == n_dofs,
	  ExcDimensionMismatch (sparsity.n_cols(), n_dofs));

  std::vector<unsigned int> dofs_on_this_cell;
  dofs_on_this_cell.reserve (max_dofs_per_cell(dof));
  typename DH::active_cell_iterator cell = dof.begin_active(),
				    endc = dof.end();

				   // In case we work with a distributed
				   // sparsity pattern of Trilinos type, we
				   // only have to do the work if the
				   // current cell is owned by the calling
				   // processor. Otherwise, just continue.
  for (; cell!=endc; ++cell)
    if (((subdomain_id == types::invalid_subdomain_id)
	 ||
	 (subdomain_id == cell->subdomain_id()))
	&&
	!cell->is_artificial()
	&&
	!cell->is_ghost())
      {
	const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
	dofs_on_this_cell.resize (dofs_per_cell);
	cell->get_dof_indices (dofs_on_this_cell);

					 // make sparsity pattern for this
					 // cell. if no constraints pattern was
					 // given, then the following call acts
					 // as if simply no constraints existed
	constraints.add_entries_local_to_global (dofs_on_this_cell,
						 sparsity,
						 keep_constrained_dofs);
      }
}



template <class DH, class SparsityPattern>
void
DoFTools::make_sparsity_pattern (const DH                &dof,
				 const Table<2,Coupling> &couplings,
				 SparsityPattern         &sparsity,
				 const ConstraintMatrix  &constraints,
				 const bool               keep_constrained_dofs,
				 const types::subdomain_id_t subdomain_id)
{
  const unsigned int n_dofs = dof.n_dofs();

  Assert (sparsity.n_rows() == n_dofs,
	  ExcDimensionMismatch (sparsity.n_rows(), n_dofs));
  Assert (sparsity.n_cols() == n_dofs,
	  ExcDimensionMismatch (sparsity.n_cols(), n_dofs));
  Assert (couplings.n_rows() == dof.get_fe().n_components(),
	  ExcDimensionMismatch(couplings.n_rows(), dof.get_fe().n_components()));
  Assert (couplings.n_cols() == dof.get_fe().n_components(),
	  ExcDimensionMismatch(couplings.n_cols(), dof.get_fe().n_components()));

  const hp::FECollection<DH::dimension,DH::space_dimension> fe_collection (dof.get_fe());

				   // first, for each finite element, build a
				   // mask for each dof, not like the one
				   // given which represents components. make
				   // sure we do the right thing also with
				   // respect to non-primitive shape
				   // functions, which takes some additional
				   // thought
  std::vector<Table<2,bool> > dof_mask(fe_collection.size());

				   // check whether the table of couplings
				   // contains only true arguments, i.e., we
				   // do not exclude any index. that is the
				   // easy case, since we don't have to set
				   // up the tables
  bool need_dof_mask = false;
  for (unsigned int i=0; i<couplings.n_rows(); ++i)
    for (unsigned int j=0; j<couplings.n_cols(); ++j)
      if (couplings(i,j) == none)
	need_dof_mask = true;

  if (need_dof_mask == true)
    for (unsigned int f=0; f<fe_collection.size(); ++f)
      {
	const unsigned int dofs_per_cell = fe_collection[f].dofs_per_cell;

	dof_mask[f].reinit (dofs_per_cell, dofs_per_cell);

	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  for (unsigned int j=0; j<dofs_per_cell; ++j)
	    if (fe_collection[f].is_primitive(i) &&
		fe_collection[f].is_primitive(j))
	      dof_mask[f](i,j)
		= (couplings(fe_collection[f].system_to_component_index(i).first,
			     fe_collection[f].system_to_component_index(j).first) != none);
	    else
	      {
		const unsigned int first_nonzero_comp_i
		  = (std::find (fe_collection[f].get_nonzero_components(i).begin(),
				fe_collection[f].get_nonzero_components(i).end(),
				true)
		     -
		     fe_collection[f].get_nonzero_components(i).begin());
		const unsigned int first_nonzero_comp_j
		  = (std::find (fe_collection[f].get_nonzero_components(j).begin(),
				fe_collection[f].get_nonzero_components(j).end(),
				true)
		     -
		     fe_collection[f].get_nonzero_components(j).begin());
		Assert (first_nonzero_comp_i < fe_collection[f].n_components(),
			ExcInternalError());
		Assert (first_nonzero_comp_j < fe_collection[f].n_components(),
			ExcInternalError());

		dof_mask[f](i,j)
		  = (couplings(first_nonzero_comp_i,first_nonzero_comp_j) != none);
	      }
      }


  std::vector<unsigned int> dofs_on_this_cell(fe_collection.max_dofs_per_cell());
  typename DH::active_cell_iterator cell = dof.begin_active(),
				    endc = dof.end();

				   // In case we work with a distributed
				   // sparsity pattern of Trilinos type, we
				   // only have to do the work if the
				   // current cell is owned by the calling
				   // processor. Otherwise, just continue.
  for (; cell!=endc; ++cell)
    if (((subdomain_id == types::invalid_subdomain_id)
	 ||
	 (subdomain_id == cell->subdomain_id()))
	&&
	!cell->is_artificial()
	&&
	!cell->is_ghost())
      {
	const unsigned int fe_index = cell->active_fe_index();
	const unsigned int dofs_per_cell =fe_collection[fe_index].dofs_per_cell;

	dofs_on_this_cell.resize (dofs_per_cell);
	cell->get_dof_indices (dofs_on_this_cell);


					 // make sparsity pattern for this
					 // cell. if no constraints pattern was
					 // given, then the following call acts
					 // as if simply no constraints existed
	constraints.add_entries_local_to_global (dofs_on_this_cell,
						 sparsity,
						 keep_constrained_dofs,
						 dof_mask[fe_index]);
      }
}



template <class DH, class SparsityPattern>
void
DoFTools::make_sparsity_pattern (
  const DH        &dof_row,
  const DH        &dof_col,
  SparsityPattern &sparsity)
{
  const unsigned int n_dofs_row = dof_row.n_dofs();
  const unsigned int n_dofs_col = dof_col.n_dofs();

  Assert (sparsity.n_rows() == n_dofs_row,
          ExcDimensionMismatch (sparsity.n_rows(), n_dofs_row));
  Assert (sparsity.n_cols() == n_dofs_col,
          ExcDimensionMismatch (sparsity.n_cols(), n_dofs_col));


  const std::list<std::pair<typename DH::cell_iterator,
    typename DH::cell_iterator> >
    cell_list
    = GridTools::get_finest_common_cells (dof_row, dof_col);


  typename std::list<std::pair<typename DH::cell_iterator,
    typename DH::cell_iterator> >
    ::const_iterator
    cell_iter = cell_list.begin();

  for (; cell_iter!=cell_list.end(); ++cell_iter)
    {
      const typename DH::cell_iterator cell_row = cell_iter->first;
      const typename DH::cell_iterator cell_col = cell_iter->second;

      if (!cell_row->has_children() && !cell_col->has_children())
        {
          const unsigned int dofs_per_cell_row =
	    cell_row->get_fe().dofs_per_cell;
          const unsigned int dofs_per_cell_col =
	    cell_col->get_fe().dofs_per_cell;
          std::vector<unsigned int>
	    local_dof_indices_row(dofs_per_cell_row);
          std::vector<unsigned int>
	    local_dof_indices_col(dofs_per_cell_col);
          cell_row->get_dof_indices (local_dof_indices_row);
          cell_col->get_dof_indices (local_dof_indices_col);
          for (unsigned int i=0; i<dofs_per_cell_row; ++i)
	    sparsity.add_entries (local_dof_indices_row[i],
				  local_dof_indices_col.begin(),
				  local_dof_indices_col.end());
        }
      else if (cell_row->has_children())
        {
          const std::vector<typename DH::active_cell_iterator >
	    child_cells = GridTools::get_active_child_cells<DH> (cell_row);
          for (unsigned int i=0; i<child_cells.size(); i++)
            {
	      const typename DH::active_cell_iterator
		cell_row_child = child_cells[i];
              const unsigned int dofs_per_cell_row =
		cell_row_child->get_fe().dofs_per_cell;
              const unsigned int dofs_per_cell_col =
		cell_col->get_fe().dofs_per_cell;
              std::vector<unsigned int>
		local_dof_indices_row(dofs_per_cell_row);
              std::vector<unsigned int>
		local_dof_indices_col(dofs_per_cell_col);
              cell_row_child->get_dof_indices (local_dof_indices_row);
              cell_col->get_dof_indices (local_dof_indices_col);
              for (unsigned int i=0; i<dofs_per_cell_row; ++i)
		sparsity.add_entries (local_dof_indices_row[i],
				      local_dof_indices_col.begin(),
				      local_dof_indices_col.end());
            }
        }
      else
        {
          std::vector<typename DH::active_cell_iterator>
	    child_cells = GridTools::get_active_child_cells<DH> (cell_col);
          for (unsigned int i=0; i<child_cells.size(); i++)
            {
	      const typename DH::active_cell_iterator
		cell_col_child = child_cells[i];
              const unsigned int dofs_per_cell_row =
		cell_row->get_fe().dofs_per_cell;
              const unsigned int dofs_per_cell_col =
		cell_col_child->get_fe().dofs_per_cell;
              std::vector<unsigned int>
		local_dof_indices_row(dofs_per_cell_row);
              std::vector<unsigned int>
		local_dof_indices_col(dofs_per_cell_col);
              cell_row->get_dof_indices (local_dof_indices_row);
              cell_col_child->get_dof_indices (local_dof_indices_col);
              for (unsigned int i=0; i<dofs_per_cell_row; ++i)
		sparsity.add_entries (local_dof_indices_row[i],
				      local_dof_indices_col.begin(),
				      local_dof_indices_col.end());
            }
        }
    }
}



template <class DH, class SparsityPattern>
void
DoFTools::make_boundary_sparsity_pattern (
  const DH                        &dof,
  const std::vector<unsigned int> &dof_to_boundary_mapping,
  SparsityPattern                 &sparsity)
{
  if (DH::dimension == 1)
    {
				   // there are only 2 boundary
				   // indicators in 1d, so it is no
				   // performance problem to call the
				   // other function
      typename DH::FunctionMap boundary_indicators;
      boundary_indicators[0] = 0;
      boundary_indicators[1] = 0;
      make_boundary_sparsity_pattern<DH, SparsityPattern> (dof,
							   boundary_indicators,
							   dof_to_boundary_mapping,
							   sparsity);
      return;
    }

  const unsigned int n_dofs = dof.n_dofs();

  AssertDimension (dof_to_boundary_mapping.size(), n_dofs);
  AssertDimension (sparsity.n_rows(), dof.n_boundary_dofs());
  AssertDimension (sparsity.n_cols(), dof.n_boundary_dofs());
#ifdef DEBUG
  if (sparsity.n_rows() != 0)
    {
      unsigned int max_element = 0;
      for (std::vector<unsigned int>::const_iterator i=dof_to_boundary_mapping.begin();
	   i!=dof_to_boundary_mapping.end(); ++i)
	if ((*i != DH::invalid_dof_index) &&
	    (*i > max_element))
	  max_element = *i;
      AssertDimension (max_element, sparsity.n_rows()-1);
    };
#endif

  std::vector<unsigned int> dofs_on_this_face;
  dofs_on_this_face.reserve (max_dofs_per_face(dof));

				   // loop over all faces to check
				   // whether they are at a
				   // boundary. note that we need not
				   // take special care of single
				   // lines (using
				   // @p{cell->has_boundary_lines}),
				   // since we do not support
				   // boundaries of dimension dim-2,
				   // and so every boundary line is
				   // also part of a boundary face.
  typename DH::active_cell_iterator cell = dof.begin_active(),
				    endc = dof.end();
  for (; cell!=endc; ++cell)
    for (unsigned int f=0; f<GeometryInfo<DH::dimension>::faces_per_cell; ++f)
      if (cell->at_boundary(f))
        {
          const unsigned int dofs_per_face = cell->get_fe().dofs_per_face;
          dofs_on_this_face.resize (dofs_per_face);
          cell->face(f)->get_dof_indices (dofs_on_this_face,
					  cell->active_fe_index());

                                           // make sparsity pattern for this cell
          for (unsigned int i=0; i<dofs_per_face; ++i)
            for (unsigned int j=0; j<dofs_per_face; ++j)
              sparsity.add (dof_to_boundary_mapping[dofs_on_this_face[i]],
                            dof_to_boundary_mapping[dofs_on_this_face[j]]);
	}
}



template <class DH, class SparsityPattern>
void DoFTools::make_boundary_sparsity_pattern (
  const DH                                        &dof,
  const typename FunctionMap<DH::space_dimension>::type &boundary_indicators,
  const std::vector<unsigned int>                 &dof_to_boundary_mapping,
  SparsityPattern                                 &sparsity)
{
  if (DH::dimension == 1)
    {
				   // first check left, then right
				   // boundary point
      for (unsigned int direction=0; direction<2; ++direction)
	{
				       // if this boundary is not
				       // requested, then go on with next one
	  if (boundary_indicators.find(direction) ==
	      boundary_indicators.end())
	    continue;

				       // find active cell at that
				       // boundary: first go to
				       // left/right, then to children
	  typename DH::cell_iterator cell = dof.begin(0);
	  while (!cell->at_boundary(direction))
	    cell = cell->neighbor(direction);
	  while (!cell->active())
	    cell = cell->child(direction);

	  const unsigned int dofs_per_vertex = cell->get_fe().dofs_per_vertex;
	  std::vector<unsigned int> boundary_dof_boundary_indices (dofs_per_vertex);

				       // next get boundary mapped dof
				       // indices of boundary dofs
	  for (unsigned int i=0; i<dofs_per_vertex; ++i)
	    boundary_dof_boundary_indices[i]
	      = dof_to_boundary_mapping[cell->vertex_dof_index(direction,i)];

	  for (unsigned int i=0; i<dofs_per_vertex; ++i)
	    sparsity.add_entries (boundary_dof_boundary_indices[i],
				  boundary_dof_boundary_indices.begin(),
				  boundary_dof_boundary_indices.end());
	};
      return;
    }

  const unsigned int n_dofs = dof.n_dofs();

  AssertDimension (dof_to_boundary_mapping.size(), n_dofs);
  Assert (boundary_indicators.find(255) == boundary_indicators.end(),
	  typename DH::ExcInvalidBoundaryIndicator());
  Assert (sparsity.n_rows() == dof.n_boundary_dofs (boundary_indicators),
	  ExcDimensionMismatch (sparsity.n_rows(), dof.n_boundary_dofs (boundary_indicators)));
  Assert (sparsity.n_cols() == dof.n_boundary_dofs (boundary_indicators),
	  ExcDimensionMismatch (sparsity.n_cols(), dof.n_boundary_dofs (boundary_indicators)));
#ifdef DEBUG
  if (sparsity.n_rows() != 0)
    {
      unsigned int max_element = 0;
      for (std::vector<unsigned int>::const_iterator i=dof_to_boundary_mapping.begin();
	   i!=dof_to_boundary_mapping.end(); ++i)
	if ((*i != DH::invalid_dof_index) &&
	    (*i > max_element))
	  max_element = *i;
      AssertDimension (max_element, sparsity.n_rows()-1);
    };
#endif

  std::vector<unsigned int> dofs_on_this_face;
  dofs_on_this_face.reserve (max_dofs_per_face(dof));
  typename DH::active_cell_iterator cell = dof.begin_active(),
				    endc = dof.end();
  for (; cell!=endc; ++cell)
    for (unsigned int f=0; f<GeometryInfo<DH::dimension>::faces_per_cell; ++f)
      if (boundary_indicators.find(cell->face(f)->boundary_indicator()) !=
          boundary_indicators.end())
        {
          const unsigned int dofs_per_face = cell->get_fe().dofs_per_face;
          dofs_on_this_face.resize (dofs_per_face);
          cell->face(f)->get_dof_indices (dofs_on_this_face,
					  cell->active_fe_index());

                                           // make sparsity pattern for this cell
          for (unsigned int i=0; i<dofs_per_face; ++i)
            for (unsigned int j=0; j<dofs_per_face; ++j)
              sparsity.add (dof_to_boundary_mapping[dofs_on_this_face[i]],
                            dof_to_boundary_mapping[dofs_on_this_face[j]]);
        }
}



template <class DH, class SparsityPattern>
void
DoFTools::make_flux_sparsity_pattern (const DH        &dof,
				      SparsityPattern &sparsity,
				      const ConstraintMatrix &constraints,
				      const bool              keep_constrained_dofs,
				      const types::subdomain_id_t subdomain_id)
{
  const unsigned int n_dofs = dof.n_dofs();

  AssertDimension (sparsity.n_rows(), n_dofs);
  AssertDimension (sparsity.n_cols(), n_dofs);

  std::vector<unsigned int> dofs_on_this_cell;
  std::vector<unsigned int> dofs_on_other_cell;
  dofs_on_this_cell.reserve (max_dofs_per_cell(dof));
  dofs_on_other_cell.reserve (max_dofs_per_cell(dof));
  typename DH::active_cell_iterator cell = dof.begin_active(),
    endc = dof.end();

				// TODO: in an old implementation, we used
				// user flags before to tag faces that were
				// already touched. this way, we could reduce
				// the work a little bit. now, we instead add
				// only data from one side. this should be OK,
				// but we need to actually verify it.

				// In case we work with a distributed
				// sparsity pattern of Trilinos type, we
				// only have to do the work if the
				// current cell is owned by the calling
				// processor. Otherwise, just continue.
  for (; cell!=endc; ++cell)
    if ((subdomain_id == types::invalid_subdomain_id)
        ||
        (subdomain_id == cell->subdomain_id()))
    {
      const unsigned int n_dofs_on_this_cell = cell->get_fe().dofs_per_cell;
      dofs_on_this_cell.resize (n_dofs_on_this_cell);
      cell->get_dof_indices (dofs_on_this_cell);

				// make sparsity pattern for this
				// cell. if no constraints pattern was
				// given, then the following call acts
				// as if simply no constraints existed
      constraints.add_entries_local_to_global (dofs_on_this_cell,
					       sparsity,
					       keep_constrained_dofs);

      for (unsigned int face = 0;
	   face < GeometryInfo<DH::dimension>::faces_per_cell;
	   ++face)
	{
	  typename DH::face_iterator cell_face = cell->face(face);
	  if (! cell->at_boundary(face) )
	    {
	      typename DH::cell_iterator neighbor = cell->neighbor(face);

	      if (cell_face->has_children())
		{
		  for (unsigned int sub_nr = 0;
		       sub_nr != cell_face->number_of_children();
		       ++sub_nr)
		    {
		      const typename DH::cell_iterator
			sub_neighbor
			= cell->neighbor_child_on_subface (face, sub_nr);

		      const unsigned int n_dofs_on_neighbor
			= sub_neighbor->get_fe().dofs_per_cell;
		      dofs_on_other_cell.resize (n_dofs_on_neighbor);
		      sub_neighbor->get_dof_indices (dofs_on_other_cell);

		      constraints.add_entries_local_to_global
			(dofs_on_this_cell, dofs_on_other_cell,
			 sparsity, keep_constrained_dofs);
		      constraints.add_entries_local_to_global
			(dofs_on_other_cell, dofs_on_this_cell,
			 sparsity, keep_constrained_dofs);
		    }
		}
	      else
		{
				// Refinement edges are taken care of
				// by coarser cells

				// TODO: in the distributed case, we miss out
				// the constraints when the neighbor cell is
				// coarser, but only the current cell is owned
				// locally!
		  if (cell->neighbor_is_coarser(face))
		    continue;

		  const unsigned int n_dofs_on_neighbor
		    = neighbor->get_fe().dofs_per_cell;
		  dofs_on_other_cell.resize (n_dofs_on_neighbor);

		  neighbor->get_dof_indices (dofs_on_other_cell);

		  constraints.add_entries_local_to_global
		    (dofs_on_this_cell, dofs_on_other_cell,
		     sparsity, keep_constrained_dofs);

				// only need to add these in case the neighbor
				// cell is not locally owned - otherwise, we
				// touch each face twice and hence put the
				// indices the other way around
		  if (cell->neighbor(face)->subdomain_id() !=
		      cell->subdomain_id())
		    constraints.add_entries_local_to_global
		      (dofs_on_other_cell, dofs_on_this_cell,
		       sparsity, keep_constrained_dofs);
		}
	    }
	}
    }
}


template <class DH, class SparsityPattern>
void
DoFTools::make_flux_sparsity_pattern (
  const DH        &dof,
  SparsityPattern &sparsity)
{
  ConstraintMatrix constraints;
  make_flux_sparsity_pattern (dof, sparsity, constraints);
}



template <int dim, int spacedim>
Table<2,DoFTools::Coupling>
DoFTools::dof_couplings_from_component_couplings
(const FiniteElement<dim,spacedim> &fe,
 const Table<2,Coupling> &component_couplings)
{
  Assert(component_couplings.n_rows() == fe.n_components(),
	 ExcDimensionMismatch(component_couplings.n_rows(),
			      fe.n_components()));
  Assert(component_couplings.n_cols() == fe.n_components(),
	 ExcDimensionMismatch(component_couplings.n_cols(),
			      fe.n_components()));

  const unsigned int n_dofs = fe.dofs_per_cell;

  Table<2,DoFTools::Coupling> dof_couplings (n_dofs, n_dofs);

  for (unsigned int i=0; i<n_dofs; ++i)
    {
      const unsigned int ii
	= (fe.is_primitive(i) ?
	   fe.system_to_component_index(i).first
	   :
	   (std::find (fe.get_nonzero_components(i).begin(),
		       fe.get_nonzero_components(i).end(),
		       true)
	    -
	    fe.get_nonzero_components(i).begin())
	);
      Assert (ii < fe.n_components(), ExcInternalError());

      for (unsigned int j=0; j<n_dofs; ++j)
	{
	  const unsigned int jj
	    = (fe.is_primitive(j) ?
	       fe.system_to_component_index(j).first
	       :
	       (std::find (fe.get_nonzero_components(j).begin(),
			   fe.get_nonzero_components(j).end(),
			   true)
		-
		fe.get_nonzero_components(j).begin())
	    );
	  Assert (jj < fe.n_components(), ExcInternalError());

	  dof_couplings(i,j) = component_couplings(ii,jj);
	}
    }
  return dof_couplings;
}


template <int dim, int spacedim>
std::vector<Table<2,DoFTools::Coupling> >
DoFTools::dof_couplings_from_component_couplings
(const hp::FECollection<dim,spacedim> &fe,
 const Table<2,Coupling> &component_couplings)
{
  std::vector<Table<2,DoFTools::Coupling> > return_value (fe.size());
  for (unsigned int i=0; i<fe.size(); ++i)
    return_value[i]
      = dof_couplings_from_component_couplings(fe[i], component_couplings);

  return return_value;
}

  


namespace internal
{
  namespace DoFTools
  {
				     // implementation of the same function in
				     // namespace DoFTools for non-hp
				     // DoFHandlers
    template <class DH, class SparsityPattern>
    void
    make_flux_sparsity_pattern (const DH                &dof,
				SparsityPattern         &sparsity,
				const Table<2,dealii::DoFTools::Coupling> &int_mask,
				const Table<2,dealii::DoFTools::Coupling> &flux_mask)
    {
      const FiniteElement<DH::dimension> &fe = dof.get_fe();

      std::vector<unsigned int> dofs_on_this_cell(fe.dofs_per_cell);
      std::vector<unsigned int> dofs_on_other_cell(fe.dofs_per_cell);

      const Table<2,dealii::DoFTools::Coupling>
	int_dof_mask  = dealii::DoFTools::dof_couplings_from_component_couplings(fe, int_mask),
	flux_dof_mask = dealii::DoFTools::dof_couplings_from_component_couplings(fe, flux_mask);

      Table<2,bool> support_on_face(fe.dofs_per_cell,
				    GeometryInfo<DH::dimension>::faces_per_cell);
      for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
	for (unsigned int f=0; f<GeometryInfo<DH::dimension>::faces_per_cell;++f)
	  support_on_face(i,f) = fe.has_support_on_face(i,f);

      typename DH::active_cell_iterator cell = dof.begin_active(),
					endc = dof.end();
      for (; cell!=endc; ++cell)
	{
	  cell->get_dof_indices (dofs_on_this_cell);
					   // make sparsity pattern for this cell
	  for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
	    for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
	      if (int_dof_mask(i,j) != dealii::DoFTools::none)
		sparsity.add (dofs_on_this_cell[i],
			      dofs_on_this_cell[j]);

					   // Loop over all interior neighbors
	  for (unsigned int face = 0;
	       face < GeometryInfo<DH::dimension>::faces_per_cell;
	       ++face)
	    {
	      const typename DH::face_iterator
		cell_face = cell->face(face);
	      if (cell_face->user_flag_set ())
		continue;

	      if (cell->at_boundary (face) )
		{
		  for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
		    {
		      const bool i_non_zero_i = support_on_face (i, face);
		      for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
			{
			  const bool j_non_zero_i = support_on_face (j, face);

			  if ((flux_dof_mask(i,j) == dealii::DoFTools::always)
			      ||
			      (flux_dof_mask(i,j) == dealii::DoFTools::nonzero
			       &&
			       i_non_zero_i
			       &&
			       j_non_zero_i))
			    sparsity.add (dofs_on_this_cell[i],
					  dofs_on_this_cell[j]);
			}
		    }
		}
	      else
		{
		  typename DH::cell_iterator
		    neighbor = cell->neighbor(face);
						   // Refinement edges are taken care of
						   // by coarser cells
		  if (cell->neighbor_is_coarser(face))
		    continue;

		  typename DH::face_iterator cell_face = cell->face(face);
		  const unsigned int
		    neighbor_face = cell->neighbor_of_neighbor(face);

		  if (cell_face->has_children())
		    {
		      for (unsigned int sub_nr = 0;
			   sub_nr != cell_face->n_children();
			   ++sub_nr)
			{
			  const typename DH::cell_iterator
			    sub_neighbor
			    = cell->neighbor_child_on_subface (face, sub_nr);

			  sub_neighbor->get_dof_indices (dofs_on_other_cell);
			  for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
			    {
			      const bool i_non_zero_i = support_on_face (i, face);
			      const bool i_non_zero_e = support_on_face (i, neighbor_face);
			      for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
				{
				  const bool j_non_zero_i = support_on_face (j, face);
				  const bool j_non_zero_e = support_on_face (j, neighbor_face);

				  if (flux_dof_mask(i,j) == dealii::DoFTools::always)
				    {
				      sparsity.add (dofs_on_this_cell[i],
						    dofs_on_other_cell[j]);
				      sparsity.add (dofs_on_other_cell[i],
						    dofs_on_this_cell[j]);
				      sparsity.add (dofs_on_this_cell[i],
						    dofs_on_this_cell[j]);
				      sparsity.add (dofs_on_other_cell[i],
						    dofs_on_other_cell[j]);
				    }
				  else if (flux_dof_mask(i,j) == dealii::DoFTools::nonzero)
				    {
				      if (i_non_zero_i && j_non_zero_e)
					sparsity.add (dofs_on_this_cell[i],
						      dofs_on_other_cell[j]);
				      if (i_non_zero_e && j_non_zero_i)
					sparsity.add (dofs_on_other_cell[i],
						      dofs_on_this_cell[j]);
				      if (i_non_zero_i && j_non_zero_i)
					sparsity.add (dofs_on_this_cell[i],
						      dofs_on_this_cell[j]);
				      if (i_non_zero_e && j_non_zero_e)
					sparsity.add (dofs_on_other_cell[i],
						      dofs_on_other_cell[j]);
				    }

				  if (flux_dof_mask(j,i) == dealii::DoFTools::always)
				    {
				      sparsity.add (dofs_on_this_cell[j],
						    dofs_on_other_cell[i]);
				      sparsity.add (dofs_on_other_cell[j],
						    dofs_on_this_cell[i]);
				      sparsity.add (dofs_on_this_cell[j],
						    dofs_on_this_cell[i]);
				      sparsity.add (dofs_on_other_cell[j],
						    dofs_on_other_cell[i]);
				    }
				  else if (flux_dof_mask(j,i) == dealii::DoFTools::nonzero)
				    {
				      if (j_non_zero_i && i_non_zero_e)
					sparsity.add (dofs_on_this_cell[j],
						      dofs_on_other_cell[i]);
				      if (j_non_zero_e && i_non_zero_i)
					sparsity.add (dofs_on_other_cell[j],
						      dofs_on_this_cell[i]);
				      if (j_non_zero_i && i_non_zero_i)
					sparsity.add (dofs_on_this_cell[j],
						      dofs_on_this_cell[i]);
				      if (j_non_zero_e && i_non_zero_e)
					sparsity.add (dofs_on_other_cell[j],
						      dofs_on_other_cell[i]);
				    }
				}
			    }
			  sub_neighbor->face(neighbor_face)->set_user_flag ();
			}
		    }
		  else
		    {
		      neighbor->get_dof_indices (dofs_on_other_cell);
		      for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
			{
			  const bool i_non_zero_i = support_on_face (i, face);
			  const bool i_non_zero_e = support_on_face (i, neighbor_face);
			  for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
			    {
			      const bool j_non_zero_i = support_on_face (j, face);
			      const bool j_non_zero_e = support_on_face (j, neighbor_face);
			      if (flux_dof_mask(i,j) == dealii::DoFTools::always)
				{
				  sparsity.add (dofs_on_this_cell[i],
						dofs_on_other_cell[j]);
				  sparsity.add (dofs_on_other_cell[i],
						dofs_on_this_cell[j]);
				  sparsity.add (dofs_on_this_cell[i],
						dofs_on_this_cell[j]);
				  sparsity.add (dofs_on_other_cell[i],
						dofs_on_other_cell[j]);
				}
			      if (flux_dof_mask(i,j) == dealii::DoFTools::nonzero)
				{
				  if (i_non_zero_i && j_non_zero_e)
				    sparsity.add (dofs_on_this_cell[i],
						  dofs_on_other_cell[j]);
				  if (i_non_zero_e && j_non_zero_i)
				    sparsity.add (dofs_on_other_cell[i],
						  dofs_on_this_cell[j]);
				  if (i_non_zero_i && j_non_zero_i)
				    sparsity.add (dofs_on_this_cell[i],
						  dofs_on_this_cell[j]);
				  if (i_non_zero_e && j_non_zero_e)
				    sparsity.add (dofs_on_other_cell[i],
						  dofs_on_other_cell[j]);
				}

			      if (flux_dof_mask(j,i) == dealii::DoFTools::always)
				{
				  sparsity.add (dofs_on_this_cell[j],
						dofs_on_other_cell[i]);
				  sparsity.add (dofs_on_other_cell[j],
						dofs_on_this_cell[i]);
				  sparsity.add (dofs_on_this_cell[j],
						dofs_on_this_cell[i]);
				  sparsity.add (dofs_on_other_cell[j],
						dofs_on_other_cell[i]);
				}
			      if (flux_dof_mask(j,i) == dealii::DoFTools::nonzero)
				{
				  if (j_non_zero_i && i_non_zero_e)
				    sparsity.add (dofs_on_this_cell[j],
						  dofs_on_other_cell[i]);
				  if (j_non_zero_e && i_non_zero_i)
				    sparsity.add (dofs_on_other_cell[j],
						  dofs_on_this_cell[i]);
				  if (j_non_zero_i && i_non_zero_i)
				    sparsity.add (dofs_on_this_cell[j],
						  dofs_on_this_cell[i]);
				  if (j_non_zero_e && i_non_zero_e)
				    sparsity.add (dofs_on_other_cell[j],
						  dofs_on_other_cell[i]);
				}
			    }
			}
		      neighbor->face(neighbor_face)->set_user_flag ();
		    }
		}
	    }
	}
    }


				     // implementation of the same function in
				     // namespace DoFTools for non-hp
				     // DoFHandlers
    template <int dim, int spacedim, class SparsityPattern>
    void
    make_flux_sparsity_pattern (const dealii::hp::DoFHandler<dim,spacedim> &dof,
				SparsityPattern                           &sparsity,
				const Table<2,dealii::DoFTools::Coupling> &int_mask,
				const Table<2,dealii::DoFTools::Coupling> &flux_mask)
    {
				       // while the implementation above is
				       // quite optimized and caches a lot of
				       // data (see e.g. the int/flux_dof_mask
				       // tables), this is no longer practical
				       // for the hp version since we would
				       // have to have it for all combinations
				       // of elements in the
				       // hp::FECollection. consequently, the
				       // implementation here is simpler and
				       // probably less efficient but at least
				       // readable...

      const dealii::hp::FECollection<dim,spacedim> &fe = dof.get_fe();

      std::vector<unsigned int> dofs_on_this_cell(dealii::DoFTools::max_dofs_per_cell(dof));
      std::vector<unsigned int> dofs_on_other_cell(dealii::DoFTools::max_dofs_per_cell(dof));

      const std::vector<Table<2,dealii::DoFTools::Coupling> >
	int_dof_mask
	= dealii::DoFTools::dof_couplings_from_component_couplings(fe, int_mask);

      typename dealii::hp::DoFHandler<dim,spacedim>::active_cell_iterator
	cell = dof.begin_active(),
	endc = dof.end();
      for (; cell!=endc; ++cell)
	{
	  dofs_on_this_cell.resize (cell->get_fe().dofs_per_cell);
	  cell->get_dof_indices (dofs_on_this_cell);
					   // make sparsity pattern for this cell
	  for (unsigned int i=0; i<cell->get_fe().dofs_per_cell; ++i)
	    for (unsigned int j=0; j<cell->get_fe().dofs_per_cell; ++j)
	      if (int_dof_mask[cell->active_fe_index()](i,j) != dealii::DoFTools::none)
		sparsity.add (dofs_on_this_cell[i],
			      dofs_on_this_cell[j]);

					   // Loop over all interior neighbors
	  for (unsigned int face = 0;
	       face < GeometryInfo<dim>::faces_per_cell;
	       ++face)
	    {
	      const typename dealii::hp::DoFHandler<dim,spacedim>::face_iterator
		cell_face = cell->face(face);
	      if (cell_face->user_flag_set ())
		continue;

	      if (cell->at_boundary (face) )
		{
		  for (unsigned int i=0; i<cell->get_fe().dofs_per_cell; ++i)
		    for (unsigned int j=0; j<cell->get_fe().dofs_per_cell; ++j)
		      if ((flux_mask(cell->get_fe().system_to_component_index(i).first,
				     cell->get_fe().system_to_component_index(j).first)
			   == dealii::DoFTools::always)
			  ||
			  (flux_mask(cell->get_fe().system_to_component_index(i).first,
				     cell->get_fe().system_to_component_index(j).first)
			   == dealii::DoFTools::nonzero))
			sparsity.add (dofs_on_this_cell[i],
				      dofs_on_this_cell[j]);
		}
	      else
		{
		  typename dealii::hp::DoFHandler<dim,spacedim>::cell_iterator
		    neighbor = cell->neighbor(face);
						   // Refinement edges are taken care of
						   // by coarser cells
		  if (cell->neighbor_is_coarser(face))
		    continue;

		  typename dealii::hp::DoFHandler<dim,spacedim>::face_iterator
		    cell_face = cell->face(face);
		  const unsigned int
		    neighbor_face = cell->neighbor_of_neighbor(face);

		  if (cell_face->has_children())
		    {
		      for (unsigned int sub_nr = 0;
			   sub_nr != cell_face->n_children();
			   ++sub_nr)
			{
			  const typename dealii::hp::DoFHandler<dim,spacedim>::cell_iterator
			    sub_neighbor
			    = cell->neighbor_child_on_subface (face, sub_nr);

			  dofs_on_other_cell.resize (sub_neighbor->get_fe().dofs_per_cell);
			  sub_neighbor->get_dof_indices (dofs_on_other_cell);
			  for (unsigned int i=0; i<cell->get_fe().dofs_per_cell; ++i)
			    {
			      for (unsigned int j=0; j<sub_neighbor->get_fe().dofs_per_cell;
				   ++j)
				{
				  if ((flux_mask(cell->get_fe().system_to_component_index(i).first,
						 sub_neighbor->get_fe().system_to_component_index(j).first)
				       == dealii::DoFTools::always)
				      ||
				      (flux_mask(cell->get_fe().system_to_component_index(i).first,
						 sub_neighbor->get_fe().system_to_component_index(j).first)
				       == dealii::DoFTools::nonzero))
				    {
				      sparsity.add (dofs_on_this_cell[i],
						    dofs_on_other_cell[j]);
				      sparsity.add (dofs_on_other_cell[i],
						    dofs_on_this_cell[j]);
				      sparsity.add (dofs_on_this_cell[i],
						    dofs_on_this_cell[j]);
				      sparsity.add (dofs_on_other_cell[i],
						    dofs_on_other_cell[j]);
				    }

				  if ((flux_mask(sub_neighbor->get_fe().system_to_component_index(j).first,
						 cell->get_fe().system_to_component_index(i).first)
				       == dealii::DoFTools::always)
				      ||
				      (flux_mask(sub_neighbor->get_fe().system_to_component_index(j).first,
						 cell->get_fe().system_to_component_index(i).first)
				       == dealii::DoFTools::nonzero))
				    {
				      sparsity.add (dofs_on_this_cell[j],
						    dofs_on_other_cell[i]);
				      sparsity.add (dofs_on_other_cell[j],
						    dofs_on_this_cell[i]);
				      sparsity.add (dofs_on_this_cell[j],
						    dofs_on_this_cell[i]);
				      sparsity.add (dofs_on_other_cell[j],
						    dofs_on_other_cell[i]);
				    }
				}
			    }
			  sub_neighbor->face(neighbor_face)->set_user_flag ();
			}
		    }
		  else
		    {
		      dofs_on_other_cell.resize (neighbor->get_fe().dofs_per_cell);
		      neighbor->get_dof_indices (dofs_on_other_cell);
		      for (unsigned int i=0; i<cell->get_fe().dofs_per_cell; ++i)
			{
			  for (unsigned int j=0; j<neighbor->get_fe().dofs_per_cell; ++j)
			    {
			      if ((flux_mask(cell->get_fe().system_to_component_index(i).first,
					     neighbor->get_fe().system_to_component_index(j).first)
				   == dealii::DoFTools::always)
				  ||
				  (flux_mask(cell->get_fe().system_to_component_index(i).first,
					     neighbor->get_fe().system_to_component_index(j).first)
				   == dealii::DoFTools::nonzero))
				{
				  sparsity.add (dofs_on_this_cell[i],
						dofs_on_other_cell[j]);
				  sparsity.add (dofs_on_other_cell[i],
						dofs_on_this_cell[j]);
				  sparsity.add (dofs_on_this_cell[i],
						dofs_on_this_cell[j]);
				  sparsity.add (dofs_on_other_cell[i],
						dofs_on_other_cell[j]);
				}

			      if ((flux_mask(neighbor->get_fe().system_to_component_index(j).first,
					     cell->get_fe().system_to_component_index(i).first)
				   == dealii::DoFTools::always)
				  ||
				  (flux_mask(neighbor->get_fe().system_to_component_index(j).first,
					     cell->get_fe().system_to_component_index(i).first)
				   == dealii::DoFTools::nonzero))
				{
				  sparsity.add (dofs_on_this_cell[j],
						dofs_on_other_cell[i]);
				  sparsity.add (dofs_on_other_cell[j],
						dofs_on_this_cell[i]);
				  sparsity.add (dofs_on_this_cell[j],
						dofs_on_this_cell[i]);
				  sparsity.add (dofs_on_other_cell[j],
						dofs_on_other_cell[i]);
				}
			    }
			}
		      neighbor->face(neighbor_face)->set_user_flag ();
		    }
		}
	    }
	}
    }
    
  }
}



template <class DH, class SparsityPattern>
void
DoFTools::
make_flux_sparsity_pattern (const DH                &dof,
			    SparsityPattern         &sparsity,
			    const Table<2,Coupling> &int_mask,
			    const Table<2,Coupling> &flux_mask)
{
				   // do the error checking and frame code
				   // here, and then pass on to more
				   // specialized functions in the internal
				   // namespace
  const unsigned int n_dofs = dof.n_dofs();
  const unsigned int n_comp = dof.get_fe().n_components();

  Assert (sparsity.n_rows() == n_dofs,
	  ExcDimensionMismatch (sparsity.n_rows(), n_dofs));
  Assert (sparsity.n_cols() == n_dofs,
	  ExcDimensionMismatch (sparsity.n_cols(), n_dofs));
  Assert (int_mask.n_rows() == n_comp,
	  ExcDimensionMismatch (int_mask.n_rows(), n_comp));
  Assert (int_mask.n_cols() == n_comp,
	  ExcDimensionMismatch (int_mask.n_cols(), n_comp));
  Assert (flux_mask.n_rows() == n_comp,
	  ExcDimensionMismatch (flux_mask.n_rows(), n_comp));
  Assert (flux_mask.n_cols() == n_comp,
	  ExcDimensionMismatch (flux_mask.n_cols(), n_comp));

  
				   // Clear user flags because we will
				   // need them. But first we save
				   // them and make sure that we
				   // restore them later such that at
				   // the end of this function the
				   // Triangulation will be in the
				   // same state as it was at the
				   // beginning of this function.
  std::vector<bool> user_flags;
  dof.get_tria().save_user_flags(user_flags);
  const_cast<Triangulation<DH::dimension> &>(dof.get_tria()).clear_user_flags ();

  internal::DoFTools::make_flux_sparsity_pattern (dof, sparsity,
						  int_mask, flux_mask);
  
				   // finally restore the user flags
  const_cast<Triangulation<DH::dimension> &>(dof.get_tria()).load_user_flags(user_flags);
}







namespace internal
{
  namespace DoFTools
  {
    namespace
    {
      inline
      bool
      check_master_dof_list (const FullMatrix<double> &face_interpolation_matrix,
			     const std::vector<unsigned int> &master_dof_list)
      {
	const unsigned int N = master_dof_list.size();

	FullMatrix<double> tmp (N,N);
	for (unsigned int i=0; i<N; ++i)
	  for (unsigned int j=0; j<N; ++j)
	    tmp(i,j) = face_interpolation_matrix (master_dof_list[i], j);

					 // then use the algorithm
					 // from
					 // FullMatrix::gauss_jordan
					 // on this matrix to find out
					 // whether it is
					 // singular. the algorithm
					 // there does piviting and at
					 // the end swaps rows back
					 // into their proper order --
					 // we omit this step here,
					 // since we don't care about
					 // the inverse matrix, all we
					 // care about is whether the
					 // matrix is regular or
					 // singular

					 // first get an estimate of the
					 // size of the elements of this
					 // matrix, for later checks whether
					 // the pivot element is large
					 // enough, or whether we have to
					 // fear that the matrix is not
					 // regular
	double diagonal_sum = 0;
	for (unsigned int i=0; i<N; ++i)
	  diagonal_sum += std::fabs(tmp(i,i));
	const double typical_diagonal_element = diagonal_sum/N;

					 // initialize the array that holds
					 // the permutations that we find
					 // during pivot search
	std::vector<unsigned int> p(N);
	for (unsigned int i=0; i<N; ++i)
	  p[i] = i;

	for (unsigned int j=0; j<N; ++j)
	  {
					     // pivot search: search that
					     // part of the line on and
					     // right of the diagonal for
					     // the largest element
	    double       max = std::fabs(tmp(j,j));
	    unsigned int r   = j;
	    for (unsigned int i=j+1; i<N; ++i)
	      {
		if (std::fabs(tmp(i,j)) > max)
		  {
		    max = std::fabs(tmp(i,j));
		    r = i;
		  }
	      }
					     // check whether the
					     // pivot is too small. if
					     // that is the case, then
					     // the matrix is singular
					     // and we shouldn't use
					     // this set of master
					     // dofs
	    if (max < 1.e-12*typical_diagonal_element)
	      return false;

					     // row interchange
	    if (r>j)
	      {
		for (unsigned int k=0; k<N; ++k)
		  std::swap (tmp(j,k), tmp(r,k));

		std::swap (p[j], p[r]);
	      }

					     // transformation
	    const double hr = 1./tmp(j,j);
	    tmp(j,j) = hr;
	    for (unsigned int k=0; k<N; ++k)
	      {
		if (k==j) continue;
		for (unsigned int i=0; i<N; ++i)
		  {
		    if (i==j) continue;
		    tmp(i,k) -= tmp(i,j)*tmp(j,k)*hr;
		  }
	      }
	    for (unsigned int i=0; i<N; ++i)
	      {
		tmp(i,j) *= hr;
		tmp(j,i) *= -hr;
	      }
	    tmp(j,j) = hr;
	  }

					 // everything went fine, so
					 // we can accept this set of
					 // master dofs (at least as
					 // far as they have already
					 // been collected)
	return true;
      }


				       /**
					* When restricting, on a face, the
					* degrees of freedom of fe1 to the
					* space described by fe2 (for example
					* for the complex case described in
					* the @ref hp_paper "hp paper"), we have to select
					* fe2.dofs_per_face out of the
					* fe1.dofs_per_face face DoFs as the
					* master DoFs, and the rest become
					* slave dofs. This function selects
					* which ones will be masters, and
					* which ones will be slaves.
					*
					* The function assumes that
					* master_dofs already has size
					* fe1.dofs_per_face. After the
					* function, exactly fe2.dofs_per_face
					* entries will be true.
					*
					* The function is a bit
					* complicated since it has to
					* figure out a set a DoFs so
					* that the corresponding rows
					* in the face interpolation
					* matrix are all linearly
					* independent. we have a good
					* heuristic (see the function
					* body) for selecting these
					* rows, but there are cases
					* where this fails and we have
					* to pick them
					* differently. what we do is
					* to run the heuristic and
					* then go back to determine
					* whether we have a set of
					* rows with full row rank. if
					* this isn't the case, go back
					* and select dofs differently
					*/
      template <int dim, int spacedim>
      void
      select_master_dofs_for_face_restriction (const FiniteElement<dim,spacedim> &fe1,
					       const FiniteElement<dim,spacedim> &fe2,
					       const FullMatrix<double> &face_interpolation_matrix,
					       std::vector<bool>        &master_dof_mask)
      {
	Assert (fe1.dofs_per_face >= fe2.dofs_per_face,
		ExcInternalError());
	AssertDimension (master_dof_mask.size(), fe1.dofs_per_face);

	Assert (fe2.dofs_per_vertex <= fe1.dofs_per_vertex,
		ExcInternalError());
	Assert (fe2.dofs_per_line <= fe1.dofs_per_line,
		ExcInternalError());
	Assert ((dim < 3)
		||
		(fe2.dofs_per_quad <= fe1.dofs_per_quad),
		ExcInternalError());

					 // the idea here is to designate as
					 // many DoFs in fe1 per object
					 // (vertex, line, quad) as master as
					 // there are such dofs in fe2
					 // (indices are int, because we want
					 // to avoid the 'unsigned int < 0 is
					 // always false warning for the cases
					 // at the bottom in 1d and 2d)
					 //
					 // as mentioned in the paper, it is
					 // not always easy to find a set of
					 // master dofs that produces an
					 // invertible matrix. to this end, we
					 // check in each step whether the
					 // matrix is still invertible and
					 // simply discard this dof if the
					 // matrix is not invertible anymore.
					 //
					 // the cases where we did have
					 // trouble in the past were with
					 // adding more quad dofs when Q3 and
					 // Q4 elements meet at a refined face
					 // in 3d (see the hp/crash_12 test
					 // that tests that we can do exactly
					 // this, and failed before we had
					 // code to compensate for this
					 // case). the other case are system
					 // elements: if we have say a Q1Q2 vs
					 // a Q2Q3 element, then we can't just
					 // take all master dofs on a line
					 // from a single base element, since
					 // the shape functions of that base
					 // element are independent of that of
					 // the other one. this latter case
					 // shows up when running
					 // hp/hp_constraints_q_system_06
	std::vector<unsigned int> master_dof_list;
	unsigned int index = 0;
	for (int v=0;
	     v<static_cast<signed int>(GeometryInfo<dim>::vertices_per_face);
	     ++v)
	  {
	    unsigned int dofs_added = 0;
	    unsigned int i          = 0;
	    while (dofs_added < fe2.dofs_per_vertex)
	      {
						 // make sure that we
						 // were able to find
						 // a set of master
						 // dofs and that the
						 // code down below
						 // didn't just reject
						 // all our efforts
		Assert (i < fe1.dofs_per_vertex,
			ExcInternalError());
						 // tentatively push
						 // this vertex dof
		master_dof_list.push_back (index+i);

						 // then see what
						 // happens. if it
						 // succeeds, fine
		if (check_master_dof_list (face_interpolation_matrix,
					   master_dof_list)
		    == true)
		  ++dofs_added;
		else
						   // well, it
						   // didn't. simply
						   // pop that dof
						   // from the list
						   // again and try
						   // with the next
						   // dof
		  master_dof_list.pop_back ();

						 // forward counter by
						 // one
		++i;
	      }
	    index += fe1.dofs_per_vertex;
	  }

	for (int l=0;
	     l<static_cast<signed int>(GeometryInfo<dim>::lines_per_face);
	     ++l)
	  {
					     // same algorithm as above
	    unsigned int dofs_added = 0;
	    unsigned int i          = 0;
	    while (dofs_added < fe2.dofs_per_line)
	      {
		Assert (i < fe1.dofs_per_line,
			ExcInternalError());

		master_dof_list.push_back (index+i);
		if (check_master_dof_list (face_interpolation_matrix,
					   master_dof_list)
		    == true)
		  ++dofs_added;
		else
		  master_dof_list.pop_back ();

		++i;
	      }
	    index += fe1.dofs_per_line;
	  }

	for (int q=0;
	     q<static_cast<signed int>(GeometryInfo<dim>::quads_per_face);
	     ++q)
	  {
					     // same algorithm as above
	    unsigned int dofs_added = 0;
	    unsigned int i          = 0;
	    while (dofs_added < fe2.dofs_per_quad)
	      {
		Assert (i < fe1.dofs_per_quad,
			ExcInternalError());

		master_dof_list.push_back (index+i);
		if (check_master_dof_list (face_interpolation_matrix,
					   master_dof_list)
		    == true)
		  ++dofs_added;
		else
		  master_dof_list.pop_back ();

		++i;
	      }
	    index += fe1.dofs_per_quad;
	  }

	AssertDimension (index, fe1.dofs_per_face);
	AssertDimension (master_dof_list.size(), fe2.dofs_per_face);

					 // finally copy the list into the
					 // mask
	std::fill (master_dof_mask.begin(), master_dof_mask.end(), false);
	for (std::vector<unsigned int>::const_iterator i=master_dof_list.begin();
	     i!=master_dof_list.end(); ++i)
	  master_dof_mask[*i] = true;
      }



				       /**
					* Make sure that the mask exists that
					* determines which dofs will be the
					* masters on refined faces where an
					* fe1 and a fe2 meet.
					*/
      template <int dim, int spacedim>
      void
      ensure_existence_of_master_dof_mask (const FiniteElement<dim,spacedim> &fe1,
					   const FiniteElement<dim,spacedim> &fe2,
					   const FullMatrix<double> &face_interpolation_matrix,
					   std_cxx1x::shared_ptr<std::vector<bool> > &master_dof_mask)
      {
	if (master_dof_mask == std_cxx1x::shared_ptr<std::vector<bool> >())
	  {
	    master_dof_mask = std_cxx1x::shared_ptr<std::vector<bool> >
			      (new std::vector<bool> (fe1.dofs_per_face));
	    select_master_dofs_for_face_restriction (fe1,
						     fe2,
						     face_interpolation_matrix,
						     *master_dof_mask);
	  }
      }


                                       /**
                                        * Make sure that the given @p
                                        * face_interpolation_matrix
                                        * pointer points to a valid
                                        * matrix. If the pointer is zero
                                        * beforehand, create an entry
                                        * with the correct data. If it
                                        * is nonzero, don't touch it.
                                        */
      template <int dim, int spacedim>
      void
      ensure_existence_of_face_matrix (const FiniteElement<dim,spacedim> &fe1,
                                       const FiniteElement<dim,spacedim> &fe2,
                                       std_cxx1x::shared_ptr<FullMatrix<double> > &matrix)
      {
        if (matrix == std_cxx1x::shared_ptr<FullMatrix<double> >())
          {
            matrix = std_cxx1x::shared_ptr<FullMatrix<double> >
                     (new FullMatrix<double> (fe2.dofs_per_face,
                                              fe1.dofs_per_face));
            fe1.get_face_interpolation_matrix (fe2,
                                               *matrix);
          }
      }



				       /**
					* Same, but for subface
					* interpolation matrices.
					*/
      template <int dim, int spacedim>
      void
      ensure_existence_of_subface_matrix (const FiniteElement<dim,spacedim> &fe1,
                                          const FiniteElement<dim,spacedim> &fe2,
                                          const unsigned int        subface,
                                          std_cxx1x::shared_ptr<FullMatrix<double> > &matrix)
      {
        if (matrix == std_cxx1x::shared_ptr<FullMatrix<double> >())
          {
            matrix = std_cxx1x::shared_ptr<FullMatrix<double> >
                     (new FullMatrix<double> (fe2.dofs_per_face,
                                              fe1.dofs_per_face));
            fe1.get_subface_interpolation_matrix (fe2,
                                                  subface,
                                                  *matrix);
          }
      }


				       /**
					* Given the face interpolation
					* matrix between two elements,
					* split it into its master and
					* slave parts and invert the
					* master part as explained in
					* the @ref hp_paper "hp paper".
					*/
#ifdef DEAL_II_ANON_NAMESPACE_BUG
      static
#endif
      void
      ensure_existence_of_split_face_matrix (const FullMatrix<double> &face_interpolation_matrix,
					     const std::vector<bool> &master_dof_mask,
					     std_cxx1x::shared_ptr<std::pair<FullMatrix<double>,FullMatrix<double> > > &split_matrix)
      {
	AssertDimension (master_dof_mask.size(), face_interpolation_matrix.m());
	Assert (std::count (master_dof_mask.begin(), master_dof_mask.end(), true) ==
		static_cast<signed int>(face_interpolation_matrix.n()),
		ExcInternalError());

        if (split_matrix ==
	    std_cxx1x::shared_ptr<std::pair<FullMatrix<double>,FullMatrix<double> > >())
          {
            split_matrix
	      = std_cxx1x::shared_ptr<std::pair<FullMatrix<double>,FullMatrix<double> > >
	      (new std::pair<FullMatrix<double>,FullMatrix<double> >());

	    const unsigned int n_master_dofs = face_interpolation_matrix.n();
	    const unsigned int n_dofs        = face_interpolation_matrix.m();

	    Assert (n_master_dofs <= n_dofs, ExcInternalError());

					     // copy and invert the master
					     // component, copy the slave
					     // component
	    split_matrix->first.reinit (n_master_dofs, n_master_dofs);
	    split_matrix->second.reinit (n_dofs-n_master_dofs, n_master_dofs);

	    unsigned int nth_master_dof = 0,
			 nth_slave_dof  = 0;

	    for (unsigned int i=0; i<n_dofs; ++i)
	      if (master_dof_mask[i] == true)
		{
		  for (unsigned int j=0; j<n_master_dofs; ++j)
		    split_matrix->first(nth_master_dof,j)
		      = face_interpolation_matrix(i,j);
		  ++nth_master_dof;
		}
	      else
		{
		  for (unsigned int j=0; j<n_master_dofs; ++j)
		    split_matrix->second(nth_slave_dof,j)
		      = face_interpolation_matrix(i,j);
		  ++nth_slave_dof;
		}

	    AssertDimension (nth_master_dof, n_master_dofs);
	    AssertDimension (nth_slave_dof, n_dofs-n_master_dofs);

//TODO[WB]: We should make sure very small entries are removed after inversion
	    split_matrix->first.gauss_jordan ();
	  }
      }


                                       // a template that can
                                       // determine statically whether
                                       // a given DoFHandler class
                                       // supports different finite
                                       // element elements
      template <typename>
      struct DoFHandlerSupportsDifferentFEs
      {
          static const bool value = true;
      };


      template <int dim, int spacedim>
      struct DoFHandlerSupportsDifferentFEs< dealii::DoFHandler<dim,spacedim> >
      {
          static const bool value = false;
      };


      template <int dim, int spacedim>
      struct DoFHandlerSupportsDifferentFEs< dealii::MGDoFHandler<dim,spacedim> >
      {
          static const bool value = false;
      };


                                       /**
                                        * A function that returns how
                                        * many different finite
                                        * elements a dof handler
                                        * uses. This is one for non-hp
                                        * DoFHandlers and
                                        * dof_handler.get_fe().size()
                                        * for the hp-versions.
                                        */
      template <int dim, int spacedim>
      unsigned int
      n_finite_elements (const dealii::hp::DoFHandler<dim,spacedim> &dof_handler)
      {
        return dof_handler.get_fe().size();
      }


      template <class DH>
      unsigned int
      n_finite_elements (const DH &)
      {
        return 1;
      }


                                       /**
                                        * For a given face belonging
                                        * to an active cell that
                                        * borders to a more refined
                                        * cell, return the fe_index of
                                        * the most dominating finite
                                        * element used on any of the
                                        * face's subfaces.
                                        */
      template <typename face_iterator>
      unsigned int
      get_most_dominating_subface_fe_index (const face_iterator &face)
      {
	const unsigned int dim
	  = face_iterator::AccessorType::dimension;
	const unsigned int spacedim
	  = face_iterator::AccessorType::space_dimension;

        unsigned int dominating_subface_no = 0;
        for (; dominating_subface_no<face->n_children();
             ++dominating_subface_no)
          {
                                             // each of the subfaces
                                             // can have only a single
                                             // fe_index associated
                                             // with them, since there
                                             // is no cell on the
                                             // other side
            Assert (face->child(dominating_subface_no)
                    ->n_active_fe_indices()
                    == 1,
                    ExcInternalError());

            const FiniteElement<dim,spacedim> &
              this_subface_fe = (face->child(dominating_subface_no)
                                 ->get_fe (face->child(dominating_subface_no)
                                           ->nth_active_fe_index(0)));

            FiniteElementDomination::Domination
              domination = FiniteElementDomination::either_element_can_dominate;
            for (unsigned int sf=0; sf<face->n_children(); ++sf)
              if (sf != dominating_subface_no)
                {
                  const FiniteElement<dim,spacedim> &
                    that_subface_fe = (face->child(sf)
                                       ->get_fe (face->child(sf)
                                                 ->nth_active_fe_index(0)));

                  domination = domination &
                               this_subface_fe.compare_for_face_domination(that_subface_fe);
                }

                                             // see if the element
                                             // on this subface is
                                             // able to dominate
                                             // the ones on all
                                             // other subfaces,
                                             // and if so take it
            if ((domination == FiniteElementDomination::this_element_dominates)
                ||
                (domination == FiniteElementDomination::either_element_can_dominate))
              break;
          }

                                         // check that we have
                                         // found one such subface
        Assert (dominating_subface_no < face->n_children(),
                ExcNotImplemented());

                                         // return the finite element
                                         // index used on it. note
                                         // that only a single fe can
                                         // be active on such subfaces
        return face->child (dominating_subface_no)->nth_active_fe_index(0);
      }



                                       /**
					* Copy constraints into a constraint
					* matrix object.
					*
                                        * This function removes zero
                                        * constraints and those, which
                                        * constrain a DoF which was
                                        * already eliminated in one of
                                        * the previous steps of the hp
                                        * hanging node procedure.
					*
					* It also suppresses very small
					* entries in the constraint matrix to
					* avoid making the sparsity pattern
					* fuller than necessary.
                                        */
#ifdef DEAL_II_ANON_NAMESPACE_BUG
      static
#endif
      void
      filter_constraints (const std::vector<unsigned int> &master_dofs,
			  const std::vector<unsigned int> &slave_dofs,
			  const FullMatrix<double> &face_constraints,
			  ConstraintMatrix &constraints)
      {
	Assert (face_constraints.n () == master_dofs.size (),
		ExcDimensionMismatch(master_dofs.size (),
				     face_constraints.n()));
	Assert (face_constraints.m () == slave_dofs.size (),
		ExcDimensionMismatch(slave_dofs.size (),
				     face_constraints.m()));

	const unsigned int n_master_dofs = master_dofs.size ();
	const unsigned int n_slave_dofs = slave_dofs.size ();

					 // check for a couple
					 // conditions that happened
					 // in parallel distributed
					 // mode
	for (unsigned int row=0; row!=n_slave_dofs; ++row)
	  Assert (slave_dofs[row] != numbers::invalid_unsigned_int,
		  ExcInternalError());
	for (unsigned int col=0; col!=n_master_dofs; ++col)
	  Assert (master_dofs[col] != numbers::invalid_unsigned_int,
		  ExcInternalError());


	for (unsigned int row=0; row!=n_slave_dofs; ++row)
          if (constraints.is_constrained (slave_dofs[row]) == false)
            {
              bool constraint_already_satisfied = false;

                                               // Check if we have an identity
                                               // constraint, which is already
                                               // satisfied by unification of
                                               // the corresponding global dof
                                               // indices
              for (unsigned int i=0; i<n_master_dofs; ++i)
                if (face_constraints (row,i) == 1.0)
                  if (master_dofs[i] == slave_dofs[row])
                    {
                      constraint_already_satisfied = true;
                      break;
                    }

              if (constraint_already_satisfied == false)
                {
						   // add up the absolute
						   // values of all
						   // constraints in this line
						   // to get a measure of
						   // their absolute size
		  double abs_sum = 0;
		  for (unsigned int i=0; i<n_master_dofs; ++i)
		    abs_sum += std::abs (face_constraints(row,i));

						   // then enter those
						   // constraints that are
						   // larger than
						   // 1e-14*abs_sum. everything
						   // else probably originated
						   // from inexact inversion
						   // of matrices and similar
						   // effects. having those
						   // constraints in here will
						   // only lead to problems
						   // because it makes
						   // sparsity patterns fuller
						   // than necessary without
						   // producing any
						   // significant effect
                  constraints.add_line (slave_dofs[row]);
                  for (unsigned int i=0; i<n_master_dofs; ++i)
                    if ((face_constraints(row,i) != 0)
			&&
			(std::fabs(face_constraints(row,i)) >= 1e-14*abs_sum))
		      constraints.add_entry (slave_dofs[row],
					     master_dofs[i],
					     face_constraints (row,i));
		  constraints.set_inhomogeneity (slave_dofs[row], 0.);
                }
            }
      }

    }



    static
    void
    make_hp_hanging_node_constraints (const dealii::DoFHandler<1> &,
				      ConstraintMatrix    &)
    {
				       // nothing to do for regular
				       // dof handlers in 1d
    }



    static
    void
    make_oldstyle_hanging_node_constraints (const dealii::DoFHandler<1> &,
					    ConstraintMatrix    &,
					    dealii::internal::int2type<1>)
    {
				       // nothing to do for regular
				       // dof handlers in 1d
    }


    static
    void
    make_hp_hanging_node_constraints (const dealii::MGDoFHandler<1> &,
				      ConstraintMatrix    &)
    {
				       // nothing to do for regular
				       // dof handlers in 1d
    }



    static
    void
    make_oldstyle_hanging_node_constraints (const dealii::MGDoFHandler<1> &,
					    ConstraintMatrix    &,
					    dealii::internal::int2type<1>)
    {
				       // nothing to do for regular
				       // dof handlers in 1d
    }


    static
    void
    make_hp_hanging_node_constraints (const dealii::hp::DoFHandler<1> &/*dof_handler*/,
				      ConstraintMatrix        &/*constraints*/)
    {
				       // we may have to compute
				       // constraints for
				       // vertices. gotta think about
				       // that a bit more
//TODO[WB]: think about what to do here...
    }



    static
    void
    make_oldstyle_hanging_node_constraints (const dealii::hp::DoFHandler<1> &/*dof_handler*/,
					    ConstraintMatrix        &/*constraints*/,
					    dealii::internal::int2type<1>)
    {
				       // we may have to compute
				       // constraints for
				       // vertices. gotta think about
				       // that a bit more
//TODO[WB]: think about what to do here...
    }


    static
    void
    make_hp_hanging_node_constraints (const dealii::DoFHandler<1,2> &,
				      ConstraintMatrix    &)
    {
				       // nothing to do for regular
				       // dof handlers in 1d
    }



    static
    void
    make_oldstyle_hanging_node_constraints (const dealii::DoFHandler<1,2> &,
					    ConstraintMatrix    &,
					    dealii::internal::int2type<1>)
    {
				       // nothing to do for regular
				       // dof handlers in 1d
    }


//   currently not used but may be in the future:

//     static
//     void
//     make_hp_hanging_node_constraints (const dealii::MGDoFHandler<1,2> &,
// 				      ConstraintMatrix    &)
//     {
// 				       // nothing to do for regular
// 				       // dof handlers in 1d
//     }



//     static
//     void
//     make_oldstyle_hanging_node_constraints (const dealii::MGDoFHandler<1,2> &,
// 					    ConstraintMatrix    &,
// 					    dealii::internal::int2type<1>)
//     {
// 				       // nothing to do for regular
// 				       // dof handlers in 1d
//     }


//     static
//     void
//     make_oldstyle_hanging_node_constraints (const dealii::hp::DoFHandler<1,2> &/*dof_handler*/,
// 					    ConstraintMatrix        &/*constraints*/,
// 					    dealii::internal::int2type<1>)
//     {
// 				       // we may have to compute
// 				       // constraints for
// 				       // vertices. gotta think about
// 				       // that a bit more
// //TODO[WB]: think about what to do here...
//     }
//#endif



    template <class DH>
    static
    void
    make_oldstyle_hanging_node_constraints (const DH         &dof_handler,
					    ConstraintMatrix &constraints,
					    dealii::internal::int2type<2>)
    {
      const unsigned int dim = 2;

      const unsigned int spacedim = DH::space_dimension;

      std::vector<unsigned int> dofs_on_mother;
      std::vector<unsigned int> dofs_on_children;

				       // loop over all lines; only on
				       // lines there can be constraints.
				       // We do so by looping over all
				       // active cells and checking
				       // whether any of the faces are
				       // refined which can only be from
				       // the neighboring cell because
				       // this one is active. In that
				       // case, the face is subject to
				       // constraints
				       //
				       // note that even though we may
				       // visit a face twice if the
				       // neighboring cells are equally
				       // refined, we can only visit each
				       // face with hanging nodes once
      typename DH::active_cell_iterator cell = dof_handler.begin_active(),
					endc = dof_handler.end();
      for (; cell!=endc; ++cell)
					 // artificial cells can at
					 // best neighbor ghost cells,
					 // but we're not interested
					 // in these interfaces
	if (!cell->is_artificial ())
	  for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
	    if (cell->face(face)->has_children())
	      {
						 // in any case, faces
						 // can have at most two
						 // active fe indices,
						 // but here the face
						 // can have only one
						 // (namely the same as
						 // that from the cell
						 // we're sitting on),
						 // and each of the
						 // children can have
						 // only one as
						 // well. check this
		Assert (cell->face(face)->n_active_fe_indices() == 1,
			ExcInternalError());
		Assert (cell->face(face)->fe_index_is_active(cell->active_fe_index())
			== true,
			ExcInternalError());
		for (unsigned int c=0; c<cell->face(face)->n_children(); ++c)
		  if (!cell->neighbor_child_on_subface(face,c)->is_artificial())
		    Assert (cell->face(face)->child(c)->n_active_fe_indices() == 1,
			    ExcInternalError());

						 // right now, all that
						 // is implemented is
						 // the case that both
						 // sides use the same
						 // fe
		for (unsigned int c=0; c<cell->face(face)->n_children(); ++c)
		  if (!cell->neighbor_child_on_subface(face,c)->is_artificial())
		    Assert (cell->face(face)->child(c)
			    ->fe_index_is_active(cell->active_fe_index()) == true,
			    ExcNotImplemented());

						 // ok, start up the work
		const FiniteElement<dim,spacedim> &fe       = cell->get_fe();
		const unsigned int        fe_index = cell->active_fe_index();

		const unsigned int
		  n_dofs_on_mother   = 2*fe.dofs_per_vertex + fe.dofs_per_line,
		  n_dofs_on_children = fe.dofs_per_vertex + 2*fe.dofs_per_line;

		dofs_on_mother.resize (n_dofs_on_mother);
		dofs_on_children.resize (n_dofs_on_children);

		Assert(n_dofs_on_mother == fe.constraints().n(),
		       ExcDimensionMismatch(n_dofs_on_mother,
					    fe.constraints().n()));
		Assert(n_dofs_on_children == fe.constraints().m(),
		       ExcDimensionMismatch(n_dofs_on_children,
					    fe.constraints().m()));

		const typename DH::line_iterator this_face = cell->face(face);

						 // fill the dofs indices. Use same
						 // enumeration scheme as in
						 // @p{FiniteElement::constraints()}
		unsigned int next_index = 0;
		for (unsigned int vertex=0; vertex<2; ++vertex)
		  for (unsigned int dof=0; dof!=fe.dofs_per_vertex; ++dof)
		    dofs_on_mother[next_index++] = this_face->vertex_dof_index(vertex,dof,
									       fe_index);
		for (unsigned int dof=0; dof!=fe.dofs_per_line; ++dof)
		  dofs_on_mother[next_index++] = this_face->dof_index(dof, fe_index);
		AssertDimension (next_index, dofs_on_mother.size());

		next_index = 0;
		for (unsigned int dof=0; dof!=fe.dofs_per_vertex; ++dof)
		  dofs_on_children[next_index++]
		    = this_face->child(0)->vertex_dof_index(1,dof,fe_index);
		for (unsigned int child=0; child<2; ++child)
		  for (unsigned int dof=0; dof!=fe.dofs_per_line; ++dof)
		    dofs_on_children[next_index++]
		      = this_face->child(child)->dof_index(dof, fe_index);
		AssertDimension (next_index, dofs_on_children.size());

						 // for each row in the constraint
						 // matrix for this line:
		for (unsigned int row=0; row!=dofs_on_children.size(); ++row)
		  {
		    constraints.add_line (dofs_on_children[row]);
		    for (unsigned int i=0; i!=dofs_on_mother.size(); ++i)
		      constraints.add_entry (dofs_on_children[row],
					     dofs_on_mother[i],
					     fe.constraints()(row,i));

		    constraints.set_inhomogeneity (dofs_on_children[row], 0.);
		  }
	      }
	    else
	      {
						 // this face has no
						 // children, but it
						 // could still be
						 // that it is shared
						 // by two cells that
						 // use a different fe
						 // index. check a
						 // couple of things,
						 // but ignore the
						 // case that the
						 // neighbor is an
						 // artificial cell
		if (!cell->at_boundary(face) &&
		    !cell->neighbor(face)->is_artificial())
		  {
		    Assert (cell->face(face)->n_active_fe_indices() == 1,
			    ExcNotImplemented());
		    Assert (cell->face(face)
			    ->fe_index_is_active(cell->active_fe_index()) == true,
			    ExcInternalError());
		  }
	      }
    }



    template <class DH>
    static
    void
    make_oldstyle_hanging_node_constraints (const DH         &dof_handler,
					    ConstraintMatrix &constraints,
					    dealii::internal::int2type<3>)
    {
      const unsigned int dim = 3;

      std::vector<unsigned int> dofs_on_mother;
      std::vector<unsigned int> dofs_on_children;

				       // loop over all quads; only on
				       // quads there can be constraints.
				       // We do so by looping over all
				       // active cells and checking
				       // whether any of the faces are
				       // refined which can only be from
				       // the neighboring cell because
				       // this one is active. In that
				       // case, the face is subject to
				       // constraints
				       //
				       // note that even though we may
				       // visit a face twice if the
				       // neighboring cells are equally
				       // refined, we can only visit each
				       // face with hanging nodes once
      typename DH::active_cell_iterator cell = dof_handler.begin_active(),
					endc = dof_handler.end();
      for (; cell!=endc; ++cell)
					 // artificial cells can at
					 // best neighbor ghost cells,
					 // but we're not interested
					 // in these interfaces
	if (!cell->is_artificial ())
	  for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
	    if (cell->face(face)->has_children())
	      {
						 // first of all, make sure that
						 // we treat a case which is
						 // possible, i.e. either no dofs
						 // on the face at all or no
						 // anisotropic refinement
		if (cell->get_fe().dofs_per_face == 0)
		  continue;

		Assert(cell->face(face)->refinement_case()==RefinementCase<dim-1>::isotropic_refinement,
		       ExcNotImplemented());

						 // in any case, faces
						 // can have at most two
						 // active fe indices,
						 // but here the face
						 // can have only one
						 // (namely the same as
						 // that from the cell
						 // we're sitting on),
						 // and each of the
						 // children can have
						 // only one as
						 // well. check this
		AssertDimension (cell->face(face)->n_active_fe_indices(), 1);
		Assert (cell->face(face)->fe_index_is_active(cell->active_fe_index())
			== true,
			ExcInternalError());
		for (unsigned int c=0; c<cell->face(face)->n_children(); ++c)
		  AssertDimension (cell->face(face)->child(c)->n_active_fe_indices(), 1);

						 // right now, all that
						 // is implemented is
						 // the case that both
						 // sides use the same
						 // fe, and not only
						 // that but also that
						 // all lines bounding
						 // this face and the
						 // children have the
						 // same fe
		for (unsigned int c=0; c<cell->face(face)->n_children(); ++c)
		  if (!cell->neighbor_child_on_subface(face,c)->is_artificial())
		    {
		      Assert (cell->face(face)->child(c)
			      ->fe_index_is_active(cell->active_fe_index()) == true,
			      ExcNotImplemented());
		      for (unsigned int e=0; e<4; ++e)
			{
			  Assert (cell->face(face)->child(c)->line(e)
				  ->n_active_fe_indices() == 1,
				  ExcNotImplemented());
			  Assert (cell->face(face)->child(c)->line(e)
				  ->fe_index_is_active(cell->active_fe_index()) == true,
				  ExcNotImplemented());
			}
		    }
		for (unsigned int e=0; e<4; ++e)
		  {
		    Assert (cell->face(face)->line(e)
			    ->n_active_fe_indices() == 1,
			    ExcNotImplemented());
		    Assert (cell->face(face)->line(e)
			    ->fe_index_is_active(cell->active_fe_index()) == true,
			    ExcNotImplemented());
		  }

						 // ok, start up the work
		const FiniteElement<dim> &fe       = cell->get_fe();
		const unsigned int        fe_index = cell->active_fe_index();

	      const unsigned int n_dofs_on_mother = fe.dofs_per_face;
	      const unsigned int n_dofs_on_children = (5*fe.dofs_per_vertex+
				      12*fe.dofs_per_line+
				      4*fe.dofs_per_quad);
//TODO[TL]: think about this and the following in case of anisotropic refinement

		dofs_on_mother.resize (n_dofs_on_mother);
		dofs_on_children.resize (n_dofs_on_children);

		Assert(n_dofs_on_mother == fe.constraints().n(),
		       ExcDimensionMismatch(n_dofs_on_mother,
					    fe.constraints().n()));
		Assert(n_dofs_on_children == fe.constraints().m(),
		       ExcDimensionMismatch(n_dofs_on_children,
					    fe.constraints().m()));

		const typename DH::face_iterator this_face = cell->face(face);

						 // fill the dofs indices. Use same
						 // enumeration scheme as in
						 // @p{FiniteElement::constraints()}
		unsigned int next_index = 0;
		for (unsigned int vertex=0; vertex<4; ++vertex)
		  for (unsigned int dof=0; dof!=fe.dofs_per_vertex; ++dof)
		    dofs_on_mother[next_index++] = this_face->vertex_dof_index(vertex,dof,
									       fe_index);
		for (unsigned int line=0; line<4; ++line)
		  for (unsigned int dof=0; dof!=fe.dofs_per_line; ++dof)
		    dofs_on_mother[next_index++]
		      = this_face->line(line)->dof_index(dof, fe_index);
		for (unsigned int dof=0; dof!=fe.dofs_per_quad; ++dof)
		  dofs_on_mother[next_index++] = this_face->dof_index(dof, fe_index);
		AssertDimension (next_index, dofs_on_mother.size());

		next_index = 0;

						 // assert some consistency
						 // assumptions
						 //TODO[TL]: think about this in case of anisotropic refinement
		Assert (dof_handler.get_tria().get_anisotropic_refinement_flag() ||
			((this_face->child(0)->vertex_index(3) ==
			  this_face->child(1)->vertex_index(2)) &&
			 (this_face->child(0)->vertex_index(3) ==
			  this_face->child(2)->vertex_index(1)) &&
			 (this_face->child(0)->vertex_index(3) ==
			  this_face->child(3)->vertex_index(0))),
			ExcInternalError());
		for (unsigned int dof=0; dof!=fe.dofs_per_vertex; ++dof)
		  dofs_on_children[next_index++]
		    = this_face->child(0)->vertex_dof_index(3,dof);

						 // dof numbers on the centers of
						 // the lines bounding this face
		for (unsigned int line=0; line<4; ++line)
		  for (unsigned int dof=0; dof!=fe.dofs_per_vertex; ++dof)
		    dofs_on_children[next_index++]
		      = this_face->line(line)->child(0)->vertex_dof_index(1,dof, fe_index);

						 // next the dofs on the lines interior
						 // to the face; the order of these
						 // lines is laid down in the
						 // FiniteElement class documentation
		for (unsigned int dof=0; dof<fe.dofs_per_line; ++dof)
		  dofs_on_children[next_index++]
		    = this_face->child(0)->line(1)->dof_index(dof, fe_index);
		for (unsigned int dof=0; dof<fe.dofs_per_line; ++dof)
		  dofs_on_children[next_index++]
		    = this_face->child(2)->line(1)->dof_index(dof, fe_index);
		for (unsigned int dof=0; dof<fe.dofs_per_line; ++dof)
		  dofs_on_children[next_index++]
		    = this_face->child(0)->line(3)->dof_index(dof, fe_index);
		for (unsigned int dof=0; dof<fe.dofs_per_line; ++dof)
		  dofs_on_children[next_index++]
		    = this_face->child(1)->line(3)->dof_index(dof, fe_index);

						 // dofs on the bordering lines
		for (unsigned int line=0; line<4; ++line)
		  for (unsigned int child=0; child<2; ++child)
		    for (unsigned int dof=0; dof!=fe.dofs_per_line; ++dof)
		      dofs_on_children[next_index++]
			= this_face->line(line)->child(child)->dof_index(dof, fe_index);

						 // finally, for the dofs interior
						 // to the four child faces
		for (unsigned int child=0; child<4; ++child)
		  for (unsigned int dof=0; dof!=fe.dofs_per_quad; ++dof)
		    dofs_on_children[next_index++]
		      = this_face->child(child)->dof_index(dof, fe_index);
		AssertDimension (next_index, dofs_on_children.size());

						 // for each row in the constraint
						 // matrix for this line:
		for (unsigned int row=0; row!=dofs_on_children.size(); ++row)
		  {
		    constraints.add_line (dofs_on_children[row]);
		    for (unsigned int i=0; i!=dofs_on_mother.size(); ++i)
		      constraints.add_entry (dofs_on_children[row],
					     dofs_on_mother[i],
					     fe.constraints()(row,i));

		    constraints.set_inhomogeneity(dofs_on_children[row], 0.);
		  }
	      }
	    else
	      {
						 // this face has no
						 // children, but it
						 // could still be
						 // that it is shared
						 // by two cells that
						 // use a different fe
						 // index. check a
						 // couple of things,
						 // but ignore the
						 // case that the
						 // neighbor is an
						 // artificial cell
		if (!cell->at_boundary(face) &&
		    !cell->neighbor(face)->is_artificial())
		  {
		    Assert (cell->face(face)->n_active_fe_indices() == 1,
			    ExcNotImplemented());
		    Assert (cell->face(face)
			    ->fe_index_is_active(cell->active_fe_index()) == true,
			    ExcInternalError());
		  }
	      }
    }


    template <class DH>
    static
    void
    make_hp_hanging_node_constraints (const DH         &dof_handler,
				      ConstraintMatrix &constraints)
    {
				       // note: this function is going
				       // to be hard to understand if
				       // you haven't read the hp
				       // paper. however, we try to
				       // follow the notation laid out
				       // there, so go read the paper
				       // before you try to understand
				       // what is going on here
      const unsigned int dim = DH::dimension;

      const unsigned int spacedim = DH::space_dimension;


                                       // a matrix to be used for
                                       // constraints below. declared
                                       // here and simply resized down
                                       // below to avoid permanent
                                       // re-allocation of memory
      FullMatrix<double> constraint_matrix;

				       // similarly have arrays that
				       // will hold master and slave
				       // dof numbers, as well as a
				       // scratch array needed for the
				       // complicated case below
      std::vector<unsigned int> master_dofs;
      std::vector<unsigned int> slave_dofs;
      std::vector<unsigned int> scratch_dofs;

                                       // caches for the face and
                                       // subface interpolation
                                       // matrices between different
                                       // (or the same) finite
                                       // elements. we compute them
                                       // only once, namely the first
                                       // time they are needed, and
                                       // then just reuse them
      Table<2,std_cxx1x::shared_ptr<FullMatrix<double> > >
        face_interpolation_matrices (n_finite_elements (dof_handler),
                                     n_finite_elements (dof_handler));
      Table<3,std_cxx1x::shared_ptr<FullMatrix<double> > >
        subface_interpolation_matrices (n_finite_elements (dof_handler),
                                        n_finite_elements (dof_handler),
                                        GeometryInfo<dim>::max_children_per_face);

				       // similarly have a cache for
				       // the matrices that are split
				       // into their master and slave
				       // parts, and for which the
				       // master part is
				       // inverted. these two matrices
				       // are derived from the face
				       // interpolation matrix as
				       // described in the @ref hp_paper "hp paper"
      Table<2,std_cxx1x::shared_ptr<std::pair<FullMatrix<double>,FullMatrix<double> > > >
        split_face_interpolation_matrices (n_finite_elements (dof_handler),
					   n_finite_elements (dof_handler));

				       // finally, for each pair of finite
				       // elements, have a mask that states
				       // which of the degrees of freedom on
				       // the coarse side of a refined face
				       // will act as master dofs.
      Table<2,std_cxx1x::shared_ptr<std::vector<bool> > >
	master_dof_masks (n_finite_elements (dof_handler),
			  n_finite_elements (dof_handler));

				       // loop over all faces
				       //
				       // note that even though we may
				       // visit a face twice if the
				       // neighboring cells are equally
				       // refined, we can only visit each
				       // face with hanging nodes once
      typename DH::active_cell_iterator cell = dof_handler.begin_active(),
					endc = dof_handler.end();
      for (; cell!=endc; ++cell)
					 // artificial cells can at
					 // best neighbor ghost cells,
					 // but we're not interested
					 // in these interfaces
	if (!cell->is_artificial ())
	  for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
	    if (cell->face(face)->has_children())
	      {
						 // first of all, make sure that
						 // we treat a case which is
						 // possible, i.e. either no dofs
						 // on the face at all or no
						 // anisotropic refinement
		if (cell->get_fe().dofs_per_face == 0)
		  continue;

		Assert(cell->face(face)->refinement_case()==RefinementCase<dim-1>::isotropic_refinement,
		       ExcNotImplemented());

						 // so now we've found a
						 // face of an active
						 // cell that has
						 // children. that means
						 // that there are
						 // hanging nodes here.

						 // in any case, faces
						 // can have at most two
						 // sets of active fe
						 // indices, but here
						 // the face can have
						 // only one (namely the
						 // same as that from
						 // the cell we're
						 // sitting on), and
						 // each of the children
						 // can have only one as
						 // well. check this
		Assert (cell->face(face)->n_active_fe_indices() == 1,
			ExcInternalError());
		Assert (cell->face(face)->fe_index_is_active(cell->active_fe_index())
			== true,
			ExcInternalError());
		for (unsigned int c=0; c<cell->face(face)->n_children(); ++c)
		  Assert (cell->face(face)->child(c)->n_active_fe_indices() == 1,
			  ExcInternalError());

						 // first find out
						 // whether we can
						 // constrain each of
						 // the subfaces to the
						 // mother face. in the
						 // lingo of the hp
						 // paper, this would be
						 // the simple
						 // case. note that we
						 // can short-circuit
						 // this decision if the
						 // dof_handler doesn't
						 // support hp at all
						 //
						 // ignore all
						 // interfaces with
						 // artificial cells
		FiniteElementDomination::Domination
		  mother_face_dominates = FiniteElementDomination::either_element_can_dominate;

		if (DoFHandlerSupportsDifferentFEs<DH>::value == true)
		  for (unsigned int c=0; c<cell->face(face)->number_of_children(); ++c)
		    if (!cell->neighbor_child_on_subface (face, c)->is_artificial())
		      mother_face_dominates = mother_face_dominates &
					      (cell->get_fe().compare_for_face_domination
					       (cell->neighbor_child_on_subface (face, c)->get_fe()));

		switch (mother_face_dominates)
		  {
		    case FiniteElementDomination::this_element_dominates:
		    case FiniteElementDomination::either_element_can_dominate:
		    {
						       // Case 1 (the
						       // simple case
						       // and the only
						       // case that can
						       // happen for
						       // non-hp
						       // DoFHandlers):
						       // The coarse
						       // element
						       // dominates the
						       // elements on
						       // the subfaces
						       // (or they are
						       // all the same)
						       //
						       // so we are
						       // going to
						       // constrain
						       // the DoFs on
						       // the face
						       // children
						       // against the
						       // DoFs on the
						       // face itself
		      master_dofs.resize (cell->get_fe().dofs_per_face);

		      cell->face(face)->get_dof_indices (master_dofs,
							 cell->active_fe_index ());

						       // Now create
						       // constraint matrix
						       // for the subfaces and
						       // assemble it. ignore
						       // all interfaces with
						       // artificial cells
						       // because we can only
						       // get to such
						       // interfaces if the
						       // current cell is a
						       // ghost cell. also
						       // ignore the interface
						       // if the neighboring
						       // cell is a ghost cell
						       // in 2d: what we would
						       // compute here are the
						       // constraints on the
						       // ghost cell's DoFs,
						       // but we are not
						       // interested in those:
						       // we only want
						       // constraints on
						       // *locally active*
						       // DoFs, not on
						       // *locally relevant*
						       // DoFs. However, in 3d
						       // we must still
						       // compute those
						       // constraints because
						       // it might happen that
						       // a constraint is
						       // related to an edge
						       // where the hanging
						       // node is only
						       // detected if we also
						       // look between ghosts
		      for (unsigned int c=0; c<cell->face(face)->n_children(); ++c)
			{
			  if (cell->neighbor_child_on_subface (face, c)->is_artificial()
			      ||
			      (dim == 2 && cell->neighbor_child_on_subface (face, c)->is_ghost()))
			    continue;

			  const typename DH::active_face_iterator
			    subface = cell->face(face)->child(c);

			  Assert (subface->n_active_fe_indices() == 1,
				  ExcInternalError());

			  const unsigned int
			    subface_fe_index = subface->nth_active_fe_index(0);

							   // we sometime run
							   // into the
							   // situation where
							   // for example on
							   // one big cell we
							   // have a FE_Q(1)
							   // and on the
							   // subfaces we have
							   // a mixture of
							   // FE_Q(1) and
							   // FE_Nothing. In
							   // that case, the
							   // face domination
							   // is
							   // either_element_can_dominate
							   // for the whole
							   // collection of
							   // subfaces, but on
							   // the particular
							   // subface between
							   // FE_Q(1) and
							   // FE_Nothing,
							   // there are no
							   // constraints that
							   // we need to take
							   // care of. in that
							   // case, just
							   // continue
			  if (cell->get_fe().compare_for_face_domination
			      (subface->get_fe(subface_fe_index))
			      ==
			      FiniteElementDomination::no_requirements)
			    continue;
			  
							   // Same procedure as for the
							   // mother cell. Extract the face
							   // DoFs from the cell DoFs.
			  slave_dofs.resize (subface->get_fe(subface_fe_index)
					     .dofs_per_face);
			  subface->get_dof_indices (slave_dofs, subface_fe_index);

			  for (unsigned int i=0; i<slave_dofs.size(); ++i)
			    Assert (slave_dofs[i] != numbers::invalid_unsigned_int,
				    ExcInternalError());

							   // Now create the
							   // element constraint
							   // for this subface.
							   //
							   // As a side remark,
							   // one may wonder the
							   // following:
							   // neighbor_child is
							   // clearly computed
							   // correctly,
							   // i.e. taking into
							   // account
							   // face_orientation
							   // (just look at the
							   // implementation of
							   // that
							   // function). however,
							   // we don't care about
							   // this here, when we
							   // ask for
							   // subface_interpolation
							   // on subface c. the
							   // question rather is:
							   // do we have to
							   // translate 'c' here
							   // as well?
							   //
							   // the answer is in
							   // fact 'no'. if one
							   // does that, results
							   // are wrong:
							   // constraints are
							   // added twice for the
							   // same pair of nodes
							   // but with differing
							   // weights. in
							   // addition, one can
							   // look at the
							   // deal.II/project_*_03
							   // tests that look at
							   // exactly this case:
							   // there, we have a
							   // mesh with at least
							   // one
							   // face_orientation==false
							   // and hanging nodes,
							   // and the results of
							   // those tests show
							   // that the result of
							   // projection verifies
							   // the approximation
							   // properties of a
							   // finite element onto
							   // that mesh
			  ensure_existence_of_subface_matrix
			    (cell->get_fe(),
			     subface->get_fe(subface_fe_index),
			     c,
			     subface_interpolation_matrices
			     [cell->active_fe_index()][subface_fe_index][c]);

							   // Add constraints to global constraint
							   // matrix.
			  filter_constraints (master_dofs,
					      slave_dofs,
					      *(subface_interpolation_matrices
						[cell->active_fe_index()][subface_fe_index][c]),
					      constraints);
			}

		      break;
		    }

		    case FiniteElementDomination::other_element_dominates:
		    case FiniteElementDomination::neither_element_dominates:
		    {
						       // Case 2 (the "complex"
						       // case): at least one
						       // (the neither_... case)
						       // of the finer elements
						       // or all of them (the
						       // other_... case) is
						       // dominating. See the hp
						       // paper for a way how to
						       // deal with this
						       // situation
						       //
						       // since this is
						       // something that
						       // can only
						       // happen for hp
						       // dof handlers,
						       // add a check
						       // here...
		      Assert (DoFHandlerSupportsDifferentFEs<DH>::value == true,
			      ExcInternalError());

						       // we first have
						       // to find the
						       // finite element
						       // that is able
						       // to generate a
						       // space that all
						       // the other ones
						       // can be
						       // constrained to
		      const unsigned int dominating_fe_index
			= get_most_dominating_subface_fe_index (cell->face(face));

		      const FiniteElement<dim,spacedim> &dominating_fe
			= dof_handler.get_fe()[dominating_fe_index];

						       // check also
						       // that it is
						       // able to
						       // constrain the
						       // mother
						       // face. it
						       // should be, or
						       // we wouldn't
						       // have gotten
						       // into the
						       // branch for the
						       // 'complex' case
		      Assert ((dominating_fe.compare_for_face_domination
			       (cell->face(face)->get_fe(cell->face(face)->nth_active_fe_index(0)))
			       == FiniteElementDomination::this_element_dominates)
			      ||
			      (dominating_fe.compare_for_face_domination
			       (cell->face(face)->get_fe(cell->face(face)->nth_active_fe_index(0)))
			       == FiniteElementDomination::either_element_can_dominate),
			      ExcInternalError());


						       // first get the
						       // interpolation matrix
						       // from the mother to the
						       // virtual dofs
		      Assert (dominating_fe.dofs_per_face <=
			      cell->get_fe().dofs_per_face,
			      ExcInternalError());

		      ensure_existence_of_face_matrix
			(dominating_fe,
			 cell->get_fe(),
			 face_interpolation_matrices
			 [dominating_fe_index][cell->active_fe_index()]);

						       // split this matrix into
						       // master and slave
						       // components. invert the
						       // master component
		      ensure_existence_of_master_dof_mask
			(cell->get_fe(),
			 dominating_fe,
			 (*face_interpolation_matrices
			  [dominating_fe_index]
			  [cell->active_fe_index()]),
			 master_dof_masks
			 [dominating_fe_index]
			 [cell->active_fe_index()]);

		      ensure_existence_of_split_face_matrix
			(*face_interpolation_matrices
			 [dominating_fe_index][cell->active_fe_index()],
			 (*master_dof_masks
			  [dominating_fe_index][cell->active_fe_index()]),
			 split_face_interpolation_matrices
			 [dominating_fe_index][cell->active_fe_index()]);

		      const FullMatrix<double> &restrict_mother_to_virtual_master_inv
			= (split_face_interpolation_matrices
			   [dominating_fe_index][cell->active_fe_index()]->first);

		      const FullMatrix<double> &restrict_mother_to_virtual_slave
			= (split_face_interpolation_matrices
			   [dominating_fe_index][cell->active_fe_index()]->second);

						       // now compute
						       // the constraint
						       // matrix as the
						       // product
						       // between the
						       // inverse matrix
						       // and the slave
						       // part
		      constraint_matrix.reinit (cell->get_fe().dofs_per_face -
						dominating_fe.dofs_per_face,
						dominating_fe.dofs_per_face);
		      restrict_mother_to_virtual_slave
			.mmult (constraint_matrix,
				restrict_mother_to_virtual_master_inv);

						       // then figure
						       // out the global
						       // numbers of
						       // master and
						       // slave dofs and
						       // apply
						       // constraints
		      scratch_dofs.resize (cell->get_fe().dofs_per_face);
		      cell->face(face)->get_dof_indices (scratch_dofs,
							 cell->active_fe_index ());

						       // split dofs into master
						       // and slave components
		      master_dofs.clear ();
		      slave_dofs.clear ();
		      for (unsigned int i=0; i<cell->get_fe().dofs_per_face; ++i)
			if ((*master_dof_masks
			     [dominating_fe_index][cell->active_fe_index()])[i] == true)
			  master_dofs.push_back (scratch_dofs[i]);
			else
			  slave_dofs.push_back (scratch_dofs[i]);

		      AssertDimension (master_dofs.size(), dominating_fe.dofs_per_face);
		      AssertDimension (slave_dofs.size(),
				       cell->get_fe().dofs_per_face - dominating_fe.dofs_per_face);

		      filter_constraints (master_dofs,
					  slave_dofs,
					  constraint_matrix,
					  constraints);



						       // next we have
						       // to deal with
						       // the
						       // subfaces. do
						       // as discussed
						       // in the hp
						       // paper
		      for (unsigned int sf=0;
			   sf<cell->face(face)->n_children(); ++sf)
			{
							   // ignore
							   // interfaces
							   // with
							   // artificial
							   // cells as
							   // well as
							   // interfaces
							   // between
							   // ghost
							   // cells in 2d
			  if (cell->neighbor_child_on_subface (face, sf)->is_artificial()
			      ||
			      (dim==2 && cell->is_ghost()
			       &&
			       cell->neighbor_child_on_subface (face, sf)->is_ghost()))
			    continue;

			  Assert (cell->face(face)->child(sf)
				  ->n_active_fe_indices() == 1,
				  ExcInternalError());

			  const unsigned int subface_fe_index
			    = cell->face(face)->child(sf)->nth_active_fe_index(0);
			  const FiniteElement<dim,spacedim> &subface_fe
			    = dof_handler.get_fe()[subface_fe_index];

							   // first get the
							   // interpolation
							   // matrix from the
							   // subface to the
							   // virtual dofs
			  Assert (dominating_fe.dofs_per_face <=
				  subface_fe.dofs_per_face,
				  ExcInternalError());
			  ensure_existence_of_subface_matrix
			    (dominating_fe,
			     subface_fe,
			     sf,
			     subface_interpolation_matrices
			     [dominating_fe_index][subface_fe_index][sf]);

			  const FullMatrix<double> &restrict_subface_to_virtual
			    = *(subface_interpolation_matrices
				[dominating_fe_index][subface_fe_index][sf]);

			  constraint_matrix.reinit (subface_fe.dofs_per_face,
						    dominating_fe.dofs_per_face);

			  restrict_subface_to_virtual
			    .mmult (constraint_matrix,
				    restrict_mother_to_virtual_master_inv);

			  slave_dofs.resize (subface_fe.dofs_per_face);
			  cell->face(face)->child(sf)->get_dof_indices (slave_dofs,
									subface_fe_index);

			  filter_constraints (master_dofs,
					      slave_dofs,
					      constraint_matrix,
					      constraints);
			}

		      break;
		    }

		    case FiniteElementDomination::no_requirements:
							   // there
							   // are no
							   // continuity
							   // requirements
							   // between
							   // the two
							   // elements. record
							   // no
							   // constraints
			  break;

		    default:
							   // we shouldn't get here
			  Assert (false, ExcInternalError());
		  }
	      }
	    else
	      {
						 // this face has no
						 // children, but it
						 // could still be that
						 // it is shared by two
						 // cells that use a
						 // different fe index
		Assert (cell->face(face)
			->fe_index_is_active(cell->active_fe_index()) == true,
			ExcInternalError());

						 // see if there is a
						 // neighbor that is
						 // an artificial
						 // cell. in that
						 // case, we're not
						 // interested in this
						 // interface. we test
						 // this case first
						 // since artificial
						 // cells may not have
						 // an active_fe_index
						 // set, etc
		if (!cell->at_boundary(face)
		    &&
		    cell->neighbor(face)->is_artificial())
		  continue;

						 // Only if there is
						 // a neighbor with
						 // a different
						 // active_fe_index
						 // and the same h-level,
						 // some action has
						 // to be taken.
		if ((DoFHandlerSupportsDifferentFEs<DH>::value == true)
		    &&
		    !cell->face(face)->at_boundary ()
		    &&
		    (cell->neighbor(face)->active_fe_index () !=
		     cell->active_fe_index ())
		    &&
		    (!cell->face(face)->has_children() &&
		     !cell->neighbor_is_coarser(face) ))
		  {
		    const typename DH::cell_iterator neighbor = cell->neighbor (face);

						     // see which side of the
						     // face we have to
						     // constrain
		    switch (cell->get_fe().compare_for_face_domination (neighbor->get_fe ()))
		      {
			case FiniteElementDomination::this_element_dominates:
			{
							   // Get DoFs on
							   // dominating and
							   // dominated side of
							   // the face
			  master_dofs.resize (cell->get_fe().dofs_per_face);
			  cell->face(face)->get_dof_indices (master_dofs,
							     cell->active_fe_index ());

			  slave_dofs.resize (neighbor->get_fe().dofs_per_face);
			  cell->face(face)->get_dof_indices (slave_dofs,
							     neighbor->active_fe_index ());

							   // break if the n_master_dofs == 0,
							   // because we are attempting to
							   // constrain to an element that has
							   // has no face dofs
			  if(master_dofs.size() == 0) break;

							   // make sure
							   // the element
							   // constraints
							   // for this
							   // face are
							   // available
			  ensure_existence_of_face_matrix
			    (cell->get_fe(),
			     neighbor->get_fe(),
			     face_interpolation_matrices
			     [cell->active_fe_index()][neighbor->active_fe_index()]);

							   // Add constraints to global constraint
							   // matrix.
			  filter_constraints (master_dofs,
					      slave_dofs,
					      *(face_interpolation_matrices
						[cell->active_fe_index()]
						[neighbor->active_fe_index()]),
					      constraints);

			  break;
			}

			case FiniteElementDomination::other_element_dominates:
			{
							   // we don't do anything
							   // here since we will
							   // come back to this
							   // face from the other
							   // cell, at which time
							   // we will fall into
							   // the first case
							   // clause above
			  break;
			}

			case FiniteElementDomination::either_element_can_dominate:
			{
							   // it appears as if
							   // neither element has
							   // any constraints on
							   // its neighbor. this
							   // may be because
							   // neither element has
							   // any DoFs on faces at
							   // all. or that the two
							   // elements are
							   // actually the same,
							   // although they happen
							   // to run under
							   // different fe_indices
							   // (this is what
							   // happens in
							   // hp/hp_hanging_nodes_01
							   // for example).
							   //
							   // another possibility
							   // is what happens in
							   // crash_13. there, we
							   // have
							   // FESystem(FE_Q(1),FE_DGQ(0))
							   // vs. FESystem(FE_Q(1),FE_DGQ(1)).
							   // neither of them
							   // dominates the
							   // other. the point is
							   // that it doesn't
							   // matter, since
							   // hopefully for this
							   // case, both sides'
							   // dofs should have
							   // been unified.
							   //
							   // make sure this is
							   // actually true. this
							   // actually only
							   // matters, of course,
							   // if either of the two
							   // finite elements
							   // actually do have
							   // dofs on the face
			  if ((cell->get_fe().dofs_per_face != 0)
			      ||
			      (cell->neighbor(face)->get_fe().dofs_per_face != 0))
			    {
			      Assert (cell->get_fe().dofs_per_face
				      ==
				      cell->neighbor(face)->get_fe().dofs_per_face,
				      ExcNotImplemented());

							       // (ab)use the master
							       // and slave dofs
							       // arrays for a
							       // moment here
			      master_dofs.resize (cell->get_fe().dofs_per_face);
			      cell->face(face)
				->get_dof_indices (master_dofs,
						   cell->active_fe_index ());

			      slave_dofs.resize (cell->neighbor(face)->get_fe().dofs_per_face);
			      cell->face(face)
				->get_dof_indices (slave_dofs,
						   cell->neighbor(face)->active_fe_index ());

			      for (unsigned int i=0; i<cell->get_fe().dofs_per_face; ++i)
				AssertDimension (master_dofs[i], slave_dofs[i]);
			    }

			  break;
			}

			case FiniteElementDomination::neither_element_dominates:
			{
							   // we don't presently
							   // know what exactly to
							   // do here. it isn't quite
							   // clear what exactly we
							   // would have to do
							   // here. sit tight until
							   // someone trips over the
							   // following statement
							   // and see what exactly
							   // is going on
			  Assert (false, ExcNotImplemented());
			  break;
			}

			case FiniteElementDomination::no_requirements:
			{
							   // nothing to do here
			  break;
			}
			
			default:
							       // we shouldn't get
							       // here
			      Assert (false, ExcInternalError());
		      }
		  }
	      }
    }
  }
}




template <class DH>
void
DoFTools::make_hanging_node_constraints (const DH &dof_handler,
					 ConstraintMatrix &constraints)
{
				   // Decide whether to use the
				   // new or old make_hanging_node_constraints
				   // function. If all the FiniteElement
				   // or all elements in a FECollection support
				   // the new face constraint matrix, the
				   // new code will be used.
				   // Otherwise, the old implementation is used
				   // for the moment.
  if (dof_handler.get_fe().hp_constraints_are_implemented ())
    internal::DoFTools::
      make_hp_hanging_node_constraints (dof_handler,
					constraints);
  else
    internal::DoFTools::
      make_oldstyle_hanging_node_constraints (dof_handler,
					      constraints,
					      dealii::internal::int2type<DH::dimension>());
}



namespace internal
{
				// this internal function assigns to each dof
				// the respective component of the vector
				// system. if use_blocks is set, then the
				// assignment is done by blocks, not by
				// components, as specified by
				// component_select. The additional argument
				// component_select only is used for
				// non-primitive FEs, where we need it since
				// more components couple, and no unique
				// component can be assigned. Then, we sort
				// them to the first selected component of the
				// vector system.
  template <class DH>
  inline
  void
  extract_dofs_by_component (const DH                   &dof,
			     const std::vector<bool>    &component_select,
			     const bool                  sort_by_blocks,
			     std::vector<unsigned char> &dofs_by_component)
  {
    const dealii::hp::FECollection<DH::dimension,DH::space_dimension>
      fe_collection (dof.get_fe());
    Assert (fe_collection.n_components() < 256, ExcNotImplemented());
    Assert (dofs_by_component.size() == dof.n_locally_owned_dofs(),
	    ExcDimensionMismatch(dofs_by_component.size(),
				 dof.n_locally_owned_dofs()));

                                   // next set up a table for the
                                   // degrees of freedom on each of
                                   // the cells whether it is
                                   // something interesting or not
    std::vector<std::vector<unsigned char> > local_component_association
      (fe_collection.size());
    for (unsigned int f=0; f<fe_collection.size(); ++f)
      {
	const FiniteElement<DH::dimension,DH::space_dimension> &fe =
	  fe_collection[f];
	local_component_association[f].resize(fe.dofs_per_cell);
	if (sort_by_blocks == true)
	  {
	    for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
	      local_component_association[f][i]
		= fe.system_to_block_index(i).first;
	  }
	else
	  for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
	    if (fe.is_primitive(i))
	      local_component_association[f][i] =
		fe.system_to_component_index(i).first;
	    else
					 // if this shape function is
					 // not primitive, then we have
					 // to work harder. we have to
					 // find out whether _any_ of
					 // the vector components of
					 // this element is selected or
					 // not
					 //
					 // to do so, get the a list of
					 // nonzero elements and see which are
					 // actually active
	      {
		const unsigned int first_comp =
		  (std::find(fe.get_nonzero_components(i).begin(),
			     fe.get_nonzero_components(i).end(),
			     true) -
		   fe.get_nonzero_components(i).begin());
		const unsigned int end_comp =
		  (std::find(fe.get_nonzero_components(i).begin()+first_comp,
			     fe.get_nonzero_components(i).end(),
			     false)-
		   fe.get_nonzero_components(i).begin());

					   // now check whether any of
					   // the components in between
					   // is set
		if (component_select.size() == 0 ||
		    (component_select[first_comp] == true ||
		     std::count(component_select.begin()+first_comp,
				component_select.begin()+end_comp, true) == 0))
		  local_component_association[f][i] = first_comp;
		else
		  for (unsigned int c=first_comp; c<end_comp; ++c)
		    if (component_select[c] == true)
		      {
			local_component_association[f][i] = c;
			break;
		      }
	      }
      }

                                   // then loop over all cells and do
                                   // the work
    std::vector<unsigned int> indices;
    for (typename DH::active_cell_iterator c=dof.begin_active();
	 c!=dof.end(); ++ c)
      if (!c->is_artificial() && !c->is_ghost())
	{
	  const unsigned int fe_index = c->active_fe_index();
	  const unsigned int dofs_per_cell = c->get_fe().dofs_per_cell;
	  indices.resize(dofs_per_cell);
	  c->get_dof_indices(indices);
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    if (dof.locally_owned_dofs().is_element(indices[i]))
	      dofs_by_component[dof.locally_owned_dofs().index_within_set(indices[i])]
		= local_component_association[fe_index][i];
	}
  }
}



template <class DH, typename Number>
void DoFTools::distribute_cell_to_dof_vector (
  const DH             &dof_handler,
  const Vector<Number> &cell_data,
  Vector<double>       &dof_data,
  const unsigned int    component)
{
  const Triangulation<DH::dimension> &tria = dof_handler.get_tria();

  Assert (cell_data.size()==tria.n_active_cells(),
	  ExcWrongSize (cell_data.size(), tria.n_active_cells()));
  Assert (dof_data.size()==dof_handler.n_dofs(),
	  ExcWrongSize (dof_data.size(), dof_handler.n_dofs()));
  Assert (component < n_components(dof_handler),
	  ExcInvalidComponent(component, n_components(dof_handler)));
  Assert (fe_is_primitive(dof_handler) == true,
          ExcFENotPrimitive());

				   // store a flag whether we should care
				   // about different components. this is
				   // just a simplification, we could ask
				   // for this at every single place
				   // equally well
  const bool consider_components = (n_components(dof_handler) != 1);

				   // zero out the components that we
				   // will touch
  if (consider_components == false)
    dof_data = 0;
  else
    {
      std::vector<unsigned char> component_dofs (dof_handler.n_locally_owned_dofs());
      std::vector<bool> component_mask (dof_handler.get_fe().n_components(),
					false);
      component_mask[component] = true;
      internal::extract_dofs_by_component (dof_handler, component_mask,
					   false, component_dofs);

      for (unsigned int i=0; i<dof_data.size(); ++i)
	if (component_dofs[i] == static_cast<unsigned char>(component))
	  dof_data(i) = 0;
    }

				   // count how often we have added a value
				   // in the sum for each dof
  std::vector<unsigned char> touch_count (dof_handler.n_dofs(), 0);

  typename DH::active_cell_iterator cell = dof_handler.begin_active(),
				    endc = dof_handler.end();
  std::vector<unsigned int> dof_indices;
  dof_indices.reserve (max_dofs_per_cell(dof_handler));

  for (unsigned int present_cell = 0; cell!=endc; ++cell, ++present_cell)
    {
      const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
      dof_indices.resize (dofs_per_cell);
      cell->get_dof_indices (dof_indices);

      for (unsigned int i=0; i<dofs_per_cell; ++i)
					 // consider this dof only if it
					 // is the right component. if there
					 // is only one component, short cut
					 // the test
	if (!consider_components ||
	    (cell->get_fe().system_to_component_index(i).first == component))
	  {
					     // sum up contribution of the
					     // present_cell to this dof
	    dof_data(dof_indices[i]) += cell_data(present_cell);
					     // note that we added another
					     // summand
	    ++touch_count[dof_indices[i]];
	  }
    }

				   // compute the mean value on all the
				   // dofs by dividing with the number
				   // of summands.
  for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
    {
				       // assert that each dof was used
				       // at least once. this needs not be
				       // the case if the vector has more than
				       // one component
      Assert (consider_components || (touch_count[i]!=0),
	      ExcInternalError());
      if (touch_count[i] != 0)
	dof_data(i) /=  touch_count[i];
    }
}



template <int dim, int spacedim>
void
DoFTools::extract_dofs (
  const DoFHandler<dim,spacedim>   &dof,
  const std::vector<bool> &component_select,
  std::vector<bool>       &selected_dofs,
  const bool               count_by_blocks)
{
  const FiniteElement<dim,spacedim> &fe = dof.get_fe();

  if (count_by_blocks == true)
    {
      Assert(component_select.size() == fe.n_blocks(),
	     ExcDimensionMismatch(component_select.size(), fe.n_blocks()));
    }
  else
    {
      Assert(component_select.size() == n_components(dof),
	     ExcDimensionMismatch(component_select.size(), n_components(dof)));
    }

  Assert(selected_dofs.size() == dof.n_locally_owned_dofs(),
	 ExcDimensionMismatch(selected_dofs.size(), dof.n_locally_owned_dofs()));

                                   // two special cases: no component
                                   // is selected, and all components
                                   // are selected; both rather
                                   // stupid, but easy to catch
  if (std::count (component_select.begin(), component_select.end(), true)
      == 0)
    {
      std::fill_n (selected_dofs.begin(), dof.n_locally_owned_dofs(), false);
      return;
    }
  else if (std::count (component_select.begin(), component_select.end(), true)
      == static_cast<signed int>(component_select.size()))
    {
      std::fill_n (selected_dofs.begin(), dof.n_locally_owned_dofs(), true);
      return;
    }


				   // preset all values by false
  std::fill_n (selected_dofs.begin(), dof.n_locally_owned_dofs(), false);

				// if we count by blocks, we need to extract
				// the association of blocks with local dofs,
				// and then go through all the cells and set
				// the properties according to this
				// info. Otherwise, we let the function
				// extract_dofs_by_component function do the
				// job.
  std::vector<unsigned char> dofs_by_component (dof.n_locally_owned_dofs());
  internal::extract_dofs_by_component (dof, component_select, count_by_blocks,
				       dofs_by_component);

  for (unsigned int i=0; i<dof.n_locally_owned_dofs(); ++i)
    if (component_select[dofs_by_component[i]] == true)
      selected_dofs[i] = true;
}



template <int dim, int spacedim>
void
DoFTools::extract_dofs (
  const hp::DoFHandler<dim,spacedim>   &dof,
  const std::vector<bool> &component_select,
  std::vector<bool>       &selected_dofs,
  const bool               count_by_blocks)
{
  const FiniteElement<dim,spacedim> &fe = dof.begin_active()->get_fe();

  if (count_by_blocks == true)
    {
      Assert(component_select.size() == fe.n_blocks(),
	     ExcDimensionMismatch(component_select.size(), fe.n_blocks()));
    }
  else
    {
      Assert(component_select.size() == n_components(dof),
	     ExcDimensionMismatch(component_select.size(), n_components(dof)));
    }

  Assert(selected_dofs.size() == dof.n_dofs(),
	 ExcDimensionMismatch(selected_dofs.size(), dof.n_dofs()));

                                   // two special cases: no component
                                   // is selected, and all components
                                   // are selected; both rather
                                   // stupid, but easy to catch
  if (std::count (component_select.begin(), component_select.end(), true)
      == 0)
    {
      std::fill_n (selected_dofs.begin(), dof.n_dofs(), false);
      return;
    };
  if (std::count (component_select.begin(), component_select.end(), true)
      == static_cast<signed int>(component_select.size()))
    {
      std::fill_n (selected_dofs.begin(), dof.n_dofs(), true);
      return;
    };


				   // preset all values by false
  std::fill_n (selected_dofs.begin(), dof.n_dofs(), false);

				// if we count by blocks, we need to extract
				// the association of blocks with local dofs,
				// and then go through all the cells and set
				// the properties according to this
				// info. Otherwise, we let the function
				// extract_dofs_by_component function do the
				// job.
  std::vector<unsigned char> dofs_by_component (dof.n_dofs());
  internal::extract_dofs_by_component (dof, component_select, count_by_blocks,
				       dofs_by_component);

  for (unsigned int i=0; i<dof.n_dofs(); ++i)
    if (component_select[dofs_by_component[i]] == true)
      selected_dofs[i] = true;
}



template<int dim, int spacedim>
void
DoFTools::extract_level_dofs(
  const unsigned int       level,
  const MGDoFHandler<dim,spacedim> &dof,
  const std::vector<bool> &component_select,
  std::vector<bool>       &selected_dofs,
  const bool               count_by_blocks)
{
  const FiniteElement<dim,spacedim>& fe = dof.get_fe();

  if (count_by_blocks == true)
    {
      Assert(component_select.size() == fe.n_blocks(),
	     ExcDimensionMismatch(component_select.size(), fe.n_blocks()));
    }
  else
    {
      Assert(component_select.size() == fe.n_components(),
	     ExcDimensionMismatch(component_select.size(), fe.n_components()));
    }

  Assert(selected_dofs.size() == dof.n_dofs(level),
	 ExcDimensionMismatch(selected_dofs.size(), dof.n_dofs(level)));

                                   // two special cases: no component
                                   // is selected, and all components
                                   // are selected, both rather
                                   // stupid, but easy to catch
  if (std::count (component_select.begin(), component_select.end(), true)
      == 0)
    {
      std::fill_n (selected_dofs.begin(), dof.n_dofs(level), false);
      return;
    };
  if (std::count (component_select.begin(), component_select.end(), true)
      == static_cast<signed int>(component_select.size()))
    {
      std::fill_n (selected_dofs.begin(), dof.n_dofs(level), true);
      return;
    };

    				   // preset all values by false
  std::fill_n (selected_dofs.begin(), dof.n_dofs(level), false);

                                   // next set up a table for the
                                   // degrees of freedom on each of
                                   // the cells whether it is
                                   // something interesting or not
  std::vector<bool> local_selected_dofs (fe.dofs_per_cell, false);
  for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
    if (count_by_blocks == true)
      local_selected_dofs[i]
        = component_select[fe.system_to_block_index(i).first];
    else
      if (fe.is_primitive(i))
	local_selected_dofs[i]
	  = component_select[fe.system_to_component_index(i).first];
      else
					 // if this shape function is
					 // not primitive, then we have
					 // to work harder. we have to
					 // find out whether _any_ of
					 // the vector components of
					 // this element is selected or
					 // not
					 //
					 // to do so, get the first and
					 // last vector components of
					 // the base element to which
					 // the local dof with index i
					 // belongs
	{
	  unsigned int first_comp = 0;
	  const unsigned int this_base = fe.system_to_base_index(i).first.first;
	  const unsigned int this_multiplicity
	    = fe.system_to_base_index(i).first.second;

	  for (unsigned int b=0; b<this_base; ++b)
	    first_comp += fe.base_element(b).n_components() *
			  fe.element_multiplicity(b);
	  for (unsigned int m=0; m<this_multiplicity; ++m)
	    first_comp += fe.base_element(this_base).n_components();
	  const unsigned int end_comp = first_comp +
					fe.base_element(this_base).n_components();

	  Assert (first_comp < fe.n_components(), ExcInternalError());
	  Assert (end_comp <= fe.n_components(),  ExcInternalError());

					   // now check whether any of
					   // the components in between
					   // is set
	  for (unsigned int c=first_comp; c<end_comp; ++c)
	    if (component_select[c] == true)
	      {
		local_selected_dofs[i] = true;
		break;
	      }
	}

                                   // then loop over all cells and do
                                   // work
  std::vector<unsigned int> indices(fe.dofs_per_cell);
  typename MGDoFHandler<dim,spacedim>::cell_iterator c;
  for (c = dof.begin(level) ; c != dof.end(level) ; ++ c)
    {
      c->get_mg_dof_indices(indices);
      for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
        selected_dofs[indices[i]] = local_selected_dofs[i];
    }
}



template <class DH>
void
DoFTools::extract_boundary_dofs (const DH                      &dof_handler,
				 const std::vector<bool>       &component_select,
				 std::vector<bool>             &selected_dofs,
				 const std::set<unsigned char> &boundary_indicators)
{
  Assert (component_select.size() == n_components(dof_handler),
	  ExcWrongSize (component_select.size(),
			n_components(dof_handler)));
  Assert (boundary_indicators.find (255) == boundary_indicators.end(),
	  ExcInvalidBoundaryIndicator());
  const unsigned int dim=DH::dimension;

				   // let's see whether we have to
				   // check for certain boundary
				   // indicators or whether we can
				   // accept all
  const bool check_boundary_indicator = (boundary_indicators.size() != 0);

                                   // also see whether we have to
                                   // check whether a certain vector
                                   // component is selected, or all
  const bool check_vector_component
    = (component_select != std::vector<bool>(component_select.size(),
                                             true));

				   // clear and reset array by default
				   // values
  selected_dofs.clear ();
  selected_dofs.resize (dof_handler.n_dofs(), false);
  std::vector<unsigned int> face_dof_indices;
  face_dof_indices.reserve (max_dofs_per_face(dof_handler));

				   // now loop over all cells and
				   // check whether their faces are at
				   // the boundary. note that we need
				   // not take special care of single
				   // lines being at the boundary
				   // (using
				   // @p{cell->has_boundary_lines}),
				   // since we do not support
				   // boundaries of dimension dim-2,
				   // and so every isolated boundary
				   // line is also part of a boundary
				   // face which we will be visiting
				   // sooner or later
  for (typename DH::active_cell_iterator cell=dof_handler.begin_active();
       cell!=dof_handler.end(); ++cell)
    for (unsigned int face=0;
	 face<GeometryInfo<DH::dimension>::faces_per_cell; ++face)
      if (cell->at_boundary(face))
	if (! check_boundary_indicator ||
	    (boundary_indicators.find (cell->face(face)->boundary_indicator())
	     != boundary_indicators.end()))
	  {
            const FiniteElement<DH::dimension> &fe = cell->get_fe();

            const unsigned int dofs_per_face = fe.dofs_per_face;
            face_dof_indices.resize (dofs_per_face);
	    cell->face(face)->get_dof_indices (face_dof_indices,
					       cell->active_fe_index());

 	    for (unsigned int i=0; i<fe.dofs_per_face; ++i)
 	      if (!check_vector_component)
  		selected_dofs[face_dof_indices[i]] = true;
              else
                                                 // check for
                                                 // component is
                                                 // required. somewhat
                                                 // tricky as usual
                                                 // for the case that
                                                 // the shape function
                                                 // is non-primitive,
                                                 // but use usual
                                                 // convention (see
                                                 // docs)
                {
                                                   // first get at the
                                                   // cell-global
                                                   // number of a face
                                                   // dof, to ask the
                                                   // fe certain
                                                   // questions
                  const unsigned int cell_index
                    = (dim == 1 ?
                       i
                       :
                       (dim == 2 ?
                        (i<2*fe.dofs_per_vertex ? i : i+2*fe.dofs_per_vertex)
                        :
                        (dim == 3 ?
                         (i<4*fe.dofs_per_vertex ?
                          i
                          :
                          (i<4*fe.dofs_per_vertex+4*fe.dofs_per_line ?
                           i+4*fe.dofs_per_vertex
                           :
                           i+4*fe.dofs_per_vertex+8*fe.dofs_per_line))
                         :
                         numbers::invalid_unsigned_int)));
                  if (fe.is_primitive (cell_index))
                    selected_dofs[face_dof_indices[i]]
                      = (component_select[fe.face_system_to_component_index(i).first]
                         == true);
                  else // not primitive
                    {
                      const unsigned int first_nonzero_comp
                        = (std::find (fe.get_nonzero_components(cell_index).begin(),
                                      fe.get_nonzero_components(cell_index).end(),
                                      true)
                           -
                           fe.get_nonzero_components(cell_index).begin());
                      Assert (first_nonzero_comp < fe.n_components(),
                              ExcInternalError());

                      selected_dofs[face_dof_indices[i]]
                        = (component_select[first_nonzero_comp]
                           == true);
                    }
                }
	  }
}



template <class DH>
void
DoFTools::extract_dofs_with_support_on_boundary (const DH                      &dof_handler,
						 const std::vector<bool>       &component_select,
						 std::vector<bool>             &selected_dofs,
						 const std::set<unsigned char> &boundary_indicators)
{
  Assert (component_select.size() == n_components(dof_handler),
	  ExcWrongSize (component_select.size(),
			n_components(dof_handler)));
  Assert (boundary_indicators.find (255) == boundary_indicators.end(),
	  ExcInvalidBoundaryIndicator());

				   // let's see whether we have to
				   // check for certain boundary
				   // indicators or whether we can
				   // accept all
  const bool check_boundary_indicator = (boundary_indicators.size() != 0);

                                   // also see whether we have to
                                   // check whether a certain vector
                                   // component is selected, or all
  const bool check_vector_component
    = (component_select != std::vector<bool>(component_select.size(),
                                             true));

				   // clear and reset array by default
				   // values
  selected_dofs.clear ();
  selected_dofs.resize (dof_handler.n_dofs(), false);
  std::vector<unsigned int> cell_dof_indices;
  cell_dof_indices.reserve (max_dofs_per_cell(dof_handler));

				   // now loop over all cells and
				   // check whether their faces are at
				   // the boundary. note that we need
				   // not take special care of single
				   // lines being at the boundary
				   // (using
				   // @p{cell->has_boundary_lines}),
				   // since we do not support
				   // boundaries of dimension dim-2,
				   // and so every isolated boundary
				   // line is also part of a boundary
				   // face which we will be visiting
				   // sooner or later
  for (typename DH::active_cell_iterator cell=dof_handler.begin_active();
       cell!=dof_handler.end(); ++cell)
    for (unsigned int face=0;
	 face<GeometryInfo<DH::dimension>::faces_per_cell; ++face)
      if (cell->at_boundary(face))
	if (! check_boundary_indicator ||
	    (boundary_indicators.find (cell->face(face)->boundary_indicator())
	     != boundary_indicators.end()))
	  {
            const FiniteElement<DH::dimension> &fe = cell->get_fe();

            const unsigned int dofs_per_cell = fe.dofs_per_cell;
            cell_dof_indices.resize (dofs_per_cell);
	    cell->get_dof_indices (cell_dof_indices);

 	    for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
	      if (fe.has_support_on_face(i,face))
		{
		  if (!check_vector_component)
		    selected_dofs[cell_dof_indices[i]] = true;
		  else
						     // check for
						     // component is
						     // required. somewhat
						     // tricky as usual
						     // for the case that
						     // the shape function
						     // is non-primitive,
						     // but use usual
						     // convention (see
						     // docs)
		    {
		      if (fe.is_primitive (i))
			selected_dofs[cell_dof_indices[i]]
			  = (component_select[fe.system_to_component_index(i).first]
			     == true);
		      else // not primitive
			{
			  const unsigned int first_nonzero_comp
			    = (std::find (fe.get_nonzero_components(i).begin(),
					  fe.get_nonzero_components(i).end(),
					  true)
			       -
			       fe.get_nonzero_components(i).begin());
			  Assert (first_nonzero_comp < fe.n_components(),
				  ExcInternalError());

			  selected_dofs[cell_dof_indices[i]]
			    = (component_select[first_nonzero_comp]
			       == true);
			}
		    }
		}
	  }
}




namespace internal
{
  namespace DoFTools
  {
    template <int spacedim>
    void extract_hanging_node_dofs (const dealii::DoFHandler<1,spacedim> &dof_handler,
				     std::vector<bool>           &selected_dofs)
    {
      Assert(selected_dofs.size() == dof_handler.n_dofs(),
	     ExcDimensionMismatch(selected_dofs.size(), dof_handler.n_dofs()));
				       // preset all values by false
      std::fill_n (selected_dofs.begin(), dof_handler.n_dofs(), false);

				       // there are no hanging nodes in 1d
    }


    template <int spacedim>
    void extract_hanging_node_dofs (const dealii::DoFHandler<2,spacedim> &dof_handler,
				     std::vector<bool>           &selected_dofs)
    {
      const unsigned int dim = 2;

      Assert(selected_dofs.size() == dof_handler.n_dofs(),
	     ExcDimensionMismatch(selected_dofs.size(), dof_handler.n_dofs()));
				       // preset all values by false
      std::fill_n (selected_dofs.begin(), dof_handler.n_dofs(), false);

      const FiniteElement<dim,spacedim> &fe = dof_handler.get_fe();

				       // this function is similar to the
				       // make_sparsity_pattern function,
				       // see there for more information
      typename dealii::DoFHandler<dim,spacedim>::active_cell_iterator
	cell = dof_handler.begin_active(),
	endc = dof_handler.end();
      for (; cell!=endc; ++cell)
	for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
	  if (cell->face(face)->has_children())
	    {
	      const typename dealii::DoFHandler<dim,spacedim>::line_iterator
		line = cell->face(face);

	      for (unsigned int dof=0; dof!=fe.dofs_per_vertex; ++dof)
		selected_dofs[line->child(0)->vertex_dof_index(1,dof)] = true;

	      for (unsigned int child=0; child<2; ++child)
		for (unsigned int dof=0; dof!=fe.dofs_per_line; ++dof)
		  selected_dofs[line->child(child)->dof_index(dof)] = true;
	    }
    }


    template <int spacedim>
    void extract_hanging_node_dofs (const dealii::DoFHandler<3,spacedim> &dof_handler,
				     std::vector<bool>           &selected_dofs)
    {
      const unsigned int dim = 3;

      Assert(selected_dofs.size() == dof_handler.n_dofs(),
	     ExcDimensionMismatch(selected_dofs.size(), dof_handler.n_dofs()));
				       // preset all values by false
      std::fill_n (selected_dofs.begin(), dof_handler.n_dofs(), false);

      const FiniteElement<dim,spacedim> &fe = dof_handler.get_fe();

				       // this function is similar to the
				       // make_sparsity_pattern function,
				       // see there for more information

      typename dealii::DoFHandler<dim,spacedim>::active_cell_iterator
	cell = dof_handler.begin_active(),
	endc = dof_handler.end();
      for (; cell!=endc; ++cell)
	for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
	  if (cell->face(f)->has_children())
	    {
	      const typename dealii::DoFHandler<dim,spacedim>::face_iterator
		face = cell->face(f);

	      for (unsigned int dof=0; dof!=fe.dofs_per_vertex; ++dof)
		selected_dofs[face->child(0)->vertex_dof_index(2,dof)] = true;

					       // dof numbers on the centers of
					       // the lines bounding this face
	      for (unsigned int line=0; line<4; ++line)
		for (unsigned int dof=0; dof!=fe.dofs_per_vertex; ++dof)
		  selected_dofs[face->line(line)->child(0)->vertex_dof_index(1,dof)] = true;

					       // next the dofs on the lines interior
					       // to the face; the order of these
					       // lines is laid down in the
					       // FiniteElement class documentation
	      for (unsigned int dof=0; dof<fe.dofs_per_line; ++dof)
		selected_dofs[face->child(0)->line(1)->dof_index(dof)] = true;
	      for (unsigned int dof=0; dof<fe.dofs_per_line; ++dof)
		selected_dofs[face->child(1)->line(2)->dof_index(dof)] = true;
	      for (unsigned int dof=0; dof<fe.dofs_per_line; ++dof)
		selected_dofs[face->child(2)->line(3)->dof_index(dof)] = true;
	      for (unsigned int dof=0; dof<fe.dofs_per_line; ++dof)
		selected_dofs[face->child(3)->line(0)->dof_index(dof)] = true;

					       // dofs on the bordering lines
	      for (unsigned int line=0; line<4; ++line)
		for (unsigned int child=0; child<2; ++child)
		  for (unsigned int dof=0; dof!=fe.dofs_per_line; ++dof)
		    selected_dofs[face->line(line)->child(child)->dof_index(dof)] = true;

					       // finally, for the dofs interior
					       // to the four child faces
	      for (unsigned int child=0; child<4; ++child)
		for (unsigned int dof=0; dof!=fe.dofs_per_quad; ++dof)
		  selected_dofs[face->child(child)->dof_index(dof)] = true;
	    }
    }
  }
}



template <int dim, int spacedim>
void
DoFTools::
extract_hanging_node_dofs (const DoFHandler<dim,spacedim> &dof_handler,
			   std::vector<bool>              &selected_dofs)
{
  internal::DoFTools::extract_hanging_node_dofs (dof_handler,
						 selected_dofs);
}



template <class DH>
void
DoFTools::extract_subdomain_dofs (const DH                   &dof_handler,
				  const types::subdomain_id_t subdomain_id,
				  std::vector<bool>          &selected_dofs)
{
  Assert(selected_dofs.size() == dof_handler.n_dofs(),
	 ExcDimensionMismatch(selected_dofs.size(), dof_handler.n_dofs()));

                                   // preset all values by false
  std::fill_n (selected_dofs.begin(), dof_handler.n_dofs(), false);

  std::vector<unsigned int> local_dof_indices;
  local_dof_indices.reserve (max_dofs_per_cell(dof_handler));

				   // this function is similar to the
				   // make_sparsity_pattern function,
				   // see there for more information
  typename DH::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    if (cell->subdomain_id() == subdomain_id)
      {
        const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
        local_dof_indices.resize (dofs_per_cell);
	cell->get_dof_indices (local_dof_indices);
	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  selected_dofs[local_dof_indices[i]] = true;
      };
}



template <class DH>
void
DoFTools::extract_locally_owned_dofs (const DH & dof_handler,
				      IndexSet & dof_set)
{
				   // collect all the locally owned dofs
  dof_set = dof_handler.locally_owned_dofs();
  dof_set.compress ();
}



template <class DH>
void
DoFTools::extract_locally_active_dofs (const DH & dof_handler,
				       IndexSet & dof_set)
{
				   // collect all the locally owned dofs
  dof_set = dof_handler.locally_owned_dofs();

				   // add the DoF on the adjacent ghost cells
				   // to the IndexSet, cache them in a
				   // set. need to check each dof manually
				   // because we can't be sure that the dof
				   // range of locally_owned_dofs is really
				   // contiguous.
  std::vector<unsigned int> dof_indices;
  std::set<unsigned int> global_dof_indices;

  typename DH::active_cell_iterator cell = dof_handler.begin_active(),
				    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    if (!cell->is_ghost() && !cell->is_artificial())
      {
	dof_indices.resize(cell->get_fe().dofs_per_cell);
	cell->get_dof_indices(dof_indices);

	for (std::vector<unsigned int>::iterator it=dof_indices.begin();
	     it!=dof_indices.end();
	     ++it)
	  if (!dof_set.is_element(*it))
	    global_dof_indices.insert(*it);
      }

  dof_set.add_indices(global_dof_indices.begin(), global_dof_indices.end());

  dof_set.compress();
}



template <class DH>
void
DoFTools::extract_locally_relevant_dofs (const DH & dof_handler,
					 IndexSet & dof_set)
{
				   // collect all the locally owned dofs
  dof_set = dof_handler.locally_owned_dofs();

				   // add the DoF on the adjacent ghost cells
				   // to the IndexSet, cache them in a
				   // set. need to check each dof manually
				   // because we can't be sure that the dof
				   // range of locally_owned_dofs is really
				   // contiguous.
  std::vector<unsigned int> dof_indices;
  std::set<unsigned int> global_dof_indices;

  typename DH::active_cell_iterator cell = dof_handler.begin_active(),
				    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    if (cell->is_ghost())
      {
	dof_indices.resize(cell->get_fe().dofs_per_cell);
	cell->get_dof_indices(dof_indices);

	for (std::vector<unsigned int>::iterator it=dof_indices.begin();
	     it!=dof_indices.end();
	     ++it)
	  if (!dof_set.is_element(*it))
	    global_dof_indices.insert(*it);
      }

  dof_set.add_indices(global_dof_indices.begin(), global_dof_indices.end());

  dof_set.compress();
}



template <class DH>
void
DoFTools::extract_constant_modes (const DH                        &dof_handler,
				  const std::vector<bool>         &component_select,
				  std::vector<std::vector<bool> > &constant_modes)
{
  const unsigned int n_components = dof_handler.get_fe().n_components();
  Assert (n_components == component_select.size(),
	  ExcDimensionMismatch(n_components,
			       component_select.size()));
  std::vector<unsigned int> localized_component (n_components,
						 numbers::invalid_unsigned_int);
  unsigned int n_components_selected = 0;
  for (unsigned int i=0; i<n_components; ++i)
    if (component_select[i] == true)
      localized_component[i] = n_components_selected++;

  std::vector<unsigned char> dofs_by_component (dof_handler.n_locally_owned_dofs());
  internal::extract_dofs_by_component (dof_handler, component_select, false,
				       dofs_by_component);
  unsigned int n_selected_dofs = 0;
  for (unsigned int i=0; i<n_components; ++i)
    if (component_select[i] == true)
      n_selected_dofs += std::count (dofs_by_component.begin(),
				     dofs_by_component.end(), i);

				 // First count the number of dofs
				 // in the current component.
  constant_modes.resize (n_components_selected, std::vector<bool>(n_selected_dofs,
								  false));
  std::vector<unsigned int> component_list (n_components, 0);
  for (unsigned int d=0; d<n_components; ++d)
    component_list[d] = component_select[d];

  unsigned int counter = 0;
  for (unsigned int i=0; i<dof_handler.n_locally_owned_dofs(); ++i)
    if (component_select[dofs_by_component[i]])
      {
	constant_modes[localized_component[dofs_by_component[i]]][counter] = true;
	++counter;
      }
}



template <class DH>
void
DoFTools::get_active_fe_indices (const DH                  &dof_handler,
				 std::vector<unsigned int> &active_fe_indices)
{
  Assert (active_fe_indices.size() == dof_handler.get_tria().n_active_cells(),
	  ExcWrongSize (active_fe_indices.size(),
			dof_handler.get_tria().n_active_cells()));

  typename DH::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (unsigned int index=0; cell!=endc; ++cell, ++index)
    active_fe_indices[index] = cell->active_fe_index();
}



template <class DH>
void
DoFTools::get_subdomain_association (const DH                  &dof_handler,
				     std::vector<types::subdomain_id_t> &subdomain_association)
{
				   // if the Triangulation is distributed, the
				   // only thing we can usefully ask is for
				   // its locally owned subdomain
  Assert ((dynamic_cast<const parallel::distributed::
	   Triangulation<DH::dimension,DH::space_dimension>*>
	   (&dof_handler.get_tria()) == 0),
	  ExcMessage ("For parallel::distributed::Triangulation objects and "
		      "associated DoF handler objects, asking for any subdomain other "
		      "than the locally owned one does not make sense."));

  Assert(subdomain_association.size() == dof_handler.n_dofs(),
	 ExcDimensionMismatch(subdomain_association.size(),
                              dof_handler.n_dofs()));

                                   // preset all values by an invalid value
  std::fill_n (subdomain_association.begin(), dof_handler.n_dofs(),
               types::invalid_subdomain_id);

  std::vector<unsigned int> local_dof_indices;
  local_dof_indices.reserve (max_dofs_per_cell(dof_handler));

				   // pseudo-randomly assign variables
				   // which lie on the interface
				   // between subdomains to each of
				   // the two or more
  bool coin_flip = true;

				   // loop over all cells and record
				   // which subdomain a DoF belongs
				   // to. toss a coin in case it is on
				   // an interface
  typename DH::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      Assert (cell->is_artificial() == false,
	      ExcMessage ("You can't call this function for meshes that "
			  "have artificial cells."));

      const types::subdomain_id_t subdomain_id = cell->subdomain_id();
      const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
      local_dof_indices.resize (dofs_per_cell);
      cell->get_dof_indices (local_dof_indices);

                                       // set subdomain ids. if dofs
                                       // already have their values
                                       // set then they must be on
                                       // partition interfaces. in
                                       // that case randomly assign
                                       // them to either the previous
                                       // association or the current
                                       // one, where we take "random"
                                       // to be "once this way once
                                       // that way"
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	if (subdomain_association[local_dof_indices[i]] ==
	    numbers::invalid_unsigned_int)
	  subdomain_association[local_dof_indices[i]] = subdomain_id;
	else
	  {
	    if (coin_flip == true)
	      subdomain_association[local_dof_indices[i]] = subdomain_id;
	    coin_flip = !coin_flip;
	  }
    }

  Assert (std::find (subdomain_association.begin(),
                     subdomain_association.end(),
                     types::invalid_subdomain_id)
          == subdomain_association.end(),
          ExcInternalError());
}



template <class DH>
unsigned int
DoFTools::count_dofs_with_subdomain_association (const DH           &dof_handler,
						 const types::subdomain_id_t subdomain)
{
  std::vector<types::subdomain_id_t> subdomain_association (dof_handler.n_dofs());
  get_subdomain_association (dof_handler, subdomain_association);

  return std::count (subdomain_association.begin(),
                     subdomain_association.end(),
                     subdomain);
}



template <class DH>
IndexSet
DoFTools::dof_indices_with_subdomain_association (const DH           &dof_handler,
						  const types::subdomain_id_t subdomain)
{
#ifdef DEAL_II_USE_P4EST
				   // if the DoFHandler is distributed, the
				   // only thing we can usefully ask is for
				   // its locally owned subdomain
  Assert ((dynamic_cast<const parallel::distributed::
	   Triangulation<DH::dimension,DH::space_dimension>*>
	   (&dof_handler.get_tria()) == 0)
	  ||
	  (subdomain ==
	   dynamic_cast<const parallel::distributed::
	   Triangulation<DH::dimension,DH::space_dimension>*>
	   (&dof_handler.get_tria())->locally_owned_subdomain()),
	  ExcMessage ("For parallel::distributed::Triangulation objects and "
		      "associated DoF handler objects, asking for any subdomain other "
		      "than the locally owned one does not make sense."));
#endif

  IndexSet index_set (dof_handler.n_dofs());

  std::vector<unsigned int> local_dof_indices;
  local_dof_indices.reserve (max_dofs_per_cell(dof_handler));

  typename DH::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    if ((cell->is_artificial() == false)
	&&
	(cell->subdomain_id() == subdomain))
      {
	const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
	local_dof_indices.resize (dofs_per_cell);
	cell->get_dof_indices (local_dof_indices);
	index_set.add_indices (local_dof_indices.begin(),
			       local_dof_indices.end());
      }
  index_set.compress ();

  return index_set;
}



template <class DH>
void
DoFTools::count_dofs_with_subdomain_association (const DH           &dof_handler,
						 const types::subdomain_id_t subdomain,
						 std::vector<unsigned int> &n_dofs_on_subdomain)
{
  Assert (n_dofs_on_subdomain.size() == dof_handler.get_fe().n_components(),
	  ExcDimensionMismatch (n_dofs_on_subdomain.size(),
				dof_handler.get_fe().n_components()));
  std::fill (n_dofs_on_subdomain.begin(), n_dofs_on_subdomain.end(), 0);

                                   // in debug mode, make sure that there are
                                   // some cells at least with this subdomain
                                   // id
#ifdef DEBUG
  {
    bool found = false;
    for (typename Triangulation<DH::dimension,DH::space_dimension>::active_cell_iterator
           cell=dof_handler.get_tria().begin_active();
         cell!=dof_handler.get_tria().end(); ++cell)
      if (cell->subdomain_id() == subdomain)
        {
          found = true;
          break;
        }
    Assert (found == true,
            ExcMessage ("There are no cells for the given subdomain!"));
  }
#endif

  std::vector<types::subdomain_id_t> subdomain_association (dof_handler.n_dofs());
  get_subdomain_association (dof_handler, subdomain_association);

  std::vector<unsigned char> component_association (dof_handler.n_dofs());
  internal::extract_dofs_by_component (dof_handler, std::vector<bool>(), false,
				       component_association);

  for (unsigned int c=0; c<dof_handler.get_fe().n_components(); ++c)
    {
      for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
	if ((subdomain_association[i] == subdomain) &&
	    (component_association[i] == static_cast<unsigned char>(c)))
	  ++n_dofs_on_subdomain[c];
    }
}



namespace internal
{
  template <int dim, int spacedim>
  void
  resolve_components (const FiniteElement<dim,spacedim>&fe,
		      const std::vector<unsigned char> &dofs_by_component,
		      const std::vector<unsigned int>  &target_component,
		      const bool                        only_once,
		      std::vector<unsigned int>        &dofs_per_component,
		      unsigned int                     &component)
  {
    for (unsigned int b=0;b<fe.n_base_elements();++b)
      {
	const FiniteElement<dim,spacedim>& base = fe.base_element(b);
				       // Dimension of base element
	unsigned int d = base.n_components();

	for (unsigned int m=0;m<fe.element_multiplicity(b);++m)
	  {
	    if (base.n_base_elements() > 1)
	      resolve_components(base, dofs_by_component, target_component,
				 only_once, dofs_per_component, component);
	    else
	      {
		for (unsigned int dd=0;dd<d;++dd,++component)
		  dofs_per_component[target_component[component]]
		    += std::count(dofs_by_component.begin(),
				  dofs_by_component.end(),
				  component);

				// if we have non-primitive FEs and want all
				// components to show the number of dofs, need
				// to copy the result to those components
		if (!base.is_primitive() && !only_once)
		  for (unsigned int dd=1;dd<d;++dd)
		    dofs_per_component[target_component[component-d+dd]] =
		      dofs_per_component[target_component[component-d]];
	      }
	  }
      }
  }
}


template <int dim, int spacedim>
void
DoFTools::count_dofs_per_component (
  const DoFHandler<dim,spacedim>&     dof_handler,
  std::vector<unsigned int>& dofs_per_component,
  bool only_once,
  std::vector<unsigned int>  target_component)
{
  const FiniteElement<dim,spacedim>& fe = dof_handler.get_fe();

  std::fill (dofs_per_component.begin(), dofs_per_component.end(), 0U);

				   // If the empty vector was given as
				   // default argument, set up this
				   // vector as identity.
  if (target_component.size()==0)
    {
      target_component.resize(fe.n_components());
      for (unsigned int i=0; i<fe.n_components(); ++i)
	target_component[i] = i;
    }
  else
    Assert (target_component.size()==fe.n_components(),
	    ExcDimensionMismatch(target_component.size(),
				 fe.n_components()));


  const unsigned int max_component
    = *std::max_element (target_component.begin(),
			 target_component.end());
  const unsigned int n_target_components = max_component + 1;
  const unsigned int n_components = fe.n_components();

  AssertDimension (dofs_per_component.size(), n_target_components);

				   // special case for only one
				   // component. treat this first
				   // since it does not require any
				   // computations
  if (n_components == 1)
    {
      dofs_per_component[0] = dof_handler.n_locally_owned_dofs();
      return;
    }


				   // otherwise determine the number
				   // of dofs in each component
				   // separately. do so in parallel
  std::vector<unsigned char> dofs_by_component (dof_handler.n_locally_owned_dofs());
  internal::extract_dofs_by_component (dof_handler, std::vector<bool>(), false,
				       dofs_by_component);

				   // next count what we got
  unsigned int component = 0;
  internal::resolve_components(fe, dofs_by_component, target_component,
			       only_once, dofs_per_component, component);
  Assert (n_components == component, ExcInternalError());

				   // finally sanity check. this is
				   // only valid if the finite element
				   // is actually primitive, so
				   // exclude other elements from this
  Assert (!dof_handler.get_fe().is_primitive()
          ||
          (std::accumulate (dofs_per_component.begin(),
                            dofs_per_component.end(), 0U)
           == dof_handler.n_locally_owned_dofs()),
	  ExcInternalError());

				     // reduce information from all CPUs
#ifdef DEAL_II_USE_P4EST
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
  if (const parallel::distributed::Triangulation<dim> * tria
      = (dynamic_cast<const parallel::distributed::Triangulation<dim>*>
	 (&dof_handler.get_tria())))
    {
      std::vector<unsigned int> local_dof_count = dofs_per_component;

      MPI_Allreduce ( &local_dof_count[0], &dofs_per_component[0], n_target_components,
		      MPI_UNSIGNED, MPI_SUM, tria->get_communicator());
    }
#endif
#endif
}



template <int dim, int spacedim>
void
DoFTools::
count_dofs_per_block (const DoFHandler<dim,spacedim>& dof_handler,
		      std::vector<unsigned int> &dofs_per_block,
		      std::vector<unsigned int>  target_block)
{
  const FiniteElement<dim,spacedim>& fe = dof_handler.get_fe();

  std::fill (dofs_per_block.begin(), dofs_per_block.end(), 0U);

				   // If the empty vector was given as
				   // default argument, set up this
				   // vector as identity.
  if (target_block.size()==0)
    {
      target_block.resize(fe.n_blocks());
      for (unsigned int i=0; i<fe.n_blocks(); ++i)
	target_block[i] = i;
    }
  else
    Assert (target_block.size()==fe.n_blocks(),
	    ExcDimensionMismatch(target_block.size(),
				 fe.n_blocks()));



  const unsigned int max_block
    = *std::max_element (target_block.begin(),
			 target_block.end());
  const unsigned int n_target_blocks = max_block + 1;
  const unsigned int n_blocks = fe.n_blocks();

  AssertDimension (dofs_per_block.size(), n_target_blocks);

				   // special case for only one
				   // block. treat this first
				   // since it does not require any
				   // computations
  if (n_blocks == 1)
    {
      dofs_per_block[0] = dof_handler.n_dofs();
      return;
    }
				   // otherwise determine the number
				   // of dofs in each block
				   // separately.
  std::vector<unsigned char> dofs_by_block (dof_handler.n_locally_owned_dofs());
  internal::extract_dofs_by_component (dof_handler, std::vector<bool>(),
				       true, dofs_by_block);

				   // next count what we got
  for (unsigned int block=0; block<fe.n_blocks(); ++block)
    dofs_per_block[target_block[block]]
      += std::count(dofs_by_block.begin(), dofs_by_block.end(),
		    block);

#ifdef DEAL_II_USE_P4EST
#if DEAL_II_COMPILER_SUPPORTS_MPI
				   // if we are working on a parallel
				   // mesh, we now need to collect
				   // this information from all
				   // processors
  if (const parallel::distributed::Triangulation<dim> * tria
      = (dynamic_cast<const parallel::distributed::Triangulation<dim>*>
	 (&dof_handler.get_tria())))
    {
      std::vector<unsigned int> local_dof_count = dofs_per_block;
      MPI_Allreduce ( &local_dof_count[0], &dofs_per_block[0], n_target_blocks,
		      MPI_UNSIGNED, MPI_SUM, tria->get_communicator());
    }
#endif
#endif
}



template <int dim, int spacedim>
void
DoFTools::
count_dofs_per_component (const DoFHandler<dim,spacedim> &dof_handler,
			  std::vector<unsigned int>      &dofs_per_component,
			  std::vector<unsigned int>       target_component)
{
  count_dofs_per_component (dof_handler, dofs_per_component,
			    false, target_component);
}




namespace internal
{
  namespace
  {
				     /**
				      * This is a helper function that
				      * is used in the computation of
				      * integrid constraints. See the
				      * function for a thorough
				      * description of how it works.
				      */
    template <int dim, int spacedim>
    unsigned int
    compute_intergrid_weights_1 (
      const dealii::DoFHandler<dim,spacedim>              &coarse_grid,
      const unsigned int                  coarse_component,
      const dealii::DoFHandler<dim,spacedim>              &fine_grid,
      const unsigned int                  fine_component,
      const InterGridMap<dealii::DoFHandler<dim,spacedim> > &coarse_to_fine_grid_map,
      std::vector<std::map<unsigned int, float> > &weights,
      std::vector<int>                   &weight_mapping)
    {
				       // aliases to the finite elements
				       // used by the dof handlers:
      const FiniteElement<dim,spacedim> &coarse_fe = coarse_grid.get_fe(),
					&fine_fe   = fine_grid.get_fe();

				       // global numbers of dofs
      const unsigned int n_coarse_dofs = coarse_grid.n_dofs(),
			 n_fine_dofs   = fine_grid.n_dofs();

				       // local numbers of dofs
      const unsigned int fine_dofs_per_cell   = fine_fe.dofs_per_cell;

				       // alias the number of dofs per
				       // cell belonging to the
				       // coarse_component which is to be
				       // the restriction of the fine
				       // grid:
      const unsigned int coarse_dofs_per_cell_component
	= coarse_fe.base_element(coarse_fe.component_to_base_index(coarse_component).first).dofs_per_cell;


				       // Try to find out whether the
				       // grids stem from the same coarse
				       // grid. This is a rather crude
				       // test, but better than nothing
      Assert (coarse_grid.get_tria().n_cells(0) == fine_grid.get_tria().n_cells(0),
	      dealii::DoFTools::ExcGridsDontMatch());

				       // check whether the map correlates
				       // the right objects
      Assert (&coarse_to_fine_grid_map.get_source_grid() == &coarse_grid,
	      dealii::DoFTools::ExcGridsDontMatch ());
      Assert (&coarse_to_fine_grid_map.get_destination_grid() == &fine_grid,
	      dealii::DoFTools::ExcGridsDontMatch ());


				       // check whether component numbers
				       // are valid
      Assert (coarse_component < coarse_fe.n_components(),
	      dealii::DoFTools::ExcInvalidComponent (coarse_component, coarse_fe.n_components()));
      Assert (fine_component < fine_fe.n_components(),
	      dealii::DoFTools::ExcInvalidComponent (fine_component, fine_fe.n_components()));
				       // check whether respective finite
				       // elements are equal
      Assert (coarse_fe.base_element (coarse_fe.component_to_base_index(coarse_component).first)
	      ==
	      fine_fe.base_element (fine_fe.component_to_base_index(fine_component).first),
	      dealii::DoFTools::ExcFiniteElementsDontMatch());

#ifdef DEBUG
				       // if in debug mode, check whether
				       // the coarse grid is indeed
				       // coarser everywhere than the fine
				       // grid
      for (typename dealii::DoFHandler<dim,spacedim>::active_cell_iterator
	     cell=coarse_grid.begin_active();
	   cell != coarse_grid.end(); ++cell)
	Assert (cell->level() <= coarse_to_fine_grid_map[cell]->level(),
		dealii::DoFTools::ExcGridNotCoarser());
#endif



/*
 * From here on: the term `parameter' refers to the selected component
 * on the coarse grid and its analogon on the fine grid. The naming of
 * variables containing this term is due to the fact that
 * `selected_component' is longer, but also due to the fact that the
 * code of this function was initially written for a program where the
 * component which we wanted to match between grids was actually the
 * `parameter' variable.
 *
 * Likewise, the terms `parameter grid' and `state grid' refer to the
 * coarse and fine grids, respectively.
 *
 * Changing the names of variables would in principle be a good idea,
 * but would not make things simpler and would be another source of
 * errors. If anyone feels like doing so: patches would be welcome!
 */



				       // set up vectors of cell-local
				       // data; each vector represents one
				       // degree of freedom of the
				       // coarse-grid variable in the
				       // fine-grid element
      std::vector<dealii::Vector<double> >
	parameter_dofs (coarse_dofs_per_cell_component,
			dealii::Vector<double>(fine_dofs_per_cell));
				       // for each coarse dof: find its
				       // position within the fine element
				       // and set this value to one in the
				       // respective vector (all other values
				       // are zero by construction)
      for (unsigned int local_coarse_dof=0;
	   local_coarse_dof<coarse_dofs_per_cell_component;
	   ++local_coarse_dof)
	for (unsigned int fine_dof=0; fine_dof<fine_fe.dofs_per_cell; ++fine_dof)
	  if (fine_fe.system_to_component_index(fine_dof)
	      ==
	      std::make_pair (fine_component, local_coarse_dof))
	    {
	      parameter_dofs[local_coarse_dof](fine_dof) = 1.;
	      break;
	    };


				       // find out how many DoFs there are
				       // on the grids belonging to the
				       // components we want to match
      unsigned int n_parameters_on_fine_grid=0;
      if (true)
	{
					   // have a flag for each dof on
					   // the fine grid and set it
					   // to true if this is an
					   // interesting dof. finally count
					   // how many true's there
	  std::vector<bool> dof_is_interesting (fine_grid.n_dofs(), false);
	  std::vector<unsigned int>  local_dof_indices (fine_fe.dofs_per_cell);

	  for (typename dealii::DoFHandler<dim,spacedim>::active_cell_iterator
		 cell=fine_grid.begin_active();
	       cell!=fine_grid.end(); ++cell)
	    {
	      cell->get_dof_indices (local_dof_indices);
	      for (unsigned int i=0; i<fine_fe.dofs_per_cell; ++i)
		if (fine_fe.system_to_component_index(i).first == fine_component)
		  dof_is_interesting[local_dof_indices[i]] = true;
	    };

	  n_parameters_on_fine_grid = std::count (dof_is_interesting.begin(),
						  dof_is_interesting.end(),
						  true);
	};


				       // set up the weights mapping
      weights.clear ();
      weights.resize (n_coarse_dofs);

      weight_mapping.clear ();
      weight_mapping.resize (n_fine_dofs, -1);

      if (true)
	{
	  std::vector<unsigned int> local_dof_indices(fine_fe.dofs_per_cell);
	  unsigned int next_free_index=0;
	  for (typename dealii::DoFHandler<dim,spacedim>::active_cell_iterator
		 cell=fine_grid.begin_active();
	       cell != fine_grid.end(); ++cell)
	    {
	      cell->get_dof_indices (local_dof_indices);
	      for (unsigned int i=0; i<fine_fe.dofs_per_cell; ++i)
						 // if this DoF is a
						 // parameter dof and has
						 // not yet been numbered,
						 // then do so
		if ((fine_fe.system_to_component_index(i).first == fine_component) &&
		    (weight_mapping[local_dof_indices[i]] == -1))
		  {
		    weight_mapping[local_dof_indices[i]] = next_free_index;
		    ++next_free_index;
		  };
	    };

	  Assert (next_free_index == n_parameters_on_fine_grid,
		  ExcInternalError());
	};


				       // for each cell on the parameter grid:
				       // find out which degrees of freedom on the
				       // fine grid correspond in which way to
				       // the degrees of freedom on the parameter
				       // grid
				       //
				       // do this in a separate function
				       // to allow for multithreading
				       // there. see this function also if
				       // you want to read more
				       // information on the algorithm
				       // used.
      compute_intergrid_weights_2 (coarse_grid, coarse_component,
				   coarse_to_fine_grid_map, parameter_dofs,
				   weight_mapping, weights);


				       // ok, now we have all weights for each
				       // dof on the fine grid. if in debug
				       // mode lets see if everything went smooth,
				       // i.e. each dof has sum of weights one
				       //
				       // in other words this means that
				       // if the sum of all shape
				       // functions on the parameter grid
				       // is one (which is always the
				       // case), then the representation
				       // on the state grid should be as
				       // well (division of unity)
				       //
				       // if the parameter grid has more
				       // than one component, then the
				       // respective dofs of the other
				       // components have sum of weights
				       // zero, of course. we do not
				       // explicitly ask which component
				       // a dof belongs to, but this at
				       // least tests some errors
#ifdef DEBUG
      for (unsigned int col=0; col<n_parameters_on_fine_grid; ++col)
	{
	  double sum=0;
	  for (unsigned int row=0; row<n_coarse_dofs; ++row)
	    if (weights[row].find(col) != weights[row].end())
	      sum += weights[row][col];
	  Assert ((std::fabs(sum-1) < 1.e-12) ||
		  ((coarse_fe.n_components()>1) && (sum==0)), ExcInternalError());
	};
#endif


      return n_parameters_on_fine_grid;
    }


				     /**
				      * This is a function that is
				      * called by the _2 function and
				      * that operates on a range of
				      * cells only. It is used to
				      * split up the whole range of
				      * cells into chunks which are
				      * then worked on in parallel, if
				      * multithreading is available.
				      */
    template <int dim, int spacedim>
    void
    compute_intergrid_weights_3 (
      const dealii::DoFHandler<dim,spacedim>              &coarse_grid,
      const unsigned int                  coarse_component,
      const InterGridMap<dealii::DoFHandler<dim,spacedim> > &coarse_to_fine_grid_map,
      const std::vector<dealii::Vector<double> > &parameter_dofs,
      const std::vector<int>             &weight_mapping,
      std::vector<std::map<unsigned int, float> > &weights,
      const typename dealii::DoFHandler<dim,spacedim>::active_cell_iterator &begin,
      const typename dealii::DoFHandler<dim,spacedim>::active_cell_iterator &end)
    {
				       // aliases to the finite elements
				       // used by the dof handlers:
      const FiniteElement<dim,spacedim> &coarse_fe = coarse_grid.get_fe();

				       // for each cell on the parameter grid:
				       // find out which degrees of freedom on the
				       // fine grid correspond in which way to
				       // the degrees of freedom on the parameter
				       // grid
				       //
				       // since for continuous FEs some
				       // dofs exist on more than one
				       // cell, we have to track which
				       // ones were already visited. the
				       // problem is that if we visit a
				       // dof first on one cell and
				       // compute its weight with respect
				       // to some global dofs to be
				       // non-zero, and later visit the
				       // dof again on another cell and
				       // (since we are on another cell)
				       // recompute the weights with
				       // respect to the same dofs as
				       // above to be zero now, we have to
				       // preserve them. we therefore
				       // overwrite all weights if they
				       // are nonzero and do not enforce
				       // zero weights since that might be
				       // only due to the fact that we are
				       // on another cell.
				       //
				       // example:
				       // coarse grid
				       //  |     |     |
				       //  *-----*-----*
				       //  | cell|cell |
				       //  |  1  |  2  |
				       //  |     |     |
				       //  0-----1-----*
				       //
				       // fine grid
				       //  |  |  |  |  |
				       //  *--*--*--*--*
				       //  |  |  |  |  |
				       //  *--*--*--*--*
				       //  |  |  |  |  |
				       //  *--x--y--*--*
				       //
				       // when on cell 1, we compute the
				       // weights of dof 'x' to be 1/2
				       // from parameter dofs 0 and 1,
				       // respectively. however, when
				       // later we are on cell 2, we again
				       // compute the prolongation of
				       // shape function 1 restricted to
				       // cell 2 to the globla grid and
				       // find that the weight of global
				       // dof 'x' now is zero. however, we
				       // should not overwrite the old
				       // value.
				       //
				       // we therefore always only set
				       // nonzero values. why adding up is
				       // not useful: dof 'y' would get
				       // weight 1 from parameter dof 1 on
				       // both cells 1 and 2, but the
				       // correct weight is nevertheless
				       // only 1.

				       // vector to hold the representation of
				       // a single degree of freedom on the
				       // coarse grid (for the selected fe)
				       // on the fine grid
      const unsigned int n_fine_dofs = weight_mapping.size();
      dealii::Vector<double> global_parameter_representation (n_fine_dofs);

      typename dealii::DoFHandler<dim,spacedim>::active_cell_iterator cell;
      std::vector<unsigned int> parameter_dof_indices (coarse_fe.dofs_per_cell);

      for (cell=begin; cell!=end; ++cell)
	{
					   // get the global indices of the
					   // parameter dofs on this parameter
					   // grid cell
	  cell->get_dof_indices (parameter_dof_indices);

					   // loop over all dofs on this
					   // cell and check whether they
					   // are interesting for us
	  for (unsigned int local_dof=0;
	       local_dof<coarse_fe.dofs_per_cell;
	       ++local_dof)
	    if (coarse_fe.system_to_component_index(local_dof).first
		==
		coarse_component)
	      {
						 // the how-many-th
						 // parameter is this on
						 // this cell?
		const unsigned int local_parameter_dof
		  = coarse_fe.system_to_component_index(local_dof).second;

		global_parameter_representation = 0;

						 // distribute the representation of
						 // @p{local_parameter_dof} on the
						 // parameter grid cell @p{cell} to
						 // the global data space
		coarse_to_fine_grid_map[cell]->
		  set_dof_values_by_interpolation (parameter_dofs[local_parameter_dof],
						   global_parameter_representation);
						 // now that we've got the global
						 // representation of each parameter
						 // dof, we've only got to clobber the
						 // non-zero entries in that vector and
						 // store the result
						 //
						 // what we have learned: if entry @p{i}
						 // of the global vector holds the value
						 // @p{v[i]}, then this is the weight with
						 // which the present dof contributes
						 // to @p{i}. there may be several such
						 // @p{i}s and their weights' sum should
						 // be one. Then, @p{v[i]} should be
						 // equal to @p{\sum_j w_{ij} p[j]} with
						 // @p{p[j]} be the values of the degrees
						 // of freedom on the coarse grid. we
						 // can thus compute constraints which
						 // link the degrees of freedom @p{v[i]}
						 // on the fine grid to those on the
						 // coarse grid, @p{p[j]}. Now to use
						 // these as real constraints, rather
						 // than as additional equations, we
						 // have to identify representants
						 // among the @p{i} for each @p{j}. this will
						 // be done by simply taking the first
						 // @p{i} for which @p{w_{ij}==1}.
						 //
						 // guard modification of
						 // the weights array by a
						 // Mutex. since it should
						 // happen rather rarely
						 // that there are several
						 // threads operating on
						 // different intergrid
						 // weights, have only one
						 // mutex for all of them
		static Threads::ThreadMutex mutex;
		Threads::ThreadMutex::ScopedLock lock (mutex);
		for (unsigned int i=0; i<global_parameter_representation.size(); ++i)
						   // set this weight if it belongs
						   // to a parameter dof.
		  if (weight_mapping[i] != -1)
		    {
						       // only overwrite old
						       // value if not by
						       // zero
		      if (global_parameter_representation(i) != 0)
			{
			  const unsigned int wi = parameter_dof_indices[local_dof],
					     wj = weight_mapping[i];
			  weights[wi][wj] = global_parameter_representation(i);
			};
		    }
		  else
		    Assert (global_parameter_representation(i) == 0,
			    ExcInternalError());
	      }
	}
    }


				     /**
				      * This is a helper function that
				      * is used in the computation of
				      * integrid constraints. See the
				      * function for a thorough
				      * description of how it works.
				      */
    template <int dim, int spacedim>
    void
    compute_intergrid_weights_2 (
      const dealii::DoFHandler<dim,spacedim>              &coarse_grid,
      const unsigned int                  coarse_component,
      const InterGridMap<dealii::DoFHandler<dim,spacedim> > &coarse_to_fine_grid_map,
      const std::vector<dealii::Vector<double> > &parameter_dofs,
      const std::vector<int>             &weight_mapping,
      std::vector<std::map<unsigned int,float> > &weights)
    {
				       // simply distribute the range of
				       // cells to different threads
      typedef typename dealii::DoFHandler<dim,spacedim>::active_cell_iterator active_cell_iterator;
      std::vector<std::pair<active_cell_iterator,active_cell_iterator> >
	cell_intervals = Threads::split_range<active_cell_iterator> (coarse_grid.begin_active(),
								     coarse_grid.end(),
								     multithread_info.n_default_threads);

//TODO: use WorkStream here
      Threads::TaskGroup<> tasks;
      void (*fun_ptr) (const dealii::DoFHandler<dim,spacedim>              &,
		       const unsigned int                  ,
		       const InterGridMap<dealii::DoFHandler<dim,spacedim> > &,
		       const std::vector<dealii::Vector<double> > &,
		       const std::vector<int>             &,
		       std::vector<std::map<unsigned int, float> > &,
		       const typename dealii::DoFHandler<dim,spacedim>::active_cell_iterator &,
		       const typename dealii::DoFHandler<dim,spacedim>::active_cell_iterator &)
	= &compute_intergrid_weights_3<dim>;
      for (unsigned int i=0; i<multithread_info.n_default_threads; ++i)
	tasks += Threads::new_task (fun_ptr,
				    coarse_grid, coarse_component,
				    coarse_to_fine_grid_map, parameter_dofs,
				    weight_mapping, weights,
				    cell_intervals[i].first,
				    cell_intervals[i].second);

				       // wait for the tasks to finish
      tasks.join_all ();
    }
  }
}



template <int dim, int spacedim>
void
DoFTools::compute_intergrid_constraints (
  const DoFHandler<dim,spacedim>              &coarse_grid,
  const unsigned int                  coarse_component,
  const DoFHandler<dim,spacedim>              &fine_grid,
  const unsigned int                  fine_component,
  const InterGridMap<DoFHandler<dim,spacedim> > &coarse_to_fine_grid_map,
  ConstraintMatrix                   &constraints)
{
				   // store the weights with which a dof
				   // on the parameter grid contributes
				   // to a dof on the fine grid. see the
				   // long doc below for more info
				   //
				   // allocate as many rows as there are
				   // parameter dofs on the coarse grid
				   // and as many columns as there are
				   // parameter dofs on the fine grid.
				   //
				   // weight_mapping is used to map the
				   // global (fine grid) parameter dof
				   // indices to the columns
				   //
				   // in the original implementation,
				   // the weights array was actually
				   // of FullMatrix<double> type. this
				   // wasted huge amounts of memory,
				   // but was fast. nonetheless, since
				   // the memory consumption was
				   // quadratic in the number of
				   // degrees of freedom, this was not
				   // very practical, so we now use a
				   // vector of rows of the matrix,
				   // and in each row a vector of
				   // pairs (colnum,value). this seems
				   // like the best tradeoff between
				   // memory and speed, as it is now
				   // linear in memory and still fast
				   // enough.
				   //
				   // to save some memory and since
				   // the weights are usually
				   // (negative) powers of 2, we
				   // choose the value type of the
				   // matrix to be @p{float} rather
				   // than @p{double}.
  std::vector<std::map<unsigned int, float> > weights;

				   // this is this mapping. there is one
				   // entry for each dof on the fine grid;
				   // if it is a parameter dof, then its
				   // value is the column in weights for
				   // that parameter dof, if it is any
				   // other dof, then its value is -1,
				   // indicating an error
  std::vector<int> weight_mapping;

  const unsigned int n_parameters_on_fine_grid
    = internal::compute_intergrid_weights_1 (coarse_grid, coarse_component,
					     fine_grid, fine_component,
					     coarse_to_fine_grid_map,
					     weights, weight_mapping);

				   // global numbers of dofs
  const unsigned int n_coarse_dofs = coarse_grid.n_dofs(),
		     n_fine_dofs   = fine_grid.n_dofs();


				   // get an array in which we store
				   // which dof on the coarse grid is
				   // a parameter and which is not
  std::vector<bool> coarse_dof_is_parameter (coarse_grid.n_dofs());
  if (true)
    {
      std::vector<bool> mask (coarse_grid.get_fe().n_components(),
			      false);
      mask[coarse_component] = true;
      extract_dofs (coarse_grid, mask, coarse_dof_is_parameter);
    };

				   // now we know that the weights in
				   // each row constitute a
				   // constraint. enter this into the
				   // constraints object
				   //
				   // first task: for each parameter
				   // dof on the parameter grid, find
				   // a representant on the fine,
				   // global grid. this is possible
				   // since we use conforming finite
				   // element. we take this
				   // representant to be the first
				   // element in this row with weight
				   // identical to one. the
				   // representant will become an
				   // unconstrained degree of freedom,
				   // while all others will be
				   // constrained to this dof (and
				   // possibly others)
  std::vector<int> representants(n_coarse_dofs, -1);
  for (unsigned int parameter_dof=0; parameter_dof<n_coarse_dofs;
       ++parameter_dof)
    if (coarse_dof_is_parameter[parameter_dof] == true)
      {
					 // if this is the line of a
					 // parameter dof on the
					 // coarse grid, then it
					 // should have at least one
					 // dependent node on the fine
					 // grid
	Assert (weights[parameter_dof].size() > 0, ExcInternalError());

					 // find the column where the
					 // representant is mentioned
	std::map<unsigned int,float>::const_iterator i = weights[parameter_dof].begin();
	for (; i!=weights[parameter_dof].end(); ++i)
	  if (i->second == 1)
	    break;
	Assert (i!=weights[parameter_dof].end(), ExcInternalError());
	const unsigned int column = i->first;

					 // now we know in which column of
					 // weights the representant is, but
					 // we don't know its global index. get
					 // it using the inverse operation of
					 // the weight_mapping
	unsigned int global_dof=0;
	for (; global_dof<weight_mapping.size(); ++global_dof)
	  if (weight_mapping[global_dof] == static_cast<int>(column))
	    break;
	Assert (global_dof < weight_mapping.size(), ExcInternalError());

					 // now enter the representants global
					 // index into our list
	representants[parameter_dof] = global_dof;
      }
    else
      {
					 // consistency check: if this
					 // is no parameter dof on the
					 // coarse grid, then the
					 // respective row must be
					 // empty!
	Assert (weights[parameter_dof].size() == 0, ExcInternalError());
      };



				   // note for people that want to
				   // optimize this function: the
				   // largest part of the computing
				   // time is spent in the following,
				   // rather innocent block of
				   // code. basically, it must be the
				   // ConstraintMatrix::add_entry call
				   // which takes the bulk of the
				   // time, but it is not known to the
				   // author how to make it faster...
  std::vector<std::pair<unsigned int,double> > constraint_line;
  for (unsigned int global_dof=0; global_dof<n_fine_dofs; ++global_dof)
    if (weight_mapping[global_dof] != -1)
				       // this global dof is a parameter
				       // dof, so it may carry a constraint
				       // note that for each global dof,
				       // the sum of weights shall be one,
				       // so we can find out whether this
				       // dof is constrained in the following
				       // way: if the only weight in this row
				       // is a one, and the representant for
				       // the parameter dof of the line in
				       // which this one is is the present
				       // dof, then we consider this dof
				       // to be unconstrained. otherwise,
				       // all other dofs are constrained
      {
	const unsigned int col = weight_mapping[global_dof];
	Assert (col < n_parameters_on_fine_grid, ExcInternalError());

	unsigned int first_used_row=0;

	{
	  Assert (weights.size() > 0, ExcInternalError());
	  std::map<unsigned int,float>::const_iterator
	    col_entry = weights[0].end();
	  for (; first_used_row<n_coarse_dofs; ++first_used_row)
	    {
	      col_entry = weights[first_used_row].find(col);
	      if (col_entry != weights[first_used_row].end())
		break;
	    }

	  Assert (col_entry != weights[first_used_row].end(), ExcInternalError());

	  if ((col_entry->second == 1) &&
	      (representants[first_used_row] == static_cast<int>(global_dof)))
					     // dof unconstrained or
					     // constrained to itself
					     // (in case this cell is
					     // mapped to itself, rather
					     // than to children of
					     // itself)
	    continue;
	}


					 // otherwise enter all constraints
	constraints.add_line (global_dof);

	constraint_line.clear ();
	for (unsigned int row=first_used_row; row<n_coarse_dofs; ++row)
	  {
	    const std::map<unsigned int,float>::const_iterator
	      j = weights[row].find(col);
	    if ((j != weights[row].end()) && (j->second != 0))
	      constraint_line.push_back (std::make_pair(representants[row],
							j->second));
	  };

	constraints.add_entries (global_dof, constraint_line);
      };
}



template <int dim, int spacedim>
void
DoFTools::
compute_intergrid_transfer_representation (
  const DoFHandler<dim,spacedim>              &coarse_grid,
  const unsigned int                  coarse_component,
  const DoFHandler<dim,spacedim>              &fine_grid,
  const unsigned int                  fine_component,
  const InterGridMap<DoFHandler<dim,spacedim> > &coarse_to_fine_grid_map,
  std::vector<std::map<unsigned int, float> > &transfer_representation)
{
				   // store the weights with which a dof
				   // on the parameter grid contributes
				   // to a dof on the fine grid. see the
				   // long doc below for more info
				   //
				   // allocate as many rows as there are
				   // parameter dofs on the coarse grid
				   // and as many columns as there are
				   // parameter dofs on the fine grid.
				   //
				   // weight_mapping is used to map the
				   // global (fine grid) parameter dof
				   // indices to the columns
				   //
				   // in the original implementation,
				   // the weights array was actually
				   // of FullMatrix<double> type. this
				   // wasted huge amounts of memory,
				   // but was fast. nonetheless, since
				   // the memory consumption was
				   // quadratic in the number of
				   // degrees of freedom, this was not
				   // very practical, so we now use a
				   // vector of rows of the matrix,
				   // and in each row a vector of
				   // pairs (colnum,value). this seems
				   // like the best tradeoff between
				   // memory and speed, as it is now
				   // linear in memory and still fast
				   // enough.
				   //
				   // to save some memory and since
				   // the weights are usually
				   // (negative) powers of 2, we
				   // choose the value type of the
				   // matrix to be @p{float} rather
				   // than @p{double}.
  std::vector<std::map<unsigned int, float> > weights;

				   // this is this mapping. there is one
				   // entry for each dof on the fine grid;
				   // if it is a parameter dof, then its
				   // value is the column in weights for
				   // that parameter dof, if it is any
				   // other dof, then its value is -1,
				   // indicating an error
  std::vector<int> weight_mapping;

  internal::compute_intergrid_weights_1 (coarse_grid, coarse_component,
					 fine_grid, fine_component,
					 coarse_to_fine_grid_map,
					 weights, weight_mapping);

				   // now compute the requested
				   // representation
  const unsigned int n_global_parm_dofs
    = std::count_if (weight_mapping.begin(), weight_mapping.end(),
		     std::bind2nd (std::not_equal_to<int> (), -1));

				   // first construct the inverse
				   // mapping of weight_mapping
  std::vector<unsigned int> inverse_weight_mapping (n_global_parm_dofs,
						    DoFHandler<dim,spacedim>::invalid_dof_index);
  for (unsigned int i=0; i<weight_mapping.size(); ++i)
    {
      const unsigned int parameter_dof = weight_mapping[i];
				       // if this global dof is a
				       // parameter
      if (parameter_dof != numbers::invalid_unsigned_int)
	{
	  Assert (parameter_dof < n_global_parm_dofs, ExcInternalError());
	  Assert ((inverse_weight_mapping[parameter_dof] == DoFHandler<dim,spacedim>::invalid_dof_index),
		  ExcInternalError());

	  inverse_weight_mapping[parameter_dof] = i;
	};
    };

				   // next copy over weights array
				   // and replace respective
				   // numbers
  const unsigned int n_rows = weight_mapping.size();

  transfer_representation.clear ();
  transfer_representation.resize (n_rows);

  const unsigned int n_coarse_dofs = coarse_grid.n_dofs();
  for (unsigned int i=0; i<n_coarse_dofs; ++i)
    {
      std::map<unsigned int, float>::const_iterator j = weights[i].begin();
      for (; j!=weights[i].end(); ++j)
	{
	  const unsigned int p = inverse_weight_mapping[j->first];
	  Assert (p<n_rows, ExcInternalError());

	  transfer_representation[p][i] = j->second;
	};
    };
}



template <class DH>
void
DoFTools::map_dof_to_boundary_indices (const DH                  &dof_handler,
                                       std::vector<unsigned int> &mapping)
{
  Assert (&dof_handler.get_fe() != 0, ExcNoFESelected());

  mapping.clear ();
  mapping.insert (mapping.end(), dof_handler.n_dofs(),
		  DH::invalid_dof_index);

  std::vector<unsigned int> dofs_on_face;
  dofs_on_face.reserve (max_dofs_per_face(dof_handler));
  unsigned int next_boundary_index = 0;

				   // now loop over all cells and
				   // check whether their faces are at
				   // the boundary. note that we need
				   // not take special care of single
				   // lines being at the boundary
				   // (using
				   // @p{cell->has_boundary_lines}),
				   // since we do not support
				   // boundaries of dimension dim-2,
				   // and so every isolated boundary
				   // line is also part of a boundary
				   // face which we will be visiting
				   // sooner or later
  typename DH::active_cell_iterator cell = dof_handler.begin_active(),
				    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    for (unsigned int f=0; f<GeometryInfo<DH::dimension>::faces_per_cell; ++f)
      if (cell->at_boundary(f))
        {
          const unsigned int dofs_per_face = cell->get_fe().dofs_per_face;
          dofs_on_face.resize (dofs_per_face);
          cell->face(f)->get_dof_indices (dofs_on_face,
					  cell->active_fe_index());
          for (unsigned int i=0; i<dofs_per_face; ++i)
            if (mapping[dofs_on_face[i]] == DH::invalid_dof_index)
              mapping[dofs_on_face[i]] = next_boundary_index++;
        }

  AssertDimension (next_boundary_index, dof_handler.n_boundary_dofs());
}



template <class DH>
void DoFTools::map_dof_to_boundary_indices (
  const DH                      &dof_handler,
  const std::set<unsigned char> &boundary_indicators,
  std::vector<unsigned int>     &mapping)
{
  Assert (&dof_handler.get_fe() != 0, ExcNoFESelected());
  Assert (boundary_indicators.find (255) == boundary_indicators.end(),
	  ExcInvalidBoundaryIndicator());

  mapping.clear ();
  mapping.insert (mapping.end(), dof_handler.n_dofs(),
		  DH::invalid_dof_index);

				   // return if there is nothing to do
  if (boundary_indicators.size() == 0)
    return;

  std::vector<unsigned int> dofs_on_face;
  dofs_on_face.reserve (max_dofs_per_face(dof_handler));
  unsigned int next_boundary_index = 0;

  typename DH::active_cell_iterator cell = dof_handler.begin_active(),
				    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    for (unsigned int f=0; f<GeometryInfo<DH::dimension>::faces_per_cell; ++f)
      if (boundary_indicators.find (cell->face(f)->boundary_indicator()) !=
          boundary_indicators.end())
	{
	  const unsigned int dofs_per_face = cell->get_fe().dofs_per_face;
	  dofs_on_face.resize (dofs_per_face);
	  cell->face(f)->get_dof_indices (dofs_on_face, cell->active_fe_index());
	  for (unsigned int i=0; i<dofs_per_face; ++i)
	    if (mapping[dofs_on_face[i]] == DH::invalid_dof_index)
	      mapping[dofs_on_face[i]] = next_boundary_index++;
	}

  AssertDimension (next_boundary_index, dof_handler.n_boundary_dofs (boundary_indicators));
}



template <int dim, int spacedim>
void
DoFTools::map_dofs_to_support_points (const Mapping<dim,spacedim>       &mapping,
				      const DoFHandler<dim,spacedim>    &dof_handler,
				      std::vector<Point<spacedim> > &support_points)
{
  const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;

				   // check whether fe has support
				   // points
  Assert (dof_handler.get_fe().has_support_points(),
	  ExcFEHasNoSupportPoints());
  Assert (support_points.size() == dof_handler.n_dofs(),
	  ExcWrongSize (support_points.size(), dof_handler.n_dofs()));

				   // now loop over all cells and
				   // enquire the support points on
				   // each of these. use a dummy
				   // quadrature formula where the
				   // quadrature points are located at
				   // the unit support points to
				   // enquire the location of the
				   // support points in real space
				   //
				   // the weights of the quadrature
				   // rule are set to invalid values
				   // by the used constructor.
  Quadrature<dim> q_dummy(dof_handler.get_fe().get_unit_support_points());
  FEValues<dim,spacedim> fe_values (mapping, dof_handler.get_fe(),
				    q_dummy, update_quadrature_points);
  typename DoFHandler<dim,spacedim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);
  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      cell->get_dof_indices (local_dof_indices);
      const std::vector<Point<spacedim> > & points
  	= fe_values.get_quadrature_points ();
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	support_points[local_dof_indices[i]] = points[i];
    };
}


template<int dim, int spacedim>
void
DoFTools::convert_couplings_to_blocks (
  const DoFHandler<dim,spacedim>& dof_handler,
  const Table<2, Coupling>& table,
  std::vector<Table<2,Coupling> >& tables_by_block)
{
  const FiniteElement<dim,spacedim>& fe = dof_handler.get_fe();
  const unsigned int nb = fe.n_blocks();

  tables_by_block.resize(1);
  tables_by_block[0].reinit(nb, nb);
  tables_by_block[0].fill(none);

  for (unsigned int i=0;i<fe.n_components();++i)
    {
      const unsigned int ib = fe.component_to_block_index(i);
      for (unsigned int j=0;j<fe.n_components();++j)
	{
	  const unsigned int jb = fe.component_to_block_index(j);
	  tables_by_block[0](ib,jb) |= table(i,j);
	}
    }
}


template<int dim, int spacedim>
void
DoFTools::convert_couplings_to_blocks (
  const hp::DoFHandler<dim,spacedim>& dof_handler,
  const Table<2, Coupling>& table,
  std::vector<Table<2,Coupling> >& tables_by_block)
{
  const hp::FECollection<dim>& fe_collection = dof_handler.get_fe();
  tables_by_block.resize(fe_collection.size());

  for (unsigned int f=0;f<fe_collection.size();++f)
    {
      const FiniteElement<dim,spacedim>& fe = fe_collection[f];

      const unsigned int nb = fe.n_blocks();
      tables_by_block[f].reinit(nb, nb);
      tables_by_block[f].fill(none);
      for (unsigned int i=0;i<fe.n_components();++i)
	{
	  const unsigned int ib = fe.component_to_block_index(i);
	  for (unsigned int j=0;j<fe.n_components();++j)
	    {
	      const unsigned int jb = fe.component_to_block_index(j);
	      tables_by_block[f](ib,jb) |= table(i,j);
	    }
	}
    }
}



template <int dim, int spacedim, template <int,int> class DH>
void
DoFTools::make_zero_boundary_constraints (const DH<dim, spacedim> &dof,
 					  ConstraintMatrix        &zero_boundary_constraints,
 					  const std::vector<bool> &component_mask_)
{
  Assert ((component_mask_.size() == 0) ||
 	  (component_mask_.size() == dof.get_fe().n_components()),
 	  ExcMessage ("The number of components in the mask has to be either "
 		      "zero or equal to the number of components in the finite "
 		      "element."));

  const unsigned int        n_components = DoFTools::n_components(dof);

				   // set the component mask to either
				   // the original value or a vector
				   // of trues
  const std::vector<bool> component_mask ((component_mask_.size() == 0) ?
					  std::vector<bool> (n_components, true) :
					  component_mask_);
  Assert (std::count(component_mask.begin(), component_mask.end(), true) > 0,
	  VectorTools::ExcNoComponentSelected());

				   // a field to store the indices
  std::vector<unsigned int> face_dofs;
  face_dofs.reserve (DoFTools::max_dofs_per_face(dof));

  typename DH<dim,spacedim>::active_cell_iterator
    cell = dof.begin_active(),
    endc = dof.end();
  for (; cell!=endc; ++cell)
    for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell;
	 ++face_no)
      {
        const FiniteElement<dim,spacedim> &fe = cell->get_fe();

	typename DH<dim,spacedim>::face_iterator face = cell->face(face_no);
	if (face->boundary_indicator () == 0)
	                           // face is of the right component
	  {
				   // get indices and physical
				   // location on this face
 	    face_dofs.resize (fe.dofs_per_face);
 	    face->get_dof_indices (face_dofs, cell->active_fe_index());

					     // enter those dofs into the list
					     // that match the component
					     // signature.
	    for (unsigned int i=0; i<face_dofs.size(); ++i)
	      {
						 // Find out if a dof
						 // has a contribution
						 // in this component,
						 // and if so, add it
						 // to the list
		const std::vector<bool> &nonzero_component_array
		  = cell->get_fe().get_nonzero_components (i);
		bool nonzero = false;
		for (unsigned int c=0; c<n_components; ++c)
		  if (nonzero_component_array[c] && component_mask[c])
		    {
		      nonzero = true;
		      break;
		    }

		if (nonzero)
		  zero_boundary_constraints.add_line (face_dofs[i]);
	      }
	  }
      }
}


// explicit instantiations

#include "dof_tools.inst"



DEAL_II_NAMESPACE_CLOSE
