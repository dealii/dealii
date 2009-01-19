//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009 by the deal.II authors
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
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <grid/intergrid_map.h>
#include <grid/grid_tools.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_constraints.h>
#include <fe/fe.h>
#include <fe/fe_values.h>
#include <fe/fe_tools.h>
#include <hp/fe_collection.h>
#include <dofs/dof_tools.h>
#include <lac/sparsity_pattern.h>
#include <lac/compressed_sparsity_pattern.h>
#include <lac/compressed_set_sparsity_pattern.h>
#include <lac/compressed_simple_sparsity_pattern.h>
#include <lac/trilinos_sparsity_pattern.h>
#include <lac/block_sparsity_pattern.h>
#include <lac/vector.h>
#include <numerics/vectors.h>

#include <multigrid/mg_dof_handler.h>
#include <multigrid/mg_dof_accessor.h>

#include <algorithm>
#include <numeric>

DEAL_II_NAMESPACE_OPEN

//#define WOLFGANG


template <class DH, class SparsityPattern>
void
DoFTools::make_sparsity_pattern (const DH               &dof,
				 SparsityPattern        &sparsity,
				 const ConstraintMatrix &constraints,
				 const bool              keep_constrained_dofs)
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
#ifdef DEAL_II_USE_TRILINOS
    if ((types_are_equal<SparsityPattern,TrilinosWrappers::SparsityPattern>::value
	 ||
	 types_are_equal<SparsityPattern,TrilinosWrappers::BlockSparsityPattern>::value)
	&&
	cell->subdomain_id() != 
	Utilities::Trilinos::get_this_mpi_process(Utilities::Trilinos::comm_world()))
      continue;
    else
	 
#endif
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
DoFTools::make_sparsity_pattern (
  const DH                &dof,
  const Table<2,Coupling> &couplings,
  SparsityPattern         &sparsity,
  const ConstraintMatrix  &constraints,
  const bool               keep_constrained_dofs)
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
  for (unsigned int f=0; f<fe_collection.size(); ++f)
    {
      const unsigned int dofs_per_cell = fe_collection[f].dofs_per_cell;

      dof_mask[f].reinit (dofs_per_cell, dofs_per_cell);
      
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=0; j<dofs_per_cell; ++j)
	  if (fe_collection[f].is_primitive(i) &&
	      fe_collection[f].is_primitive(j))
	    dof_mask[f][i][j]
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
	      
	      dof_mask[f][i][j]
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
#ifdef DEAL_II_USE_TRILINOS
    if ((types_are_equal<SparsityPattern,TrilinosWrappers::SparsityPattern>::value
	 ||
	 types_are_equal<SparsityPattern,TrilinosWrappers::BlockSparsityPattern>::value)
	&&
	cell->subdomain_id() != 
	Utilities::Trilinos::get_this_mpi_process(Utilities::Trilinos::comm_world()))
      continue;
    else
	 
#endif
    {
      const unsigned int fe_index
	= cell->active_fe_index();

      const unsigned int dofs_per_cell
	= fe_collection[fe_index].dofs_per_cell;

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
            for (unsigned int j=0; j<dofs_per_cell_col; ++j)
              sparsity.add (local_dof_indices_row[i],
                            local_dof_indices_col[j]);
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
                for (unsigned int j=0; j<dofs_per_cell_col; ++j)
                  sparsity.add (local_dof_indices_row[i],
                                local_dof_indices_col[j]);
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
                for (unsigned int j=0; j<dofs_per_cell_col; ++j)
                  sparsity.add (local_dof_indices_row[i],
                                local_dof_indices_col[j]);
            }
        }
    }
}



#if deal_II_dimension == 1

template <class DH, class SparsityPattern>
void
DoFTools::make_boundary_sparsity_pattern (
  const DH                               &dof_handler,
  const typename FunctionMap<DH::space_dimension>::type &function_map,
  const std::vector<unsigned int>        &dof_to_boundary_mapping,
  SparsityPattern                        &sparsity)
{
				   // first check left, the right
				   // boundary point
  for (unsigned int direction=0; direction<2; ++direction)
    {
				       // if this boundary is not
				       // requested, then go on with next one
      if (function_map.find(direction) ==
	  function_map.end())
	continue;

				       // find active cell at that
				       // boundary: first go to
				       // left/right, then to children
      typename DH::cell_iterator cell = dof_handler.begin(0);
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
	for (unsigned int j=0; j<dofs_per_vertex; ++j)
	  sparsity.add (boundary_dof_boundary_indices[i],
			boundary_dof_boundary_indices[j]);
    };
}



template <class DH, class SparsityPattern>
void DoFTools::make_boundary_sparsity_pattern (
  const DH                        &dof_handler,
  const std::vector<unsigned int> &dof_to_boundary_mapping,
  SparsityPattern                 &sparsity)
{
				   // there are only 2 boundary
				   // indicators in 1d, so it is no
				   // performance problem to call the
				   // other function
  typename DH::FunctionMap boundary_indicators;
  boundary_indicators[0] = 0;
  boundary_indicators[1] = 0;
  make_boundary_sparsity_pattern<DH, SparsityPattern> (dof_handler,
						       boundary_indicators,
						       dof_to_boundary_mapping,
						       sparsity);
}


#else


template <class DH, class SparsityPattern>
void
DoFTools::make_boundary_sparsity_pattern (
  const DH                        &dof,
  const std::vector<unsigned int> &dof_to_boundary_mapping,
  SparsityPattern                 &sparsity)
{
  const unsigned int n_dofs = dof.n_dofs();

  Assert (dof_to_boundary_mapping.size() == n_dofs, ExcInternalError());
  Assert (sparsity.n_rows() == dof.n_boundary_dofs(),
	  ExcDimensionMismatch (sparsity.n_rows(), dof.n_boundary_dofs()));
  Assert (sparsity.n_cols() == dof.n_boundary_dofs(),
	  ExcDimensionMismatch (sparsity.n_cols(), dof.n_boundary_dofs()));
#ifdef DEBUG
  if (sparsity.n_rows() != 0)
    {
      unsigned int max_element = 0;
      for (std::vector<unsigned int>::const_iterator i=dof_to_boundary_mapping.begin();
	   i!=dof_to_boundary_mapping.end(); ++i)
	if ((*i != DH::invalid_dof_index) &&
	    (*i > max_element))
	  max_element = *i;
      Assert (max_element  == sparsity.n_rows()-1,
	      ExcInternalError());
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
  const unsigned int n_dofs = dof.n_dofs();

  Assert (dof_to_boundary_mapping.size() == n_dofs, ExcInternalError());
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
      Assert (max_element  == sparsity.n_rows()-1,
	      ExcInternalError());
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

#endif


#if deal_II_dimension != 1

template <class DH, class SparsityPattern>
void
DoFTools::make_flux_sparsity_pattern (
  const DH        &dof,
  SparsityPattern &sparsity)
{
  const unsigned int n_dofs = dof.n_dofs();
  
  Assert (sparsity.n_rows() == n_dofs,
	  ExcDimensionMismatch (sparsity.n_rows(), n_dofs));
  Assert (sparsity.n_cols() == n_dofs,
	  ExcDimensionMismatch (sparsity.n_cols(), n_dofs));

  std::vector<unsigned int> dofs_on_this_cell;
  std::vector<unsigned int> dofs_on_other_cell;
  dofs_on_this_cell.reserve (max_dofs_per_cell(dof));
  dofs_on_other_cell.reserve (max_dofs_per_cell(dof));
  typename DH::active_cell_iterator cell = dof.begin_active(),
				    endc = dof.end();

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
  
  for (; cell!=endc; ++cell)
    {
      const unsigned int n_dofs_on_this_cell = cell->get_fe().dofs_per_cell;
      dofs_on_this_cell.resize (n_dofs_on_this_cell);
      cell->get_dof_indices (dofs_on_this_cell);
				       // make sparsity pattern for this cell
      for (unsigned int i=0; i<n_dofs_on_this_cell; ++i)
	for (unsigned int j=0; j<n_dofs_on_this_cell; ++j)
	  sparsity.add (dofs_on_this_cell[i],
			dofs_on_this_cell[j]);

				       // Loop over all interior neighbors
      for (unsigned int face = 0;
	   face < GeometryInfo<DH::dimension>::faces_per_cell;
	   ++face)
	{
	  typename DH::face_iterator cell_face = cell->face(face);
	  if (cell_face->user_flag_set ())
	    continue;

	  if (! cell_face->at_boundary() )
	    {
	      typename DH::cell_iterator neighbor = cell->neighbor(face);

	      const unsigned int neighbor_face
                = cell->neighbor_face_no(face);

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

                      for (unsigned int i=0; i<n_dofs_on_this_cell; ++i)
                        for (unsigned int j=0; j<n_dofs_on_neighbor; ++j)
                          {
                            sparsity.add (dofs_on_this_cell[i],
                                          dofs_on_other_cell[j]);
                            sparsity.add (dofs_on_other_cell[j],
                                          dofs_on_this_cell[i]);
                          }
		      sub_neighbor->face(neighbor_face)->set_user_flag ();
		    }
		}
              else
                {
						   // Refinement edges are
						   // taken care of by
						   // coarser cells
		  if (cell->neighbor_is_coarser(face))
		    continue;
		  
		  const unsigned int n_dofs_on_neighbor
                    = neighbor->get_fe().dofs_per_cell;
                  dofs_on_other_cell.resize (n_dofs_on_neighbor);

                  neighbor->get_dof_indices (dofs_on_other_cell);
		  for (unsigned int i=0; i<n_dofs_on_this_cell; ++i)
                    for (unsigned int j=0; j<n_dofs_on_neighbor; ++j)
                      {
                        sparsity.add (dofs_on_this_cell[i],
                                      dofs_on_other_cell[j]);
                        sparsity.add (dofs_on_other_cell[j],
                                      dofs_on_this_cell[i]);
                      }
		  neighbor->face(neighbor_face)->set_user_flag (); 
		}
	    } 
	}
    }

				   // finally restore the user flags
  const_cast<Triangulation<DH::dimension> &>(dof.get_tria()).load_user_flags(user_flags);
}

#else // deal_II_dimension == 1


template <class DH, class SparsityPattern>
void
DoFTools::make_flux_sparsity_pattern (
  const DH        &dof,
  SparsityPattern &sparsity)
{
  typedef typename DH::cell_iterator        cell_iterator;
  typedef typename DH::active_cell_iterator active_cell_iterator;

  std::vector<unsigned int> local_dof_indices;
  std::vector<unsigned int> neighbor_dof_indices;
  local_dof_indices.reserve (max_dofs_per_cell(dof));
  neighbor_dof_indices.reserve (max_dofs_per_cell(dof));
  
  active_cell_iterator cell = dof.begin_active(),
		       endc = dof.end();
  for (; cell!=endc; ++cell)
    {
				       // first do couplings of dofs
				       // locally on this cell
      const unsigned int n_dofs_on_this_cell = cell->get_fe().dofs_per_cell;
      local_dof_indices.resize (n_dofs_on_this_cell);
      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<n_dofs_on_this_cell; ++i)
	for (unsigned int j=0; j<n_dofs_on_this_cell; ++j)
	  sparsity.add (local_dof_indices[i], local_dof_indices[j]);

				       // then do the same for the up
				       // to 2 neighbors
      for (unsigned int nb=0; nb<2; ++nb)
	if (! cell->at_boundary(nb))
	  {
					     // find active neighbor
	    cell_iterator neighbor = cell->neighbor(nb);
	    while (neighbor->has_children())
	      neighbor = neighbor->child(nb==0 ? 1 : 0);

					     // get dofs on it
            const unsigned int n_dofs_on_neighbor
              = neighbor->get_fe().dofs_per_cell;
            neighbor_dof_indices.resize (n_dofs_on_neighbor);
	    neighbor->get_dof_indices (neighbor_dof_indices);

					     // compute couplings
	    for (unsigned int i=0; i<n_dofs_on_this_cell; ++i)
	      for (unsigned int j=0; j<n_dofs_on_neighbor; ++j)
		sparsity.add (local_dof_indices[i], neighbor_dof_indices[j]);
	  };
    };
}

#endif


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
      Assert (ii < fe.n_components(),
	      ExcInternalError());

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
	  Assert (jj < fe.n_components(),
		  ExcInternalError());          

	  dof_couplings(i,j) = component_couplings(ii,jj);
	}
    }
  return dof_couplings;
}



template <class DH, class SparsityPattern>
void
DoFTools::make_flux_sparsity_pattern (
  const DH                &dof,
  SparsityPattern         &sparsity,
  const Table<2,Coupling> &int_mask,
  const Table<2,Coupling> &flux_mask)
{
  const unsigned int n_dofs = dof.n_dofs();
  const FiniteElement<DH::dimension> &fe = dof.get_fe();
  const unsigned int n_comp = fe.n_components();
  
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
  
  const unsigned int total_dofs = fe.dofs_per_cell;
  std::vector<unsigned int> dofs_on_this_cell(total_dofs);
  std::vector<unsigned int> dofs_on_other_cell(total_dofs);
  Table<2,bool> support_on_face(
    total_dofs, GeometryInfo<DH::dimension>::faces_per_cell);
  
  typename DH::active_cell_iterator cell = dof.begin_active(),
				    endc = dof.end();
  
  const Table<2,Coupling>
    int_dof_mask  = dof_couplings_from_component_couplings(fe, int_mask),
    flux_dof_mask = dof_couplings_from_component_couplings(fe, flux_mask);
  
  for (unsigned int i=0; i<total_dofs; ++i)
    for (unsigned int f=0; f<GeometryInfo<DH::dimension>::faces_per_cell;++f)
      support_on_face(i,f) = fe.has_support_on_face(i,f);
  
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
  
  for (; cell!=endc; ++cell)
    {
      cell->get_dof_indices (dofs_on_this_cell);
				       // make sparsity pattern for this cell
      for (unsigned int i=0; i<total_dofs; ++i)
	for (unsigned int j=0; j<total_dofs; ++j)
	  if (int_dof_mask(i,j) != none)
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
	      for (unsigned int i=0; i<total_dofs; ++i)
		{
		  const bool i_non_zero_i = support_on_face (i, face);
		  for (unsigned int j=0; j<total_dofs; ++j)
		    {
		      const bool j_non_zero_i = support_on_face (j, face);
		      
		      if (flux_dof_mask(i,j) == always)
                        sparsity.add (dofs_on_this_cell[i],
                                      dofs_on_this_cell[j]);
		      if (flux_dof_mask(i,j) == nonzero
			  && i_non_zero_i && j_non_zero_i)
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
		      for (unsigned int i=0; i<total_dofs; ++i)
			{
			  const bool i_non_zero_i = support_on_face (i, face);
			  const bool i_non_zero_e = support_on_face (i, neighbor_face);
			  for (unsigned int j=0; j<total_dofs; ++j)
			    {
			      const bool j_non_zero_i = support_on_face (j, face);
			      const bool j_non_zero_e  =support_on_face (j, neighbor_face);
			      if (flux_dof_mask(i,j) == always)
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
			      if (flux_dof_mask(i,j) == nonzero)
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
			      
			      if (flux_dof_mask(j,i) == always)
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
			      if (flux_dof_mask(j,i) == nonzero)
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
		  for (unsigned int i=0; i<total_dofs; ++i)
		    {
		      const bool i_non_zero_i = support_on_face (i, face);
		      const bool i_non_zero_e = support_on_face (i, neighbor_face);
		      for (unsigned int j=0; j<total_dofs; ++j)
			{
			  const bool j_non_zero_i = support_on_face (j, face);
			  const bool j_non_zero_e = support_on_face (j, neighbor_face);
			  if (flux_dof_mask(i,j) == always)
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
			  if (flux_dof_mask(i,j) == nonzero)
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

			  if (flux_dof_mask(j,i) == always)
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
			  if (flux_dof_mask(j,i) == nonzero)
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
	Assert (master_dof_mask.size() == fe1.dofs_per_face,
		ExcInternalError());

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

	Assert (index == fe1.dofs_per_face, ExcInternalError());
	Assert (master_dof_list.size() == fe2.dofs_per_face,
		ExcInternalError());
	
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
					   std_cxx0x::shared_ptr<std::vector<bool> > &master_dof_mask)
      {
	if (master_dof_mask == 0)
	  {
	    master_dof_mask = std_cxx0x::shared_ptr<std::vector<bool> >
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
                                       std_cxx0x::shared_ptr<FullMatrix<double> > &matrix)
      {
        if (matrix == 0)
          {
            matrix = std_cxx0x::shared_ptr<FullMatrix<double> >
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
                                          std_cxx0x::shared_ptr<FullMatrix<double> > &matrix)
      {
        if (matrix == 0)
          {
            matrix = std_cxx0x::shared_ptr<FullMatrix<double> >
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
					     std_cxx0x::shared_ptr<std::pair<FullMatrix<double>,FullMatrix<double> > > &split_matrix)
      {
	Assert (master_dof_mask.size() == face_interpolation_matrix.m(),
		ExcInternalError());
	Assert (std::count (master_dof_mask.begin(), master_dof_mask.end(), true) ==
		static_cast<signed int>(face_interpolation_matrix.n()),
		ExcInternalError());
	
        if (split_matrix == 0)
          {
            split_matrix
	      = std_cxx0x::shared_ptr<std::pair<FullMatrix<double>,FullMatrix<double> > >
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

	    Assert (nth_master_dof == n_master_dofs, ExcInternalError());
	    Assert (nth_slave_dof == n_dofs-n_master_dofs, ExcInternalError());	    

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
      template <int dim, typename face_iterator>
      unsigned int
      get_most_dominating_subface_fe_index (const face_iterator &face)
      {

	const unsigned int spacedim=dim;
 
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
                      {
#ifdef WOLFGANG
                        std::cout << "   " << slave_dofs[row]
                                  << " -> " << face_constraints (row,i) << " * "
                                  << master_dofs[i]
                                  << std::endl;
#endif
                        constraints.add_entry (slave_dofs[row],
                                               master_dofs[i],
                                               face_constraints (row,i));
                      }
                }
            }      
      }

    }
    
    
#if deal_II_dimension == 1
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
#endif



#if deal_II_dimension == 2
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
	for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
	  if (cell->face(face)->has_children()) 
	    {
					       // so now we've found a
					       // face of an active
					       // cell that has
					       // children. that means
					       // that there are
					       // hanging nodes here.

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
		Assert (cell->face(face)->child(c)->n_active_fe_indices() == 1,
			ExcInternalError());

					       // right now, all that
					       // is implemented is
					       // the case that both
					       // sides use the same
					       // fe
	      for (unsigned int c=0; c<cell->face(face)->n_children(); ++c)
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
	      Assert (next_index == dofs_on_mother.size(),
		      ExcInternalError());
	  
	      next_index = 0;
	      for (unsigned int dof=0; dof!=fe.dofs_per_vertex; ++dof)
		dofs_on_children[next_index++]
		  = this_face->child(0)->vertex_dof_index(1,dof,fe_index);
	      for (unsigned int child=0; child<2; ++child)
		for (unsigned int dof=0; dof!=fe.dofs_per_line; ++dof)
		  dofs_on_children[next_index++]
		    = this_face->child(child)->dof_index(dof, fe_index);
	      Assert (next_index == dofs_on_children.size(),
		      ExcInternalError());
	  
					       // for each row in the constraint
					       // matrix for this line:
	      for (unsigned int row=0; row!=dofs_on_children.size(); ++row) 
		{
		  constraints.add_line (dofs_on_children[row]);
		  for (unsigned int i=0; i!=dofs_on_mother.size(); ++i)
		    constraints.add_entry (dofs_on_children[row],
					   dofs_on_mother[i],
					   fe.constraints()(row,i));
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
	      Assert (cell->face(face)->n_active_fe_indices() == 1,
		      ExcNotImplemented());
	      Assert (cell->face(face)
		      ->fe_index_is_active(cell->active_fe_index()) == true,
		      ExcInternalError());
	    } 
    }
#endif      


#if deal_II_dimension == 3
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
		Assert (cell->face(face)->child(c)->n_active_fe_indices() == 1,
			ExcInternalError());

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

	      const unsigned int
		n_dofs_on_mother   = (4*fe.dofs_per_vertex+
				      4*fe.dofs_per_line+
				      fe.dofs_per_quad),
		n_dofs_on_children = (5*fe.dofs_per_vertex+
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
	      Assert (next_index == dofs_on_mother.size(),
		      ExcInternalError());
	  
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
	      Assert (next_index == dofs_on_children.size(),
		      ExcInternalError());
	  
					       // for each row in the constraint
					       // matrix for this line:
	      for (unsigned int row=0; row!=dofs_on_children.size(); ++row) 
		{
		  constraints.add_line (dofs_on_children[row]);
		  for (unsigned int i=0; i!=dofs_on_mother.size(); ++i)
		    constraints.add_entry (dofs_on_children[row],
					   dofs_on_mother[i],
					   fe.constraints()(row,i));
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
	      Assert (cell->face(face)->n_active_fe_indices() == 1,
		      ExcNotImplemented());
	      Assert (cell->face(face)
		      ->fe_index_is_active(cell->active_fe_index()) == true,
		      ExcInternalError());
	    } 
    }
#endif

    
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
      Table<2,std_cxx0x::shared_ptr<FullMatrix<double> > >
        face_interpolation_matrices (n_finite_elements (dof_handler),
                                     n_finite_elements (dof_handler));
      Table<3,std_cxx0x::shared_ptr<FullMatrix<double> > >
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
      Table<2,std_cxx0x::shared_ptr<std::pair<FullMatrix<double>,FullMatrix<double> > > >
        split_face_interpolation_matrices (n_finite_elements (dof_handler),
					   n_finite_elements (dof_handler));

				       // finally, for each pair of finite
				       // elements, have a mask that states
				       // which of the degrees of freedom on
				       // the coarse side of a refined face
				       // will act as master dofs.
      Table<2,std_cxx0x::shared_ptr<std::vector<bool> > >
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
	      FiniteElementDomination::Domination
		mother_face_dominates = FiniteElementDomination::either_element_can_dominate;

              if (DoFHandlerSupportsDifferentFEs<DH>::value == true)
                for (unsigned int c=0; c<cell->face(face)->number_of_children(); ++c)
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
		    master_dofs.resize (cell->get_fe().dofs_per_face);

		    cell->face(face)->get_dof_indices (master_dofs,
                                                       cell->active_fe_index ());
		  
						     // Now create constraint matrix for
						     // the subfaces and assemble it.
		    for (unsigned int c=0; c<cell->face(face)->n_children(); ++c)
		      {
			const typename DH::active_face_iterator
			  subface = cell->face(face)->child(c);

			Assert (subface->n_active_fe_indices() == 1,
				ExcInternalError());

			const unsigned int
			  subface_fe_index = subface->nth_active_fe_index(0);

							 // Same procedure as for the
							 // mother cell. Extract the face
							 // DoFs from the cell DoFs.
			slave_dofs.resize (subface->get_fe(subface_fe_index)
					   .dofs_per_face);
			subface->get_dof_indices (slave_dofs, subface_fe_index);
					      
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
#ifdef WOLFGANG
			std::cout << "Constraints for cell=" << cell
				  << ", face=" << face
				  << ", subface=" << c
				  << std::endl;
#endif
			
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
                      = get_most_dominating_subface_fe_index<dim> (cell->face(face));
                    
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
		    
		    Assert (master_dofs.size() == dominating_fe.dofs_per_face,
			    ExcInternalError());
		    Assert (slave_dofs.size() ==
			    cell->get_fe().dofs_per_face - dominating_fe.dofs_per_face,
			    ExcInternalError());
		    
#ifdef WOLFGANG
		    std::cout << "Constraints for cell=" << cell
			      << ", face=" << face
			      << " (complicated case, mother)"
			      << std::endl;
#endif
		    filter_constraints (master_dofs,
					slave_dofs,
					constraint_matrix,
					constraints);



						     // next we have to
						     // deal with the
						     // subfaces. do as
						     // discussed in the
						     // paper
		    for (unsigned int sf=0;
			 sf<cell->face(face)->n_children(); ++sf)
		      {
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

#ifdef WOLFGANG
			std::cout << "Constraints for cell=" << cell
				  << ", face=" << face
				  << ", subface=" << sf
				  << " (complicated case, children)"
				  << std::endl;
#endif
			filter_constraints (master_dofs,
					    slave_dofs,
					    constraint_matrix,
					    constraints);
		      }
			
		    break;
		  }

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
#ifdef WOLFGANG		      
			std::cout << "p-constraints for cell=" << cell
				  << ", face=" << face << std::endl;
#endif
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
			      Assert (master_dofs[i] == slave_dofs[i],
				      ExcInternalError());
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
      std::vector<bool> component_dofs (dof_handler.n_dofs());
      std::vector<bool> component_mask (dof_handler.get_fe().n_components(),
					false);
      component_mask[component] = true;
      DoFTools::extract_dofs (dof_handler, component_mask, component_dofs);

      for (unsigned int i=0; i<dof_data.size(); ++i)
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
                                   // the work
  std::vector<unsigned int> indices(fe.dofs_per_cell);
  typename DoFHandler<dim,spacedim>::active_cell_iterator c;
  for (c=dof.begin_active(); c!=dof.end(); ++ c)
    {
      c->get_dof_indices(indices);
      for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
        selected_dofs[indices[i]] = local_selected_dofs[i];
    }
}



template <int dim, int spacedim>
void
DoFTools::extract_dofs (
  const hp::DoFHandler<dim,spacedim>   &/*dof*/,
  const std::vector<bool> &/*component_select*/,
  std::vector<bool>       &/*selected_dofs*/,
  const bool               /*count_by_blocks*/)
{
//TODO[WB]: implement this function  
  Assert (false, ExcNotImplemented());
/*  
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
                                   // the work
  std::vector<unsigned int> indices(fe.dofs_per_cell);
  typename DoFHandler<dim,spacedim>::active_cell_iterator c;
  for (c=dof.begin_active(); c!=dof.end(); ++ c)
    {
      c->get_dof_indices(indices);
      for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
        selected_dofs[indices[i]] = local_selected_dofs[i];
    }
*/
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


#if deal_II_dimension != 1

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


#else  // 1d


template <class DH>
void
DoFTools::extract_boundary_dofs (const DH                 &dof_handler,
				 const std::vector<bool>  &component_select,
				 std::vector<bool>        &selected_dofs,
				 const std::set<unsigned char> &boundary_indicators)
{
  Assert (component_select.size() == n_components(dof_handler),
	  ExcWrongSize (component_select.size(),
			n_components(dof_handler)));
	  
				   // clear and reset array by default
				   // values
  selected_dofs.clear ();
  selected_dofs.resize (dof_handler.n_dofs(), false);

				   // let's see whether we have to
				   // check for certain boundary
				   // indicators or whether we can
				   // accept all
  const bool check_left_vertex  = ((boundary_indicators.size() == 0) ||
				   (boundary_indicators.find(0) !=
				    boundary_indicators.end()));
  const bool check_right_vertex = ((boundary_indicators.size() == 0) ||
				   (boundary_indicators.find(1) !=
				    boundary_indicators.end()));

                                   // see whether we have to check
                                   // whether a certain vector
                                   // component is selected, or all
  const bool check_vector_component
    = (component_select != std::vector<bool>(component_select.size(),
                                             true));
  
				   // loop over coarse grid cells
  for (typename DH::cell_iterator cell=dof_handler.begin(0);
       cell!=dof_handler.end(0); ++cell)
    {
      const FiniteElement<1> &fe = cell->get_fe();
            
				       // check left-most vertex
      if (check_left_vertex)
	if (cell->neighbor(0) == dof_handler.end())
          {
					     // In 1D the number of DoFs
					     // on the faces should be
					     // equal to the number of DoFs
					     // on the vertices.
            Assert (fe.dofs_per_face == 
		    fe.dofs_per_vertex,
                    ExcInternalError());
            
            for (unsigned int i=0; i<fe.dofs_per_face; ++i)
              if (!check_vector_component)
                selected_dofs[cell->vertex_dof_index(0,i)] = true;
              else
                                                 // check
                                                 // component. make sure
                                                 // we don't ask the
                                                 // wrong question
                                                 // (leading to an
                                                 // exception) in case
                                                 // the shape function
                                                 // is non-primitive. note
                                                 // that the face dof
                                                 // index i is also the
                                                 // cell dof index of a
                                                 // corresponding dof in 1d
                {
                  const unsigned int component =
                    (fe.is_primitive(i) ?
                     fe.face_system_to_component_index(i).first :
                     (std::find (fe.get_nonzero_components(i).begin(),
                                 fe.get_nonzero_components(i).end(),
                                 true)
                      -
                      fe.get_nonzero_components(i).begin()));
                  Assert (component < fe.n_components(),
                          ExcInternalError());
                 
                  if (component_select[component] == true)
                    selected_dofs[cell->vertex_dof_index(0,i)] = true;
                }
          }
      
				       // check right-most
				       // vertex. same procedure here
				       // as above
      if (check_right_vertex)
	if (cell->neighbor(1) == dof_handler.end())
          {
            Assert (fe.dofs_per_face ==
                    fe.dofs_per_vertex,
                    ExcInternalError());
            
            for (unsigned int i=0; i<fe.dofs_per_face; ++i)
              if (!check_vector_component)
                selected_dofs[cell->vertex_dof_index(1,i)] = true;
              else
                {
                  const unsigned int component =
                    (fe.is_primitive(i) ?
                     fe.face_system_to_component_index(i).first :
                     (std::find (fe.get_nonzero_components(i).begin(),
                                 fe.get_nonzero_components(i).end(),
                                 true)
                      -
                      fe.get_nonzero_components(i).begin()));
                  Assert (component < fe.n_components(),
                          ExcInternalError());
                 
                  if (component_select[component] == true)
                    selected_dofs[cell->vertex_dof_index(1,i)] = true;
                }
          }
    }
}


#endif


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
DoFTools::extract_subdomain_dofs (const DH           &dof_handler,
				  const unsigned int  subdomain_id,
				  std::vector<bool>  &selected_dofs)
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
DoFTools::extract_constant_modes (const DH                        &dof_handler,
				  const std::vector<bool>         &component_select,
				  std::vector<std::vector<bool> > &constant_modes)
{
  const unsigned int n_components = dof_handler.get_fe().n_components();
  Assert (n_components == component_select.size(),
	  ExcDimensionMismatch(n_components,
			       component_select.size()));

				 // First count the number of dofs
				 // in the current component.
  unsigned int n_components_selected = 0;
  std::vector<unsigned int> component_list (n_components, 0);
  for (unsigned int d=0; d<n_components; ++d)
    {
      component_list[d] = component_select[d];
      n_components_selected += component_select[d];
    }

  std::vector<unsigned int> dofs_per_block(2);
  count_dofs_per_block(dof_handler, dofs_per_block, component_list);

  const unsigned int n_u = dofs_per_block[1];
  std::vector<bool> selection_dof_list (dof_handler.n_dofs(), false);
  std::vector<bool> temporary_dof_list (dof_handler.n_dofs(), false);
  extract_dofs (dof_handler, component_select, selection_dof_list);

  constant_modes.resize (n_components_selected, std::vector<bool>(n_u, false));
  
  for (unsigned int component=0, component_used=0; 
       component < n_components; ++component, ++component_used)
    if (component_select[component])
      {
	std::vector<bool> selection_mask (n_components, false);
	selection_mask[component] = true;
	extract_dofs (dof_handler, selection_mask, temporary_dof_list);
  
	unsigned int counter = 0;
	for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
	{
	  if (selection_dof_list[i])
	    {
	      if (temporary_dof_list[i])
		constant_modes [component][counter] = true;
	      else
		constant_modes [component][counter] = false;

	      ++counter;
	    }
	}
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
				     std::vector<unsigned int> &subdomain_association)
{
  Assert(subdomain_association.size() == dof_handler.n_dofs(),
	 ExcDimensionMismatch(subdomain_association.size(),
                              dof_handler.n_dofs()));

                                   // preset all values by an invalid value
  std::fill_n (subdomain_association.begin(), dof_handler.n_dofs(),
               numbers::invalid_unsigned_int);

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
      const unsigned int subdomain_id = cell->subdomain_id();
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
                     numbers::invalid_unsigned_int)
          == subdomain_association.end(),
          ExcInternalError());
}



template <class DH>
unsigned int
DoFTools::count_dofs_with_subdomain_association (const DH           &dof_handler,
						 const unsigned int  subdomain)
{
                                   // in debug mode, make sure that there are
                                   // some cells at least with this subdomain
                                   // id
#ifdef DEBUG
  {
    bool found = false;
    for (typename Triangulation<DH::dimension, DH::space_dimension>::active_cell_iterator
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

  std::vector<unsigned int> subdomain_association (dof_handler.n_dofs());
  get_subdomain_association (dof_handler, subdomain_association);

  return std::count (subdomain_association.begin(),
                     subdomain_association.end(),
                     subdomain);
}



template <class DH>
void
DoFTools::count_dofs_with_subdomain_association (const DH           &dof_handler,
						 const unsigned int  subdomain,
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

  std::vector<unsigned int> subdomain_association (dof_handler.n_dofs());
  get_subdomain_association (dof_handler, subdomain_association);

  std::vector<bool> component_association (dof_handler.n_dofs());
  for (unsigned int c=0; c<dof_handler.get_fe().n_components(); ++c)
    {
      std::vector<bool> component_mask (dof_handler.get_fe().n_components(), false);
      component_mask[c] = true;
      extract_dofs (dof_handler, component_mask, component_association);

      for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
	if ((subdomain_association[i] == subdomain) &&
	    (component_association[i] == true))
	  ++n_dofs_on_subdomain[c];
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
  const unsigned int n_components = fe.n_components();
  dofs_per_component.resize (n_components);
  std::fill (dofs_per_component.begin(), dofs_per_component.end(), 0U);
  
				   // If the empty vector was given as
				   // default argument, set up this
				   // vector as identity.
  if (target_component.size()==0)
    {
      target_component.resize(n_components);
      for (unsigned int i=0;i<n_components;++i)
	target_component[i] = i;
    }
  
  Assert(target_component.size()==n_components,
	 ExcDimensionMismatch(target_component.size(),n_components));

				   // special case for only one
				   // component. treat this first
				   // since it does not require any
				   // computations
  if (n_components == 1)
    {
      dofs_per_component[0] = dof_handler.n_dofs();
      return;
    }

      
				   // otherwise determine the number
				   // of dofs in each component
				   // separately. do so in parallel
  std::vector<std::vector<bool> >
    dofs_in_component (n_components,
                       std::vector<bool>(dof_handler.n_dofs(), false));
  std::vector<std::vector<bool> >
    component_select (n_components,
                      std::vector<bool>(n_components, false));
  Threads::ThreadGroup<> threads;
  for (unsigned int i=0; i<n_components; ++i)
    {
      void (*fun_ptr) (const DoFHandler<dim,spacedim>   &,
		       const std::vector<bool> &,
		       std::vector<bool>       &,
		       bool)
        = &DoFTools::template extract_dofs<dim>;
      component_select[i][i] = true;
      threads += Threads::spawn (fun_ptr)(dof_handler, component_select[i],
                                          dofs_in_component[i], false);
    };
  threads.join_all ();

				   // next count what we got
  unsigned int component = 0;
  for (unsigned int b=0;b<fe.n_base_elements();++b)
    {
      const FiniteElement<dim,spacedim>& base = fe.base_element(b);
				       // Dimension of base element
      unsigned int d = base.n_components();
      
      for (unsigned int m=0;m<fe.element_multiplicity(b);++m)
	{
	  for (unsigned int dd=0;dd<d;++dd)
	    {
	      if (base.is_primitive() || (!only_once || dd==0))
		dofs_per_component[target_component[component]]
		  += std::count(dofs_in_component[component].begin(),
				dofs_in_component[component].end(),
				true);
	      ++component;
	    }
	}
    }
  
				   // finally sanity check. this is
				   // only valid if the finite element
				   // is actually primitive, so
				   // exclude other elements from this
  Assert (!dof_handler.get_fe().is_primitive()
          ||
          (std::accumulate (dofs_per_component.begin(),
                            dofs_per_component.end(), 0U)
           == dof_handler.n_dofs()),
	  ExcInternalError());
}


template <int dim, int spacedim>
void
DoFTools::count_dofs_per_block (
  const DoFHandler<dim,spacedim>&     dof_handler,
  std::vector<unsigned int>& dofs_per_block,
  std::vector<unsigned int>  target_block)
{
  const FiniteElement<dim,spacedim>& fe = dof_handler.get_fe();
  const unsigned int n_blocks = fe.n_blocks();
  dofs_per_block.resize (n_blocks);
  std::fill (dofs_per_block.begin(), dofs_per_block.end(), 0U);
  
				   // If the empty vector was given as
				   // default argument, set up this
				   // vector as identity.
  if (target_block.size()==0)
    {
      target_block.resize(n_blocks);
      for (unsigned int i=0;i<n_blocks;++i)
	target_block[i] = i;
    }
  
  Assert(target_block.size()==n_blocks,
	 ExcDimensionMismatch(target_block.size(),n_blocks));

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
				   // separately. do so in parallel
  std::vector<std::vector<bool> >
    dofs_in_block (n_blocks, std::vector<bool>(dof_handler.n_dofs(), false));
  std::vector<std::vector<bool> >
    block_select (n_blocks, std::vector<bool>(n_blocks, false));
  Threads::ThreadGroup<> threads;
  for (unsigned int i=0; i<n_blocks; ++i)
    {
      void (*fun_ptr) (const DoFHandler<dim,spacedim>   &,
		       const std::vector<bool> &,
		       std::vector<bool>       &,
		       bool)
        = &DoFTools::template extract_dofs<dim>;
      block_select[i][i] = true;
      threads += Threads::spawn (fun_ptr)(dof_handler, block_select[i],
                                          dofs_in_block[i], true);
    };
  threads.join_all ();

				   // next count what we got
  for (unsigned int block=0;block<fe.n_blocks();++block)
    dofs_per_block[target_block[block]]
      += std::count(dofs_in_block[block].begin(),
		    dofs_in_block[block].end(),
		    true);
}


template <int dim, int spacedim>
void
DoFTools::count_dofs_per_component (
  const DoFHandler<dim,spacedim>&     dof_handler,
  std::vector<unsigned int>& dofs_per_component,
  std::vector<unsigned int>  target_component)
{
  count_dofs_per_component (dof_handler, dofs_per_component, false, target_component);
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
      std::vector<Vector<double> > parameter_dofs (coarse_dofs_per_cell_component,
						   Vector<double>(fine_dofs_per_cell));
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
      const std::vector<Vector<double> > &parameter_dofs,
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
      Vector<double> global_parameter_representation (n_fine_dofs);
  
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
      const std::vector<Vector<double> > &parameter_dofs,
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

      Threads::ThreadGroup<> threads;
      void (*fun_ptr) (const dealii::DoFHandler<dim,spacedim>              &,
		       const unsigned int                  ,
		       const InterGridMap<dealii::DoFHandler<dim,spacedim> > &,
		       const std::vector<Vector<double> > &,
		       const std::vector<int>             &,
		       std::vector<std::map<unsigned int, float> > &,
		       const typename dealii::DoFHandler<dim,spacedim>::active_cell_iterator &,
		       const typename dealii::DoFHandler<dim,spacedim>::active_cell_iterator &)
	= &compute_intergrid_weights_3<dim>;
      for (unsigned int i=0; i<multithread_info.n_default_threads; ++i)
	threads += Threads::spawn (fun_ptr)(coarse_grid, coarse_component,
					    coarse_to_fine_grid_map, parameter_dofs,
					    weight_mapping, weights,
					    cell_intervals[i].first,
					    cell_intervals[i].second);

				       // wait for the threads to finish
      threads.join_all ();
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



#if deal_II_dimension == 1


template <class DH>
void DoFTools::map_dof_to_boundary_indices (
  const DH                      &dof_handler,
  const std::set<unsigned char> &boundary_indicators,
  std::vector<unsigned int> &mapping)
{
  Assert (&dof_handler.get_fe() != 0, ExcNoFESelected());

  mapping.clear ();
  mapping.insert (mapping.end(), dof_handler.n_dofs(),
		  DH::invalid_dof_index);

  unsigned int next_free_index = 0;
  
				   // first check left, the right
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
      typename DH::cell_iterator cell = dof_handler.begin(0);
      while (!cell->at_boundary(direction))
	cell = cell->neighbor(direction);
      while (!cell->active())
	cell = cell->child(direction);

				       // next enumerate these degrees
				       // of freedom
      for (unsigned int i=0; i<cell->get_fe().dofs_per_vertex; ++i)
	mapping[cell->vertex_dof_index(direction,i)] = next_free_index++;
    };
}



template <>
void
DoFTools::map_dof_to_boundary_indices (const DoFHandler<1>       &dof_handler,
                                       std::vector<unsigned int> &mapping)
{
  Assert (&dof_handler.get_fe() != 0, ExcNoFESelected());

				   // in 1d, there are only 2 boundary
				   // indicators, so enumerate them
				   // and pass on to the other
				   // function
  std::set<unsigned char> boundary_indicators;
  boundary_indicators.insert (0U);
  boundary_indicators.insert (1U);

  map_dof_to_boundary_indices<DoFHandler<1> > (dof_handler, boundary_indicators, mapping);
}

#else


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

  Assert (next_boundary_index == dof_handler.n_boundary_dofs(),
	  ExcInternalError());
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

  Assert (next_boundary_index == dof_handler.n_boundary_dofs (boundary_indicators),
	  ExcInternalError());
}

#endif



template <int dim, int spacedim>
void
DoFTools::map_dofs_to_support_points (const Mapping<dim,spacedim>       &mapping,
				      const DoFHandler<dim,spacedim>    &dof_handler,
				      std::vector<Point<dim> > &support_points)
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
  FEValues<dim> fe_values (mapping, dof_handler.get_fe(),
			   q_dummy, update_quadrature_points);
  typename DoFHandler<dim,spacedim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);
  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      cell->get_dof_indices (local_dof_indices);
      const std::vector<Point<dim> > & points
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

#if deal_II_dimension == 1

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

                                   // check for boundary cells in both
                                   // directions to secure indicators
                                   // for the entire boundary, i.e. 0
                                   // is left boundary and 1 is right
                                   // boundary
  for (unsigned int direction = 0; direction < 2; ++direction)
    {

                                   // first find the outermost active
                                   // cell by traversing the coarse
                                   // grid to its end and then looking
                                   // for the children
      typename DH<dim,spacedim>::cell_iterator outermost_cell = dof.begin(0);
      while (outermost_cell->neighbor(direction).state() == IteratorState::valid)
	outermost_cell = outermost_cell->neighbor(direction);
      
      while (outermost_cell->has_children())
	outermost_cell = outermost_cell->child(direction);

                                   // then get the FE corresponding to
                                   // this cell
      const FiniteElement<dim,spacedim> &fe = outermost_cell->get_fe();

				   // set the component mask to either
				   // the original value or a vector
				   // of trues
      const std::vector<bool> component_mask ((component_mask_.size() == 0) ?
					      std::vector<bool> (fe.n_components(), true) :
					      component_mask_);
      Assert (std::count(component_mask.begin(), component_mask.end(), true) > 0,
	      VectorTools::ExcNoComponentSelected());

                                   // cast zero boundary constraints
                                   // onto a matrix if component_mask
                                   // == true
      for (unsigned int i=0; i<fe.dofs_per_vertex; ++i)
	if (component_mask[fe.face_system_to_component_index(i).first])
	  zero_boundary_constraints.add_line (outermost_cell->vertex_dof_index (direction, i));

    }
}

#else

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
  const bool                fe_is_system = (n_components != 1);

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

				   // we can presently deal only with
				   // primitive elements for boundary
				   // values. make sure that all shape
				   // functions that are non-zero for
				   // the components we are interested
				   // in, are in fact primitive
	for (unsigned int i=0; i<cell->get_fe().dofs_per_cell; ++i)
	  {
	    const std::vector<bool> &nonzero_component_array
	      = cell->get_fe().get_nonzero_components (i);
	    for (unsigned int c=0; c<n_components; ++c)
	      if ((nonzero_component_array[c] == true)
		  &&
		  (component_mask[c] == true))
		Assert (cell->get_fe().is_primitive (i),
			ExcMessage ("This function can only deal with requested boundary "
				    "values that correspond to primitive (scalar) base "
				    "elements"));
	  }

	typename DH<dim,spacedim>::face_iterator face = cell->face(face_no);
	if (face->boundary_indicator () == 0)
	                           // face is of the right component
	  {
				   // get indices and physical
				   // location on this face
 	    face_dofs.resize (fe.dofs_per_face);
 	    face->get_dof_indices (face_dofs, cell->active_fe_index());

	    if (fe_is_system)
	      {
				    // enter those dofs into the list
				    // that match the component
				    // signature.
		for (unsigned int i=0; i<face_dofs.size(); ++i)
                  {
                    unsigned int component;
                    if (fe.is_primitive())
                      component = fe.face_system_to_component_index(i).first;
                    else
                      {
                                    // non-primitive case. make sure
                                    // that this particular shape
                                    // function _is_ primitive, and
                                    // get at it's component. use
                                    // usual trick to transfer face
                                    // dof index to cell dof index
                        const unsigned int cell_i
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
                        Assert (cell_i < fe.dofs_per_cell, ExcInternalError());

                                    // make sure that if this is not a
                                    // primitive shape function, then
                                    // all the corresponding
                                    // components in the mask are not
                                    // set
                        if (!fe.is_primitive(cell_i))
                          for (unsigned int c=0; c<n_components; ++c)
                            if (fe.get_nonzero_components(cell_i)[c])
                              Assert (component_mask[c] == false,
                                      FETools::ExcFENotPrimitive());

                                    // pick the first possibly of more
                                    // than one non-zero component. if
                                    // the shape function is
                                    // non-primitive, then we will
                                    // ignore the result in the
                                    // following anyway, otherwise
                                    // there's only one non-zero
                                    // component which we will use
                        component = (std::find (fe.get_nonzero_components(cell_i).begin(),
                                                fe.get_nonzero_components(cell_i).end(),
                                                true)
                                     -
                                     fe.get_nonzero_components(cell_i).begin());
                      }

				    // cast zero boundary constraints onto
				    // a matrix
		    for (unsigned int i=0; i<fe.dofs_per_vertex; ++i)
		      if (component_mask[fe.face_system_to_component_index(i).first])
			zero_boundary_constraints.add_line (face_dofs[i]);

                  } 
	      }
	    else 
	      {
				    // get the one component that this
				    // function has and cast zero
				    // boundary constraints onto a
				    // matrix
		for (unsigned int i=0; i<fe.dofs_per_face; ++i)
		  if (component_mask[fe.face_system_to_component_index(i).first])
		    zero_boundary_constraints.add_line (face_dofs[i]);

	      }

	  }
      } 
}
    
#endif

// explicit instantiations

#include "dof_tools.inst"

// #if deal_II_dimension > 1
template void
DoFTools::make_boundary_sparsity_pattern<DoFHandler<deal_II_dimension>,SparsityPattern>
(const DoFHandler<deal_II_dimension>& dof,
 const std::vector<unsigned int>  &,
 SparsityPattern    &);
template void
DoFTools::make_boundary_sparsity_pattern<DoFHandler<deal_II_dimension>,CompressedSparsityPattern>
(const DoFHandler<deal_II_dimension>& dof,
 const std::vector<unsigned int>  &,
 CompressedSparsityPattern    &);
template void
DoFTools::make_boundary_sparsity_pattern<DoFHandler<deal_II_dimension>,CompressedSetSparsityPattern>
(const DoFHandler<deal_II_dimension>& dof,
 const std::vector<unsigned int>  &,
 CompressedSetSparsityPattern    &);
template void
DoFTools::make_boundary_sparsity_pattern<DoFHandler<deal_II_dimension>,CompressedSimpleSparsityPattern>
(const DoFHandler<deal_II_dimension>& dof,
 const std::vector<unsigned int>  &,
 CompressedSimpleSparsityPattern    &);
template void
DoFTools::make_boundary_sparsity_pattern<DoFHandler<deal_II_dimension>,BlockSparsityPattern>
(const DoFHandler<deal_II_dimension>& dof,
 const std::vector<unsigned int>  &,
 BlockSparsityPattern    &);
template void
DoFTools::make_boundary_sparsity_pattern<DoFHandler<deal_II_dimension>,BlockCompressedSparsityPattern>
(const DoFHandler<deal_II_dimension>& dof,
 const std::vector<unsigned int>  &,
 BlockCompressedSparsityPattern    &);
template void
DoFTools::make_boundary_sparsity_pattern<DoFHandler<deal_II_dimension>,BlockCompressedSetSparsityPattern>
(const DoFHandler<deal_II_dimension>& dof,
 const std::vector<unsigned int>  &,
 BlockCompressedSetSparsityPattern    &);
template void
DoFTools::make_boundary_sparsity_pattern<DoFHandler<deal_II_dimension>,BlockCompressedSimpleSparsityPattern>
(const DoFHandler<deal_II_dimension>& dof,
 const std::vector<unsigned int>  &,
 BlockCompressedSimpleSparsityPattern    &);

template void
DoFTools::make_boundary_sparsity_pattern<hp::DoFHandler<deal_II_dimension>,SparsityPattern>
(const hp::DoFHandler<deal_II_dimension>& dof,
 const std::vector<unsigned int>  &,
 SparsityPattern    &);
template void
DoFTools::make_boundary_sparsity_pattern<hp::DoFHandler<deal_II_dimension>,CompressedSparsityPattern>
(const hp::DoFHandler<deal_II_dimension>& dof,
 const std::vector<unsigned int>  &,
 CompressedSparsityPattern    &);
template void
DoFTools::make_boundary_sparsity_pattern<hp::DoFHandler<deal_II_dimension>,CompressedSetSparsityPattern>
(const hp::DoFHandler<deal_II_dimension>& dof,
 const std::vector<unsigned int>  &,
 CompressedSetSparsityPattern    &);
template void
DoFTools::make_boundary_sparsity_pattern<hp::DoFHandler<deal_II_dimension>,CompressedSimpleSparsityPattern>
(const hp::DoFHandler<deal_II_dimension>& dof,
 const std::vector<unsigned int>  &,
 CompressedSimpleSparsityPattern    &);
template void
DoFTools::make_boundary_sparsity_pattern<hp::DoFHandler<deal_II_dimension>,BlockSparsityPattern>
(const hp::DoFHandler<deal_II_dimension>& dof,
 const std::vector<unsigned int>  &,
 BlockSparsityPattern    &);
template void
DoFTools::make_boundary_sparsity_pattern<hp::DoFHandler<deal_II_dimension>,BlockCompressedSparsityPattern>
(const hp::DoFHandler<deal_II_dimension>& dof,
 const std::vector<unsigned int>  &,
 BlockCompressedSparsityPattern    &);
template void
DoFTools::make_boundary_sparsity_pattern<hp::DoFHandler<deal_II_dimension>,BlockCompressedSetSparsityPattern>
(const hp::DoFHandler<deal_II_dimension>& dof,
 const std::vector<unsigned int>  &,
 BlockCompressedSetSparsityPattern    &);
template void
DoFTools::make_boundary_sparsity_pattern<hp::DoFHandler<deal_II_dimension>,BlockCompressedSimpleSparsityPattern>
(const hp::DoFHandler<deal_II_dimension>& dof,
 const std::vector<unsigned int>  &,
 BlockCompressedSimpleSparsityPattern    &);


template void
DoFTools::make_boundary_sparsity_pattern<DoFHandler<deal_II_dimension>,SparsityPattern>
(const DoFHandler<deal_II_dimension>& dof,
 const FunctionMap<deal_II_dimension>::type  &boundary_indicators,
 const std::vector<unsigned int>  &dof_to_boundary_mapping,
 SparsityPattern    &sparsity);
template void
DoFTools::make_boundary_sparsity_pattern<DoFHandler<deal_II_dimension>,CompressedSparsityPattern>
(const DoFHandler<deal_II_dimension>& dof,
 const FunctionMap<deal_II_dimension>::type  &boundary_indicators,
 const std::vector<unsigned int>  &dof_to_boundary_mapping,
 CompressedSparsityPattern    &sparsity);
template void
DoFTools::make_boundary_sparsity_pattern<DoFHandler<deal_II_dimension>,CompressedSetSparsityPattern>
(const DoFHandler<deal_II_dimension>& dof,
 const FunctionMap<deal_II_dimension>::type  &boundary_indicators,
 const std::vector<unsigned int>  &dof_to_boundary_mapping,
 CompressedSetSparsityPattern    &sparsity);
template void
DoFTools::make_boundary_sparsity_pattern<DoFHandler<deal_II_dimension>,CompressedSimpleSparsityPattern>
(const DoFHandler<deal_II_dimension>& dof,
 const FunctionMap<deal_II_dimension>::type  &boundary_indicators,
 const std::vector<unsigned int>  &dof_to_boundary_mapping,
 CompressedSimpleSparsityPattern    &sparsity);
template void
DoFTools::make_boundary_sparsity_pattern<DoFHandler<deal_II_dimension>,BlockSparsityPattern>
(const DoFHandler<deal_II_dimension>& dof,
 const FunctionMap<deal_II_dimension>::type  &boundary_indicators,
 const std::vector<unsigned int>  &dof_to_boundary_mapping,
 BlockSparsityPattern    &sparsity);
template void
DoFTools::make_boundary_sparsity_pattern<DoFHandler<deal_II_dimension>,BlockCompressedSparsityPattern>
(const DoFHandler<deal_II_dimension>& dof,
 const FunctionMap<deal_II_dimension>::type  &boundary_indicators,
 const std::vector<unsigned int>  &dof_to_boundary_mapping,
 BlockCompressedSparsityPattern    &sparsity);
template void
DoFTools::make_boundary_sparsity_pattern<DoFHandler<deal_II_dimension>,BlockCompressedSetSparsityPattern>
(const DoFHandler<deal_II_dimension>& dof,
 const FunctionMap<deal_II_dimension>::type  &boundary_indicators,
 const std::vector<unsigned int>  &dof_to_boundary_mapping,
 BlockCompressedSetSparsityPattern    &sparsity);
template void
DoFTools::make_boundary_sparsity_pattern<DoFHandler<deal_II_dimension>,BlockCompressedSimpleSparsityPattern>
(const DoFHandler<deal_II_dimension>& dof,
 const FunctionMap<deal_II_dimension>::type  &boundary_indicators,
 const std::vector<unsigned int>  &dof_to_boundary_mapping,
 BlockCompressedSimpleSparsityPattern    &sparsity);

template void
DoFTools::make_boundary_sparsity_pattern<hp::DoFHandler<deal_II_dimension>,SparsityPattern>
(const hp::DoFHandler<deal_II_dimension>& dof,
 const FunctionMap<deal_II_dimension>::type  &boundary_indicators,
 const std::vector<unsigned int>  &dof_to_boundary_mapping,
 SparsityPattern    &sparsity);
template void
DoFTools::make_boundary_sparsity_pattern<hp::DoFHandler<deal_II_dimension>,CompressedSparsityPattern>
(const hp::DoFHandler<deal_II_dimension>& dof,
 const FunctionMap<deal_II_dimension>::type  &boundary_indicators,
 const std::vector<unsigned int>  &dof_to_boundary_mapping,
 CompressedSparsityPattern    &sparsity);
template void
DoFTools::make_boundary_sparsity_pattern<hp::DoFHandler<deal_II_dimension>,CompressedSetSparsityPattern>
(const hp::DoFHandler<deal_II_dimension>& dof,
 const FunctionMap<deal_II_dimension>::type  &boundary_indicators,
 const std::vector<unsigned int>  &dof_to_boundary_mapping,
 CompressedSetSparsityPattern    &sparsity);
template void
DoFTools::make_boundary_sparsity_pattern<hp::DoFHandler<deal_II_dimension>,CompressedSimpleSparsityPattern>
(const hp::DoFHandler<deal_II_dimension>& dof,
 const FunctionMap<deal_II_dimension>::type  &boundary_indicators,
 const std::vector<unsigned int>  &dof_to_boundary_mapping,
 CompressedSimpleSparsityPattern    &sparsity);
template void
DoFTools::make_boundary_sparsity_pattern<hp::DoFHandler<deal_II_dimension>,BlockSparsityPattern>
(const hp::DoFHandler<deal_II_dimension>& dof,
 const FunctionMap<deal_II_dimension>::type  &boundary_indicators,
 const std::vector<unsigned int>  &dof_to_boundary_mapping,
 BlockSparsityPattern    &sparsity);
template void
DoFTools::make_boundary_sparsity_pattern<hp::DoFHandler<deal_II_dimension>,BlockCompressedSparsityPattern>
(const hp::DoFHandler<deal_II_dimension>& dof,
 const FunctionMap<deal_II_dimension>::type  &boundary_indicators,
 const std::vector<unsigned int>  &dof_to_boundary_mapping,
 BlockCompressedSparsityPattern    &sparsity);
template void
DoFTools::make_boundary_sparsity_pattern<hp::DoFHandler<deal_II_dimension>,BlockCompressedSetSparsityPattern>
(const hp::DoFHandler<deal_II_dimension>& dof,
 const FunctionMap<deal_II_dimension>::type  &boundary_indicators,
 const std::vector<unsigned int>  &dof_to_boundary_mapping,
 BlockCompressedSetSparsityPattern    &sparsity);
template void
DoFTools::make_boundary_sparsity_pattern<hp::DoFHandler<deal_II_dimension>,BlockCompressedSimpleSparsityPattern>
(const hp::DoFHandler<deal_II_dimension>& dof,
 const FunctionMap<deal_II_dimension>::type  &boundary_indicators,
 const std::vector<unsigned int>  &dof_to_boundary_mapping,
 BlockCompressedSimpleSparsityPattern    &sparsity);


template void
DoFTools::make_flux_sparsity_pattern<DoFHandler<deal_II_dimension>,SparsityPattern>
(const DoFHandler<deal_II_dimension> &dof,
 SparsityPattern    &sparsity);
template void
DoFTools::make_flux_sparsity_pattern<DoFHandler<deal_II_dimension>,CompressedSparsityPattern>
(const DoFHandler<deal_II_dimension> &dof,
 CompressedSparsityPattern    &sparsity);
template void
DoFTools::make_flux_sparsity_pattern<DoFHandler<deal_II_dimension>,CompressedSetSparsityPattern>
(const DoFHandler<deal_II_dimension> &dof,
 CompressedSetSparsityPattern    &sparsity);
template void
DoFTools::make_flux_sparsity_pattern<DoFHandler<deal_II_dimension>,CompressedSimpleSparsityPattern>
(const DoFHandler<deal_II_dimension> &dof,
 CompressedSimpleSparsityPattern    &sparsity);
template void
DoFTools::make_flux_sparsity_pattern<DoFHandler<deal_II_dimension>,BlockSparsityPattern>
(const DoFHandler<deal_II_dimension> &dof,
 BlockSparsityPattern    &sparsity);
template void
DoFTools::make_flux_sparsity_pattern<DoFHandler<deal_II_dimension>,BlockCompressedSparsityPattern>
(const DoFHandler<deal_II_dimension> &dof,
 BlockCompressedSparsityPattern    &sparsity);
template void
DoFTools::make_flux_sparsity_pattern<DoFHandler<deal_II_dimension>,BlockCompressedSetSparsityPattern>
(const DoFHandler<deal_II_dimension> &dof,
 BlockCompressedSetSparsityPattern    &sparsity);
template void
DoFTools::make_flux_sparsity_pattern<DoFHandler<deal_II_dimension>,BlockCompressedSimpleSparsityPattern>
(const DoFHandler<deal_II_dimension> &dof,
 BlockCompressedSimpleSparsityPattern    &sparsity);

template void
DoFTools::make_flux_sparsity_pattern<hp::DoFHandler<deal_II_dimension>,SparsityPattern>
(const hp::DoFHandler<deal_II_dimension> &dof,
 SparsityPattern    &sparsity);
template void
DoFTools::make_flux_sparsity_pattern<hp::DoFHandler<deal_II_dimension>,CompressedSparsityPattern>
(const hp::DoFHandler<deal_II_dimension> &dof,
 CompressedSparsityPattern    &sparsity);
template void
DoFTools::make_flux_sparsity_pattern<hp::DoFHandler<deal_II_dimension>,CompressedSetSparsityPattern>
(const hp::DoFHandler<deal_II_dimension> &dof,
 CompressedSetSparsityPattern    &sparsity);
template void
DoFTools::make_flux_sparsity_pattern<hp::DoFHandler<deal_II_dimension>,CompressedSimpleSparsityPattern>
(const hp::DoFHandler<deal_II_dimension> &dof,
 CompressedSimpleSparsityPattern    &sparsity);
template void
DoFTools::make_flux_sparsity_pattern<hp::DoFHandler<deal_II_dimension>,BlockSparsityPattern>
(const hp::DoFHandler<deal_II_dimension> &dof,
 BlockSparsityPattern    &sparsity);
template void
DoFTools::make_flux_sparsity_pattern<hp::DoFHandler<deal_II_dimension>,BlockCompressedSparsityPattern>
(const hp::DoFHandler<deal_II_dimension> &dof,
 BlockCompressedSparsityPattern    &sparsity);
template void
DoFTools::make_flux_sparsity_pattern<hp::DoFHandler<deal_II_dimension>,BlockCompressedSetSparsityPattern>
(const hp::DoFHandler<deal_II_dimension> &dof,
 BlockCompressedSetSparsityPattern    &sparsity);
template void
DoFTools::make_flux_sparsity_pattern<hp::DoFHandler<deal_II_dimension>,BlockCompressedSimpleSparsityPattern>
(const hp::DoFHandler<deal_II_dimension> &dof,
 BlockCompressedSimpleSparsityPattern    &sparsity);


#if deal_II_dimension > 1
template void
DoFTools::make_flux_sparsity_pattern<DoFHandler<deal_II_dimension>,SparsityPattern>
(const DoFHandler<deal_II_dimension> &dof,
 SparsityPattern    &,
 const Table<2,Coupling>&,
 const Table<2,Coupling>&);
template void
DoFTools::make_flux_sparsity_pattern<DoFHandler<deal_II_dimension>,CompressedSparsityPattern>
(const DoFHandler<deal_II_dimension> &dof,
 CompressedSparsityPattern    &,
 const Table<2,Coupling>&,
 const Table<2,Coupling>&);
template void
DoFTools::make_flux_sparsity_pattern<DoFHandler<deal_II_dimension>,CompressedSetSparsityPattern>
(const DoFHandler<deal_II_dimension> &dof,
 CompressedSetSparsityPattern    &,
 const Table<2,Coupling>&,
 const Table<2,Coupling>&);
template void
DoFTools::make_flux_sparsity_pattern<DoFHandler<deal_II_dimension>,CompressedSimpleSparsityPattern>
(const DoFHandler<deal_II_dimension> &dof,
 CompressedSimpleSparsityPattern    &,
 const Table<2,Coupling>&,
 const Table<2,Coupling>&);
template void
DoFTools::make_flux_sparsity_pattern<DoFHandler<deal_II_dimension>,BlockSparsityPattern>
(const DoFHandler<deal_II_dimension> &dof,
 BlockSparsityPattern    &,
 const Table<2,Coupling>&,
 const Table<2,Coupling>&);
template void
DoFTools::make_flux_sparsity_pattern<DoFHandler<deal_II_dimension>,BlockCompressedSparsityPattern>
(const DoFHandler<deal_II_dimension> &dof,
 BlockCompressedSparsityPattern    &,
 const Table<2,Coupling>&,
 const Table<2,Coupling>&);
template void
DoFTools::make_flux_sparsity_pattern<DoFHandler<deal_II_dimension>,BlockCompressedSetSparsityPattern>
(const DoFHandler<deal_II_dimension> &dof,
 BlockCompressedSetSparsityPattern    &,
 const Table<2,Coupling>&,
 const Table<2,Coupling>&);
template void
DoFTools::make_flux_sparsity_pattern<DoFHandler<deal_II_dimension>,BlockCompressedSimpleSparsityPattern>
(const DoFHandler<deal_II_dimension> &dof,
 BlockCompressedSimpleSparsityPattern    &,
 const Table<2,Coupling>&,
 const Table<2,Coupling>&);
#endif



template
void
DoFTools::distribute_cell_to_dof_vector<DoFHandler<deal_II_dimension> >
(const DoFHandler<deal_II_dimension> &dof_handler,
 const Vector<float> &cell_data,
 Vector<double>      &dof_data,
 const unsigned int   component);
template
void
DoFTools::distribute_cell_to_dof_vector<DoFHandler<deal_II_dimension> >
(const DoFHandler<deal_II_dimension> &dof_handler,
 const Vector<double> &cell_data,
 Vector<double>       &dof_data,
 const unsigned int    component);

template
void
DoFTools::distribute_cell_to_dof_vector<hp::DoFHandler<deal_II_dimension> >
(const hp::DoFHandler<deal_II_dimension> &dof_handler,
 const Vector<float> &cell_data,
 Vector<double>      &dof_data,
 const unsigned int   component);
template
void
DoFTools::distribute_cell_to_dof_vector<hp::DoFHandler<deal_II_dimension> >
(const hp::DoFHandler<deal_II_dimension> &dof_handler,
 const Vector<double> &cell_data,
 Vector<double>       &dof_data,
 const unsigned int    component);


template void DoFTools::extract_dofs<deal_II_dimension>
(const DoFHandler<deal_II_dimension>&,
 const std::vector<bool>&, std::vector<bool>&, bool);

template void DoFTools::extract_dofs<deal_II_dimension>
(const hp::DoFHandler<deal_II_dimension>&,
 const std::vector<bool>&, std::vector<bool>&, bool);

template void DoFTools::extract_level_dofs<deal_II_dimension>
(const unsigned int level, const MGDoFHandler<deal_II_dimension>&,
 const std::vector<bool>&, std::vector<bool>&, bool);

template
void
DoFTools::extract_boundary_dofs<DoFHandler<deal_II_dimension> >
(const DoFHandler<deal_II_dimension> &,
 const std::vector<bool>                  &,
 std::vector<bool>                        &,
 const std::set<unsigned char> &);
template
void
DoFTools::extract_boundary_dofs<hp::DoFHandler<deal_II_dimension> >
(const hp::DoFHandler<deal_II_dimension> &,
 const std::vector<bool>                  &,
 std::vector<bool>                        &,
 const std::set<unsigned char> &);

template
void
DoFTools::extract_hanging_node_dofs
(const DoFHandler<deal_II_dimension> &dof_handler,
 std::vector<bool>     &selected_dofs);

template
void
DoFTools::extract_subdomain_dofs<DoFHandler<deal_II_dimension> >
(const DoFHandler<deal_II_dimension> &dof_handler,
 const unsigned int     subdomain_id,
 std::vector<bool>     &selected_dofs);
template
void
DoFTools::extract_subdomain_dofs<hp::DoFHandler<deal_II_dimension> >
(const hp::DoFHandler<deal_II_dimension> &dof_handler,
 const unsigned int     subdomain_id,
 std::vector<bool>     &selected_dofs);

template
void
DoFTools::extract_constant_modes<DoFHandler<deal_II_dimension> >
(const DoFHandler<deal_II_dimension> &dof_handler,
 const std::vector<bool> &selected_components,
 std::vector<std::vector<bool> > &constant_modes);

template
void
DoFTools::get_active_fe_indices<DoFHandler<deal_II_dimension> >
(const DoFHandler<deal_II_dimension> &dof_handler,
 std::vector<unsigned int> &active_fe_indices);

template
void
DoFTools::get_active_fe_indices<hp::DoFHandler<deal_II_dimension> >
(const hp::DoFHandler<deal_II_dimension> &dof_handler,
 std::vector<unsigned int> &active_fe_indices);

template
void
DoFTools::get_subdomain_association<DoFHandler<deal_II_dimension> >
(const DoFHandler<deal_II_dimension> &dof_handler,
 std::vector<unsigned int>           &subdomain_association);
template
void
DoFTools::get_subdomain_association<hp::DoFHandler<deal_II_dimension> >
(const hp::DoFHandler<deal_II_dimension> &dof_handler,
 std::vector<unsigned int>           &subdomain_association);


template
unsigned int
DoFTools::count_dofs_with_subdomain_association<DoFHandler<deal_II_dimension> >
(const DoFHandler<deal_II_dimension> &,
 const unsigned int);
template
void
DoFTools::count_dofs_with_subdomain_association<DoFHandler<deal_II_dimension> >
(const DoFHandler<deal_II_dimension> &,
 const unsigned int,
 std::vector<unsigned int> &);
template
unsigned int
DoFTools::count_dofs_with_subdomain_association<hp::DoFHandler<deal_II_dimension> >
(const hp::DoFHandler<deal_II_dimension> &,
 const unsigned int);
template
void
DoFTools::count_dofs_with_subdomain_association<hp::DoFHandler<deal_II_dimension> >
(const hp::DoFHandler<deal_II_dimension> &,
 const unsigned int,
 std::vector<unsigned int> &);
template
unsigned int
DoFTools::count_dofs_with_subdomain_association<MGDoFHandler<deal_II_dimension> >
(const MGDoFHandler<deal_II_dimension> &,
 const unsigned int);
template
void
DoFTools::count_dofs_with_subdomain_association<MGDoFHandler<deal_II_dimension> >
(const MGDoFHandler<deal_II_dimension> &,
 const unsigned int,
 std::vector<unsigned int> &);


template
void
DoFTools::count_dofs_per_component<deal_II_dimension> (
  const DoFHandler<deal_II_dimension>&,
  std::vector<unsigned int>&, bool, std::vector<unsigned int>);

template
void
DoFTools::count_dofs_per_block<deal_II_dimension> (
  const DoFHandler<deal_II_dimension>&,
  std::vector<unsigned int>&, std::vector<unsigned int>);

template
void
DoFTools::count_dofs_per_component<deal_II_dimension> (
  const DoFHandler<deal_II_dimension>&,
  std::vector<unsigned int>&, std::vector<unsigned int>);

template
void
DoFTools::compute_intergrid_constraints<deal_II_dimension> (
  const DoFHandler<deal_II_dimension> &, const unsigned int,
  const DoFHandler<deal_II_dimension> &, const unsigned int,
  const InterGridMap<DoFHandler<deal_II_dimension> > &,
  ConstraintMatrix&);

template
void
DoFTools::compute_intergrid_transfer_representation<deal_II_dimension>
(const DoFHandler<deal_II_dimension> &, const unsigned int,
 const DoFHandler<deal_II_dimension> &, const unsigned int,
 const InterGridMap<DoFHandler<deal_II_dimension> > &,
 std::vector<std::map<unsigned int, float> > &);


template
void
DoFTools::map_dof_to_boundary_indices<DoFHandler<deal_II_dimension> >
(const DoFHandler<deal_II_dimension> &,
 const std::set<unsigned char> &,
 std::vector<unsigned int> &);

#if deal_II_dimension != 1

template
void
DoFTools::map_dof_to_boundary_indices<DoFHandler<deal_II_dimension> >
(const DoFHandler<deal_II_dimension> &,
 std::vector<unsigned int> &);

#endif



template
void
DoFTools::map_dof_to_boundary_indices<hp::DoFHandler<deal_II_dimension> >
(const hp::DoFHandler<deal_II_dimension> &,
 const std::set<unsigned char> &,
 std::vector<unsigned int> &);

#if deal_II_dimension != 1

template
void
DoFTools::map_dof_to_boundary_indices<hp::DoFHandler<deal_II_dimension> >
(const hp::DoFHandler<deal_II_dimension> &,
 std::vector<unsigned int> &);

#endif




template
void
DoFTools::map_dofs_to_support_points<deal_II_dimension>
(const Mapping<deal_II_dimension,deal_II_dimension>&,
 const DoFHandler<deal_II_dimension>&,
 std::vector<Point<deal_II_dimension> >&);

template
void
DoFTools::convert_couplings_to_blocks (
  const DoFHandler<deal_II_dimension>&, const Table<2, Coupling>&,
  std::vector<Table<2,Coupling> >&);

template
void
DoFTools::convert_couplings_to_blocks (
  const hp::DoFHandler<deal_II_dimension>&, const Table<2, Coupling>&,
  std::vector<Table<2,Coupling> >&);

template 
void 
DoFTools::make_zero_boundary_constraints
(const DoFHandler<deal_II_dimension> &, 
 ConstraintMatrix                    &,
 const std::vector<bool>             &);


#if deal_II_dimension < 3
template void
DoFTools::make_boundary_sparsity_pattern<DoFHandler<deal_II_dimension,deal_II_dimension+1>,SparsityPattern>
(const DoFHandler<deal_II_dimension,deal_II_dimension+1>& dof,
 const std::vector<unsigned int>  &,
 SparsityPattern    &);
template void
DoFTools::make_boundary_sparsity_pattern<DoFHandler<deal_II_dimension,deal_II_dimension+1>,CompressedSparsityPattern>
(const DoFHandler<deal_II_dimension,deal_II_dimension+1>& dof,
 const std::vector<unsigned int>  &,
 CompressedSparsityPattern    &);
template void
DoFTools::make_boundary_sparsity_pattern<DoFHandler<deal_II_dimension,deal_II_dimension+1>,CompressedSetSparsityPattern>
(const DoFHandler<deal_II_dimension,deal_II_dimension+1>& dof,
 const std::vector<unsigned int>  &,
 CompressedSetSparsityPattern    &);
template void
DoFTools::make_boundary_sparsity_pattern<DoFHandler<deal_II_dimension,deal_II_dimension+1>,BlockSparsityPattern>
(const DoFHandler<deal_II_dimension,deal_II_dimension+1>& dof,
 const std::vector<unsigned int>  &,
 BlockSparsityPattern    &);
template void
DoFTools::make_boundary_sparsity_pattern<DoFHandler<deal_II_dimension,deal_II_dimension+1>,BlockCompressedSparsityPattern>
(const DoFHandler<deal_II_dimension,deal_II_dimension+1>& dof,
 const std::vector<unsigned int>  &,
 BlockCompressedSparsityPattern    &);
template void
DoFTools::make_boundary_sparsity_pattern<DoFHandler<deal_II_dimension,deal_II_dimension+1>,BlockCompressedSetSparsityPattern>
(const DoFHandler<deal_II_dimension,deal_II_dimension+1>& dof,
 const std::vector<unsigned int>  &,
 BlockCompressedSetSparsityPattern    &);

// template void
// DoFTools::make_boundary_sparsity_pattern<hp::DoFHandler<deal_II_dimension,deal_II_dimension+1>,SparsityPattern>
// (const hp::DoFHandler<deal_II_dimension,deal_II_dimension+1>& dof,
//  const std::vector<unsigned int>  &,
//  SparsityPattern    &);
// template void
// DoFTools::make_boundary_sparsity_pattern<hp::DoFHandler<deal_II_dimension,deal_II_dimension+1>,CompressedSparsityPattern>
// (const hp::DoFHandler<deal_II_dimension,deal_II_dimension+1>& dof,
//  const std::vector<unsigned int>  &,
//  CompressedSparsityPattern    &);
// template void
// DoFTools::make_boundary_sparsity_pattern<hp::DoFHandler<deal_II_dimension,deal_II_dimension+1>,CompressedSetSparsityPattern>
// (const hp::DoFHandler<deal_II_dimension,deal_II_dimension+1>& dof,
//  const std::vector<unsigned int>  &,
//  CompressedSetSparsityPattern    &);
// template void
// DoFTools::make_boundary_sparsity_pattern<hp::DoFHandler<deal_II_dimension,deal_II_dimension+1>,BlockSparsityPattern>
// (const hp::DoFHandler<deal_II_dimension,deal_II_dimension+1>& dof,
//  const std::vector<unsigned int>  &,
//  BlockSparsityPattern    &);
// template void
// DoFTools::make_boundary_sparsity_pattern<hp::DoFHandler<deal_II_dimension,deal_II_dimension+1>,BlockCompressedSparsityPattern>
// (const hp::DoFHandler<deal_II_dimension,deal_II_dimension+1>& dof,
//  const std::vector<unsigned int>  &,
//  BlockCompressedSparsityPattern    &);
// template void
// DoFTools::make_boundary_sparsity_pattern<hp::DoFHandler<deal_II_dimension,deal_II_dimension+1>,BlockCompressedSetSparsityPattern>
// (const hp::DoFHandler<deal_II_dimension,deal_II_dimension+1>& dof,
//  const std::vector<unsigned int>  &,
//  BlockCompressedSetSparsityPattern    &);


template void
DoFTools::make_boundary_sparsity_pattern<DoFHandler<deal_II_dimension,deal_II_dimension+1>,SparsityPattern>
(const DoFHandler<deal_II_dimension,deal_II_dimension+1>& dof,
 const FunctionMap<deal_II_dimension+1>::type  &boundary_indicators,
 const std::vector<unsigned int>  &dof_to_boundary_mapping,
 SparsityPattern    &sparsity);
template void
DoFTools::make_boundary_sparsity_pattern<DoFHandler<deal_II_dimension,deal_II_dimension+1>,CompressedSparsityPattern>
(const DoFHandler<deal_II_dimension,deal_II_dimension+1>& dof,
 const FunctionMap<deal_II_dimension+1>::type  &boundary_indicators,
 const std::vector<unsigned int>  &dof_to_boundary_mapping,
 CompressedSparsityPattern    &sparsity);
template void
DoFTools::make_boundary_sparsity_pattern<DoFHandler<deal_II_dimension,deal_II_dimension+1>,CompressedSetSparsityPattern>
(const DoFHandler<deal_II_dimension,deal_II_dimension+1>& dof,
 const FunctionMap<deal_II_dimension+1>::type  &boundary_indicators,
 const std::vector<unsigned int>  &dof_to_boundary_mapping,
 CompressedSetSparsityPattern    &sparsity);
template void
DoFTools::make_boundary_sparsity_pattern<DoFHandler<deal_II_dimension,deal_II_dimension+1>,BlockSparsityPattern>
(const DoFHandler<deal_II_dimension,deal_II_dimension+1>& dof,
 const FunctionMap<deal_II_dimension+1>::type  &boundary_indicators,
 const std::vector<unsigned int>  &dof_to_boundary_mapping,
 BlockSparsityPattern    &sparsity);
template void
DoFTools::make_boundary_sparsity_pattern<DoFHandler<deal_II_dimension,deal_II_dimension+1>,BlockCompressedSparsityPattern>
(const DoFHandler<deal_II_dimension,deal_II_dimension+1>& dof,
 const FunctionMap<deal_II_dimension+1>::type  &boundary_indicators,
 const std::vector<unsigned int>  &dof_to_boundary_mapping,
 BlockCompressedSparsityPattern    &sparsity);
template void
DoFTools::make_boundary_sparsity_pattern<DoFHandler<deal_II_dimension,deal_II_dimension+1>,BlockCompressedSetSparsityPattern>
(const DoFHandler<deal_II_dimension,deal_II_dimension+1>& dof,
 const FunctionMap<deal_II_dimension+1>::type  &boundary_indicators,
 const std::vector<unsigned int>  &dof_to_boundary_mapping,
 BlockCompressedSetSparsityPattern    &sparsity);

// template void
// DoFTools::make_boundary_sparsity_pattern<hp::DoFHandler<deal_II_dimension,deal_II_dimension+1>,SparsityPattern>
// (const hp::DoFHandler<deal_II_dimension,deal_II_dimension+1>& dof,
//  const FunctionMap<deal_II_dimension+1>::type  &boundary_indicators,
//  const std::vector<unsigned int>  &dof_to_boundary_mapping,
//  SparsityPattern    &sparsity);
// template void
// DoFTools::make_boundary_sparsity_pattern<hp::DoFHandler<deal_II_dimension,deal_II_dimension+1>,CompressedSparsityPattern>
// (const hp::DoFHandler<deal_II_dimension,deal_II_dimension+1>& dof,
//  const FunctionMap<deal_II_dimension+1>::type  &boundary_indicators,
//  const std::vector<unsigned int>  &dof_to_boundary_mapping,
//  CompressedSparsityPattern    &sparsity);
// template void
// DoFTools::make_boundary_sparsity_pattern<hp::DoFHandler<deal_II_dimension,deal_II_dimension+1>,CompressedSetSparsityPattern>
// (const hp::DoFHandler<deal_II_dimension,deal_II_dimension+1>& dof,
//  const FunctionMap<deal_II_dimension+1>::type  &boundary_indicators,
//  const std::vector<unsigned int>  &dof_to_boundary_mapping,
//  CompressedSetSparsityPattern    &sparsity);
// template void
// DoFTools::make_boundary_sparsity_pattern<hp::DoFHandler<deal_II_dimension,deal_II_dimension+1>,BlockSparsityPattern>
// (const hp::DoFHandler<deal_II_dimension,deal_II_dimension+1>& dof,
//  const FunctionMap<deal_II_dimension+1>::type  &boundary_indicators,
//  const std::vector<unsigned int>  &dof_to_boundary_mapping,
//  BlockSparsityPattern    &sparsity);
// template void
// DoFTools::make_boundary_sparsity_pattern<hp::DoFHandler<deal_II_dimension,deal_II_dimension+1>,BlockCompressedSparsityPattern>
// (const hp::DoFHandler<deal_II_dimension,deal_II_dimension+1>& dof,
//  const FunctionMap<deal_II_dimension+1>::type  &boundary_indicators,
//  const std::vector<unsigned int>  &dof_to_boundary_mapping,
//  BlockCompressedSparsityPattern    &sparsity);
// template void
// DoFTools::make_boundary_sparsity_pattern<hp::DoFHandler<deal_II_dimension,deal_II_dimension+1>,BlockCompressedSetSparsityPattern>
// (const hp::DoFHandler<deal_II_dimension,deal_II_dimension+1>& dof,
//  const FunctionMap<deal_II_dimension+1>::type  &boundary_indicators,
//  const std::vector<unsigned int>  &dof_to_boundary_mapping,
//  BlockCompressedSetSparsityPattern    &sparsity);


template
void
DoFTools::map_dof_to_boundary_indices<DoFHandler<deal_II_dimension,deal_II_dimension+1> >
(const DoFHandler<deal_II_dimension,deal_II_dimension+1> &,
 const std::set<unsigned char> &,
 std::vector<unsigned int> &);

#if deal_II_dimension != 1

template
void
DoFTools::map_dof_to_boundary_indices<DoFHandler<deal_II_dimension,deal_II_dimension+1> >
(const DoFHandler<deal_II_dimension,deal_II_dimension+1> &,
 std::vector<unsigned int> &);

#endif


template
void
DoFTools::extract_hanging_node_dofs
(const DoFHandler<deal_II_dimension,deal_II_dimension+1> &dof_handler,
 std::vector<bool>     &selected_dofs);

// template
// void
// DoFTools::map_dof_to_boundary_indices<hp::DoFHandler<deal_II_dimension,deal_II_dimension+1> >
// (const hp::DoFHandler<deal_II_dimension,deal_II_dimension+1> &,
//  const std::set<unsigned char> &,
//  std::vector<unsigned int> &);

// #if deal_II_dimension != 1

// template
// void
// DoFTools::map_dof_to_boundary_indices<hp::DoFHandler<deal_II_dimension,deal_II_dimension+1> >
// (const hp::DoFHandler<deal_II_dimension,deal_II_dimension+1> &,
//  std::vector<unsigned int> &);

// #endif



#endif




DEAL_II_NAMESPACE_CLOSE
