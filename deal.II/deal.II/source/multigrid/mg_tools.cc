//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <base/multithread_info.h>
#include <base/logstream.h>
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
#include <multigrid/mg_tools.h>
#include <multigrid/mg_base.h>
#include <multigrid/mg_level_object.h>
#include <dofs/dof_tools.h>
#include <fe/fe.h>

#include <vector>
#include <algorithm>
#include <numeric>

DEAL_II_NAMESPACE_OPEN


#if deal_II_dimension == 1

template <int dim, int spacedim>
void
MGTools::compute_row_length_vector(
  const MGDoFHandler<dim,spacedim>&,
  const unsigned int,
  std::vector<unsigned int>&,
  const DoFTools::Coupling)
{
  Assert(false, ExcNotImplemented());
}


template <int dim, int spacedim>
void
MGTools::compute_row_length_vector(
  const MGDoFHandler<dim,spacedim>&,
  const unsigned int,
  std::vector<unsigned int>&,
  const Table<2,DoFTools::Coupling>&,
  const Table<2,DoFTools::Coupling>&)
{
  Assert(false, ExcNotImplemented());
}

#else

// Template for 2D and 3D. For 1D see specialization above
template <int dim, int spacedim>
void
MGTools::compute_row_length_vector(
  const MGDoFHandler<dim,spacedim>& dofs,
  const unsigned int level,
  std::vector<unsigned int>& row_lengths,
  const DoFTools::Coupling             flux_coupling)
{
  Assert (row_lengths.size() == dofs.n_dofs(),
	  ExcDimensionMismatch(row_lengths.size(), dofs.n_dofs()));

				   // Function starts here by
				   // resetting the counters.
  std::fill(row_lengths.begin(), row_lengths.end(), 0);
				   // We need the user flags, so we
				   // save them for later restoration
  std::vector<bool> old_flags;
				   // We need a non-constant
				   // triangulation for the user
				   // flags. Since we restore them in
				   // the end, this cast is safe.
  Triangulation<dim,spacedim>& user_flags_triangulation =
    const_cast<Triangulation<dim,spacedim>&> (dofs.get_tria());
  user_flags_triangulation.save_user_flags(old_flags);
  user_flags_triangulation.clear_user_flags();

  const typename MGDoFHandler<dim,spacedim>::cell_iterator end = dofs.end(level);
  typename MGDoFHandler<dim,spacedim>::active_cell_iterator cell;
  std::vector<unsigned int> cell_indices;
  std::vector<unsigned int> neighbor_indices;

				   // We loop over cells and go from
				   // cells to lower dimensional
				   // objects. This is the only way to
				   // cope with the fact, that an
				   // unknown number of cells may
				   // share an object of dimension
				   // smaller than dim-1.
  for (cell = dofs.begin(level); cell != end; ++cell)
    {
      const FiniteElement<dim>& fe = cell->get_fe();
      cell_indices.resize(fe.dofs_per_cell);
      cell->get_mg_dof_indices(cell_indices);
      unsigned int i = 0;
				       // First, dofs on
				       // vertices. We assume that
				       // each vertex dof couples
				       // with all dofs on
				       // adjacent grid cells.

				       // Adding all dofs of the cells
				       // will add dofs of the faces
				       // of the cell adjacent to the
				       // vertex twice. Therefore, we
				       // subtract these here and add
				       // them in a loop over the
				       // faces below.

				       // in 1d, faces and vertices
				       // are identical. Nevertheless,
				       // this will only work if
				       // dofs_per_face is zero and
				       // dofs_per_vertex is
				       // arbitrary, not the other way
				       // round.
//TODO: This assumes that even in hp context, the dofs per face coincide!
      unsigned int increment = fe.dofs_per_cell - dim * fe.dofs_per_face;
      while (i < fe.first_line_index)
	row_lengths[cell_indices[i++]] += increment;
				       // From now on, if an object is
				       // a cell, its dofs only couple
				       // inside the cell. Since the
				       // faces are handled below, we
				       // have to subtract ALL faces
				       // in this case.

				       // In all other cases we
				       // subtract adjacent faces to be
				       // added in the loop below.
      increment = (dim>1)
		  ? fe.dofs_per_cell - (dim-1) * fe.dofs_per_face
		  : fe.dofs_per_cell - GeometryInfo<dim>::faces_per_cell * fe.dofs_per_face;
      while (i < fe.first_quad_index)
	row_lengths[cell_indices[i++]] += increment;

				       // Now quads in 2D and 3D
      increment = (dim>2)
		  ? fe.dofs_per_cell - (dim-2) * fe.dofs_per_face
		  : fe.dofs_per_cell - GeometryInfo<dim>::faces_per_cell * fe.dofs_per_face;
      while (i < fe.first_hex_index)
	row_lengths[cell_indices[i++]] += increment;
				       // Finally, cells in 3D
      increment = fe.dofs_per_cell - GeometryInfo<dim>::faces_per_cell * fe.dofs_per_face;
      while (i < fe.dofs_per_cell)
	row_lengths[cell_indices[i++]] += increment;

				   // At this point, we have
				   // counted all dofs
				   // contributiong from cells
				   // coupled topologically to the
				   // adjacent cells, but we
				   // subtracted some faces.

				   // Now, let's go by the faces
				   // and add the missing
				   // contribution as well as the
				   // flux contributions.
      for (unsigned int iface=0;iface<GeometryInfo<dim>::faces_per_cell;++iface)
	{
	  bool level_boundary = cell->at_boundary(iface);
	  typename MGDoFHandler<dim,spacedim>::cell_iterator neighbor;
	  if (!level_boundary)
	    {
	      neighbor = cell->neighbor(iface);
	      if (static_cast<unsigned int>(neighbor->level()) != level)
		level_boundary = true;
	    }

	  if (level_boundary)
	    {
	      for (unsigned int i=0;i<fe.dofs_per_cell;++i)
		row_lengths[cell_indices[i]] += fe.dofs_per_face;
	      continue;
	    }

	  const FiniteElement<dim>& nfe = neighbor->get_fe();
	  typename MGDoFHandler<dim,spacedim>::face_iterator face = cell->face(iface);

					   // Flux couplings are
					   // computed from both sides
					   // for simplicity.

					   // The dofs on the common face
					   // will be handled below,
					   // therefore, we subtract them
					   // here.
	  if (flux_coupling != DoFTools::none)
	    {
	      unsigned int increment = nfe.dofs_per_cell - nfe.dofs_per_face;
	      for (unsigned int i=0;i<fe.dofs_per_cell;++i)
		row_lengths[cell_indices[i]] += increment;
	    }

					   // Do this only once per
					   // face.
	  if (face->user_flag_set())
	    continue;
	  face->set_user_flag();
					   // At this point, we assume
					   // that each cell added its
					   // dofs minus the face to
					   // the couplings of the
					   // face dofs. Since we
					   // subtracted two faces, we
					   // have to re-add one.

					   // If one side of the face
					   // is refined, all the fine
					   // face dofs couple with
					   // the coarse one.
	  neighbor_indices.resize(nfe.dofs_per_cell);
	  neighbor->get_mg_dof_indices(neighbor_indices);
	  for (unsigned int i=0;i<fe.dofs_per_cell;++i)
	    row_lengths[cell_indices[i]] += nfe.dofs_per_face;
	  for (unsigned int i=0;i<nfe.dofs_per_cell;++i)
	    row_lengths[neighbor_indices[i]] += fe.dofs_per_face;
	}
    }
  user_flags_triangulation.load_user_flags(old_flags);
}


// This is the template for 2D and 3D. See version for 1D above
template <int dim, int spacedim>
void
MGTools::compute_row_length_vector(
  const MGDoFHandler<dim,spacedim>& dofs,
  const unsigned int level,
  std::vector<unsigned int>& row_lengths,
  const Table<2,DoFTools::Coupling>& couplings,
  const Table<2,DoFTools::Coupling>& flux_couplings)
{
  Assert (row_lengths.size() == dofs.n_dofs(),
	  ExcDimensionMismatch(row_lengths.size(), dofs.n_dofs()));

				   // Function starts here by
				   // resetting the counters.
  std::fill(row_lengths.begin(), row_lengths.end(), 0);
				   // We need the user flags, so we
				   // save them for later restoration
  std::vector<bool> old_flags;
				   // We need a non-constant
				   // triangulation for the user
				   // flags. Since we restore them in
				   // the end, this cast is safe.
  Triangulation<dim,spacedim>& user_flags_triangulation =
    const_cast<Triangulation<dim,spacedim>&> (dofs.get_tria());
  user_flags_triangulation.save_user_flags(old_flags);
  user_flags_triangulation.clear_user_flags();

  const typename MGDoFHandler<dim,spacedim>::cell_iterator end = dofs.end(level);
  typename MGDoFHandler<dim,spacedim>::active_cell_iterator cell;
  std::vector<unsigned int> cell_indices;
  std::vector<unsigned int> neighbor_indices;

				   // We have to translate the
				   // couplings from components to
				   // blocks, so this works for
				   // nonprimitive elements as well.
  std::vector<Table<2, DoFTools::Coupling> > couple_cell;
  std::vector<Table<2, DoFTools::Coupling> > couple_face;
  DoFTools::convert_couplings_to_blocks(dofs, couplings, couple_cell);
  DoFTools::convert_couplings_to_blocks(dofs, flux_couplings, couple_face);

				   // We loop over cells and go from
				   // cells to lower dimensional
				   // objects. This is the only way to
				   // cope withthe fact, that an
				   // unknown number of cells may
				   // share an object of dimension
				   // smaller than dim-1.
  for (cell = dofs.begin_active(); cell != end; ++cell)
    {
      const FiniteElement<dim>& fe = cell->get_fe();
      const unsigned int fe_index = cell->active_fe_index();

      Assert (couplings.n_rows()==fe.n_components(),
	      ExcDimensionMismatch(couplings.n_rows(), fe.n_components()));
      Assert (couplings.n_cols()==fe.n_components(),
	      ExcDimensionMismatch(couplings.n_cols(), fe.n_components()));
      Assert (flux_couplings.n_rows()==fe.n_components(),
	      ExcDimensionMismatch(flux_couplings.n_rows(), fe.n_components()));
      Assert (flux_couplings.n_cols()==fe.n_components(),
	      ExcDimensionMismatch(flux_couplings.n_cols(), fe.n_components()));

      cell_indices.resize(fe.dofs_per_cell);
      cell->get_mg_dof_indices(cell_indices);
      unsigned int i = 0;
				       // First, dofs on
				       // vertices. We assume that
				       // each vertex dof couples
				       // with all dofs on
				       // adjacent grid cells.

				       // Adding all dofs of the cells
				       // will add dofs of the faces
				       // of the cell adjacent to the
				       // vertex twice. Therefore, we
				       // subtract these here and add
				       // them in a loop over the
				       // faces below.

				       // in 1d, faces and vertices
				       // are identical. Nevertheless,
				       // this will only work if
				       // dofs_per_face is zero and
				       // dofs_per_vertex is
				       // arbitrary, not the other way
				       // round.
      unsigned int increment;
      while (i < fe.first_line_index)
	{
	  for (unsigned int base=0;base<fe.n_base_elements();++base)
	    for (unsigned int mult=0;mult<fe.element_multiplicity(base);++mult)
	      if (couple_cell[fe_index](fe.system_to_block_index(i).first,
					fe.first_block_of_base(base) + mult) != DoFTools::none)
		{
		  increment = fe.base_element(base).dofs_per_cell
			      - dim * fe.base_element(base).dofs_per_face;
		  row_lengths[cell_indices[i]] += increment;
		}
	  ++i;
	}
				       // From now on, if an object is
				       // a cell, its dofs only couple
				       // inside the cell. Since the
				       // faces are handled below, we
				       // have to subtract ALL faces
				       // in this case.

				       // In all other cases we
				       // subtract adjacent faces to be
				       // added in the loop below.
      while (i < fe.first_quad_index)
	{
	  for (unsigned int base=0;base<fe.n_base_elements();++base)
	    for (unsigned int mult=0;mult<fe.element_multiplicity(base);++mult)
	      if (couple_cell[fe_index](fe.system_to_block_index(i).first,
					fe.first_block_of_base(base) + mult) != DoFTools::none)
		{
		  increment = fe.base_element(base).dofs_per_cell
			      - ((dim>1)
				 ? (dim-1)
				 : GeometryInfo<dim>::faces_per_cell)
			      * fe.base_element(base).dofs_per_face;
		row_lengths[cell_indices[i]] += increment;
		}
	  ++i;
	}

				       // Now quads in 2D and 3D
      while (i < fe.first_hex_index)
	{
	  for (unsigned int base=0;base<fe.n_base_elements();++base)
	    for (unsigned int mult=0;mult<fe.element_multiplicity(base);++mult)
	      if (couple_cell[fe_index](fe.system_to_block_index(i).first,
					fe.first_block_of_base(base) + mult) != DoFTools::none)
		{
		  increment = fe.base_element(base).dofs_per_cell
			      - ((dim>2)
				 ? (dim-2)
				 : GeometryInfo<dim>::faces_per_cell)
			      * fe.base_element(base).dofs_per_face;
		row_lengths[cell_indices[i]] += increment;
	      }
	  ++i;
	}

				       // Finally, cells in 3D
      while (i < fe.dofs_per_cell)
	{
	  for (unsigned int base=0;base<fe.n_base_elements();++base)
	    for (unsigned int mult=0;mult<fe.element_multiplicity(base);++mult)
	      if (couple_cell[fe_index](fe.system_to_block_index(i).first,
					fe.first_block_of_base(base) + mult) != DoFTools::none)
		{
		  increment = fe.base_element(base).dofs_per_cell
			      - GeometryInfo<dim>::faces_per_cell
			      * fe.base_element(base).dofs_per_face;
		  row_lengths[cell_indices[i]] += increment;
		}
	  ++i;
	}

				   // At this point, we have
				   // counted all dofs
				   // contributiong from cells
				   // coupled topologically to the
				   // adjacent cells, but we
				   // subtracted some faces.

				   // Now, let's go by the faces
				   // and add the missing
				   // contribution as well as the
				   // flux contributions.
      for (unsigned int iface=0;iface<GeometryInfo<dim>::faces_per_cell;++iface)
	{
	  bool level_boundary = cell->at_boundary(iface);
	  typename MGDoFHandler<dim,spacedim>::cell_iterator neighbor;
	  if (!level_boundary)
	    {
	      neighbor = cell->neighbor(iface);
	      if (static_cast<unsigned int>(neighbor->level()) != level)
		level_boundary = true;
	    }

	  if (level_boundary)
	    {
	      for (unsigned int i=0;i<fe.dofs_per_cell;++i)
		row_lengths[cell_indices[i]] += fe.dofs_per_face;
	      continue;
	    }

	  const FiniteElement<dim>& nfe = neighbor->get_fe();
	  typename MGDoFHandler<dim,spacedim>::face_iterator face = cell->face(iface);

					   // Flux couplings are
					   // computed from both sides
					   // for simplicity.

					   // The dofs on the common face
					   // will be handled below,
					   // therefore, we subtract them
					   // here.
	  for (unsigned int base=0;base<nfe.n_base_elements();++base)
	    for (unsigned int mult=0;mult<nfe.element_multiplicity(base);++mult)
	      for (unsigned int i=0;i<fe.dofs_per_cell;++i)
		if (couple_face[fe_index](fe.system_to_block_index(i).first,
					  nfe.first_block_of_base(base) + mult) != DoFTools::none)
		{
		  unsigned int increment = nfe.base_element(base).dofs_per_cell
					   - nfe.base_element(base).dofs_per_face;
		  row_lengths[cell_indices[i]] += increment;
		}

					   // Do this only once per
					   // face and not on the
					   // hanging faces.
	  if (face->user_flag_set())
	    continue;
	  face->set_user_flag();
					   // At this point, we assume
					   // that each cell added its
					   // dofs minus the face to
					   // the couplings of the
					   // face dofs. Since we
					   // subtracted two faces, we
					   // have to re-add one.

					   // If one side of the face
					   // is refined, all the fine
					   // face dofs couple with
					   // the coarse one.

					   // Wolfgang, do they couple
					   // with each other by
					   // constraints?

					   // This will not work with
					   // different couplings on
					   // different cells.
	  neighbor_indices.resize(nfe.dofs_per_cell);
	  neighbor->get_mg_dof_indices(neighbor_indices);
	  for (unsigned int base=0;base<nfe.n_base_elements();++base)
	    for (unsigned int mult=0;mult<nfe.element_multiplicity(base);++mult)
	      for (unsigned int i=0;i<fe.dofs_per_cell;++i)
		if (couple_cell[fe_index](fe.system_to_component_index(i).first,
					  nfe.first_block_of_base(base) + mult) != DoFTools::none)
		  row_lengths[cell_indices[i]]
		    += nfe.base_element(base).dofs_per_face;
	  for (unsigned int base=0;base<fe.n_base_elements();++base)
	    for (unsigned int mult=0;mult<fe.element_multiplicity(base);++mult)
	      for (unsigned int i=0;i<nfe.dofs_per_cell;++i)
		if (couple_cell[fe_index](nfe.system_to_component_index(i).first,
					  fe.first_block_of_base(base) + mult) != DoFTools::none)
		  row_lengths[neighbor_indices[i]]
		    += fe.base_element(base).dofs_per_face;
	}
    }
  user_flags_triangulation.load_user_flags(old_flags);
}
#endif


template <int dim, class SparsityPattern, int spacedim>
void MGTools::make_sparsity_pattern (
  const MGDoFHandler<dim,spacedim> &dof,
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
  typename MGDoFHandler<dim,spacedim>::cell_iterator cell = dof.begin(level),
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



template <int dim, class SparsityPattern, int spacedim>
void
MGTools::make_flux_sparsity_pattern (
  const MGDoFHandler<dim,spacedim> &dof,
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
  typename MGDoFHandler<dim,spacedim>::cell_iterator cell = dof.begin(level),
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
	      typename MGDoFHandler<dim,spacedim>::cell_iterator
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



template <int dim, class SparsityPattern, int spacedim>
void
MGTools::make_flux_sparsity_pattern_edge (
  const MGDoFHandler<dim,spacedim> &dof,
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
  typename MGDoFHandler<dim,spacedim>::cell_iterator cell = dof.begin(level),
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
	      typename MGDoFHandler<dim,spacedim>::cell_iterator
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






template <int dim, class SparsityPattern, int spacedim>
void
MGTools::make_flux_sparsity_pattern (
  const MGDoFHandler<dim,spacedim> &dof,
  SparsityPattern       &sparsity,
  const unsigned int level,
  const Table<2,DoFTools::Coupling> &int_mask,
  const Table<2,DoFTools::Coupling> &flux_mask)
{
  const FiniteElement<dim>& fe = dof.get_fe();
  const unsigned int n_dofs = dof.n_dofs(level);
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
  Table<2,bool> support_on_face(total_dofs, GeometryInfo<dim>::faces_per_cell);

  typename MGDoFHandler<dim,spacedim>::cell_iterator cell = dof.begin(level),
					    endc = dof.end(level);

  const Table<2,DoFTools::Coupling>
    int_dof_mask  = DoFTools::dof_couplings_from_component_couplings(fe, int_mask),
    flux_dof_mask = DoFTools::dof_couplings_from_component_couplings(fe, flux_mask);

  for (unsigned int i=0; i<total_dofs; ++i)
    for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell;++f)
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
  const_cast<Triangulation<dim,spacedim> &>(dof.get_tria()).clear_user_flags ();

  for (; cell!=endc; ++cell)
    {
      cell->get_mg_dof_indices (dofs_on_this_cell);
				       // make sparsity pattern for this cell
      for (unsigned int i=0; i<total_dofs; ++i)
	for (unsigned int j=0; j<total_dofs; ++j)
	  if (int_dof_mask[i][j] != DoFTools::none)
	    sparsity.add (dofs_on_this_cell[i],
			  dofs_on_this_cell[j]);

				       // Loop over all interior neighbors
      for (unsigned int face = 0;
	   face < GeometryInfo<dim>::faces_per_cell;
	   ++face)
	{
	  typename MGDoFHandler<dim,spacedim>::face_iterator cell_face = cell->face(face);
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

		      if (flux_dof_mask(i,j) == DoFTools::always)
                        sparsity.add (dofs_on_this_cell[i],
                                      dofs_on_this_cell[j]);
		      if (flux_dof_mask(i,j) == DoFTools::nonzero
			  && i_non_zero_i && j_non_zero_i)
			sparsity.add (dofs_on_this_cell[i],
				      dofs_on_this_cell[j]);
		    }
		}
	    }
	  else
	    {
	      typename MGDoFHandler<dim,spacedim>::cell_iterator
		neighbor = cell->neighbor(face);

	      if (neighbor->level() < cell->level())
		continue;

	      unsigned int neighbor_face = cell->neighbor_of_neighbor(face);

	      neighbor->get_mg_dof_indices (dofs_on_other_cell);
	      for (unsigned int i=0; i<total_dofs; ++i)
		{
		  const bool i_non_zero_i = support_on_face (i, face);
		  const bool i_non_zero_e = support_on_face (i, neighbor_face);
		  for (unsigned int j=0; j<total_dofs; ++j)
		    {
		      const bool j_non_zero_i = support_on_face (j, face);
		      const bool j_non_zero_e = support_on_face (j, neighbor_face);
		      if (flux_dof_mask(i,j) == DoFTools::always)
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
		      if (flux_dof_mask(i,j) == DoFTools::nonzero)
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

		      if (flux_dof_mask(j,i) == DoFTools::always)
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
		      if (flux_dof_mask(j,i) == DoFTools::nonzero)
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

				   // finally restore the user flags
  const_cast<Triangulation<dim,spacedim> &>(dof.get_tria()).load_user_flags(user_flags);
}



template <int dim, class SparsityPattern, int spacedim>
void
MGTools::make_flux_sparsity_pattern_edge (
  const MGDoFHandler<dim,spacedim> &dof,
  SparsityPattern       &sparsity,
  const unsigned int level,
  const Table<2,DoFTools::Coupling> &flux_mask)
{
  const FiniteElement<dim>& fe = dof.get_fe();
  const unsigned int n_comp = fe.n_components();

  Assert ((level>=1) && (level<dof.get_tria().n_levels()),
	  ExcIndexRange(level, 1, dof.get_tria().n_levels()));

  const unsigned int fine_dofs = dof.n_dofs(level);
  const unsigned int coarse_dofs = dof.n_dofs(level-1);

  // Matrix maps from fine level to coarse level

  Assert (sparsity.n_rows() == coarse_dofs,
	  ExcDimensionMismatch (sparsity.n_rows(), coarse_dofs));
  Assert (sparsity.n_cols() == fine_dofs,
	  ExcDimensionMismatch (sparsity.n_cols(), fine_dofs));
  Assert (flux_mask.n_rows() == n_comp,
	  ExcDimensionMismatch (flux_mask.n_rows(), n_comp));
  Assert (flux_mask.n_cols() == n_comp,
	  ExcDimensionMismatch (flux_mask.n_cols(), n_comp));

  const unsigned int dofs_per_cell = dof.get_fe().dofs_per_cell;
  std::vector<unsigned int> dofs_on_this_cell(dofs_per_cell);
  std::vector<unsigned int> dofs_on_other_cell(dofs_per_cell);
  Table<2,bool> support_on_face(dofs_per_cell, GeometryInfo<dim>::faces_per_cell);

  typename MGDoFHandler<dim,spacedim>::cell_iterator cell = dof.begin(level),
						     endc = dof.end(level);

  const Table<2,DoFTools::Coupling> flux_dof_mask
    = DoFTools::dof_couplings_from_component_couplings(fe, flux_mask);

  for (unsigned int i=0; i<dofs_per_cell; ++i)
    for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell;++f)
      support_on_face(i,f) = fe.has_support_on_face(i,f);

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
	      typename MGDoFHandler<dim,spacedim>::cell_iterator
		neighbor = cell->neighbor(face);
	      neighbor->get_mg_dof_indices (dofs_on_other_cell);

	      for (unsigned int i=0; i<dofs_per_cell; ++i)
		{
		  for (unsigned int j=0; j<dofs_per_cell; ++j)
		    {
		      if (flux_dof_mask(i,j) != DoFTools::none)
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
}


template <int dim, int spacedim>
void
MGTools::
count_dofs_per_component (const MGDoFHandler<dim,spacedim> &dof_handler,
			  std::vector<std::vector<unsigned int> > &result,
			  bool                              only_once,
			  std::vector<unsigned int>         target_component)
{
  const FiniteElement<dim>& fe = dof_handler.get_fe();
  const unsigned int n_components = fe.n_components();
  const unsigned int nlevels = dof_handler.get_tria().n_levels();

  Assert (result.size() == nlevels,
	  ExcDimensionMismatch(result.size(), nlevels));

  if (target_component.size() == 0)
    {
      target_component.resize(n_components);
      for (unsigned int i=0;i<n_components;++i)
	target_component[i] = i;
    }

  Assert(target_component.size() == n_components,
	 ExcDimensionMismatch(target_component.size(), n_components));

  for (unsigned int l=0;l<nlevels;++l)
    {
      result[l].resize (n_components);
      std::fill (result[l].begin(),result[l].end(), 0U);

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
	  Threads::TaskGroup<> tasks;
	  for (unsigned int i=0; i<n_components; ++i)
	    {
	      void (*fun_ptr) (const unsigned int       level,
			       const MGDoFHandler<dim,spacedim>    &,
			       const std::vector<bool>    &,
			       std::vector<bool>          &,
			       bool)
		= &DoFTools::template extract_level_dofs<dim>;
	      component_select[i][i] = true;
	      tasks += Threads::new_task (fun_ptr,
					  l, dof_handler,
					  component_select[i],
					  dofs_in_component[i], false);
	    }
	  tasks.join_all();

					   // next count what we got
	  unsigned int component = 0;
	  for (unsigned int b=0;b<fe.n_base_elements();++b)
	    {
	      const FiniteElement<dim>& base = fe.base_element(b);
					       // Dimension of base element
	      unsigned int d = base.n_components();

	      for (unsigned int m=0;m<fe.element_multiplicity(b);++m)
		{
		  for (unsigned int dd=0;dd<d;++dd)
		    {
		      if (base.is_primitive() || (!only_once || dd==0))
			result[l][target_component[component]]
			  += std::count(dofs_in_component[component].begin(),
					dofs_in_component[component].end(),
					true);
		      ++component;
		    }
		}
	    }
					   // finally sanity check
	  Assert (!dof_handler.get_fe().is_primitive()
		  ||
		  std::accumulate (result[l].begin(),
				   result[l].end(), 0U)
		  ==
		  dof_handler.n_dofs(l),
		  ExcInternalError());
	}
    }
}



template <int dim, int spacedim>
void
MGTools::
count_dofs_per_component (const MGDoFHandler<dim,spacedim>        &dof_handler,
			  std::vector<std::vector<unsigned int> > &result,
			  std::vector<unsigned int>            target_component)
{
  count_dofs_per_component (dof_handler, result,
			    false, target_component);
}


template <int dim, int spacedim>
void
MGTools::count_dofs_per_block (
  const MGDoFHandler<dim,spacedim>&     dof_handler,
  std::vector<std::vector<unsigned int> >& dofs_per_block,
  std::vector<unsigned int>  target_block)
{
  const FiniteElement<dim>& fe = dof_handler.get_fe();
  const unsigned int n_blocks = fe.n_blocks();
  const unsigned int n_levels = dof_handler.get_tria().n_levels();

  AssertDimension (dofs_per_block.size(), n_levels);

  for (unsigned int l=0;l<n_levels;++l)
    std::fill (dofs_per_block[l].begin(), dofs_per_block[l].end(), 0U);
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

  const unsigned int max_block
    = *std::max_element (target_block.begin(),
			 target_block.end());
  const unsigned int n_target_blocks = max_block + 1;

  for (unsigned int l=0;l<n_levels;++l)
    AssertDimension (dofs_per_block[l].size(), n_target_blocks);

				   // special case for only one
				   // block. treat this first
				   // since it does not require any
				   // computations
  if (n_blocks == 1)
    {
      for (unsigned int l=0;l<n_levels;++l)
	dofs_per_block[l][0] = dof_handler.n_dofs(l);
      return;
    }
				   // otherwise determine the number
				   // of dofs in each block
				   // separately. do so in parallel
  for (unsigned int l=0;l<n_levels;++l)
    {
      std::vector<std::vector<bool> >
	dofs_in_block (n_blocks, std::vector<bool>(dof_handler.n_dofs(l), false));
      std::vector<std::vector<bool> >
	block_select (n_blocks, std::vector<bool>(n_blocks, false));
      Threads::TaskGroup<> tasks;
      for (unsigned int i=0; i<n_blocks; ++i)
	{
	  void (*fun_ptr) (const unsigned int level,
			   const MGDoFHandler<dim,spacedim>&,
			   const std::vector<bool>&,
			   std::vector<bool>&,
			   bool)
	    = &DoFTools::template extract_level_dofs<dim>;
	  block_select[i][i] = true;
	  tasks += Threads::new_task (fun_ptr,
				      l, dof_handler, block_select[i],
				      dofs_in_block[i], true);
	};
      tasks.join_all ();

				       // next count what we got
      for (unsigned int block=0;block<fe.n_blocks();++block)
	dofs_per_block[l][target_block[block]]
	  += std::count(dofs_in_block[block].begin(),
			dofs_in_block[block].end(),
			true);
    }
}


#if deal_II_dimension == 1

template <int dim, int spacedim>
void
MGTools::make_boundary_list(
  const MGDoFHandler<dim,spacedim>&,
  const typename FunctionMap<dim>::type&,
  std::vector<std::set<unsigned int> >&,
  const std::vector<bool>&)
{
  Assert(false, ExcNotImplemented());
}


#else

template <int dim, int spacedim>
void
MGTools::make_boundary_list(
  const MGDoFHandler<dim,spacedim>& dof,
  const typename FunctionMap<dim>::type& function_map,
  std::vector<std::set<unsigned int> >& boundary_indices,
  const std::vector<bool>& component_mask_)
{
                                   // if for whatever reason we were
                                   // passed an empty map, return
                                   // immediately
  if (function_map.size() == 0)
    return;

  const unsigned int n_levels = dof.get_tria().n_levels();

  
  
  const unsigned int n_components = DoFTools::n_components(dof);
  const bool          fe_is_system = (n_components != 1);

  AssertDimension (boundary_indices.size(), n_levels);

  std::vector<unsigned int> local_dofs;
  
				   // First, deal with the simpler
				   // case when we have to identify
				   // all boundary dofs
  if (component_mask_.size() == 0)
    {
      typename MGDoFHandler<dim,spacedim>::cell_iterator
	cell = dof.begin(),
	endc = dof.end();
      for (; cell!=endc; ++cell)
	{
	  const FiniteElement<dim> &fe = cell->get_fe();
	  const unsigned int level = cell->level();
	  local_dofs.resize(fe.dofs_per_face);
	  
	  for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell;
	       ++face_no)
	    {
	      const typename MGDoFHandler<dim,spacedim>::face_iterator
		face = cell->face(face_no);
	      face->get_mg_dof_indices(level, local_dofs);
	      const unsigned char bi = face->boundary_indicator();
					       // Face is listed in
					       // boundary map
	      if (function_map.find(bi) != function_map.end())
		{
		  for (unsigned int i=0;i<fe.dofs_per_face;++i)
		    boundary_indices[level].insert(i);
		}
	    }
	}
    }
  else
    {
      
//    for (typename FunctionMap<dim>::type::const_iterator i=function_map.begin();
//         i!=function_map.end(); ++i)
//      Assert (n_components == i->second->n_components,
//           ExcInvalidFE());
      
				       // set the component mask to either
				       // the original value or a vector
				       // of @p{true}s
      const std::vector<bool> component_mask ((component_mask_.size() == 0) ?
					      std::vector<bool> (n_components, true) :
					      component_mask_);
//   Assert (std::count(component_mask.begin(), component_mask.end(), true) > 0,
//        ExcComponentMismatch());
      
				       // field to store the indices
      std::vector<unsigned int> local_dofs;
      local_dofs.reserve (DoFTools::max_dofs_per_face(dof));
      std::fill (local_dofs.begin (), local_dofs.end (),
		 DoFHandler<dim,spacedim>::invalid_dof_index);
      
      typename MGDoFHandler<dim,spacedim>::cell_iterator
	cell = dof.begin(),
	endc = dof.end();
      for (; cell!=endc; ++cell)
	for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell;
	     ++face_no)
	  {
	    const FiniteElement<dim> &fe = cell->get_fe();
	    const unsigned int level = cell->level();
	    
					     // we can presently deal only with
					     // primitive elements for boundary
					     // values. this does not preclude
					     // us using non-primitive elements
					     // in components that we aren't
					     // interested in, however. make
					     // sure that all shape functions
					     // that are non-zero for the
					     // components we are interested in,
					     // are in fact primitive
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
	    
	    typename MGDoFHandler<dim,spacedim>::face_iterator face = cell->face(face_no);
	    const unsigned char boundary_component = face->boundary_indicator();
	    if (function_map.find(boundary_component) != function_map.end())
					       // face is of the right component
	      {
						 // get indices, physical location and
						 // boundary values of dofs on this
						 // face
		local_dofs.resize (fe.dofs_per_face);
		face->get_mg_dof_indices (level, local_dofs);
		if (fe_is_system)
		  {
						     // enter those dofs
						     // into the list that
						     // match the
						     // component
						     // signature. avoid
						     // the usual
						     // complication that
						     // we can't just use
						     // *_system_to_component_index
						     // for non-primitive
						     // FEs
		    for (unsigned int i=0; i<local_dofs.size(); ++i)
		      {
			unsigned int component;
			if (fe.is_primitive())
			  component = fe.face_system_to_component_index(i).first;
			else
			  {
							     // non-primitive
							     // case. make
							     // sure that
							     // this
							     // particular
							     // shape
							     // function
							     // _is_
							     // primitive,
							     // and get at
							     // it's
							     // component. use
							     // usual
							     // trick to
							     // transfer
							     // face dof
							     // index to
							     // cell dof
			    
							     // index
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

							     // make sure
							     // that if
							     // this is
							     // not a
							     // primitive
							     // shape function,
							     // then all
							     // the
							     // corresponding
							     // components
							     // in the
							     // mask are
							     // not set
//                         if (!fe.is_primitive(cell_i))
//                           for (unsigned int c=0; c<n_components; ++c)
//                             if (fe.get_nonzero_components(cell_i)[c])
//                               Assert (component_mask[c] == false,
//                                       ExcFENotPrimitive());

// let's pick the first of possibly more than one non-zero
// components. if shape function is non-primitive, then we will ignore
// the result in the following anyway, otherwise there's only one
// non-zero component which we will use
			    component = (std::find (fe.get_nonzero_components(cell_i).begin(),
						    fe.get_nonzero_components(cell_i).end(),
						    true)
					 -
					 fe.get_nonzero_components(cell_i).begin());
			  }

			if (component_mask[component] == true)
			  boundary_indices[level].insert(local_dofs[i]);
		      }
		  }
		else
		  for (unsigned int i=0; i<local_dofs.size(); ++i)
		    boundary_indices[level].insert(local_dofs[i]);
	      }
	  }
    }
}

#endif


template <int dim, int spacedim>
void
MGTools::
make_boundary_list(const MGDoFHandler<dim,spacedim>& dof,
		   const typename FunctionMap<dim>::type& function_map,
		   std::vector<IndexSet>& boundary_indices,
		   const std::vector<bool>& component_mask)
{
  Assert (boundary_indices.size() == dof.get_tria().n_levels(),
	  ExcDimensionMismatch (boundary_indices.size(),
				dof.get_tria().n_levels()));

  std::vector<std::set<unsigned int> >
    my_boundary_indices (dof.get_tria().n_levels());
  make_boundary_list (dof, function_map, my_boundary_indices, component_mask);
  for (unsigned int i=0; i<dof.get_tria().n_levels(); ++i)
    {
      boundary_indices[i] = IndexSet (dof.n_dofs(i));
      boundary_indices[i].add_indices (my_boundary_indices[i].begin(),
				       my_boundary_indices[i].end());
    }
}


#if deal_II_dimension == 1

template <int dim, int spacedim>
void
MGTools::
extract_inner_interface_dofs (
  const MGDoFHandler<dim,spacedim>&,
  std::vector<std::vector<bool> > &,
  std::vector<std::vector<bool> > &)
{
  Assert(false, ExcNotImplemented());
}

template <int dim, int spacedim>
void
MGTools::
extract_inner_interface_dofs (
  const MGDoFHandler<dim,spacedim>&,
  std::vector<std::vector<bool> > &)
{
  Assert(false, ExcNotImplemented());
}


#else

template <int dim, int spacedim>
void
MGTools::
extract_inner_interface_dofs (const MGDoFHandler<dim,spacedim> &mg_dof_handler,
			      std::vector<std::vector<bool> >  &interface_dofs)
{
  Assert (interface_dofs.size() == mg_dof_handler.get_tria().n_levels(),
	  ExcDimensionMismatch (interface_dofs.size(),
				mg_dof_handler.get_tria().n_levels()));

  for (unsigned int l=0; l<mg_dof_handler.get_tria().n_levels(); ++l)
    {
      Assert (interface_dofs[l].size() == mg_dof_handler.n_dofs(l),
	      ExcDimensionMismatch (interface_dofs[l].size(),
				    mg_dof_handler.n_dofs(l)));

      std::fill (interface_dofs[l].begin(),
		 interface_dofs[l].end(),
		 false);
    }

  const FiniteElement<dim,spacedim> &fe = mg_dof_handler.get_fe();

  const unsigned int   dofs_per_cell   = fe.dofs_per_cell;
  const unsigned int   dofs_per_face   = fe.dofs_per_face;

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);
  std::vector<bool> cell_dofs(dofs_per_cell, false);

  typename MGDoFHandler<dim>::cell_iterator cell = mg_dof_handler.begin(),
					    endc = mg_dof_handler.end();

  for (; cell!=endc; ++cell)
    {
      std::fill (cell_dofs.begin(), cell_dofs.end(), false);

      for (unsigned int face_nr=0; face_nr<GeometryInfo<dim>::faces_per_cell; ++face_nr)
	{
	  const typename DoFHandler<dim,spacedim>::face_iterator face = cell->face(face_nr);
	  if (!face->at_boundary())
	    {
					       //interior face
	      const typename MGDoFHandler<dim>::cell_iterator
		neighbor = cell->neighbor(face_nr);

	      if (neighbor->level() < cell->level())
		{
		  for (unsigned int j=0; j<dofs_per_face; ++j)
		    cell_dofs[fe.face_to_cell_index(j,face_nr)] = true;

		}
	    }
	}

      const unsigned int level = cell->level();
      cell->get_mg_dof_indices (local_dof_indices);

      for(unsigned int i=0; i<dofs_per_cell; ++i)
	  if (cell_dofs[i])
	    interface_dofs[level][local_dof_indices[i]] = true;
    }
}


template <int dim, int spacedim>
void
MGTools::
extract_inner_interface_dofs (const MGDoFHandler<dim,spacedim> &mg_dof_handler,
			      std::vector<std::vector<bool> >  &interface_dofs,
			      std::vector<std::vector<bool> >  &boundary_interface_dofs)
{
  Assert (interface_dofs.size() == mg_dof_handler.get_tria().n_levels(),
	  ExcDimensionMismatch (interface_dofs.size(),
				mg_dof_handler.get_tria().n_levels()));
  Assert (boundary_interface_dofs.size() == mg_dof_handler.get_tria().n_levels(),
	  ExcDimensionMismatch (boundary_interface_dofs.size(),
				mg_dof_handler.get_tria().n_levels()));

  for (unsigned int l=0; l<mg_dof_handler.get_tria().n_levels(); ++l)
    {
      Assert (interface_dofs[l].size() == mg_dof_handler.n_dofs(l),
	      ExcDimensionMismatch (interface_dofs[l].size(),
				    mg_dof_handler.n_dofs(l)));
      Assert (boundary_interface_dofs[l].size() == mg_dof_handler.n_dofs(l),
	      ExcDimensionMismatch (boundary_interface_dofs[l].size(),
				    mg_dof_handler.n_dofs(l)));

      std::fill (interface_dofs[l].begin(),
		 interface_dofs[l].end(),
		 false);
      std::fill (boundary_interface_dofs[l].begin(),
		 boundary_interface_dofs[l].end(),
		 false);
    }

  const FiniteElement<dim,spacedim> &fe = mg_dof_handler.get_fe();

  const unsigned int   dofs_per_cell   = fe.dofs_per_cell;
  const unsigned int   dofs_per_face   = fe.dofs_per_face;

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);
  std::vector<unsigned int> face_dof_indices (dofs_per_face);

  std::vector<bool> cell_dofs(dofs_per_cell, false);
  std::vector<bool> boundary_cell_dofs(dofs_per_cell, false);

  typename MGDoFHandler<dim>::cell_iterator cell = mg_dof_handler.begin(),
					    endc = mg_dof_handler.end();

  for (; cell!=endc; ++cell)
    {
      bool has_coarser_neighbor = false;

      std::fill (cell_dofs.begin(), cell_dofs.end(), false);
      std::fill (boundary_cell_dofs.begin(), boundary_cell_dofs.end(), false);

      for (unsigned int face_nr=0; face_nr<GeometryInfo<dim>::faces_per_cell; ++face_nr)
	{
	  const typename DoFHandler<dim,spacedim>::face_iterator face = cell->face(face_nr);
	  if (!face->at_boundary())
	    {
					       //interior face
	      const typename MGDoFHandler<dim>::cell_iterator
		neighbor = cell->neighbor(face_nr);

					       // Do refinement face
					       // from the coarse side
	      if (neighbor->level() < cell->level())
		{
		  for (unsigned int j=0; j<dofs_per_face; ++j)
		    cell_dofs[fe.face_to_cell_index(j,face_nr)] = true;

		  has_coarser_neighbor = true;
		}
	    }
	}

      if (has_coarser_neighbor == true)
	for (unsigned int face_nr=0; face_nr<GeometryInfo<dim>::faces_per_cell; ++face_nr)
	  if(cell->at_boundary(face_nr))
	    for(unsigned int j=0; j<dofs_per_face; ++j)
	      if (cell_dofs[fe.face_to_cell_index(j,face_nr)] == true) //is this necessary?
		boundary_cell_dofs[fe.face_to_cell_index(j,face_nr)] = true;


      const unsigned int level = cell->level();
      cell->get_mg_dof_indices (local_dof_indices);

      for(unsigned int i=0; i<dofs_per_cell; ++i)
	{
	  if (cell_dofs[i])
	    interface_dofs[level][local_dof_indices[i]] = true;

	  if (boundary_cell_dofs[i])
	    boundary_interface_dofs[level][local_dof_indices[i]] = true;
	}
    }
}

#endif


DEAL_II_NAMESPACE_CLOSE



// explicit instantiations

// Functions for building sparsity patterns are in a driver file. We
// call it for all different patterns known.
#define PATTERN SparsityPattern
#include "mg_tools.pattern.in.h"
#undef PATTERN

#define PATTERN CompressedSparsityPattern
#include "mg_tools.pattern.in.h"
#undef PATTERN

#define PATTERN CompressedSetSparsityPattern
#include "mg_tools.pattern.in.h"
#undef PATTERN

#define PATTERN CompressedSimpleSparsityPattern
#include "mg_tools.pattern.in.h"
#undef PATTERN

#define PATTERN BlockSparsityPattern
#include "mg_tools.pattern.in.h"
#undef PATTERN

#define PATTERN BlockCompressedSparsityPattern
#include "mg_tools.pattern.in.h"
#undef PATTERN

#define PATTERN BlockCompressedSetSparsityPattern
#include "mg_tools.pattern.in.h"
#undef PATTERN

#define PATTERN BlockCompressedSimpleSparsityPattern
#include "mg_tools.pattern.in.h"
#undef PATTERN


DEAL_II_NAMESPACE_OPEN

template void
MGTools::compute_row_length_vector(
  const MGDoFHandler<deal_II_dimension>&, unsigned int,
  std::vector<unsigned int>&, const DoFTools::Coupling);
template void
MGTools::compute_row_length_vector(
  const MGDoFHandler<deal_II_dimension>&, unsigned int,
  std::vector<unsigned int>&,
  const Table<2,DoFTools::Coupling>&, const Table<2,DoFTools::Coupling>&);

template void MGTools::count_dofs_per_component<deal_II_dimension> (
  const MGDoFHandler<deal_II_dimension>&, std::vector<std::vector<unsigned int> >&,
  bool, std::vector<unsigned int>);
template void MGTools::count_dofs_per_component<deal_II_dimension> (
  const MGDoFHandler<deal_II_dimension>&, std::vector<std::vector<unsigned int> >&,
  std::vector<unsigned int>);
template void MGTools::count_dofs_per_block<deal_II_dimension> (
  const MGDoFHandler<deal_II_dimension>&, std::vector<std::vector<unsigned int> >&,
  std::vector<unsigned int>);

template void MGTools::make_boundary_list(
  const MGDoFHandler<deal_II_dimension>&,
  const FunctionMap<deal_II_dimension>::type&,
  std::vector<std::set<unsigned int> >&,
  const std::vector<bool>&);

template void MGTools::make_boundary_list(
  const MGDoFHandler<deal_II_dimension>&,
  const FunctionMap<deal_II_dimension>::type&,
  std::vector<IndexSet>&,
  const std::vector<bool>&);

template
void
MGTools::
extract_inner_interface_dofs (const MGDoFHandler<deal_II_dimension> &mg_dof_handler,
			      std::vector<std::vector<bool> >  &interface_dofs,
			      std::vector<std::vector<bool> >  &boundary_interface_dofs);
template
void
MGTools::
extract_inner_interface_dofs (const MGDoFHandler<deal_II_dimension> &mg_dof_handler,
			      std::vector<std::vector<bool> >  &interface_dofs);

DEAL_II_NAMESPACE_CLOSE
