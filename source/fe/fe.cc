// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2014 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#include <deal.II/base/memory_consumption.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/grid/tria_boundary.h>

#include <algorithm>
#include <functional>
#include <numeric>

DEAL_II_NAMESPACE_OPEN


/*------------------------------- FiniteElement ----------------------*/


template <int dim, int spacedim>
const double FiniteElement<dim,spacedim>::fd_step_length = 1.0e-6;


template <int dim, int spacedim>
void
FiniteElement<dim,spacedim>::
InternalDataBase::initialize_2nd (const FiniteElement<dim,spacedim> *element,
                                  const Mapping<dim,spacedim>       &mapping,
                                  const Quadrature<dim>    &quadrature)
{
  // if we shall compute second
  // derivatives, then we do so by
  // finite differencing the
  // gradients. that we do by
  // evaluating the gradients of
  // shape values at points shifted
  // star-like a little in each
  // coordinate direction around each
  // quadrature point.
  //
  // therefore generate 2*dim (the
  // number of evaluation points)
  // FEValues objects with slightly
  // shifted positions
  std::vector<Point<dim> > diff_points (quadrature.size());

  differences.resize(2*dim);
  for (unsigned int d=0; d<dim; ++d)
    {
      Point<dim> shift;
      shift (d) = fd_step_length;

      // generate points and FEValues
      // objects shifted in
      // plus-direction. note that
      // they only need to compute
      // gradients, not more
      for (unsigned int i=0; i<diff_points.size(); ++i)
        diff_points[i] = quadrature.point(i) + shift;
      const Quadrature<dim> plus_quad (diff_points, quadrature.get_weights());
      differences[d] = new FEValues<dim,spacedim> (mapping, *element,
                                                   plus_quad, update_gradients);

      // now same in minus-direction
      for (unsigned int i=0; i<diff_points.size(); ++i)
        diff_points[i] = quadrature.point(i) - shift;
      const Quadrature<dim> minus_quad (diff_points, quadrature.get_weights());
      differences[d+dim] = new FEValues<dim,spacedim> (mapping, *element,
                                                       minus_quad, update_gradients);
    }
}




template <int dim, int spacedim>
FiniteElement<dim,spacedim>::InternalDataBase::~InternalDataBase ()
{
  for (unsigned int i=0; i<differences.size (); ++i)
    if (differences[i] != 0)
      {
        // delete pointer and set it
        // to zero to avoid
        // inadvertent use
        delete differences[i];
        differences[i] = 0;
      };
}




template <int dim, int spacedim>
FiniteElement<dim,spacedim>::FiniteElement (
  const FiniteElementData<dim> &fe_data,
  const std::vector<bool> &r_i_a_f,
  const std::vector<ComponentMask> &nonzero_c)
  :
  FiniteElementData<dim> (fe_data),
  adjust_quad_dof_index_for_face_orientation_table (dim == 3 ?
                                                    this->dofs_per_quad : 0 ,
                                                    dim==3 ? 8 : 0),
  adjust_line_dof_index_for_line_orientation_table (dim == 3 ?
                                                    this->dofs_per_line : 0),
  system_to_base_table(this->dofs_per_cell),
  face_system_to_base_table(this->dofs_per_face),
  component_to_base_table (this->components,
                           std::make_pair(std::make_pair(0U, 0U), 0U)),
  restriction_is_additive_flags(r_i_a_f),
  nonzero_components (nonzero_c)
{
  // Special handling of vectors of
  // length one: in this case, we
  // assume that all entries were
  // supposed to be equal.

//TODO: Do the following in a better way by expanding these arrays in the member initializer section above
  // Normally, we should be careful
  // with const_cast, but since this
  // is the constructor and we do it
  // here only, we are fine.
  unsigned int ndofs = this->dofs_per_cell;
  if (restriction_is_additive_flags.size() == 1 && ndofs > 1)
    {
      std::vector<bool> &aux
        = const_cast<std::vector<bool>&> (restriction_is_additive_flags);
      aux.resize(ndofs, restriction_is_additive_flags[0]);
    }

  if (nonzero_components.size() == 1 && ndofs > 1)
    {
      std::vector<ComponentMask> &aux
        = const_cast<std::vector<ComponentMask>&> (nonzero_components);
      aux.resize(ndofs, nonzero_components[0]);
    }

  // These used to be initialized in
  // the constructor, but here we
  // have the possibly corrected
  // nonzero_components vector.
  const_cast<std::vector<unsigned int>&>
  (n_nonzero_components_table) = compute_n_nonzero_components(nonzero_components);
  this->set_primitivity(std::find_if (n_nonzero_components_table.begin(),
                                      n_nonzero_components_table.end(),
                                      std::bind2nd(std::not_equal_to<unsigned int>(),
                                                   1U))
                        == n_nonzero_components_table.end());


  Assert (restriction_is_additive_flags.size() == this->dofs_per_cell,
          ExcDimensionMismatch(restriction_is_additive_flags.size(),
                               this->dofs_per_cell));
  AssertDimension (nonzero_components.size(), this->dofs_per_cell);
  for (unsigned int i=0; i<nonzero_components.size(); ++i)
    {
      Assert (nonzero_components[i].size() == this->n_components(),
              ExcInternalError());
      Assert (nonzero_components[i].n_selected_components ()
              >= 1,
              ExcInternalError());
      Assert (n_nonzero_components_table[i] >= 1,
              ExcInternalError());
      Assert (n_nonzero_components_table[i] <= this->n_components(),
              ExcInternalError());
    };

  // initialize some tables in the default way, i.e. if there is only one
  // (vector-)component; if the element is not primitive, leave these tables
  // empty.
  if (this->is_primitive())
    {
      system_to_component_table.resize(this->dofs_per_cell);
      face_system_to_component_table.resize(this->dofs_per_face);
      for (unsigned int j=0 ; j<this->dofs_per_cell ; ++j)
        system_to_component_table[j] = std::pair<unsigned,unsigned>(0,j);
      for (unsigned int j=0 ; j<this->dofs_per_face ; ++j)
        face_system_to_component_table[j] = std::pair<unsigned,unsigned>(0,j);
    }

  for (unsigned int j=0 ; j<this->dofs_per_cell ; ++j)
    system_to_base_table[j] = std::make_pair(std::make_pair(0U,0U),j);
  for (unsigned int j=0 ; j<this->dofs_per_face ; ++j)
    face_system_to_base_table[j] = std::make_pair(std::make_pair(0U,0U),j);

  // Fill with default value; may be changed by constructor of derived class.
  base_to_block_indices.reinit(1,1);

  // initialize the restriction and prolongation matrices. the default
  // contructur of FullMatrix<dim> initializes them with size zero
  prolongation.resize(RefinementCase<dim>::isotropic_refinement);
  restriction.resize(RefinementCase<dim>::isotropic_refinement);
  for (unsigned int ref=RefinementCase<dim>::cut_x;
       ref<RefinementCase<dim>::isotropic_refinement+1; ++ref)
    {
      prolongation[ref-1].resize (GeometryInfo<dim>::
                                  n_children(RefinementCase<dim>(ref)),
                                  FullMatrix<double>());
      restriction[ref-1].resize (GeometryInfo<dim>::
                                 n_children(RefinementCase<dim>(ref)),
                                 FullMatrix<double>());
    }

  adjust_quad_dof_index_for_face_orientation_table.fill(0);
}



template <int dim, int spacedim>
FiniteElement<dim,spacedim>::~FiniteElement ()
{}




template <int dim, int spacedim>
double
FiniteElement<dim,spacedim>::shape_value (const unsigned int,
                                          const Point<dim> &) const
{
  AssertThrow(false, ExcUnitShapeValuesDoNotExist());
  return 0.;
}



template <int dim, int spacedim>
double
FiniteElement<dim,spacedim>::shape_value_component (const unsigned int,
                                                    const Point<dim> &,
                                                    const unsigned int) const
{
  AssertThrow(false, ExcUnitShapeValuesDoNotExist());
  return 0.;
}



template <int dim, int spacedim>
Tensor<1,dim>
FiniteElement<dim,spacedim>::shape_grad (const unsigned int,
                                         const Point<dim> &) const
{
  AssertThrow(false, ExcUnitShapeValuesDoNotExist());
  return Tensor<1,dim> ();
}



template <int dim, int spacedim>
Tensor<1,dim>
FiniteElement<dim,spacedim>::shape_grad_component (const unsigned int,
                                                   const Point<dim> &,
                                                   const unsigned int) const
{
  AssertThrow(false, ExcUnitShapeValuesDoNotExist());
  return Tensor<1,dim> ();
}



template <int dim, int spacedim>
Tensor<2,dim>
FiniteElement<dim,spacedim>::shape_grad_grad (const unsigned int,
                                              const Point<dim> &) const
{
  AssertThrow(false, ExcUnitShapeValuesDoNotExist());
  return Tensor<2,dim> ();
}



template <int dim, int spacedim>
Tensor<2,dim>
FiniteElement<dim,spacedim>::shape_grad_grad_component (const unsigned int,
                                                        const Point<dim> &,
                                                        const unsigned int) const
{
  AssertThrow(false, ExcUnitShapeValuesDoNotExist());
  return Tensor<2,dim> ();
}


template <int dim, int spacedim>
void
FiniteElement<dim,spacedim>::reinit_restriction_and_prolongation_matrices (
  const bool isotropic_restriction_only,
  const bool isotropic_prolongation_only)
{
  for (unsigned int ref_case=RefinementCase<dim>::cut_x;
       ref_case <= RefinementCase<dim>::isotropic_refinement; ++ref_case)
    {
      const unsigned int nc = GeometryInfo<dim>::n_children(RefinementCase<dim>(ref_case));

      for (unsigned int i=0; i<nc; ++i)
        {
          if (this->restriction[ref_case-1][i].m() != this->dofs_per_cell
              &&
              (!isotropic_restriction_only || ref_case==RefinementCase<dim>::isotropic_refinement))
            this->restriction[ref_case-1][i].reinit (this->dofs_per_cell,
                                                     this->dofs_per_cell);
          if (this->prolongation[ref_case-1][i].m() != this->dofs_per_cell
              &&
              (!isotropic_prolongation_only || ref_case==RefinementCase<dim>::isotropic_refinement))
            this->prolongation[ref_case-1][i].reinit (this->dofs_per_cell,
                                                      this->dofs_per_cell);
        }
    }
}


template <int dim, int spacedim>
const FullMatrix<double> &
FiniteElement<dim,spacedim>::get_restriction_matrix (const unsigned int child,
                                                     const RefinementCase<dim> &refinement_case) const
{
  Assert (refinement_case<RefinementCase<dim>::isotropic_refinement+1,
          ExcIndexRange(refinement_case,0,RefinementCase<dim>::isotropic_refinement+1));
  Assert (refinement_case!=RefinementCase<dim>::no_refinement,
          ExcMessage("Restriction matrices are only available for refined cells!"));
  Assert (child<GeometryInfo<dim>::n_children(RefinementCase<dim>(refinement_case)),
          ExcIndexRange(child,0,GeometryInfo<dim>::n_children(RefinementCase<dim>(refinement_case))));
  // we use refinement_case-1 here. the -1 takes care of the origin of the
  // vector, as for RefinementCase<dim>::no_refinement (=0) there is no data
  // available and so the vector indices are shifted
  Assert (restriction[refinement_case-1][child].n() == this->dofs_per_cell, ExcProjectionVoid());
  return restriction[refinement_case-1][child];
}



template <int dim, int spacedim>
const FullMatrix<double> &
FiniteElement<dim,spacedim>::get_prolongation_matrix (const unsigned int child,
                                                      const RefinementCase<dim> &refinement_case) const
{
  Assert (refinement_case<RefinementCase<dim>::isotropic_refinement+1,
          ExcIndexRange(refinement_case,0,RefinementCase<dim>::isotropic_refinement+1));
  Assert (refinement_case!=RefinementCase<dim>::no_refinement,
          ExcMessage("Prolongation matrices are only available for refined cells!"));
  Assert (child<GeometryInfo<dim>::n_children(RefinementCase<dim>(refinement_case)),
          ExcIndexRange(child,0,GeometryInfo<dim>::n_children(RefinementCase<dim>(refinement_case))));
  // we use refinement_case-1 here. the -1 takes care
  // of the origin of the vector, as for
  // RefinementCase::no_refinement (=0) there is no
  // data available and so the vector indices
  // are shifted
  Assert (prolongation[refinement_case-1][child].n() == this->dofs_per_cell, ExcEmbeddingVoid());
  return prolongation[refinement_case-1][child];
}


//TODO:[GK] This is probably not the most efficient way of doing this.
template <int dim, int spacedim>
unsigned int
FiniteElement<dim,spacedim>::component_to_block_index (const unsigned int index) const
{
  Assert (index < this->n_components(),
          ExcIndexRange(index, 0, this->n_components()));

  return first_block_of_base(component_to_base_table[index].first.first)
         + component_to_base_table[index].second;
}


template <int dim, int spacedim>
ComponentMask
FiniteElement<dim,spacedim>::
component_mask (const FEValuesExtractors::Scalar &scalar) const
{
  AssertIndexRange(scalar.component, this->n_components());

//TODO: it would be nice to verify that it is indeed possible
// to select this scalar component, i.e., that it is not part
// of a non-primitive element. unfortunately, there is no simple
// way to write such a condition...

  std::vector<bool> mask (this->n_components(), false);
  mask[scalar.component] = true;
  return mask;
}


template <int dim, int spacedim>
ComponentMask
FiniteElement<dim,spacedim>::
component_mask (const FEValuesExtractors::Vector &vector) const
{
  AssertIndexRange(vector.first_vector_component+dim-1, this->n_components());

  //TODO: it would be nice to verify that it is indeed possible
  // to select these vector components, i.e., that they don't span
  // beyond the beginning or end of anon-primitive element.
  // unfortunately, there is no simple way to write such a condition...

  std::vector<bool> mask (this->n_components(), false);
  for (unsigned int c=vector.first_vector_component; c<vector.first_vector_component+dim; ++c)
    mask[c] = true;
  return mask;
}


template <int dim, int spacedim>
ComponentMask
FiniteElement<dim,spacedim>::
component_mask (const FEValuesExtractors::SymmetricTensor<2> &sym_tensor) const
{
  AssertIndexRange((sym_tensor.first_tensor_component +
                    SymmetricTensor<2,dim>::n_independent_components-1),
                   this->n_components());

  //TODO: it would be nice to verify that it is indeed possible
  // to select these vector components, i.e., that they don't span
  // beyond the beginning or end of anon-primitive element.
  // unfortunately, there is no simple way to write such a condition...

  std::vector<bool> mask (this->n_components(), false);
  for (unsigned int c=sym_tensor.first_tensor_component;
       c<sym_tensor.first_tensor_component+SymmetricTensor<2,dim>::n_independent_components; ++c)
    mask[c] = true;
  return mask;
}



template <int dim, int spacedim>
ComponentMask
FiniteElement<dim,spacedim>::
component_mask (const BlockMask &block_mask) const
{
  // if we get a block mask that represents all blocks, then
  // do the same for the returned component mask
  if (block_mask.represents_the_all_selected_mask())
    return ComponentMask();

  AssertDimension(block_mask.size(), this->n_blocks());

  std::vector<bool> component_mask (this->n_components(), false);
  for (unsigned int c=0; c<this->n_components(); ++c)
    if (block_mask[component_to_block_index(c)] == true)
      component_mask[c] = true;

  return component_mask;
}



template <int dim, int spacedim>
BlockMask
FiniteElement<dim,spacedim>::
block_mask (const FEValuesExtractors::Scalar &scalar) const
{
  // simply create the corresponding component mask (a simpler
  // process) and then convert it to a block mask
  return block_mask(component_mask(scalar));
}


template <int dim, int spacedim>
BlockMask
FiniteElement<dim,spacedim>::
block_mask (const FEValuesExtractors::Vector &vector) const
{
  // simply create the corresponding component mask (a simpler
  // process) and then convert it to a block mask
  return block_mask(component_mask(vector));
}


template <int dim, int spacedim>
BlockMask
FiniteElement<dim,spacedim>::
block_mask (const FEValuesExtractors::SymmetricTensor<2> &sym_tensor) const
{
  // simply create the corresponding component mask (a simpler
  // process) and then convert it to a block mask
  return block_mask(component_mask(sym_tensor));
}



template <int dim, int spacedim>
BlockMask
FiniteElement<dim,spacedim>::
block_mask (const ComponentMask &component_mask) const
{
  // if we get a component mask that represents all component, then
  // do the same for the returned block mask
  if (component_mask.represents_the_all_selected_mask())
    return BlockMask();

  AssertDimension(component_mask.size(), this->n_components());

  // walk over all of the components
  // of this finite element and see
  // if we need to set the
  // corresponding block. inside the
  // block, walk over all the
  // components that correspond to
  // this block and make sure the
  // component mask is set for all of
  // them
  std::vector<bool> block_mask (this->n_blocks(), false);
  for (unsigned int c=0; c<this->n_components();)
    {
      const unsigned int block = component_to_block_index(c);
      if (component_mask[c] == true)
        block_mask[block] = true;

      // now check all of the other
      // components that correspond
      // to this block
      ++c;
      while ((c<this->n_components())
             &&
             (component_to_block_index(c) == block))
        {
          Assert (component_mask[c] == block_mask[block],
                  ExcMessage ("The component mask argument given to this function "
                              "is not a mask where the individual components belonging "
                              "to one block of the finite element are either all "
                              "selected or not selected. You can't call this function "
                              "with a component mask that splits blocks."));
          ++c;
        }
    }


  return block_mask;
}



template <int dim, int spacedim>
unsigned int
FiniteElement<dim,spacedim>::
face_to_cell_index (const unsigned int face_index,
                    const unsigned int face,
                    const bool face_orientation,
                    const bool face_flip,
                    const bool face_rotation) const
{
  Assert (face_index < this->dofs_per_face,
          ExcIndexRange(face_index, 0, this->dofs_per_face));
  Assert (face < GeometryInfo<dim>::faces_per_cell,
          ExcIndexRange(face, 0, GeometryInfo<dim>::faces_per_cell));

//TODO: we could presumably solve the 3d case below using the
// adjust_quad_dof_index_for_face_orientation_table field. for the
// 2d case, we can't use adjust_line_dof_index_for_line_orientation_table
// since that array is empty (presumably because we thought that
// there are no flipped edges in 2d, but these can happen in
// DoFTools::make_periodicity_constraints, for example). so we
// would need to either fill this field, or rely on derived classes
// implementing this function, as we currently do

  // see the function's documentation for an explanation of this
  // assertion -- in essence, derived classes have to implement
  // an overloaded version of this function if we are to use any
  // other than standard orientation
  if ((face_orientation != true) || (face_flip != false) || (face_rotation != false))
    Assert ((this->dofs_per_line <= 1) && (this->dofs_per_quad <= 1),
            ExcMessage ("The function in this base class can not handle this case. "
                        "Rather, the derived class you are using must provide "
                        "an overloaded version but apparently hasn't done so. See "
                        "the documentation of this function for more information."));

  // we need to distinguish between DoFs on vertices, lines and in 3d quads.
  // do so in a sequence of if-else statements
  if (face_index < this->first_face_line_index)
    // DoF is on a vertex
    {
      // get the number of the vertex on the face that corresponds to this DoF,
      // along with the number of the DoF on this vertex
      const unsigned int face_vertex         = face_index / this->dofs_per_vertex;
      const unsigned int dof_index_on_vertex = face_index % this->dofs_per_vertex;

      // then get the number of this vertex on the cell and translate
      // this to a DoF number on the cell
      return (GeometryInfo<dim>::face_to_cell_vertices(face, face_vertex,
                                                       face_orientation,
                                                       face_flip,
                                                       face_rotation)
              * this->dofs_per_vertex
              +
              dof_index_on_vertex);
    }
  else if (face_index < this->first_face_quad_index)
    // DoF is on a face
    {
      // do the same kind of translation as before. we need to only consider
      // DoFs on the lines, i.e., ignoring those on the vertices
      const unsigned int index = face_index - this->first_face_line_index;

      const unsigned int face_line         = index / this->dofs_per_line;
      const unsigned int dof_index_on_line = index % this->dofs_per_line;

      return (this->first_line_index
              + GeometryInfo<dim>::face_to_cell_lines(face, face_line,
                                                      face_orientation,
                                                      face_flip,
                                                      face_rotation)
              * this->dofs_per_line
              +
              dof_index_on_line);
    }
  else
    // DoF is on a quad
    {
      Assert (dim >= 3, ExcInternalError());

      // ignore vertex and line dofs
      const unsigned int index = face_index - this->first_face_quad_index;

      return (this->first_quad_index
              + face * this->dofs_per_quad
              + index);
    }
}




template <int dim, int spacedim>
unsigned int
FiniteElement<dim,spacedim>::adjust_quad_dof_index_for_face_orientation (const unsigned int index,
    const bool face_orientation,
    const bool face_flip,
    const bool face_rotation) const
{
  // general template for 1D and 2D: not
  // implemented. in fact, the function
  // shouldn't even be called unless we are
  // in 3d, so throw an internal error
  Assert (dim==3, ExcInternalError());
  if (dim < 3)
    return index;

  // adjust dofs on 3d faces if the face is
  // flipped. note that we query a table that
  // derived elements need to have set up
  // front. the exception are discontinuous
  // elements for which there should be no
  // face dofs anyway (i.e. dofs_per_quad==0
  // in 3d), so we don't need the table, but
  // the function should also not have been
  // called
  Assert (index<this->dofs_per_quad, ExcIndexRange(index,0,this->dofs_per_quad));
  Assert (adjust_quad_dof_index_for_face_orientation_table.n_elements()==8*this->dofs_per_quad,
          ExcInternalError());
  return index+adjust_quad_dof_index_for_face_orientation_table(index,4*face_orientation+2*face_flip+face_rotation);
}



template <int dim, int spacedim>
unsigned int
FiniteElement<dim,spacedim>::adjust_line_dof_index_for_line_orientation (const unsigned int index,
    const bool line_orientation) const
{
  // general template for 1D and 2D: do
  // nothing. Do not throw an Assertion,
  // however, in order to allow to call this
  // function in 2D as well
  if (dim<3)
    return index;

  Assert (index<this->dofs_per_line, ExcIndexRange(index,0,this->dofs_per_line));
  Assert (adjust_line_dof_index_for_line_orientation_table.size()==this->dofs_per_line,
          ExcInternalError());
  if (line_orientation)
    return index;
  else
    return index+adjust_line_dof_index_for_line_orientation_table[index];
}



template <int dim, int spacedim>
bool
FiniteElement<dim,spacedim>::prolongation_is_implemented () const
{
  for (unsigned int ref_case=RefinementCase<dim>::cut_x;
       ref_case<RefinementCase<dim>::isotropic_refinement+1; ++ref_case)
    for (unsigned int c=0;
         c<GeometryInfo<dim>::n_children(RefinementCase<dim>(ref_case)); ++c)
      {
        // make sure also the lazily initialized matrices are created
        get_prolongation_matrix(c, RefinementCase<dim>(ref_case));
        Assert ((prolongation[ref_case-1][c].m() == this->dofs_per_cell) ||
                (prolongation[ref_case-1][c].m() == 0),
                ExcInternalError());
        Assert ((prolongation[ref_case-1][c].n() == this->dofs_per_cell) ||
                (prolongation[ref_case-1][c].n() == 0),
                ExcInternalError());
        if ((prolongation[ref_case-1][c].m() == 0) ||
            (prolongation[ref_case-1][c].n() == 0))
          return false;
      }
  return true;
}



template <int dim, int spacedim>
bool
FiniteElement<dim,spacedim>::restriction_is_implemented () const
{
  for (unsigned int ref_case=RefinementCase<dim>::cut_x;
       ref_case<RefinementCase<dim>::isotropic_refinement+1; ++ref_case)
    for (unsigned int c=0;
         c<GeometryInfo<dim>::n_children(RefinementCase<dim>(ref_case)); ++c)
      {
        // make sure also the lazily initialized matrices are created
        get_restriction_matrix(c, RefinementCase<dim>(ref_case));
        Assert ((restriction[ref_case-1][c].m() == this->dofs_per_cell) ||
                (restriction[ref_case-1][c].m() == 0),
                ExcInternalError());
        Assert ((restriction[ref_case-1][c].n() == this->dofs_per_cell) ||
                (restriction[ref_case-1][c].n() == 0),
                ExcInternalError());
        if ((restriction[ref_case-1][c].m() == 0) ||
            (restriction[ref_case-1][c].n() == 0))
          return false;
      }
  return true;
}



template <int dim, int spacedim>
bool
FiniteElement<dim,spacedim>::isotropic_prolongation_is_implemented () const
{
  const RefinementCase<dim> ref_case=RefinementCase<dim>::isotropic_refinement;

  for (unsigned int c=0;
       c<GeometryInfo<dim>::n_children(RefinementCase<dim>(ref_case)); ++c)
    {
      // make sure also the lazily initialized matrices are created
      get_prolongation_matrix(c, RefinementCase<dim>(ref_case));
      Assert ((prolongation[ref_case-1][c].m() == this->dofs_per_cell) ||
              (prolongation[ref_case-1][c].m() == 0),
              ExcInternalError());
      Assert ((prolongation[ref_case-1][c].n() == this->dofs_per_cell) ||
              (prolongation[ref_case-1][c].n() == 0),
              ExcInternalError());
      if ((prolongation[ref_case-1][c].m() == 0) ||
          (prolongation[ref_case-1][c].n() == 0))
        return false;
    }
  return true;
}



template <int dim, int spacedim>
bool
FiniteElement<dim,spacedim>::isotropic_restriction_is_implemented () const
{
  const RefinementCase<dim> ref_case = RefinementCase<dim>::isotropic_refinement;

  for (unsigned int c=0;
       c<GeometryInfo<dim>::n_children(RefinementCase<dim>(ref_case)); ++c)
    {
      // make sure also the lazily initialized matrices are created
      get_restriction_matrix(c, RefinementCase<dim>(ref_case));
      Assert ((restriction[ref_case-1][c].m() == this->dofs_per_cell) ||
              (restriction[ref_case-1][c].m() == 0),
              ExcInternalError());
      Assert ((restriction[ref_case-1][c].n() == this->dofs_per_cell) ||
              (restriction[ref_case-1][c].n() == 0),
              ExcInternalError());
      if ((restriction[ref_case-1][c].m() == 0) ||
          (restriction[ref_case-1][c].n() == 0))
        return false;
    }
  return true;
}



template <int dim, int spacedim>
bool
FiniteElement<dim,spacedim>::constraints_are_implemented (const internal::SubfaceCase<dim> &subface_case) const
{
  if (subface_case==internal::SubfaceCase<dim>::case_isotropic)
    return (this->dofs_per_face  == 0) || (interface_constraints.m() != 0);
  else
    return false;
}



template <int dim, int spacedim>
bool
FiniteElement<dim,spacedim>::hp_constraints_are_implemented () const
{
  return false;
}



template <int dim, int spacedim>
const FullMatrix<double> &
FiniteElement<dim,spacedim>::constraints (const internal::SubfaceCase<dim> &subface_case) const
{
  Assert (subface_case==internal::SubfaceCase<dim>::case_isotropic, ExcConstraintsVoid());
  Assert ((this->dofs_per_face  == 0) || (interface_constraints.m() != 0),
          ExcConstraintsVoid());

  if (dim==1)
    Assert ((interface_constraints.m()==0) && (interface_constraints.n()==0),
            ExcWrongInterfaceMatrixSize(interface_constraints.m(),
                                        interface_constraints.n()));

  return interface_constraints;
}



template <int dim, int spacedim>
TableIndices<2>
FiniteElement<dim,spacedim>::interface_constraints_size () const
{
  switch (dim)
    {
    case 1:
      return TableIndices<2> (0U, 0U);
    case 2:
      return TableIndices<2> (this->dofs_per_vertex +
                              2*this->dofs_per_line,
                              this->dofs_per_face);
    case 3:
      return TableIndices<2> (5*this->dofs_per_vertex +
                              12*this->dofs_per_line  +
                              4*this->dofs_per_quad,
                              this->dofs_per_face);
    default:
      Assert (false, ExcNotImplemented());
    };
  return TableIndices<2> (numbers::invalid_unsigned_int,
                          numbers::invalid_unsigned_int);
}



template <int dim, int spacedim>
void
FiniteElement<dim,spacedim>::
get_interpolation_matrix (const FiniteElement<dim,spacedim> &,
                          FullMatrix<double> &) const
{
  // by default, no interpolation
  // implemented. so throw exception,
  // as documentation says
  typedef FiniteElement<dim,spacedim> FEE;
  AssertThrow (false,
               typename FEE::
               ExcInterpolationNotImplemented());
}



template <int dim, int spacedim>
void
FiniteElement<dim,spacedim>::
get_face_interpolation_matrix (const FiniteElement<dim,spacedim> &,
                               FullMatrix<double> &) const
{
  // by default, no interpolation
  // implemented. so throw exception,
  // as documentation says
  typedef    FiniteElement<dim,spacedim> FEE;
  AssertThrow (false,
               typename FEE::
               ExcInterpolationNotImplemented());
}



template <int dim, int spacedim>
void
FiniteElement<dim,spacedim>::
get_subface_interpolation_matrix (const FiniteElement<dim,spacedim> &,
                                  const unsigned int,
                                  FullMatrix<double> &) const
{
  // by default, no interpolation
  // implemented. so throw exception,
  // as documentation says
  typedef    FiniteElement<dim,spacedim> FEE;
  AssertThrow (false,
               typename FEE::ExcInterpolationNotImplemented());
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int> >
FiniteElement<dim,spacedim>::
hp_vertex_dof_identities (const FiniteElement<dim,spacedim> &) const
{
  Assert (false, ExcNotImplemented());
  return std::vector<std::pair<unsigned int, unsigned int> > ();
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int> >
FiniteElement<dim,spacedim>::
hp_line_dof_identities (const FiniteElement<dim,spacedim> &) const
{
  Assert (false, ExcNotImplemented());
  return std::vector<std::pair<unsigned int, unsigned int> > ();
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int> >
FiniteElement<dim,spacedim>::
hp_quad_dof_identities (const FiniteElement<dim,spacedim> &) const
{
  Assert (false, ExcNotImplemented());
  return std::vector<std::pair<unsigned int, unsigned int> > ();
}



template <int dim, int spacedim>
FiniteElementDomination::Domination
FiniteElement<dim,spacedim>::
compare_for_face_domination (const FiniteElement<dim,spacedim> &) const
{
  Assert (false, ExcNotImplemented());
  return FiniteElementDomination::neither_element_dominates;
}



template <int dim, int spacedim>
bool
FiniteElement<dim,spacedim>::operator == (const FiniteElement<dim,spacedim> &f) const
{
  return ((static_cast<const FiniteElementData<dim>&>(*this) ==
           static_cast<const FiniteElementData<dim>&>(f)) &&
          (interface_constraints == f.interface_constraints));
}



template <int dim, int spacedim>
const std::vector<Point<dim> > &
FiniteElement<dim,spacedim>::get_unit_support_points () const
{
  // a finite element may define
  // support points, but only if
  // there are as many as there are
  // degrees of freedom
  Assert ((unit_support_points.size() == 0) ||
          (unit_support_points.size() == this->dofs_per_cell),
          ExcInternalError());
  return unit_support_points;
}



template <int dim, int spacedim>
bool
FiniteElement<dim,spacedim>::has_support_points () const
{
  return (unit_support_points.size() != 0);
}



template <int dim, int spacedim>
const std::vector<Point<dim> > &
FiniteElement<dim,spacedim>::get_generalized_support_points () const
{
  // a finite element may define
  // support points, but only if
  // there are as many as there are
  // degrees of freedom
  return ((generalized_support_points.size() == 0)
          ? unit_support_points
          : generalized_support_points);
}



template <int dim, int spacedim>
bool
FiniteElement<dim,spacedim>::has_generalized_support_points () const
{
  return (get_generalized_support_points().size() != 0);
}



template <int dim, int spacedim>
Point<dim>
FiniteElement<dim,spacedim>::unit_support_point (const unsigned int index) const
{
  Assert (index < this->dofs_per_cell,
          ExcIndexRange (index, 0, this->dofs_per_cell));
  Assert (unit_support_points.size() == this->dofs_per_cell,
          ExcFEHasNoSupportPoints ());
  return unit_support_points[index];
}



template <int dim, int spacedim>
const std::vector<Point<dim-1> > &
FiniteElement<dim,spacedim>::get_unit_face_support_points () const
{
  // a finite element may define
  // support points, but only if
  // there are as many as there are
  // degrees of freedom on a face
  Assert ((unit_face_support_points.size() == 0) ||
          (unit_face_support_points.size() == this->dofs_per_face),
          ExcInternalError());
  return unit_face_support_points;
}



template <int dim, int spacedim>
bool
FiniteElement<dim,spacedim>::has_face_support_points () const
{
  return (unit_face_support_points.size() != 0);
}



template <int dim, int spacedim>
const std::vector<Point<dim-1> > &
FiniteElement<dim,spacedim>::get_generalized_face_support_points () const
{
  // a finite element may define
  // support points, but only if
  // there are as many as there are
  // degrees of freedom on a face
  return ((generalized_face_support_points.size() == 0)
          ? unit_face_support_points
          : generalized_face_support_points);
}



template <int dim, int spacedim>
bool
FiniteElement<dim,spacedim>::has_generalized_face_support_points () const
{
  return (generalized_face_support_points.size() != 0);
}



template <int dim, int spacedim>
Point<dim-1>
FiniteElement<dim,spacedim>::unit_face_support_point (const unsigned int index) const
{
  Assert (index < this->dofs_per_face,
          ExcIndexRange (index, 0, this->dofs_per_face));
  Assert (unit_face_support_points.size() == this->dofs_per_face,
          ExcFEHasNoSupportPoints ());
  return unit_face_support_points[index];
}



template <int dim, int spacedim>
bool
FiniteElement<dim,spacedim>::has_support_on_face (
  const unsigned int,
  const unsigned int) const
{
  return true;
}



template <int dim, int spacedim>
std::pair<Table<2,bool>, std::vector<unsigned int> >
FiniteElement<dim,spacedim>::get_constant_modes () const
{
  Assert (false, ExcNotImplemented());
  return std::pair<Table<2,bool>, std::vector<unsigned int> >
         (Table<2,bool>(this->n_components(), this->dofs_per_cell),
          std::vector<unsigned int>(this->n_components()));
}



template <int dim, int spacedim>
void
FiniteElement<dim,spacedim>::interpolate(
  std::vector<double>       &local_dofs,
  const std::vector<double> &values) const
{
  Assert (has_support_points(), ExcFEHasNoSupportPoints());
  Assert (values.size() == unit_support_points.size(),
          ExcDimensionMismatch(values.size(), unit_support_points.size()));
  Assert (local_dofs.size() == this->dofs_per_cell,
          ExcDimensionMismatch(local_dofs.size(),this->dofs_per_cell));
  Assert (this->n_components() == 1,
          ExcDimensionMismatch(this->n_components(), 1));

  std::copy(values.begin(), values.end(), local_dofs.begin());
}




template <int dim, int spacedim>
void
FiniteElement<dim,spacedim>::interpolate(
  std::vector<double>    &local_dofs,
  const std::vector<Vector<double> > &values,
  unsigned int offset) const
{
  Assert (has_support_points(), ExcFEHasNoSupportPoints());
  Assert (values.size() == unit_support_points.size(),
          ExcDimensionMismatch(values.size(), unit_support_points.size()));
  Assert (local_dofs.size() == this->dofs_per_cell,
          ExcDimensionMismatch(local_dofs.size(),this->dofs_per_cell));
  Assert (values[0].size() >= offset+this->n_components(),
          ExcDimensionMismatch(values[0].size(),offset+this->n_components()));

  for (unsigned int i=0; i<this->dofs_per_cell; ++i)
    {
      const std::pair<unsigned int, unsigned int> index
        = this->system_to_component_index(i);
      local_dofs[i] = values[i](offset+index.first);
    }
}




template <int dim, int spacedim>
void
FiniteElement<dim,spacedim>::interpolate(
  std::vector<double> &local_dofs,
  const VectorSlice<const std::vector<std::vector<double> > > &values) const
{
  Assert (has_support_points(), ExcFEHasNoSupportPoints());
  Assert (values[0].size() == unit_support_points.size(),
          ExcDimensionMismatch(values.size(), unit_support_points.size()));
  Assert (local_dofs.size() == this->dofs_per_cell,
          ExcDimensionMismatch(local_dofs.size(),this->dofs_per_cell));
  Assert (values.size() == this->n_components(),
          ExcDimensionMismatch(values.size(), this->n_components()));

  for (unsigned int i=0; i<this->dofs_per_cell; ++i)
    {
      const std::pair<unsigned int, unsigned int> index
        = this->system_to_component_index(i);
      local_dofs[i] = values[index.first][i];
    }
}



template <int dim, int spacedim>
std::size_t
FiniteElement<dim,spacedim>::memory_consumption () const
{
  return (sizeof(FiniteElementData<dim>) +
          MemoryConsumption::memory_consumption (restriction)+
          MemoryConsumption::memory_consumption (prolongation) +
          MemoryConsumption::memory_consumption (interface_constraints) +
          MemoryConsumption::memory_consumption (system_to_component_table) +
          MemoryConsumption::memory_consumption (face_system_to_component_table) +
          MemoryConsumption::memory_consumption (system_to_base_table) +
          MemoryConsumption::memory_consumption (face_system_to_base_table) +
          MemoryConsumption::memory_consumption (component_to_base_table) +
          MemoryConsumption::memory_consumption (restriction_is_additive_flags) +
          MemoryConsumption::memory_consumption (nonzero_components) +
          MemoryConsumption::memory_consumption (n_nonzero_components_table));
}



template<>
void
FiniteElement<1,2>::compute_2nd (
  const Mapping<1,2> &,
  const Triangulation<1,2>::cell_iterator &,
  const unsigned int,
  Mapping<1,2>::InternalDataBase &,
  InternalDataBase &,
  FEValuesData<1,2> &) const
{

  Assert(false, ExcNotImplemented());
}


template<>
void
FiniteElement<1,3>::compute_2nd (
  const Mapping<1,3> &,
  const Triangulation<1,3>::cell_iterator &,
  const unsigned int,
  Mapping<1,3>::InternalDataBase &,
  InternalDataBase &,
  FEValuesData<1,3> &) const
{

  Assert(false, ExcNotImplemented());
}



template<>
void
FiniteElement<2,3>::compute_2nd (
  const Mapping<2,3> &,
  const Triangulation<2,3>::cell_iterator &,
  const unsigned int,
  Mapping<2,3>::InternalDataBase &,
  InternalDataBase &,
  FEValuesData<2,3> &) const
{

  Assert(false, ExcNotImplemented());
}



template <int dim, int spacedim>
void
FiniteElement<dim,spacedim>::compute_2nd (
  const Mapping<dim,spacedim>                   &mapping,
  const typename Triangulation<dim,spacedim>::cell_iterator &cell,
  const unsigned int offset,
  typename Mapping<dim,spacedim>::InternalDataBase &mapping_internal,
  InternalDataBase                     &fe_internal,
  FEValuesData<dim,spacedim>                    &data) const
{
  Assert ((fe_internal.update_each | fe_internal.update_once)
          & update_hessians,
          ExcInternalError());

  // there is nothing to do if there are no degrees of freedom (e.g., in an
  // FE_Nothing)
  if (this->dofs_per_cell == 0)
    return;

// make sure we have as many entries as there are nonzero components
//  Assert (data.shape_hessians.size() ==
//          std::accumulate (n_nonzero_components_table.begin(),
//                        n_nonzero_components_table.end(),
//                        0U),
//        ExcInternalError());
  // Number of quadrature points
  const unsigned int n_q_points = data.shape_hessians[0].size();

  // first reinit the fe_values
  // objects used for the finite
  // differencing stuff
  for (unsigned int d=0; d<dim; ++d)
    {
      fe_internal.differences[d]->reinit(cell);
      fe_internal.differences[d+dim]->reinit(cell);
      Assert(offset <= fe_internal.differences[d]->n_quadrature_points - n_q_points,
             ExcIndexRange(offset, 0, fe_internal.differences[d]->n_quadrature_points
                           - n_q_points));
    }

  // collection of difference
  // quotients of gradients in each
  // direction (first index) and at
  // all q-points (second index)
  std::vector<std::vector<Tensor<1,dim> > >
  diff_quot (spacedim, std::vector<Tensor<1,dim> > (n_q_points));
  std::vector<Tensor<1,spacedim> > diff_quot2 (n_q_points);

  // for all nonzero components of
  // all shape functions at all
  // quadrature points and difference
  // quotients in all directions:
  unsigned int total_index = 0;
  for (unsigned int shape_index=0; shape_index<this->dofs_per_cell; ++shape_index)
    for (unsigned int n=0; n<n_nonzero_components(shape_index); ++n, ++total_index)
      {
        for (unsigned int d1=0; d1<dim; ++d1)
          for (unsigned int q=0; q<n_q_points; ++q)
            {
              // get gradient at points
              // shifted slightly to
              // the right and to the
              // left in the present
              // coordinate direction
              //
              // note that things
              // might be more
              // difficult if the
              // shape function has
              // more than one
              // non-zero component,
              // so find out about
              // the actual component
              // if necessary
              Tensor<1,spacedim> right, left;
              if (is_primitive(shape_index))
                {
                  right = fe_internal.differences[d1]->shape_grad(shape_index, q+offset);
                  left  = fe_internal.differences[d1+dim]->shape_grad(shape_index, q+offset);
                }
              else
                {
                  // get the
                  // component index
                  // of the n-th
                  // nonzero
                  // compoment
                  unsigned int component=0;
                  for (unsigned int nonzero_comp=0; component<this->n_components();
                       ++component)
                    if (nonzero_components[shape_index][component] == true)
                      {
                        ++nonzero_comp;
                        // check
                        // whether we
                        // have found
                        // the
                        // component
                        // we are
                        // looking
                        // for. note
                        // that
                        // nonzero_comp
                        // is 1-based
                        // by the way
                        // we compute
                        // it
                        if (nonzero_comp == n+1)
                          break;
                      }
                  Assert (component < this->n_components(),
                          ExcInternalError());

                  right = fe_internal.differences[d1]
                          ->shape_grad_component(shape_index, q+offset, component);
                  left  = fe_internal.differences[d1+dim]
                          ->shape_grad_component(shape_index, q+offset, component);
                };

              // compute the second
              // derivative from a
              // symmetric difference
              // approximation
              for (unsigned int d=0; d<spacedim; ++d)
                diff_quot[d][q][d1] = 1./(2*fd_step_length) * (right[d]-left[d]);
            }

        // up to now we still have
        // difference quotients on the
        // unit cell, so transform it
        // to something on the real
        // cell
        for (unsigned int d=0; d<spacedim; ++d)
          {
            mapping.transform (diff_quot[d], diff_quot2, mapping_internal, mapping_covariant);

            for (unsigned int q=0; q<n_q_points; ++q)
              for (unsigned int d1=0; d1<spacedim; ++d1)
                data.shape_hessians[total_index][q][d][d1]
                  = diff_quot2[q][d1];
          }
      }
}



template <int dim, int spacedim>
std::vector<unsigned int>
FiniteElement<dim,spacedim>::compute_n_nonzero_components (
  const std::vector<ComponentMask> &nonzero_components)
{
  std::vector<unsigned int> retval (nonzero_components.size());
  for (unsigned int i=0; i<nonzero_components.size(); ++i)
    retval[i] = nonzero_components[i].n_selected_components();
  return retval;
}



/*------------------------------- FiniteElement ----------------------*/

template <int dim, int spacedim>
typename Mapping<dim,spacedim>::InternalDataBase *
FiniteElement<dim,spacedim>::get_face_data (const UpdateFlags       flags,
                                            const Mapping<dim,spacedim>      &mapping,
                                            const Quadrature<dim-1> &quadrature) const
{
  return get_data (flags, mapping,
                   QProjector<dim>::project_to_all_faces(quadrature));
}



template <int dim, int spacedim>
typename Mapping<dim,spacedim>::InternalDataBase *
FiniteElement<dim,spacedim>::get_subface_data (const UpdateFlags        flags,
                                               const Mapping<dim,spacedim>      &mapping,
                                               const Quadrature<dim-1> &quadrature) const
{
  return get_data (flags, mapping,
                   QProjector<dim>::project_to_all_subfaces(quadrature));
}



template <int dim, int spacedim>
const FiniteElement<dim,spacedim> &
FiniteElement<dim,spacedim>::base_element(const unsigned int index) const
{
  Assert (index==0, ExcIndexRange(index,0,1));
  // This function should not be
  // called for a system element
  Assert (base_to_block_indices.size() == 1, ExcInternalError());
  return *this;
}



/*------------------------------- Explicit Instantiations -------------*/
#include "fe.inst"


DEAL_II_NAMESPACE_CLOSE
