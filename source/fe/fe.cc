// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/quadrature.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <algorithm>
#include <functional>
#include <numeric>
#include <typeinfo>

DEAL_II_NAMESPACE_OPEN


/*------------------------------- FiniteElement ----------------------*/


template <int dim, int spacedim>
FiniteElement<dim, spacedim>::InternalDataBase::InternalDataBase()
  : update_each(update_default)
{}



template <int dim, int spacedim>
std::size_t
FiniteElement<dim, spacedim>::InternalDataBase::memory_consumption() const
{
  return sizeof(*this);
}



template <int dim, int spacedim>
FiniteElement<dim, spacedim>::FiniteElement(
  const FiniteElementData<dim> &    fe_data,
  const std::vector<bool> &         r_i_a_f,
  const std::vector<ComponentMask> &nonzero_c)
  : FiniteElementData<dim>(fe_data)
  , adjust_quad_dof_index_for_face_orientation_table(dim == 3 ?
                                                       this->dofs_per_quad :
                                                       0,
                                                     dim == 3 ? 8 : 0)
  , adjust_line_dof_index_for_line_orientation_table(
      dim == 3 ? this->dofs_per_line : 0)
  , system_to_base_table(this->dofs_per_cell)
  , face_system_to_base_table(this->dofs_per_face)
  , component_to_base_table(this->components,
                            std::make_pair(std::make_pair(0U, 0U), 0U))
  ,

  // Special handling of vectors of length one: in this case, we
  // assume that all entries were supposed to be equal
  restriction_is_additive_flags(
    r_i_a_f.size() == 1 ? std::vector<bool>(fe_data.dofs_per_cell, r_i_a_f[0]) :
                          r_i_a_f)
  , nonzero_components(
      nonzero_c.size() == 1 ?
        std::vector<ComponentMask>(fe_data.dofs_per_cell, nonzero_c[0]) :
        nonzero_c)
  , n_nonzero_components_table(compute_n_nonzero_components(nonzero_components))
  , cached_primitivity(std::find_if(n_nonzero_components_table.begin(),
                                    n_nonzero_components_table.end(),
                                    [](const unsigned int n_components) {
                                      return n_components != 1U;
                                    }) == n_nonzero_components_table.end())
{
  Assert(restriction_is_additive_flags.size() == this->dofs_per_cell,
         ExcDimensionMismatch(restriction_is_additive_flags.size(),
                              this->dofs_per_cell));
  AssertDimension(nonzero_components.size(), this->dofs_per_cell);
  for (unsigned int i = 0; i < nonzero_components.size(); ++i)
    {
      Assert(nonzero_components[i].size() == this->n_components(),
             ExcInternalError());
      Assert(nonzero_components[i].n_selected_components() >= 1,
             ExcInternalError());
      Assert(n_nonzero_components_table[i] >= 1, ExcInternalError());
      Assert(n_nonzero_components_table[i] <= this->n_components(),
             ExcInternalError());
    }

  // initialize some tables in the default way, i.e. if there is only one
  // (vector-)component; if the element is not primitive, leave these tables
  // empty.
  if (this->is_primitive())
    {
      system_to_component_table.resize(this->dofs_per_cell);
      face_system_to_component_table.resize(this->dofs_per_face);
      for (unsigned int j = 0; j < this->dofs_per_cell; ++j)
        system_to_component_table[j] = std::pair<unsigned, unsigned>(0, j);
      for (unsigned int j = 0; j < this->dofs_per_face; ++j)
        face_system_to_component_table[j] = std::pair<unsigned, unsigned>(0, j);
    }

  for (unsigned int j = 0; j < this->dofs_per_cell; ++j)
    system_to_base_table[j] = std::make_pair(std::make_pair(0U, 0U), j);
  for (unsigned int j = 0; j < this->dofs_per_face; ++j)
    face_system_to_base_table[j] = std::make_pair(std::make_pair(0U, 0U), j);

  // Fill with default value; may be changed by constructor of derived class.
  base_to_block_indices.reinit(1, 1);

  // initialize the restriction and prolongation matrices. the default
  // constructor of FullMatrix<dim> initializes them with size zero
  prolongation.resize(RefinementCase<dim>::isotropic_refinement);
  restriction.resize(RefinementCase<dim>::isotropic_refinement);
  for (unsigned int ref = RefinementCase<dim>::cut_x;
       ref < RefinementCase<dim>::isotropic_refinement + 1;
       ++ref)
    {
      prolongation[ref - 1].resize(GeometryInfo<dim>::n_children(
                                     RefinementCase<dim>(ref)),
                                   FullMatrix<double>());
      restriction[ref - 1].resize(GeometryInfo<dim>::n_children(
                                    RefinementCase<dim>(ref)),
                                  FullMatrix<double>());
    }

  adjust_quad_dof_index_for_face_orientation_table.fill(0);
}



template <int dim, int spacedim>
std::pair<std::unique_ptr<FiniteElement<dim, spacedim>>, unsigned int>
FiniteElement<dim, spacedim>::operator^(const unsigned int multiplicity) const
{
  return {this->clone(), multiplicity};
}



template <int dim, int spacedim>
double
FiniteElement<dim, spacedim>::shape_value(const unsigned int,
                                          const Point<dim> &) const
{
  AssertThrow(false, ExcUnitShapeValuesDoNotExist());
  return 0.;
}



template <int dim, int spacedim>
double
FiniteElement<dim, spacedim>::shape_value_component(const unsigned int,
                                                    const Point<dim> &,
                                                    const unsigned int) const
{
  AssertThrow(false, ExcUnitShapeValuesDoNotExist());
  return 0.;
}



template <int dim, int spacedim>
Tensor<1, dim>
FiniteElement<dim, spacedim>::shape_grad(const unsigned int,
                                         const Point<dim> &) const
{
  AssertThrow(false, ExcUnitShapeValuesDoNotExist());
  return Tensor<1, dim>();
}



template <int dim, int spacedim>
Tensor<1, dim>
FiniteElement<dim, spacedim>::shape_grad_component(const unsigned int,
                                                   const Point<dim> &,
                                                   const unsigned int) const
{
  AssertThrow(false, ExcUnitShapeValuesDoNotExist());
  return Tensor<1, dim>();
}



template <int dim, int spacedim>
Tensor<2, dim>
FiniteElement<dim, spacedim>::shape_grad_grad(const unsigned int,
                                              const Point<dim> &) const
{
  AssertThrow(false, ExcUnitShapeValuesDoNotExist());
  return Tensor<2, dim>();
}



template <int dim, int spacedim>
Tensor<2, dim>
FiniteElement<dim, spacedim>::shape_grad_grad_component(
  const unsigned int,
  const Point<dim> &,
  const unsigned int) const
{
  AssertThrow(false, ExcUnitShapeValuesDoNotExist());
  return Tensor<2, dim>();
}



template <int dim, int spacedim>
Tensor<3, dim>
FiniteElement<dim, spacedim>::shape_3rd_derivative(const unsigned int,
                                                   const Point<dim> &) const
{
  AssertThrow(false, ExcUnitShapeValuesDoNotExist());
  return Tensor<3, dim>();
}



template <int dim, int spacedim>
Tensor<3, dim>
FiniteElement<dim, spacedim>::shape_3rd_derivative_component(
  const unsigned int,
  const Point<dim> &,
  const unsigned int) const
{
  AssertThrow(false, ExcUnitShapeValuesDoNotExist());
  return Tensor<3, dim>();
}



template <int dim, int spacedim>
Tensor<4, dim>
FiniteElement<dim, spacedim>::shape_4th_derivative(const unsigned int,
                                                   const Point<dim> &) const
{
  AssertThrow(false, ExcUnitShapeValuesDoNotExist());
  return Tensor<4, dim>();
}



template <int dim, int spacedim>
Tensor<4, dim>
FiniteElement<dim, spacedim>::shape_4th_derivative_component(
  const unsigned int,
  const Point<dim> &,
  const unsigned int) const
{
  AssertThrow(false, ExcUnitShapeValuesDoNotExist());
  return Tensor<4, dim>();
}


template <int dim, int spacedim>
void
FiniteElement<dim, spacedim>::reinit_restriction_and_prolongation_matrices(
  const bool isotropic_restriction_only,
  const bool isotropic_prolongation_only)
{
  for (unsigned int ref_case = RefinementCase<dim>::cut_x;
       ref_case <= RefinementCase<dim>::isotropic_refinement;
       ++ref_case)
    {
      const unsigned int nc =
        GeometryInfo<dim>::n_children(RefinementCase<dim>(ref_case));

      for (unsigned int i = 0; i < nc; ++i)
        {
          if (this->restriction[ref_case - 1][i].m() != this->dofs_per_cell &&
              (!isotropic_restriction_only ||
               ref_case == RefinementCase<dim>::isotropic_refinement))
            this->restriction[ref_case - 1][i].reinit(this->dofs_per_cell,
                                                      this->dofs_per_cell);
          if (this->prolongation[ref_case - 1][i].m() != this->dofs_per_cell &&
              (!isotropic_prolongation_only ||
               ref_case == RefinementCase<dim>::isotropic_refinement))
            this->prolongation[ref_case - 1][i].reinit(this->dofs_per_cell,
                                                       this->dofs_per_cell);
        }
    }
}


template <int dim, int spacedim>
const FullMatrix<double> &
FiniteElement<dim, spacedim>::get_restriction_matrix(
  const unsigned int         child,
  const RefinementCase<dim> &refinement_case) const
{
  AssertIndexRange(refinement_case,
                   RefinementCase<dim>::isotropic_refinement + 1);
  Assert(refinement_case != RefinementCase<dim>::no_refinement,
         ExcMessage(
           "Restriction matrices are only available for refined cells!"));
  AssertIndexRange(
    child, GeometryInfo<dim>::n_children(RefinementCase<dim>(refinement_case)));
  // we use refinement_case-1 here. the -1 takes care of the origin of the
  // vector, as for RefinementCase<dim>::no_refinement (=0) there is no data
  // available and so the vector indices are shifted
  Assert(restriction[refinement_case - 1][child].n() == this->dofs_per_cell,
         ExcProjectionVoid());
  return restriction[refinement_case - 1][child];
}



template <int dim, int spacedim>
const FullMatrix<double> &
FiniteElement<dim, spacedim>::get_prolongation_matrix(
  const unsigned int         child,
  const RefinementCase<dim> &refinement_case) const
{
  AssertIndexRange(refinement_case,
                   RefinementCase<dim>::isotropic_refinement + 1);
  Assert(refinement_case != RefinementCase<dim>::no_refinement,
         ExcMessage(
           "Prolongation matrices are only available for refined cells!"));
  AssertIndexRange(
    child, GeometryInfo<dim>::n_children(RefinementCase<dim>(refinement_case)));
  // we use refinement_case-1 here. the -1 takes care
  // of the origin of the vector, as for
  // RefinementCase::no_refinement (=0) there is no
  // data available and so the vector indices
  // are shifted
  Assert(prolongation[refinement_case - 1][child].n() == this->dofs_per_cell,
         ExcEmbeddingVoid());
  return prolongation[refinement_case - 1][child];
}


// TODO:[GK] This is probably not the most efficient way of doing this.
template <int dim, int spacedim>
unsigned int
FiniteElement<dim, spacedim>::component_to_block_index(
  const unsigned int index) const
{
  AssertIndexRange(index, this->n_components());

  return first_block_of_base(component_to_base_table[index].first.first) +
         component_to_base_table[index].second;
}


template <int dim, int spacedim>
ComponentMask
FiniteElement<dim, spacedim>::component_mask(
  const FEValuesExtractors::Scalar &scalar) const
{
  AssertIndexRange(scalar.component, this->n_components());

  // TODO: it would be nice to verify that it is indeed possible
  // to select this scalar component, i.e., that it is not part
  // of a non-primitive element. unfortunately, there is no simple
  // way to write such a condition...

  std::vector<bool> mask(this->n_components(), false);
  mask[scalar.component] = true;
  return mask;
}


template <int dim, int spacedim>
ComponentMask
FiniteElement<dim, spacedim>::component_mask(
  const FEValuesExtractors::Vector &vector) const
{
  AssertIndexRange(vector.first_vector_component + dim - 1,
                   this->n_components());

  // TODO: it would be nice to verify that it is indeed possible
  // to select these vector components, i.e., that they don't span
  // beyond the beginning or end of anon-primitive element.
  // unfortunately, there is no simple way to write such a condition...

  std::vector<bool> mask(this->n_components(), false);
  for (unsigned int c = vector.first_vector_component;
       c < vector.first_vector_component + dim;
       ++c)
    mask[c] = true;
  return mask;
}


template <int dim, int spacedim>
ComponentMask
FiniteElement<dim, spacedim>::component_mask(
  const FEValuesExtractors::SymmetricTensor<2> &sym_tensor) const
{
  AssertIndexRange((sym_tensor.first_tensor_component +
                    SymmetricTensor<2, dim>::n_independent_components - 1),
                   this->n_components());

  // TODO: it would be nice to verify that it is indeed possible
  // to select these vector components, i.e., that they don't span
  // beyond the beginning or end of anon-primitive element.
  // unfortunately, there is no simple way to write such a condition...

  std::vector<bool> mask(this->n_components(), false);
  for (unsigned int c = sym_tensor.first_tensor_component;
       c < sym_tensor.first_tensor_component +
             SymmetricTensor<2, dim>::n_independent_components;
       ++c)
    mask[c] = true;
  return mask;
}



template <int dim, int spacedim>
ComponentMask
FiniteElement<dim, spacedim>::component_mask(const BlockMask &block_mask) const
{
  // if we get a block mask that represents all blocks, then
  // do the same for the returned component mask
  if (block_mask.represents_the_all_selected_mask())
    return {};

  AssertDimension(block_mask.size(), this->n_blocks());

  std::vector<bool> component_mask(this->n_components(), false);
  for (unsigned int c = 0; c < this->n_components(); ++c)
    if (block_mask[component_to_block_index(c)] == true)
      component_mask[c] = true;

  return component_mask;
}



template <int dim, int spacedim>
BlockMask
FiniteElement<dim, spacedim>::block_mask(
  const FEValuesExtractors::Scalar &scalar) const
{
  // simply create the corresponding component mask (a simpler
  // process) and then convert it to a block mask
  return block_mask(component_mask(scalar));
}


template <int dim, int spacedim>
BlockMask
FiniteElement<dim, spacedim>::block_mask(
  const FEValuesExtractors::Vector &vector) const
{
  // simply create the corresponding component mask (a simpler
  // process) and then convert it to a block mask
  return block_mask(component_mask(vector));
}


template <int dim, int spacedim>
BlockMask
FiniteElement<dim, spacedim>::block_mask(
  const FEValuesExtractors::SymmetricTensor<2> &sym_tensor) const
{
  // simply create the corresponding component mask (a simpler
  // process) and then convert it to a block mask
  return block_mask(component_mask(sym_tensor));
}



template <int dim, int spacedim>
BlockMask
FiniteElement<dim, spacedim>::block_mask(
  const ComponentMask &component_mask) const
{
  // if we get a component mask that represents all component, then
  // do the same for the returned block mask
  if (component_mask.represents_the_all_selected_mask())
    return {};

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
  std::vector<bool> block_mask(this->n_blocks(), false);
  for (unsigned int c = 0; c < this->n_components();)
    {
      const unsigned int block = component_to_block_index(c);
      if (component_mask[c] == true)
        block_mask[block] = true;

      // now check all of the other
      // components that correspond
      // to this block
      ++c;
      while ((c < this->n_components()) &&
             (component_to_block_index(c) == block))
        {
          Assert(component_mask[c] == block_mask[block],
                 ExcMessage(
                   "The component mask argument given to this function "
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
FiniteElement<dim, spacedim>::face_to_cell_index(const unsigned int face_index,
                                                 const unsigned int face,
                                                 const bool face_orientation,
                                                 const bool face_flip,
                                                 const bool face_rotation) const
{
  AssertIndexRange(face_index, this->dofs_per_face);
  AssertIndexRange(face, GeometryInfo<dim>::faces_per_cell);

  // TODO: we could presumably solve the 3d case below using the
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
  if ((face_orientation != true) || (face_flip != false) ||
      (face_rotation != false))
    Assert((this->dofs_per_line <= 1) && (this->dofs_per_quad <= 1),
           ExcMessage(
             "The function in this base class can not handle this case. "
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
      const unsigned int face_vertex = face_index / this->dofs_per_vertex;
      const unsigned int dof_index_on_vertex =
        face_index % this->dofs_per_vertex;

      // then get the number of this vertex on the cell and translate
      // this to a DoF number on the cell
      return (GeometryInfo<dim>::face_to_cell_vertices(
                face, face_vertex, face_orientation, face_flip, face_rotation) *
                this->dofs_per_vertex +
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

      return (this->first_line_index +
              GeometryInfo<dim>::face_to_cell_lines(
                face, face_line, face_orientation, face_flip, face_rotation) *
                this->dofs_per_line +
              dof_index_on_line);
    }
  else
    // DoF is on a quad
    {
      Assert(dim >= 3, ExcInternalError());

      // ignore vertex and line dofs
      const unsigned int index = face_index - this->first_face_quad_index;

      return (this->first_quad_index + face * this->dofs_per_quad + index);
    }
}



template <int dim, int spacedim>
unsigned int
FiniteElement<dim, spacedim>::adjust_quad_dof_index_for_face_orientation(
  const unsigned int index,
  const bool         face_orientation,
  const bool         face_flip,
  const bool         face_rotation) const
{
  // general template for 1D and 2D: not
  // implemented. in fact, the function
  // shouldn't even be called unless we are
  // in 3d, so throw an internal error
  Assert(dim == 3, ExcInternalError());
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
  AssertIndexRange(index, this->dofs_per_quad);
  Assert(adjust_quad_dof_index_for_face_orientation_table.n_elements() ==
           8 * this->dofs_per_quad,
         ExcInternalError());
  return index + adjust_quad_dof_index_for_face_orientation_table(
                   index, 4 * face_orientation + 2 * face_flip + face_rotation);
}



template <int dim, int spacedim>
unsigned int
FiniteElement<dim, spacedim>::adjust_line_dof_index_for_line_orientation(
  const unsigned int index,
  const bool         line_orientation) const
{
  // general template for 1D and 2D: do
  // nothing. Do not throw an Assertion,
  // however, in order to allow to call this
  // function in 2D as well
  if (dim < 3)
    return index;

  AssertIndexRange(index, this->dofs_per_line);
  Assert(adjust_line_dof_index_for_line_orientation_table.size() ==
           this->dofs_per_line,
         ExcInternalError());
  if (line_orientation)
    return index;
  else
    return index + adjust_line_dof_index_for_line_orientation_table[index];
}



template <int dim, int spacedim>
bool
FiniteElement<dim, spacedim>::prolongation_is_implemented() const
{
  for (unsigned int ref_case = RefinementCase<dim>::cut_x;
       ref_case < RefinementCase<dim>::isotropic_refinement + 1;
       ++ref_case)
    for (unsigned int c = 0;
         c < GeometryInfo<dim>::n_children(RefinementCase<dim>(ref_case));
         ++c)
      {
        // make sure also the lazily initialized matrices are created
        get_prolongation_matrix(c, RefinementCase<dim>(ref_case));
        Assert((prolongation[ref_case - 1][c].m() == this->dofs_per_cell) ||
                 (prolongation[ref_case - 1][c].m() == 0),
               ExcInternalError());
        Assert((prolongation[ref_case - 1][c].n() == this->dofs_per_cell) ||
                 (prolongation[ref_case - 1][c].n() == 0),
               ExcInternalError());
        if ((prolongation[ref_case - 1][c].m() == 0) ||
            (prolongation[ref_case - 1][c].n() == 0))
          return false;
      }
  return true;
}



template <int dim, int spacedim>
bool
FiniteElement<dim, spacedim>::restriction_is_implemented() const
{
  for (unsigned int ref_case = RefinementCase<dim>::cut_x;
       ref_case < RefinementCase<dim>::isotropic_refinement + 1;
       ++ref_case)
    for (unsigned int c = 0;
         c < GeometryInfo<dim>::n_children(RefinementCase<dim>(ref_case));
         ++c)
      {
        // make sure also the lazily initialized matrices are created
        get_restriction_matrix(c, RefinementCase<dim>(ref_case));
        Assert((restriction[ref_case - 1][c].m() == this->dofs_per_cell) ||
                 (restriction[ref_case - 1][c].m() == 0),
               ExcInternalError());
        Assert((restriction[ref_case - 1][c].n() == this->dofs_per_cell) ||
                 (restriction[ref_case - 1][c].n() == 0),
               ExcInternalError());
        if ((restriction[ref_case - 1][c].m() == 0) ||
            (restriction[ref_case - 1][c].n() == 0))
          return false;
      }
  return true;
}



template <int dim, int spacedim>
bool
FiniteElement<dim, spacedim>::isotropic_prolongation_is_implemented() const
{
  const RefinementCase<dim> ref_case =
    RefinementCase<dim>::isotropic_refinement;

  for (unsigned int c = 0;
       c < GeometryInfo<dim>::n_children(RefinementCase<dim>(ref_case));
       ++c)
    {
      // make sure also the lazily initialized matrices are created
      get_prolongation_matrix(c, RefinementCase<dim>(ref_case));
      Assert((prolongation[ref_case - 1][c].m() == this->dofs_per_cell) ||
               (prolongation[ref_case - 1][c].m() == 0),
             ExcInternalError());
      Assert((prolongation[ref_case - 1][c].n() == this->dofs_per_cell) ||
               (prolongation[ref_case - 1][c].n() == 0),
             ExcInternalError());
      if ((prolongation[ref_case - 1][c].m() == 0) ||
          (prolongation[ref_case - 1][c].n() == 0))
        return false;
    }
  return true;
}



template <int dim, int spacedim>
bool
FiniteElement<dim, spacedim>::isotropic_restriction_is_implemented() const
{
  const RefinementCase<dim> ref_case =
    RefinementCase<dim>::isotropic_refinement;

  for (unsigned int c = 0;
       c < GeometryInfo<dim>::n_children(RefinementCase<dim>(ref_case));
       ++c)
    {
      // make sure also the lazily initialized matrices are created
      get_restriction_matrix(c, RefinementCase<dim>(ref_case));
      Assert((restriction[ref_case - 1][c].m() == this->dofs_per_cell) ||
               (restriction[ref_case - 1][c].m() == 0),
             ExcInternalError());
      Assert((restriction[ref_case - 1][c].n() == this->dofs_per_cell) ||
               (restriction[ref_case - 1][c].n() == 0),
             ExcInternalError());
      if ((restriction[ref_case - 1][c].m() == 0) ||
          (restriction[ref_case - 1][c].n() == 0))
        return false;
    }
  return true;
}



template <int dim, int spacedim>
bool
FiniteElement<dim, spacedim>::constraints_are_implemented(
  const internal::SubfaceCase<dim> &subface_case) const
{
  if (subface_case == internal::SubfaceCase<dim>::case_isotropic)
    return (this->dofs_per_face == 0) || (interface_constraints.m() != 0);
  else
    return false;
}



template <int dim, int spacedim>
bool
FiniteElement<dim, spacedim>::hp_constraints_are_implemented() const
{
  return false;
}



template <int dim, int spacedim>
const FullMatrix<double> &
FiniteElement<dim, spacedim>::constraints(
  const internal::SubfaceCase<dim> &subface_case) const
{
  (void)subface_case;
  Assert(subface_case == internal::SubfaceCase<dim>::case_isotropic,
         ExcMessage("Constraints for this element are only implemented "
                    "for the case that faces are refined isotropically "
                    "(which is always the case in 2d, and in 3d requires "
                    "that the neighboring cell of a coarse cell presents "
                    "exactly four children on the common face)."));
  Assert((this->dofs_per_face == 0) || (interface_constraints.m() != 0),
         ExcMessage("The finite element for which you try to obtain "
                    "hanging node constraints does not appear to "
                    "implement them."));

  if (dim == 1)
    Assert((interface_constraints.m() == 0) && (interface_constraints.n() == 0),
           ExcWrongInterfaceMatrixSize(interface_constraints.m(),
                                       interface_constraints.n()));

  return interface_constraints;
}



template <int dim, int spacedim>
TableIndices<2>
FiniteElement<dim, spacedim>::interface_constraints_size() const
{
  switch (dim)
    {
      case 1:
        return {0U, 0U};
      case 2:
        return {this->dofs_per_vertex + 2 * this->dofs_per_line,
                this->dofs_per_face};
      case 3:
        return {5 * this->dofs_per_vertex + 12 * this->dofs_per_line +
                  4 * this->dofs_per_quad,
                this->dofs_per_face};
      default:
        Assert(false, ExcNotImplemented());
    }
  return {numbers::invalid_unsigned_int, numbers::invalid_unsigned_int};
}



template <int dim, int spacedim>
void
FiniteElement<dim, spacedim>::get_interpolation_matrix(
  const FiniteElement<dim, spacedim> &,
  FullMatrix<double> &) const
{
  // by default, no interpolation
  // implemented. so throw exception,
  // as documentation says
  AssertThrow(
    false,
    (typename FiniteElement<dim, spacedim>::ExcInterpolationNotImplemented()));
}



template <int dim, int spacedim>
void
FiniteElement<dim, spacedim>::get_face_interpolation_matrix(
  const FiniteElement<dim, spacedim> &,
  FullMatrix<double> &) const
{
  // by default, no interpolation
  // implemented. so throw exception,
  // as documentation says
  AssertThrow(
    false,
    (typename FiniteElement<dim, spacedim>::ExcInterpolationNotImplemented()));
}



template <int dim, int spacedim>
void
FiniteElement<dim, spacedim>::get_subface_interpolation_matrix(
  const FiniteElement<dim, spacedim> &,
  const unsigned int,
  FullMatrix<double> &) const
{
  // by default, no interpolation
  // implemented. so throw exception,
  // as documentation says
  AssertThrow(
    false,
    (typename FiniteElement<dim, spacedim>::ExcInterpolationNotImplemented()));
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FiniteElement<dim, spacedim>::hp_vertex_dof_identities(
  const FiniteElement<dim, spacedim> &) const
{
  Assert(false, ExcNotImplemented());
  return std::vector<std::pair<unsigned int, unsigned int>>();
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FiniteElement<dim, spacedim>::hp_line_dof_identities(
  const FiniteElement<dim, spacedim> &) const
{
  Assert(false, ExcNotImplemented());
  return std::vector<std::pair<unsigned int, unsigned int>>();
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FiniteElement<dim, spacedim>::hp_quad_dof_identities(
  const FiniteElement<dim, spacedim> &) const
{
  Assert(false, ExcNotImplemented());
  return std::vector<std::pair<unsigned int, unsigned int>>();
}



template <int dim, int spacedim>
FiniteElementDomination::Domination
FiniteElement<dim, spacedim>::compare_for_face_domination(
  const FiniteElement<dim, spacedim> &fe_other) const
{
  return this->compare_for_domination(fe_other, 1);
}



template <int dim, int spacedim>
FiniteElementDomination::Domination
FiniteElement<dim, spacedim>::compare_for_domination(
  const FiniteElement<dim, spacedim> &,
  const unsigned int) const
{
  Assert(false, ExcNotImplemented());
  return FiniteElementDomination::neither_element_dominates;
}



template <int dim, int spacedim>
bool
FiniteElement<dim, spacedim>::
operator==(const FiniteElement<dim, spacedim> &f) const
{
  // Compare fields in roughly increasing order of how expensive the
  // comparison is
  return ((typeid(*this) == typeid(f)) && (this->get_name() == f.get_name()) &&
          (static_cast<const FiniteElementData<dim> &>(*this) ==
           static_cast<const FiniteElementData<dim> &>(f)) &&
          (interface_constraints == f.interface_constraints));
}



template <int dim, int spacedim>
bool
FiniteElement<dim, spacedim>::
operator!=(const FiniteElement<dim, spacedim> &f) const
{
  return !(*this == f);
}



template <int dim, int spacedim>
const std::vector<Point<dim>> &
FiniteElement<dim, spacedim>::get_unit_support_points() const
{
  // a finite element may define
  // support points, but only if
  // there are as many as there are
  // degrees of freedom
  Assert((unit_support_points.size() == 0) ||
           (unit_support_points.size() == this->dofs_per_cell),
         ExcInternalError());
  return unit_support_points;
}



template <int dim, int spacedim>
bool
FiniteElement<dim, spacedim>::has_support_points() const
{
  return (unit_support_points.size() != 0);
}



template <int dim, int spacedim>
const std::vector<Point<dim>> &
FiniteElement<dim, spacedim>::get_generalized_support_points() const
{
  // If the finite element implements generalized support points, return
  // those. Otherwise fall back to unit support points.
  return ((generalized_support_points.size() == 0) ?
            unit_support_points :
            generalized_support_points);
}



template <int dim, int spacedim>
bool
FiniteElement<dim, spacedim>::has_generalized_support_points() const
{
  return (get_generalized_support_points().size() != 0);
}



template <int dim, int spacedim>
Point<dim>
FiniteElement<dim, spacedim>::unit_support_point(const unsigned int index) const
{
  AssertIndexRange(index, this->dofs_per_cell);
  Assert(unit_support_points.size() == this->dofs_per_cell,
         ExcFEHasNoSupportPoints());
  return unit_support_points[index];
}



template <int dim, int spacedim>
const std::vector<Point<dim - 1>> &
FiniteElement<dim, spacedim>::get_unit_face_support_points() const
{
  // a finite element may define
  // support points, but only if
  // there are as many as there are
  // degrees of freedom on a face
  Assert((unit_face_support_points.size() == 0) ||
           (unit_face_support_points.size() == this->dofs_per_face),
         ExcInternalError());
  return unit_face_support_points;
}



template <int dim, int spacedim>
bool
FiniteElement<dim, spacedim>::has_face_support_points() const
{
  return (unit_face_support_points.size() != 0);
}



template <int dim, int spacedim>
Point<dim - 1>
FiniteElement<dim, spacedim>::unit_face_support_point(
  const unsigned int index) const
{
  AssertIndexRange(index, this->dofs_per_face);
  Assert(unit_face_support_points.size() == this->dofs_per_face,
         ExcFEHasNoSupportPoints());
  return unit_face_support_points[index];
}



template <int dim, int spacedim>
bool
FiniteElement<dim, spacedim>::has_support_on_face(const unsigned int,
                                                  const unsigned int) const
{
  return true;
}



template <int dim, int spacedim>
const FiniteElement<dim, spacedim> &
FiniteElement<dim, spacedim>::get_sub_fe(const ComponentMask &mask) const
{
  // Translate the ComponentMask into first_selected and n_components after
  // some error checking:
  const unsigned int n_total_components = this->n_components();
  Assert((n_total_components == mask.size()) || (mask.size() == 0),
         ExcMessage("The given ComponentMask has the wrong size."));

  const unsigned int n_selected =
    mask.n_selected_components(n_total_components);
  Assert(n_selected > 0,
         ExcMessage("You need at least one selected component."));

  const unsigned int first_selected =
    mask.first_selected_component(n_total_components);

#ifdef DEBUG
  // check that it is contiguous:
  for (unsigned int c = 0; c < n_total_components; ++c)
    Assert((c < first_selected && (!mask[c])) ||
             (c >= first_selected && c < first_selected + n_selected &&
              mask[c]) ||
             (c >= first_selected + n_selected && !mask[c]),
           ExcMessage("Error: the given ComponentMask is not contiguous!"));
#endif

  return get_sub_fe(first_selected, n_selected);
}



template <int dim, int spacedim>
const FiniteElement<dim, spacedim> &
FiniteElement<dim, spacedim>::get_sub_fe(
  const unsigned int first_component,
  const unsigned int n_selected_components) const
{
  // No complicated logic is needed here, because it is overridden in
  // FESystem<dim,spacedim>. Just make sure that what the user chose is valid:
  Assert(first_component == 0 && n_selected_components == this->n_components(),
         ExcMessage(
           "You can only select a whole FiniteElement, not a part of one."));

  (void)first_component;
  (void)n_selected_components;
  return *this;
}



template <int dim, int spacedim>
std::pair<Table<2, bool>, std::vector<unsigned int>>
FiniteElement<dim, spacedim>::get_constant_modes() const
{
  Assert(false, ExcNotImplemented());
  return std::pair<Table<2, bool>, std::vector<unsigned int>>(
    Table<2, bool>(this->n_components(), this->dofs_per_cell),
    std::vector<unsigned int>(this->n_components()));
}



template <int dim, int spacedim>
void
FiniteElement<dim, spacedim>::
  convert_generalized_support_point_values_to_dof_values(
    const std::vector<Vector<double>> &,
    std::vector<double> &) const
{
  Assert(has_generalized_support_points(),
         ExcMessage("The element for which you are calling the current "
                    "function does not have generalized support points (see "
                    "the glossary for a definition of generalized support "
                    "points). Consequently, the current function can not "
                    "be defined and is not implemented by the element."));
  Assert(false, ExcNotImplemented());
}



template <int dim, int spacedim>
std::size_t
FiniteElement<dim, spacedim>::memory_consumption() const
{
  return (
    sizeof(FiniteElementData<dim>) +
    MemoryConsumption::memory_consumption(restriction) +
    MemoryConsumption::memory_consumption(prolongation) +
    MemoryConsumption::memory_consumption(interface_constraints) +
    MemoryConsumption::memory_consumption(system_to_component_table) +
    MemoryConsumption::memory_consumption(face_system_to_component_table) +
    MemoryConsumption::memory_consumption(system_to_base_table) +
    MemoryConsumption::memory_consumption(face_system_to_base_table) +
    MemoryConsumption::memory_consumption(component_to_base_table) +
    MemoryConsumption::memory_consumption(restriction_is_additive_flags) +
    MemoryConsumption::memory_consumption(nonzero_components) +
    MemoryConsumption::memory_consumption(n_nonzero_components_table));
}



template <int dim, int spacedim>
std::vector<unsigned int>
FiniteElement<dim, spacedim>::compute_n_nonzero_components(
  const std::vector<ComponentMask> &nonzero_components)
{
  std::vector<unsigned int> retval(nonzero_components.size());
  for (unsigned int i = 0; i < nonzero_components.size(); ++i)
    retval[i] = nonzero_components[i].n_selected_components();
  return retval;
}



/*------------------------------- FiniteElement ----------------------*/

template <int dim, int spacedim>
std::unique_ptr<typename FiniteElement<dim, spacedim>::InternalDataBase>
FiniteElement<dim, spacedim>::get_face_data(
  const UpdateFlags             flags,
  const Mapping<dim, spacedim> &mapping,
  const Quadrature<dim - 1> &   quadrature,
  dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                     spacedim>
    &output_data) const
{
  return get_data(flags,
                  mapping,
                  QProjector<dim>::project_to_all_faces(quadrature),
                  output_data);
}



template <int dim, int spacedim>
std::unique_ptr<typename FiniteElement<dim, spacedim>::InternalDataBase>
FiniteElement<dim, spacedim>::get_subface_data(
  const UpdateFlags             flags,
  const Mapping<dim, spacedim> &mapping,
  const Quadrature<dim - 1> &   quadrature,
  dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                     spacedim>
    &output_data) const
{
  return get_data(flags,
                  mapping,
                  QProjector<dim>::project_to_all_subfaces(quadrature),
                  output_data);
}



template <int dim, int spacedim>
const FiniteElement<dim, spacedim> &
FiniteElement<dim, spacedim>::base_element(const unsigned int index) const
{
  (void)index;
  AssertIndexRange(index, 1);
  // This function should not be
  // called for a system element
  Assert(base_to_block_indices.size() == 1, ExcInternalError());
  return *this;
}



/*------------------------------- Explicit Instantiations -------------*/
#include "fe.inst"


DEAL_II_NAMESPACE_CLOSE
