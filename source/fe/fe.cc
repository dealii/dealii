// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1998 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/quadrature.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_orientation.h>

#include <algorithm>
#include <functional>
#include <numeric>
#include <typeinfo>

DEAL_II_NAMESPACE_OPEN


/*------------------------------- FiniteElement ----------------------*/
#ifndef DOXYGEN

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
  const FiniteElementData<dim>     &fe_data,
  const std::vector<bool>          &r_i_a_f,
  const std::vector<ComponentMask> &nonzero_c)
  : FiniteElementData<dim>(fe_data)
  , adjust_line_dof_index_for_line_orientation_table(this->n_dofs_per_line())
  , system_to_base_table(this->n_dofs_per_cell())
  , component_to_base_table(this->components,
                            std::make_pair(std::make_pair(0U, 0U), 0U))
  ,

  // Special handling of vectors of length one: in this case, we
  // assume that all entries were supposed to be equal
  restriction_is_additive_flags(
    r_i_a_f.size() == 1 ?
      std::vector<bool>(fe_data.n_dofs_per_cell(), r_i_a_f[0]) :
      r_i_a_f)
  , nonzero_components(
      nonzero_c.size() == 1 ?
        std::vector<ComponentMask>(fe_data.n_dofs_per_cell(), nonzero_c[0]) :
        nonzero_c)
  , n_nonzero_components_table(compute_n_nonzero_components(nonzero_components))
  , cached_primitivity(std::find_if(n_nonzero_components_table.begin(),
                                    n_nonzero_components_table.end(),
                                    [](const unsigned int n_components) {
                                      return n_components != 1U;
                                    }) == n_nonzero_components_table.end())
{
  Assert(restriction_is_additive_flags.size() == this->n_dofs_per_cell(),
         ExcDimensionMismatch(restriction_is_additive_flags.size(),
                              this->n_dofs_per_cell()));
  AssertDimension(nonzero_components.size(), this->n_dofs_per_cell());
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
      system_to_component_table.resize(this->n_dofs_per_cell());
      for (unsigned int j = 0; j < this->n_dofs_per_cell(); ++j)
        system_to_component_table[j] = std::pair<unsigned, unsigned>(0, j);

      face_system_to_component_table.resize(this->n_unique_faces());
      for (unsigned int f = 0; f < this->n_unique_faces(); ++f)
        {
          face_system_to_component_table[f].resize(this->n_dofs_per_face(f));
          for (unsigned int j = 0; j < this->n_dofs_per_face(f); ++j)
            face_system_to_component_table[f][j] =
              std::pair<unsigned, unsigned>(0, j);
        }
    }

  for (unsigned int j = 0; j < this->n_dofs_per_cell(); ++j)
    system_to_base_table[j] = std::make_pair(std::make_pair(0U, 0U), j);

  face_system_to_base_table.resize(this->n_unique_faces());
  for (unsigned int f = 0; f < this->n_unique_faces(); ++f)
    {
      face_system_to_base_table[f].resize(this->n_dofs_per_face(f));
      for (unsigned int j = 0; j < this->n_dofs_per_face(f); ++j)
        face_system_to_base_table[f][j] =
          std::make_pair(std::make_pair(0U, 0U), j);
    }

  // Fill with default value; may be changed by constructor of derived class.
  base_to_block_indices.reinit(1, 1);

  // initialize the restriction and prolongation matrices. the default
  // constructor of FullMatrix<dim> initializes them with size zero
  prolongation.resize(RefinementCase<dim>::isotropic_refinement);
  restriction.resize(RefinementCase<dim>::isotropic_refinement);
  for (const unsigned int ref_case :
       RefinementCase<dim>::all_refinement_cases())
    if (ref_case != RefinementCase<dim>::no_refinement)
      {
        prolongation[ref_case - 1].resize(this->reference_cell().n_children(
                                            RefinementCase<dim>(ref_case)),
                                          FullMatrix<double>());
        restriction[ref_case - 1].resize(this->reference_cell().n_children(
                                           RefinementCase<dim>(ref_case)),
                                         FullMatrix<double>());
      }


  if (dim == 3)
    {
      adjust_quad_dof_index_for_face_orientation_table.resize(
        this->n_unique_2d_subobjects());

      for (unsigned int f = 0; f < this->n_unique_2d_subobjects(); ++f)
        {
          adjust_quad_dof_index_for_face_orientation_table[f] =
            Table<2, int>(this->n_dofs_per_quad(f),
                          this->reference_cell().n_face_orientations(f));
          adjust_quad_dof_index_for_face_orientation_table[f].fill(0);
        }
    }

  unit_face_support_points.resize(this->n_unique_faces());
  generalized_face_support_points.resize(this->n_unique_faces());
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
  for (const unsigned int ref_case :
       RefinementCase<dim>::all_refinement_cases())
    if (ref_case != RefinementCase<dim>::no_refinement)
      {
        const unsigned int nc =
          this->reference_cell().n_children(RefinementCase<dim>(ref_case));

        for (unsigned int i = 0; i < nc; ++i)
          {
            if (this->restriction[ref_case - 1][i].m() !=
                  this->n_dofs_per_cell() &&
                (!isotropic_restriction_only ||
                 ref_case == RefinementCase<dim>::isotropic_refinement))
              this->restriction[ref_case - 1][i].reinit(
                this->n_dofs_per_cell(), this->n_dofs_per_cell());
            if (this->prolongation[ref_case - 1][i].m() !=
                  this->n_dofs_per_cell() &&
                (!isotropic_prolongation_only ||
                 ref_case == RefinementCase<dim>::isotropic_refinement))
              this->prolongation[ref_case - 1][i].reinit(
                this->n_dofs_per_cell(), this->n_dofs_per_cell());
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
  AssertIndexRange(child,
                   this->reference_cell().n_children(
                     RefinementCase<dim>(refinement_case)));
  // we use refinement_case-1 here. the -1 takes care of the origin of the
  // vector, as for RefinementCase<dim>::no_refinement (=0) there is no data
  // available and so the vector indices are shifted
  Assert(restriction[refinement_case - 1][child].n() == this->n_dofs_per_cell(),
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
  AssertIndexRange(child,
                   this->reference_cell().n_children(
                     RefinementCase<dim>(refinement_case)));
  // we use refinement_case-1 here. the -1 takes care
  // of the origin of the vector, as for
  // RefinementCase::no_refinement (=0) there is no
  // data available and so the vector indices
  // are shifted
  Assert(prolongation[refinement_case - 1][child].n() ==
           this->n_dofs_per_cell(),
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
  return ComponentMask(mask);
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
  return ComponentMask(mask);
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
  return ComponentMask(mask);
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

  return ComponentMask(component_mask);
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


  return BlockMask(block_mask);
}



template <int dim, int spacedim>
unsigned int
FiniteElement<dim, spacedim>::face_to_cell_index(
  const unsigned int                 face_index,
  const unsigned int                 face,
  const types::geometric_orientation combined_orientation) const
{
  AssertIndexRange(face_index, this->n_dofs_per_face(face));
  AssertIndexRange(face, this->reference_cell().n_faces());

  // see the function's documentation for an explanation of this
  // assertion -- in essence, derived classes have to implement
  // an overloaded version of this function if we are to use any
  // other than default (standard) orientation
  if (combined_orientation != numbers::default_geometric_orientation)
    Assert((this->n_dofs_per_line() <= 1) && (this->n_dofs_per_quad(face) <= 1),
           ExcMessage(
             "The function in this base class can not handle this case. "
             "Rather, the derived class you are using must provide "
             "an overloaded version but apparently hasn't done so. See "
             "the documentation of this function for more information."));

  // we need to distinguish between DoFs on vertices, lines and (in 3d) quads.
  // do so in a sequence of if-else statements
  if (face_index < this->get_first_face_line_index(face))
    // DoF is on a vertex
    {
      // get the number of the vertex on the face that corresponds to this DoF,
      // along with the number of the DoF on this vertex
      const unsigned int face_vertex = face_index / this->n_dofs_per_vertex();
      const unsigned int dof_index_on_vertex =
        face_index % this->n_dofs_per_vertex();

      // then get the number of this vertex on the cell and translate
      // this to a DoF number on the cell
      return (this->reference_cell().face_to_cell_vertices(
                face, face_vertex, combined_orientation) *
                this->n_dofs_per_vertex() +
              dof_index_on_vertex);
    }
  else if (face_index < this->get_first_face_quad_index(face))
    // DoF is on a face
    {
      // do the same kind of translation as before. we need to only consider
      // DoFs on the lines, i.e., ignoring those on the vertices
      const unsigned int index =
        face_index - this->get_first_face_line_index(face);

      const unsigned int face_line         = index / this->n_dofs_per_line();
      const unsigned int dof_index_on_line = index % this->n_dofs_per_line();

      return (this->get_first_line_index() +
              this->reference_cell().face_to_cell_lines(face,
                                                        face_line,
                                                        combined_orientation) *
                this->n_dofs_per_line() +
              dof_index_on_line);
    }
  else
    // DoF is on a quad
    {
      Assert(dim >= 3, ExcInternalError());

      // ignore vertex and line dofs
      const unsigned int index =
        face_index - this->get_first_face_quad_index(face);

      return (this->get_first_quad_index(face) + index);
    }
}



template <int dim, int spacedim>
unsigned int
FiniteElement<dim, spacedim>::adjust_quad_dof_index_for_face_orientation(
  const unsigned int                 index,
  const unsigned int                 face,
  const types::geometric_orientation combined_orientation) const
{
  // general template for 1d and 2d: not
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
  AssertIndexRange(index, this->n_dofs_per_quad(face));
  const auto table_n = this->n_unique_2d_subobjects() == 1 ? 0 : face;
  Assert(
    adjust_quad_dof_index_for_face_orientation_table[table_n].n_elements() ==
      (this->reference_cell().n_face_orientations(face)) *
        this->n_dofs_per_quad(face),
    ExcInternalError());
  return index + adjust_quad_dof_index_for_face_orientation_table[table_n](
                   index, combined_orientation);
}



template <int dim, int spacedim>
unsigned int
FiniteElement<dim, spacedim>::adjust_line_dof_index_for_line_orientation(
  const unsigned int                 index,
  const types::geometric_orientation combined_orientation) const
{
  Assert(combined_orientation == numbers::default_geometric_orientation ||
           combined_orientation == numbers::reverse_line_orientation,
         ExcInternalError());

  AssertIndexRange(index, this->n_dofs_per_line());
  Assert(adjust_line_dof_index_for_line_orientation_table.size() ==
           this->n_dofs_per_line(),
         ExcInternalError());
  if (combined_orientation == numbers::default_geometric_orientation)
    return index;
  else
    return index + adjust_line_dof_index_for_line_orientation_table[index];
}



template <int dim, int spacedim>
bool
FiniteElement<dim, spacedim>::prolongation_is_implemented() const
{
  for (const unsigned int ref_case :
       RefinementCase<dim>::all_refinement_cases())
    if (ref_case != RefinementCase<dim>::no_refinement)
      for (unsigned int c = 0;
           c < this->reference_cell().n_children(RefinementCase<dim>(ref_case));
           ++c)
        {
          // make sure also the lazily initialized matrices are created
          get_prolongation_matrix(c, RefinementCase<dim>(ref_case));
          Assert((prolongation[ref_case - 1][c].m() ==
                  this->n_dofs_per_cell()) ||
                   (prolongation[ref_case - 1][c].m() == 0),
                 ExcInternalError());
          Assert((prolongation[ref_case - 1][c].n() ==
                  this->n_dofs_per_cell()) ||
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
  for (const unsigned int ref_case :
       RefinementCase<dim>::all_refinement_cases())
    if (ref_case != RefinementCase<dim>::no_refinement)
      for (unsigned int c = 0;
           c < this->reference_cell().n_children(RefinementCase<dim>(ref_case));
           ++c)
        {
          // make sure also the lazily initialized matrices are created
          get_restriction_matrix(c, RefinementCase<dim>(ref_case));
          Assert((restriction[ref_case - 1][c].m() ==
                  this->n_dofs_per_cell()) ||
                   (restriction[ref_case - 1][c].m() == 0),
                 ExcInternalError());
          Assert((restriction[ref_case - 1][c].n() ==
                  this->n_dofs_per_cell()) ||
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
       c < this->reference_cell().n_children(RefinementCase<dim>(ref_case));
       ++c)
    {
      // make sure also the lazily initialized matrices are created
      get_prolongation_matrix(c, RefinementCase<dim>(ref_case));
      Assert((prolongation[ref_case - 1][c].m() == this->n_dofs_per_cell()) ||
               (prolongation[ref_case - 1][c].m() == 0),
             ExcInternalError());
      Assert((prolongation[ref_case - 1][c].n() == this->n_dofs_per_cell()) ||
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
       c < this->reference_cell().n_children(RefinementCase<dim>(ref_case));
       ++c)
    {
      // make sure also the lazily initialized matrices are created
      get_restriction_matrix(c, RefinementCase<dim>(ref_case));
      Assert((restriction[ref_case - 1][c].m() == this->n_dofs_per_cell()) ||
               (restriction[ref_case - 1][c].m() == 0),
             ExcInternalError());
      Assert((restriction[ref_case - 1][c].n() == this->n_dofs_per_cell()) ||
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
    {
      unsigned int n_dofs_on_faces = 0;

      for (const auto face_no : this->reference_cell().face_indices())
        n_dofs_on_faces += this->n_dofs_per_face(face_no);

      return (n_dofs_on_faces == 0) || (interface_constraints.m() != 0);
    }
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
  // TODO: the implementation makes the assumption that all faces have the
  // same number of dofs
  AssertDimension(this->n_unique_faces(), 1);
  [[maybe_unused]] const unsigned int face_no = 0;

  Assert(subface_case == internal::SubfaceCase<dim>::case_isotropic,
         ExcMessage("Constraints for this element are only implemented "
                    "for the case that faces are refined isotropically "
                    "(which is always the case in 2d, and in 3d requires "
                    "that the neighboring cell of a coarse cell presents "
                    "exactly four children on the common face)."));
  Assert((this->n_dofs_per_face(face_no) == 0) ||
           (interface_constraints.m() != 0),
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
  // TODO: the implementation makes the assumption that all faces have the
  // same number of dofs
  AssertDimension(this->n_unique_faces(), 1);
  const unsigned int face_no = 0;

  switch (dim)
    {
      case 1:
        return {0U, 0U};

      case 2:
        // We have to interpolate from the DoFs in the interior of the
        // the two child faces (=lines) and the one central vertex
        // to the DoFs of the parent face:
        return {this->n_dofs_per_vertex() + 2 * this->n_dofs_per_line(),
                this->n_dofs_per_face(face_no)};

      case 3:
        // We have to interpolate from the DoFs in the interior of the
        // the child faces (=quads or tris) and the vertices that are
        // not part of the parent face, to the DoFs of the parent face:
        if (this->reference_cell().face_reference_cell(face_no) ==
            ReferenceCells::Quadrilateral)
          return {
            5 * this->n_dofs_per_vertex() +  // 4 vertices at mid-edge points
                                             // + 1 at cell center
              12 * this->n_dofs_per_line() + // 4*2 children of the old edges
                                             // + 2*2 edges in the cell interior
              4 * this->n_dofs_per_quad(face_no), // 4 child faces
            this->n_dofs_per_face(face_no)};
        else if (this->reference_cell().face_reference_cell(face_no) ==
                 ReferenceCells::Triangle)
          return {
            3 * this->n_dofs_per_vertex() + // 3 vertices at mid-edge points
              9 * this->n_dofs_per_line() + // 3*2 children of the old edges
                                            // + 3 edges in the cell interior
              4 * this->n_dofs_per_quad(face_no), // 4 child faces
            this->n_dofs_per_face(face_no)};
        else
          DEAL_II_ASSERT_UNREACHABLE();

      default:
        DEAL_II_NOT_IMPLEMENTED();
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
  FullMatrix<double> &,
  const unsigned int) const
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
  FullMatrix<double> &,
  const unsigned int) const
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
  DEAL_II_NOT_IMPLEMENTED();
  return std::vector<std::pair<unsigned int, unsigned int>>();
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FiniteElement<dim, spacedim>::hp_line_dof_identities(
  const FiniteElement<dim, spacedim> &) const
{
  DEAL_II_NOT_IMPLEMENTED();
  return std::vector<std::pair<unsigned int, unsigned int>>();
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FiniteElement<dim, spacedim>::hp_quad_dof_identities(
  const FiniteElement<dim, spacedim> &,
  const unsigned int) const
{
  DEAL_II_NOT_IMPLEMENTED();
  return std::vector<std::pair<unsigned int, unsigned int>>();
}



template <int dim, int spacedim>
FiniteElementDomination::Domination
FiniteElement<dim, spacedim>::compare_for_domination(
  const FiniteElement<dim, spacedim> &,
  const unsigned int) const
{
  DEAL_II_NOT_IMPLEMENTED();
  return FiniteElementDomination::neither_element_dominates;
}



template <int dim, int spacedim>
bool
FiniteElement<dim, spacedim>::operator==(
  const FiniteElement<dim, spacedim> &f) const
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
FiniteElement<dim, spacedim>::operator!=(
  const FiniteElement<dim, spacedim> &f) const
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
  Assert((unit_support_points.empty()) ||
           (unit_support_points.size() == this->n_dofs_per_cell()),
         ExcInternalError());
  return unit_support_points;
}



template <int dim, int spacedim>
bool
FiniteElement<dim, spacedim>::has_support_points() const
{
  if (this->dofs_per_cell > 0)
    return (unit_support_points.size() != 0);
  else
    {
      // If the FE has no DoFs, we shouldn't expect the array
      // size to be anything other than zero:
      AssertDimension(unit_support_points.size(), 0);

      // A finite element without DoFs *has* support points
      // (which is then an empty array)
      return true;
    }
}



template <int dim, int spacedim>
const std::vector<Point<dim>> &
FiniteElement<dim, spacedim>::get_generalized_support_points() const
{
  // If the finite element implements generalized support points, return
  // those. Otherwise fall back to unit support points.
  return ((generalized_support_points.empty()) ? unit_support_points :
                                                 generalized_support_points);
}



template <int dim, int spacedim>
bool
FiniteElement<dim, spacedim>::has_generalized_support_points() const
{
  if (this->dofs_per_cell > 0)
    return (get_generalized_support_points().size() != 0);
  else
    {
      // If the FE has no DoFs, the array size should be zero:
      AssertDimension(get_generalized_support_points().size(), 0);

      // A finite element without DoFs *has* generalized support points
      // (which is then an empty array)
      return true;
    }
}



template <int dim, int spacedim>
Point<dim>
FiniteElement<dim, spacedim>::unit_support_point(const unsigned int index) const
{
  AssertIndexRange(index, this->n_dofs_per_cell());
  Assert(unit_support_points.size() == this->n_dofs_per_cell(),
         ExcFEHasNoSupportPoints());
  return unit_support_points[index];
}



template <int dim, int spacedim>
const std::vector<Point<dim - 1>> &
FiniteElement<dim, spacedim>::get_unit_face_support_points(
  const unsigned int face_no) const
{
  // a finite element may define
  // support points, but only if
  // there are as many as there are
  // degrees of freedom on a face
  Assert((unit_face_support_points[this->n_unique_faces() == 1 ? 0 : face_no]
            .empty()) ||
           (unit_face_support_points[this->n_unique_faces() == 1 ? 0 : face_no]
              .size() == this->n_dofs_per_face(face_no)),
         ExcInternalError());
  return unit_face_support_points[this->n_unique_faces() == 1 ? 0 : face_no];
}



template <int dim, int spacedim>
bool
FiniteElement<dim, spacedim>::has_face_support_points(
  const unsigned int face_no) const
{
  const unsigned int face_index = this->n_unique_faces() == 1 ? 0 : face_no;
  if (this->n_dofs_per_face(face_index) > 0)
    return (unit_face_support_points[face_index].size() != 0);
  else
    {
      // If the FE has no DoFs on face, the array size should be zero
      AssertDimension(unit_face_support_points[face_index].size(), 0);

      // A finite element without DoFs *has* face support points
      // (which is then an empty array)
      return true;
    }
}



template <int dim, int spacedim>
Point<dim - 1>
FiniteElement<dim, spacedim>::unit_face_support_point(
  const unsigned int index,
  const unsigned int face_no) const
{
  AssertIndexRange(index, this->n_dofs_per_face(face_no));
  Assert(unit_face_support_points[this->n_unique_faces() == 1 ? 0 : face_no]
             .size() == this->n_dofs_per_face(face_no),
         ExcFEHasNoSupportPoints());
  return unit_face_support_points[this->n_unique_faces() == 1 ? 0 : face_no]
                                 [index];
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

  if constexpr (running_in_debug_mode())
    {
      // check that it is contiguous:
      for (unsigned int c = 0; c < n_total_components; ++c)
        Assert((c < first_selected && (!mask[c])) ||
                 (c >= first_selected && c < first_selected + n_selected &&
                  mask[c]) ||
                 (c >= first_selected + n_selected && !mask[c]),
               ExcMessage("Error: the given ComponentMask is not contiguous!"));
    }

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

  return *this;
}



template <int dim, int spacedim>
std::pair<Table<2, bool>, std::vector<unsigned int>>
FiniteElement<dim, spacedim>::get_constant_modes() const
{
  DEAL_II_NOT_IMPLEMENTED();
  return std::pair<Table<2, bool>, std::vector<unsigned int>>(
    Table<2, bool>(this->n_components(), this->n_dofs_per_cell()),
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
  DEAL_II_NOT_IMPLEMENTED();
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

#  ifndef DOXYGEN
template <int dim, int spacedim>
std::unique_ptr<typename FiniteElement<dim, spacedim>::InternalDataBase>
FiniteElement<dim, spacedim>::get_face_data(
  const UpdateFlags               flags,
  const Mapping<dim, spacedim>   &mapping,
  const hp::QCollection<dim - 1> &quadrature,
  dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                     spacedim>
    &output_data) const
{
  return get_data(flags,
                  mapping,
                  QProjector<dim>::project_to_all_faces(this->reference_cell(),
                                                        quadrature),
                  output_data);
}



template <int dim, int spacedim>
std::unique_ptr<typename FiniteElement<dim, spacedim>::InternalDataBase>
FiniteElement<dim, spacedim>::get_face_data(
  const UpdateFlags             flags,
  const Mapping<dim, spacedim> &mapping,
  const Quadrature<dim - 1>    &quadrature,
  dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                     spacedim>
    &output_data) const
{
  return get_data(flags,
                  mapping,
                  QProjector<dim>::project_to_all_faces(this->reference_cell(),
                                                        quadrature),
                  output_data);
}



template <int dim, int spacedim>
inline void
FiniteElement<dim, spacedim>::fill_fe_face_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const unsigned int                                          face_no,
  const hp::QCollection<dim - 1>                             &quadrature,
  const Mapping<dim, spacedim>                               &mapping,
  const typename Mapping<dim, spacedim>::InternalDataBase    &mapping_internal,
  const internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
                                                                &mapping_data,
  const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
  dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                     spacedim>
    &output_data) const
{
  // base class version, implement overridden function in derived classes
  AssertDimension(quadrature.size(), 1);
  fill_fe_face_values(cell,
                      face_no,
                      quadrature[0],
                      mapping,
                      mapping_internal,
                      mapping_data,
                      fe_internal,
                      output_data);
}



template <int dim, int spacedim>
inline void
FiniteElement<dim, spacedim>::fill_fe_face_values(
  const typename Triangulation<dim, spacedim>::cell_iterator & /* cell */,
  const unsigned int /* face_no */,
  const Quadrature<dim - 1> & /* quadrature */,
  const Mapping<dim, spacedim> & /* mapping */,
  const typename Mapping<dim, spacedim>::InternalDataBase
    & /* mapping_internal */,
  const internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
    & /* mapping_data */,
  const typename FiniteElement<dim, spacedim>::InternalDataBase
    & /* fe_internal */,
  dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                     spacedim>
    & /* output_data */) const
{
  Assert(false,
         ExcMessage("Use of a deprecated interface, please implement "
                    "fill_fe_face_values taking a hp::QCollection argument"));
}
#  endif



template <int dim, int spacedim>
std::unique_ptr<typename FiniteElement<dim, spacedim>::InternalDataBase>
FiniteElement<dim, spacedim>::get_subface_data(
  const UpdateFlags             flags,
  const Mapping<dim, spacedim> &mapping,
  const Quadrature<dim - 1>    &quadrature,
  dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                     spacedim>
    &output_data) const
{
  return get_data(flags,
                  mapping,
                  QProjector<dim>::project_to_all_subfaces(
                    this->reference_cell(), quadrature),
                  output_data);
}



template <int dim, int spacedim>
const FiniteElement<dim, spacedim> &
FiniteElement<dim, spacedim>::base_element(const unsigned int index) const
{
  AssertIndexRange(index, 1);
  // This function should not be
  // called for a system element
  Assert(base_to_block_indices.size() == 1, ExcInternalError());
  return *this;
}


#endif
/*------------------------------- Explicit Instantiations -------------*/
#include "fe/fe.inst"


DEAL_II_NAMESPACE_CLOSE
