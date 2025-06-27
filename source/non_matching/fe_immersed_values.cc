// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/thread_management.h>

#include <deal.II/grid/tria_iterator.h>

#include <deal.II/non_matching/fe_immersed_values.h>


DEAL_II_NAMESPACE_OPEN
namespace NonMatching
{
  template <int dim>
  FEImmersedSurfaceValues<dim>::FEImmersedSurfaceValues(
    const Mapping<dim>                   &mapping,
    const FiniteElement<dim>             &element,
    const ImmersedSurfaceQuadrature<dim> &quadrature,
    const UpdateFlags                     update_flags)
    : FEValuesBase<dim, dim>(quadrature.size(),
                             element.dofs_per_cell,
                             update_default,
                             mapping,
                             element)
    , quadrature(quadrature)
  {
    initialize(update_flags);
  }



  template <int dim>
  void
  FEImmersedSurfaceValues<dim>::reinit(
    const typename Triangulation<dim>::cell_iterator &cell)
  {
    // Check that mapping and reference cell type are compatible:
    Assert(this->get_mapping().is_compatible_with(cell->reference_cell()),
           ExcMessage(
             "You are trying to call FEImmersedSurfaceValues::reinit() with "
             " a cell of type " +
             cell->reference_cell().to_string() +
             " with a Mapping that is not compatible with it."));

    // No FE in this cell, so no assertion necessary here.
    Assert(
      this->present_cell.is_initialized() == false,
      ExcMessage(
        "FEImmersedSurfaceValues::reinit() can only be used for one cell!"));

    this->present_cell = {cell};

    // This was the part of the work that is dependent on the actual data type
    // of the iterator. Now pass on to the function doing the real work.
    do_reinit();
  }



  template <int dim>
  template <bool level_dof_access>
  void
  FEImmersedSurfaceValues<dim>::reinit(
    const TriaIterator<DoFCellAccessor<dim, dim, level_dof_access>> &cell)
  {
    // Check that mapping and reference cell type are compatible:
    Assert(this->get_mapping().is_compatible_with(cell->reference_cell()),
           ExcMessage(
             "You are trying to call FEImmersedSurfaceValues::reinit() with "
             "a cell of type " +
             cell->reference_cell().to_string() +
             " with a Mapping that is not compatible with it."));

    // Assert that the finite elements passed to the constructor and used by the
    // DoFHandler used by this cell, are the same
    Assert(static_cast<const FiniteElementData<dim> &>(*this->fe) ==
             static_cast<const FiniteElementData<dim> &>(cell->get_fe()),
           (typename FEValuesBase<dim>::ExcFEDontMatch()));

    Assert(
      this->present_cell.is_initialized() == false,
      ExcMessage(
        "FEImmersedSurfaceValues::reinit() can only be used for one cell!"));

    this->present_cell = {cell};

    // This was the part of the work that is dependent on the actual data type
    // of the iterator. Now pass on to the function doing the real work.
    do_reinit();
  }



  template <int dim>
  void
  FEImmersedSurfaceValues<dim>::do_reinit()
  {
    // First call the mapping and let it generate the data specific to the
    // mapping.
    if (this->update_flags & update_mapping)
      {
        this->get_mapping().fill_fe_immersed_surface_values(
          this->present_cell,
          quadrature,
          *this->mapping_data,
          this->mapping_output);
      }

    // Call the finite element and, with the data already filled by the mapping,
    // let it compute the data for the mapped shape function values, gradients
    // etc.
    this->get_fe().fill_fe_values(this->present_cell,
                                  CellSimilarity::none,
                                  this->quadrature,
                                  this->get_mapping(),
                                  *this->mapping_data,
                                  this->mapping_output,
                                  *this->fe_data,
                                  this->finite_element_output);
  }



  template <int dim>
  Tensor<1, dim>
  FEImmersedSurfaceValues<dim>::shape_surface_grad(
    const unsigned int function_no,
    const unsigned int quadrature_point) const
  {
    const unsigned int component = 0;
    return shape_surface_grad_component(function_no,
                                        quadrature_point,
                                        component);
  }



  template <int dim>
  Tensor<1, dim>
  FEImmersedSurfaceValues<dim>::shape_surface_grad_component(
    const unsigned int function_no,
    const unsigned int quadrature_point,
    const unsigned int component) const
  {
    const Tensor<1, dim> gradient =
      this->shape_grad_component(function_no, quadrature_point, component);
    const Tensor<1, dim> &normal = this->normal_vector(quadrature_point);

    return gradient - (normal * gradient) * normal;
  }



  template <int dim>
  const ImmersedSurfaceQuadrature<dim> &
  FEImmersedSurfaceValues<dim>::get_quadrature() const
  {
    return quadrature;
  }



  template <int dim>
  inline void
  FEImmersedSurfaceValues<dim>::initialize(const UpdateFlags update_flags)
  {
    UpdateFlags flags = this->compute_update_flags(update_flags);

    if ((flags & (update_JxW_values | update_normal_vectors)) != 0u)
      flags |= update_covariant_transformation;

    // Initialize the base classes.
    if ((flags & update_mapping) != 0u)
      this->mapping_output.initialize(this->n_quadrature_points, flags);
    this->finite_element_output.initialize(this->n_quadrature_points,
                                           *this->fe,
                                           flags);

    // Then get objects into which the FE and the Mapping can store
    // intermediate data used across calls to reinit. We can do this in
    // parallel.
    Threads::Task<
      std::unique_ptr<typename FiniteElement<dim, dim>::InternalDataBase>>
      fe_get_data = Threads::new_task(&FiniteElement<dim, dim>::get_data,
                                      *this->fe,
                                      flags,
                                      *this->mapping,
                                      this->quadrature,
                                      this->finite_element_output);

    Threads::Task<std::unique_ptr<typename Mapping<dim>::InternalDataBase>>
      mapping_get_data;
    if ((flags & update_mapping) != 0u)
      mapping_get_data = Threads::new_task(&Mapping<dim>::get_data,
                                           *this->mapping,
                                           flags,
                                           this->quadrature);

    this->update_flags = flags;

    // Then collect answers from the two task above.
    this->fe_data = std::move(fe_get_data.return_value());
    if ((flags & update_mapping) != 0u)
      this->mapping_data = std::move(mapping_get_data.return_value());
    else
      this->mapping_data =
        std::make_unique<typename Mapping<dim>::InternalDataBase>();
  }


#ifndef DOXYGEN
#  include "non_matching/fe_immersed_values.inst"
#endif

} // namespace NonMatching
DEAL_II_NAMESPACE_CLOSE
