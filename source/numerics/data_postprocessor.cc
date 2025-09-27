// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2007 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/numerics/data_postprocessor.h>

DEAL_II_NAMESPACE_OPEN



namespace DataPostprocessorInputs
{
  template <int spacedim>
  CommonInputs<spacedim>::CommonInputs()
    : face_number(numbers::invalid_unsigned_int)
  {}



  template <int spacedim>
  unsigned int
  CommonInputs<spacedim>::get_face_number() const
  {
    Assert(
      face_number != numbers::invalid_unsigned_int,
      ExcMessage(
        "This function can only be called if set_cell_and_face() has "
        "previously been called. Typically, this would be by using DataOutFaces "
        "or a related class."));
    return face_number;
  }
} // namespace DataPostprocessorInputs

// -------------------------- DataPostprocessor ---------------------------

template <int dim>
void
DataPostprocessor<dim>::evaluate_scalar_field(
  const DataPostprocessorInputs::Scalar<dim> &,
  std::vector<Vector<double>> &) const
{
  AssertThrow(false, ExcPureFunctionCalled());
}



template <int dim>
void
DataPostprocessor<dim>::evaluate_vector_field(
  const DataPostprocessorInputs::Vector<dim> &,
  std::vector<Vector<double>> &) const
{
  AssertThrow(false, ExcPureFunctionCalled());
}



template <int dim>
std::vector<DataComponentInterpretation::DataComponentInterpretation>
DataPostprocessor<dim>::get_data_component_interpretation() const
{
  // default implementation assumes that all
  // components are independent scalars
  return std::vector<DataComponentInterpretation::DataComponentInterpretation>(
    get_names().size(), DataComponentInterpretation::component_is_scalar);
}


// -------------------- DataPostprocessorScalar -------------------------

template <int dim>
DataPostprocessorScalar<dim>::DataPostprocessorScalar(
  const std::string &name,
  const UpdateFlags  update_flags)
  : name(name)
  , update_flags(update_flags)
{}



template <int dim>
std::vector<std::string>
DataPostprocessorScalar<dim>::get_names() const
{
  return std::vector<std::string>(1, name);
}



template <int dim>
std::vector<DataComponentInterpretation::DataComponentInterpretation>
DataPostprocessorScalar<dim>::get_data_component_interpretation() const
{
  return std::vector<DataComponentInterpretation::DataComponentInterpretation>(
    1, DataComponentInterpretation::component_is_scalar);
}


template <int dim>
UpdateFlags
DataPostprocessorScalar<dim>::get_needed_update_flags() const
{
  return update_flags;
}



// ------------------------- DataPostprocessorVector ------------------------

template <int dim>
DataPostprocessorVector<dim>::DataPostprocessorVector(
  const std::string &name,
  const UpdateFlags  update_flags)
  : name(name)
  , update_flags(update_flags)
{}



template <int dim>
std::vector<std::string>
DataPostprocessorVector<dim>::get_names() const
{
  return std::vector<std::string>(dim, name);
}



template <int dim>
std::vector<DataComponentInterpretation::DataComponentInterpretation>
DataPostprocessorVector<dim>::get_data_component_interpretation() const
{
  return std::vector<DataComponentInterpretation::DataComponentInterpretation>(
    dim, DataComponentInterpretation::component_is_part_of_vector);
}


template <int dim>
UpdateFlags
DataPostprocessorVector<dim>::get_needed_update_flags() const
{
  return update_flags;
}



// ------------------------- DataPostprocessorTensor ------------------------

template <int dim>
DataPostprocessorTensor<dim>::DataPostprocessorTensor(
  const std::string &name,
  const UpdateFlags  update_flags)
  : name(name)
  , update_flags(update_flags)
{}



template <int dim>
std::vector<std::string>
DataPostprocessorTensor<dim>::get_names() const
{
  return std::vector<std::string>(dim * dim, name);
}



template <int dim>
std::vector<DataComponentInterpretation::DataComponentInterpretation>
DataPostprocessorTensor<dim>::get_data_component_interpretation() const
{
  return std::vector<DataComponentInterpretation::DataComponentInterpretation>(
    dim * dim, DataComponentInterpretation::component_is_part_of_tensor);
}


template <int dim>
UpdateFlags
DataPostprocessorTensor<dim>::get_needed_update_flags() const
{
  return update_flags;
}



namespace DataPostprocessors
{
  template <int dim>
  BoundaryIds<dim>::BoundaryIds()
    : DataPostprocessorScalar<dim>("boundary_id", update_default)
  {}


  template <int dim>
  void
  BoundaryIds<dim>::evaluate_scalar_field(
    const DataPostprocessorInputs::Scalar<dim> &inputs,
    std::vector<Vector<double>>                &computed_quantities) const
  {
    AssertDimension(computed_quantities.size(), inputs.solution_values.size());

    const typename DoFHandler<dim>::active_cell_iterator cell =
      inputs.template get_cell<dim>();
    const unsigned int face = inputs.get_face_number();

    for (auto &output : computed_quantities)
      {
        AssertDimension(output.size(), 1);

        // By default, DataOutFaces is only run on faces at the boundary of the
        // domain. But one can instruct it to also run on internal faces, and in
        // that case we cannot ask for the boundary id. Rather, we output -1, as
        // described in the documentation.
        if (cell->at_boundary(face))
          output(0) = cell->face(face)->boundary_id();
        else
          output(0) = -1;
      }
  }
} // namespace DataPostprocessors



// explicit instantiation
#include "numerics/data_postprocessor.inst"


DEAL_II_NAMESPACE_CLOSE
