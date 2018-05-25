// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2017 by the deal.II authors
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

#include <deal.II/numerics/data_postprocessor.h>

DEAL_II_NAMESPACE_OPEN



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


// -------------------------- DataPostprocessorScalar -------------------------

template <int dim>
DataPostprocessorScalar<dim>::DataPostprocessorScalar(
  const std::string &name,
  const UpdateFlags  update_flags) :
  name(name),
  update_flags(update_flags)
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



// -------------------------- DataPostprocessorVector -------------------------

template <int dim>
DataPostprocessorVector<dim>::DataPostprocessorVector(
  const std::string &name,
  const UpdateFlags  update_flags) :
  name(name),
  update_flags(update_flags)
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



// -------------------------- DataPostprocessorTensor -------------------------

template <int dim>
DataPostprocessorTensor<dim>::DataPostprocessorTensor(
  const std::string &name,
  const UpdateFlags  update_flags) :
  name(name),
  update_flags(update_flags)
{}



template <int dim>
std::vector<std::string>
DataPostprocessorTensor<dim>::get_names() const
{
  static_assert(dim <= 3,
                "The following variable needs to be expanded for dim>3");
  static const char suffixes[] = {'x', 'y', 'z'};

  std::vector<std::string> names;
  for (unsigned int d = 0; d < dim; ++d)
    for (unsigned int e = 0; e < dim; ++e)
      names.push_back(name + '_' + suffixes[d] + suffixes[e]);
  return names;
}



template <int dim>
std::vector<DataComponentInterpretation::DataComponentInterpretation>
DataPostprocessorTensor<dim>::get_data_component_interpretation() const
{
  return std::vector<DataComponentInterpretation::DataComponentInterpretation>(
    dim * dim, DataComponentInterpretation::component_is_scalar);
}


template <int dim>
UpdateFlags
DataPostprocessorTensor<dim>::get_needed_update_flags() const
{
  return update_flags;
}



// explicit instantiation
#include "data_postprocessor.inst"


DEAL_II_NAMESPACE_CLOSE
