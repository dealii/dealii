// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2018 by the deal.II authors
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

#include <deal.II/numerics/data_postprocessor.h>

DEAL_II_NAMESPACE_OPEN



// -------------------------- DataPostprocessor ---------------------------

template <int dim, typename Number>
void
DataPostprocessor<dim, Number>::evaluate_scalar_field(
  const DataPostprocessorInputs::Scalar<dim, Number> &,
  std::vector<Vector<Number>> &) const
{
  AssertThrow(false, ExcPureFunctionCalled());
}



template <int dim, typename Number>
void
DataPostprocessor<dim, Number>::evaluate_vector_field(
  const DataPostprocessorInputs::Vector<dim, Number> &,
  std::vector<Vector<Number>> &) const
{
  AssertThrow(false, ExcPureFunctionCalled());
}



template <int dim, typename Number>
std::vector<DataComponentInterpretation::DataComponentInterpretation>
DataPostprocessor<dim, Number>::get_data_component_interpretation() const
{
  // default implementation assumes that all
  // components are independent scalars
  return std::vector<DataComponentInterpretation::DataComponentInterpretation>(
    get_names().size(), DataComponentInterpretation::component_is_scalar);
}


// -------------------------- DataPostprocessorScalar -------------------------

template <int dim, typename Number>
DataPostprocessorScalar<dim, Number>::DataPostprocessorScalar(
  const std::string &name,
  const UpdateFlags  update_flags)
  : name(name)
  , update_flags(update_flags)
{}



template <int dim, typename Number>
std::vector<std::string>
DataPostprocessorScalar<dim, Number>::get_names() const
{
  return std::vector<std::string>(1, name);
}



template <int dim, typename Number>
std::vector<DataComponentInterpretation::DataComponentInterpretation>
DataPostprocessorScalar<dim, Number>::get_data_component_interpretation() const
{
  return std::vector<DataComponentInterpretation::DataComponentInterpretation>(
    1, DataComponentInterpretation::component_is_scalar);
}


template <int dim, typename Number>
UpdateFlags
DataPostprocessorScalar<dim, Number>::get_needed_update_flags() const
{
  return update_flags;
}



// -------------------------- DataPostprocessorVector -------------------------

template <int dim, typename Number>
DataPostprocessorVector<dim, Number>::DataPostprocessorVector(
  const std::string &name,
  const UpdateFlags  update_flags)
  : name(name)
  , update_flags(update_flags)
{}



template <int dim, typename Number>
std::vector<std::string>
DataPostprocessorVector<dim, Number>::get_names() const
{
  return std::vector<std::string>(dim, name);
}



template <int dim, typename Number>
std::vector<DataComponentInterpretation::DataComponentInterpretation>
DataPostprocessorVector<dim, Number>::get_data_component_interpretation() const
{
  return std::vector<DataComponentInterpretation::DataComponentInterpretation>(
    dim, DataComponentInterpretation::component_is_part_of_vector);
}


template <int dim, typename Number>
UpdateFlags
DataPostprocessorVector<dim, Number>::get_needed_update_flags() const
{
  return update_flags;
}



// -------------------------- DataPostprocessorTensor -------------------------

template <int dim, typename Number>
DataPostprocessorTensor<dim, Number>::DataPostprocessorTensor(
  const std::string &name,
  const UpdateFlags  update_flags)
  : name(name)
  , update_flags(update_flags)
{}



template <int dim, typename Number>
std::vector<std::string>
DataPostprocessorTensor<dim, Number>::get_names() const
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



template <int dim, typename Number>
std::vector<DataComponentInterpretation::DataComponentInterpretation>
DataPostprocessorTensor<dim, Number>::get_data_component_interpretation() const
{
  return std::vector<DataComponentInterpretation::DataComponentInterpretation>(
    dim * dim, DataComponentInterpretation::component_is_scalar);
}


template <int dim, typename Number>
UpdateFlags
DataPostprocessorTensor<dim, Number>::get_needed_update_flags() const
{
  return update_flags;
}



// explicit instantiation
#include "data_postprocessor.inst"


DEAL_II_NAMESPACE_CLOSE
