// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/numerics/data_out_points.h>


DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim>
void
DataOutPoints<dim, spacedim>::build_patches(
  const std::vector<Point<spacedim>>     &locations,
  const std::size_t                       offset,
  const std::vector<std::vector<double>> &data,
  const std::vector<std::string>         &data_component_names,
  const std::vector<DataComponentInterpretation::DataComponentInterpretation>
    &data_component_interpretations_)
{
  Assert(
    data_component_names.size() == data_component_interpretations_.size(),
    ExcMessage(
      "When calling PointDataOut::build_patches() with data component "
      "names and interpretations you need to provide as many data component "
      "names as interpretations. Provide the same name for components that "
      "belong to a single vector or tensor."));

  Assert(data.size() == 0 || locations.size(),
         ExcMessage(
           "You need to either provide no data or data for each point."));
  for (const auto &datum : data)
    Assert(datum.size() == data_component_names.size(),
           ExcMessage(
             "The data provided in each point needs to have the same number "
             "of components as names were provided."));

  // Prepend the "id" to the data fields provided by the user:
  dataset_names.clear();
  dataset_names.emplace_back("id");
  dataset_names.insert(dataset_names.end(),
                       data_component_names.begin(),
                       data_component_names.end());

  data_component_interpretations.clear();
  data_component_interpretations.emplace_back(
    DataComponentInterpretation::component_is_scalar);
  data_component_interpretations.insert(data_component_interpretations.end(),
                                        data_component_interpretations_.begin(),
                                        data_component_interpretations_.end());

  const unsigned int n_property_components = data_component_names.size();
  const unsigned int n_data_components     = dataset_names.size();

  patches.resize(locations.size());

  for (unsigned int i = 0; i < locations.size(); ++i)
    {
      patches[i].vertices[0] = locations[i];
      patches[i].patch_index = i;

      // Store id and properties given by the user:
      patches[i].data.reinit(n_data_components, 1);
      patches[i].data(0, 0) = static_cast<float>(offset + i);
      for (unsigned int property_index = 0;
           property_index < n_property_components;
           ++property_index)
        patches[i].data(property_index + 1, 0) = data[i][property_index];
    }
}



template <int dim, int spacedim>
const std::vector<DataOutBase::Patch<0, spacedim>> &
DataOutPoints<dim, spacedim>::get_patches() const
{
  return patches;
}



template <int dim, int spacedim>
std::vector<std::string>
DataOutPoints<dim, spacedim>::get_dataset_names() const
{
  return dataset_names;
}



template <int dim, int spacedim>
std::vector<
  std::tuple<unsigned int,
             unsigned int,
             std::string,
             DataComponentInterpretation::DataComponentInterpretation>>
DataOutPoints<dim, spacedim>::get_nonscalar_data_ranges() const
{
  std::vector<
    std::tuple<unsigned int,
               unsigned int,
               std::string,
               DataComponentInterpretation::DataComponentInterpretation>>
    ranges;

  // Make sure the data structures were set up correctly. Since they
  // can only be filled by build_patches() above, they should have
  // been checked already.
  Assert(dataset_names.size() == data_component_interpretations.size(),
         ExcInternalError());

  const unsigned int n_output_components =
    data_component_interpretations.size();
  unsigned int output_component = 0;
  for (unsigned int i = 0; i < n_output_components; /* i is updated below */)
    switch (data_component_interpretations[i])
      {
        case DataComponentInterpretation::component_is_scalar:
          {
            // OK. Just move component by one:
            ++i;
            ++output_component;

            break;
          }
        case DataComponentInterpretation::component_is_part_of_vector:
          {
            // Ensure that there is a continuous number of next spacedim
            // components that all deal with vectors
            Assert(
              i + spacedim <= n_output_components,
              Exceptions::DataOutImplementation::ExcInvalidVectorDeclaration(
                i, dataset_names[i]));
            for (unsigned int dd = 1; dd < spacedim; ++dd)
              Assert(
                data_component_interpretations[i + dd] ==
                  DataComponentInterpretation::component_is_part_of_vector,
                Exceptions::DataOutImplementation::ExcInvalidVectorDeclaration(
                  i, dataset_names[i]));

            // All seems right, so figure out whether there is a common
            // name to these components. if not, leave the name empty and
            // let the output format writer decide what to do here
            std::string name = dataset_names[i];
            for (unsigned int dd = 1; dd < spacedim; ++dd)
              if (name != dataset_names[i + dd])
                {
                  name = "";
                  break;
                }

            // Finally add a corresponding range.
            //
            // This sort of logic is also explained in some detail in
            //   DataOut::build_one_patch().
            ranges.emplace_back(std::forward_as_tuple(
              output_component,
              output_component + spacedim - 1,
              name,
              DataComponentInterpretation::component_is_part_of_vector));

            // Increase the 'component' counter by the appropriate amount,
            // same for 'i', since we have already dealt with all these
            // components
            output_component += spacedim;
            i += spacedim;

            break;
          }

        case DataComponentInterpretation::component_is_part_of_tensor:
          {
            const unsigned int size = spacedim * spacedim;
            // Ensure that there is a continuous number of next
            // spacedim*spacedim components that all deal with tensors
            Assert(
              i + size <= n_output_components,
              Exceptions::DataOutImplementation::ExcInvalidTensorDeclaration(
                i, dataset_names[i]));
            for (unsigned int dd = 1; dd < size; ++dd)
              Assert(
                data_component_interpretations[i + dd] ==
                  DataComponentInterpretation::component_is_part_of_tensor,
                Exceptions::DataOutImplementation::ExcInvalidTensorDeclaration(
                  i, dataset_names[i]));

            // All seems right, so figure out whether there is a common
            // name to these components. if not, leave the name empty and
            // let the output format writer decide what to do here
            std::string name = dataset_names[i];
            for (unsigned int dd = 1; dd < size; ++dd)
              if (name != dataset_names[i + dd])
                {
                  name = "";
                  break;
                }

            // Finally add a corresponding range.
            ranges.emplace_back(std::forward_as_tuple(
              output_component,
              output_component + size - 1,
              name,
              DataComponentInterpretation::component_is_part_of_tensor));

            // Increase the 'component' counter by the appropriate amount,
            // same for 'i', since we have already dealt with all these
            // components
            output_component += size;
            i += size;
            break;
          }

        default:
          Assert(false, ExcNotImplemented());
      }
  return ranges;
}



// explicit instantiations
#include "numerics/data_out_points.inst"

DEAL_II_NAMESPACE_CLOSE
