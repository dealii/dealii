// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2019 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#include <deal.II/base/config.h>

#include <deal.II/base/utilities.h>

#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/fe_values_views.h>

#include <string>

DEAL_II_NAMESPACE_OPEN

namespace FEValuesExtractors
{
  std::string
  Scalar::get_name() const
  {
    return "Scalar(" + Utilities::int_to_string(component) + ")";
  }


  std::string
  Vector::get_name() const
  {
    return "Vector(" + Utilities::int_to_string(first_vector_component) + ")";
  }


  template <int rank>
  std::string
  Tensor<rank>::get_name() const
  {
    return "Tensor<" + Utilities::int_to_string(rank) + ">(" +
           Utilities::int_to_string(first_tensor_component) + ")";
  }


  template <int rank>
  std::string
  SymmetricTensor<rank>::get_name() const
  {
    return "SymmetricTensor<" + Utilities::int_to_string(rank) + ">(" +
           Utilities::int_to_string(first_tensor_component) + ")";
  }

  // Explicit instantiations
  template struct Tensor<0>;
  template struct Tensor<1>;
  template struct Tensor<2>;
  template struct Tensor<3>;
  template struct Tensor<4>;
  template struct SymmetricTensor<2>;
  template struct SymmetricTensor<4>;



  namespace internal
  {
    template <int dim, int spacedim>
    std::vector<unsigned int>
    generate_component_order(
      const std::vector<FEValuesExtractors::AnyExtractor> &extractors,
      const unsigned int                                   fe_n_components)
    {
      // Initialize `component_order` with invalid unsigned ints representing
      // unassigned state.
      std::vector<unsigned int> component_order(fe_n_components,
                                                numbers::invalid_unsigned_int);

      // Extract the start component index and determine the number of
      // components for each extractor.
      unsigned int block_index = 0;
      for (const auto &extractor : extractors)
        {
          auto start_component_index = numbers::invalid_unsigned_int;
          auto n_components          = numbers::invalid_unsigned_int;

          if (std::holds_alternative<FEValuesExtractors::Scalar>(extractor))
            {
              start_component_index =
                std::get<FEValuesExtractors::Scalar>(extractor).component;
              n_components = 1;
            }
          else if (std::holds_alternative<FEValuesExtractors::Vector>(
                     extractor))
            {
              start_component_index =
                std::get<FEValuesExtractors::Vector>(extractor)
                  .first_vector_component;
              n_components = FEValuesViews::Vector<dim, spacedim>::value_type::
                n_independent_components;
            }
          else if (std::holds_alternative<FEValuesExtractors::Tensor<2>>(
                     extractor))
            {
              start_component_index =
                std::get<FEValuesExtractors::Tensor<2>>(extractor)
                  .first_tensor_component;
              n_components = FEValuesViews::Tensor<2, dim, spacedim>::
                value_type::n_independent_components;
            }
          else if (std::holds_alternative<
                     FEValuesExtractors::SymmetricTensor<2>>(extractor))
            {
              start_component_index =
                std::get<FEValuesExtractors::SymmetricTensor<2>>(extractor)
                  .first_tensor_component;
              n_components = FEValuesViews::SymmetricTensor<2, dim, spacedim>::
                value_type::n_independent_components;
            }
          else
            {
              AssertThrow(
                false,
                ExcNotImplemented(
                  "An unsupported ExtractorVariant was passed in the component_wise extractor_order argument."));
            }

          // Fill `component_order` vector with `n_components` starting at
          // `start_component_index`. Set the values to `block_index`.
          for (unsigned int i = start_component_index;
               i < start_component_index + n_components;
               ++i)
            {
              AssertThrow(i < component_order.size(),
                          ExcIndexRange(i, 0, component_order.size()));

              AssertThrow(
                component_order[i] == numbers::invalid_unsigned_int,
                ExcMessage(
                  "A component which has already been assigned a block "
                  "index is trying to be overwritten. This indicates that the "
                  "component_wise function is being called with an invalid set "
                  "of extractors in the extractor_order argument "
                  "that overlap in component indices."));

              component_order[i] = block_index;
            }
          // Increment block index
          block_index++;
        }

      AssertThrow(
        std::none_of(component_order.begin(),
                     component_order.end(),
                     [](unsigned int value) {
                       return value == numbers::invalid_unsigned_int;
                     }),
        ExcMessage(
          "All entries in component_order must contain valid indices, no "
          "numbers::invalid_unsigned_int must be present after processing all given "
          "extractors. This error typically occurs when the component_wise function "
          "is passed fewer extractors than needed or when the set of extractors "
          "provided in the extractor_order argument doesn't cover all required components "
          "or when passed total number of components f_n_components is incorrect."));

      return component_order;
    }
  } // namespace internal

} // namespace FEValuesExtractors

/*------------------------------- Explicit Instantiations -------------*/

#include "fe/fe_values_extractors.inst"

DEAL_II_NAMESPACE_CLOSE
