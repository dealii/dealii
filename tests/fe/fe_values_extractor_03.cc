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



// test that the FEValuesExtractors::internal::generate_component_order
// work as expected.

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>

#include "../tests.h"

template <typename T>
static std::vector<std::vector<T>>
permutations(std::vector<T> vec)
{
  std::vector<std::vector<T>> result;

  // Permute over indices.
  std::vector<std::size_t> indices(vec.size());
  for (std::size_t i = 0; i < indices.size(); ++i)
    indices[i] = i;

  // Reserve space for permutations
  bool has_more          = true;
  int  permutation_count = 0; // Just in case
  while (has_more)
    {
      std::vector<T> permutation;
      permutation.reserve(vec.size());
      for (std::size_t idx : indices)
        permutation.push_back(vec[idx]);
      result.push_back(std::move(permutation));

      has_more = std::next_permutation(indices.begin(), indices.end());
      permutation_count++;
      if (permutation_count > 500)
        {
          throw std::out_of_range("too many permutations. Check test logic.");
        }
    }

  return result;
}

template <int dim, int spacedim>
static std::string
dims_str()
{
  std::stringstream ss;
  ss << "<" << dim << ", " << spacedim << ">";
  return ss.str();
}


static std::string
args_str(const std::vector<FEValuesExtractors::AnyExtractor> &extractors)
{
  std::stringstream ss;
  ss << "{";
  bool first = true;
  for (const auto &extractor : extractors)
    {
      if (!first)
        ss << ", ";
      first = false;

      ss << std::visit(
        [](const auto &e) -> std::string { return e.get_name(); }, extractor);
    }
  ss << "}";
  return ss.str();
}


template <int dim, int spacedim>
static void
test_generate_component_order()
{
  // Define some extractors
  constexpr FEValuesExtractors::Vector velocity(
    0); // 1 component, starting at pos=0
  constexpr FEValuesExtractors::Scalar pressure(
    spacedim); // dim components, starting at pos=dim
  constexpr FEValuesExtractors::Tensor<2> stress(
    spacedim + 1); // dim*dim components, staring at position dim+1
  constexpr FEValuesExtractors::SymmetricTensor<2> symmetric_stress(
    1 + spacedim +
    (spacedim * spacedim)); // dim*dim components, staring at position dim+1


  // Calculate total number of components
  constexpr int total_components =
    1 + spacedim + (spacedim * spacedim) + (spacedim * spacedim + spacedim) / 2;

  std::vector<FEValuesExtractors::AnyExtractor> extractors;
  extractors.emplace_back(velocity);
  extractors.emplace_back(pressure);
  extractors.emplace_back(stress);
  extractors.emplace_back(symmetric_stress);


  // For each ordering permutation generate the `block_vector` and log it.
  for (const auto &extractor_permutation : permutations<>(extractors))
    {
      // Function under test
      const auto block_vector =
        FEValuesExtractors::internal::generate_component_order<dim, spacedim>(
          extractor_permutation, total_components);

      deallog.push(dims_str<dim, spacedim>());
      deallog << args_str(extractor_permutation) << " -> ";
      for (const auto val : block_vector)
        {
          deallog << val;
        }
      deallog << std::endl;
      deallog.pop();
    }
}


int
main()
{
  initlog();
  test_generate_component_order<1, 1>();
  test_generate_component_order<1, 2>();
  test_generate_component_order<2, 2>();
  test_generate_component_order<1, 3>();
  test_generate_component_order<2, 3>();
  test_generate_component_order<3, 3>();
}
