// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_segment_data_out_h
#define dealii_segment_data_out_h

#include <deal.II/base/config.h>

#include <deal.II/base/data_out_base.h>
#include <deal.II/base/point.h>


DEAL_II_NAMESPACE_OPEN

/**
 * This class generates graphical output for line segment objects.
 *
 * Line segments are one-dimensional objects embedded in a @p spacedim dimensional
 * space. This class creates output data that can be visualized as lines in
 * visualization tools. Each segment is represented as a pair of Point<spacedim>
 * objects defining the start and end points.
 */
template <int spacedim>
class SegmentDataOut : public DataOutInterface<1, spacedim>
{
public:
  /**
   * Constructor.
   */
  SegmentDataOut() = default;

  /**
   * Destructor.
   */
  ~SegmentDataOut() = default;

  /**
   * Generate patches from a vector of line segments and optionally attach
   * per-segment dataset values.
   *
   * If @p datasets and @p dataset_names are empty (default), this function
   * behaves like the old `build_patches()` and only creates the patches.
   * If @p datasets is non-empty, it must have the same size as @p segments and
   * each inner vector must have the same size as @p dataset_names. In that
   * case the data are attached to the generated patches (replaces the former
   * separate `add_datasets()` function).
   */
  void
  build_patches(
    const std::vector<std::pair<Point<spacedim>, Point<spacedim>>> &segments,
    const std::vector<std::vector<double>> &datasets            = {},
    const std::vector<std::string>         &dataset_names_param = {});

protected:
  // @copydoc
  virtual const std::vector<dealii::DataOutBase::Patch<1, spacedim>> &
  get_patches() const override;

  // @copydoc
  virtual std::vector<std::string>
  get_dataset_names() const override;

private:
  /**
   * The actual segments as patches.
   */
  std::vector<DataOutBase::Patch<1, spacedim>> patches;

  /**
   * Names of datasets.
   */
  std::vector<std::string> dataset_names;
};


// Template and inline functions
#ifndef DOXYGEN
template <int spacedim>
void
SegmentDataOut<spacedim>::build_patches(
  const std::vector<std::pair<Point<spacedim>, Point<spacedim>>> &segments,
  const std::vector<std::vector<double>>                         &datasets,
  const std::vector<std::string> &dataset_names_param)
{
  const unsigned int N = segments.size();

  dataset_names.clear();
  patches.resize(N);

  for (unsigned int i = 0; i < N; ++i)
    {
      // A line segment has 2 vertices
      patches[i].vertices[0]          = segments[i].first;
      patches[i].vertices[1]          = segments[i].second;
      patches[i].patch_index          = i;
      patches[i].n_subdivisions       = 1;
      patches[i].reference_cell       = ReferenceCells::Line;
      patches[i].points_are_available = false;
    }

  // If datasets were provided, attach them to the patches. We expect either
  // both datasets and dataset_names to be empty, or both to be non-empty and
  // with matching sizes.
  if (!datasets.empty() || !dataset_names_param.empty())
    {
      AssertDimension(datasets.size(), N);
      AssertDimension(dataset_names_param.size(), datasets[0].size());
      dataset_names = dataset_names_param; // assign names
      for (unsigned int i = 0; i < datasets.size(); ++i)
        {
          AssertDimension(datasets[i].size(), dataset_names.size());
          patches[i].data.reinit(dataset_names.size(), 2);
          for (unsigned int j = 0; j < dataset_names.size(); ++j)
            for (unsigned int k = 0; k < 2; ++k)
              patches[i].data(j, k) = datasets[i][j];
        }
    }
}


// Note: add_datasets has been merged into build_patches. No separate
// implementation remains.


template <int spacedim>
const std::vector<DataOutBase::Patch<1, spacedim>> &
SegmentDataOut<spacedim>::get_patches() const
{
  return patches;
}



template <int spacedim>
std::vector<std::string>
SegmentDataOut<spacedim>::get_dataset_names() const
{
  return dataset_names;
}

#endif

DEAL_II_NAMESPACE_CLOSE

#endif
