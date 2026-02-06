// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2026 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_data_out_points_h
#define dealii_data_out_points_h


#include <deal.II/base/config.h>

#include <deal.II/numerics/data_out_dof_data.h>

#include <cstddef>
#include <string>
#include <vector>

DEAL_II_NAMESPACE_OPEN



/**
 * Produce graphical output defined in a set of points provided by the user.
 *
 * This class provides a simplified interface to output data defined in a
 * set of points (or particles). Internally, it builds 0-dimensional patches
 * defined in the coordinates given by the user. An arbitrary number of
 * data fields of different types (scalar, vector, tensor) can be provided.
 * The class always includes an "ID" field that contains the index of the
 * point. These ids can be made unique in parallel by providing an offset
 * in build_patches().
 *
 * <h3>Interface and Usage</h3>
 *
 * The only method specific to this class is the build_patches() method.
 * Afterwards, output can be generated in various formats using the methods in
 * DataOutInterface. Examples include DataOutInterface::write_vtu() and
 * DataOutInterface::write_vtu_with_pvtu_record():
 *
 * @code
 * std::vector<Point<dim>> points;
 * // [Fill points with data]]
 * DataOutPoints<dim> data_out;
 * data_out.build_patches(points);
 * data_out.write_vtu(file_stream);
 * @endcode
 *
 * <h3>Parallel computations</h3>
 *
 * If the output is to be written in parallel with different points on different
 * MPI ranks, the build_patches() method can be used as follows:
 * @code
 * std::vector<Point<dim>> points;
 * // [Fill points with data on each rank]
 * DataOutPoints<dim> data_out;
 * const std::pair<std::size_t,std::size_t> offset_and_total =
 * Utilities::MPI::partial_and_total_sum(locations.size(), MPI_COMM_WORLD);
 * data_out.build_patches(points, offset_and_total.first);
 * data_out.write_vtu_in_parallel("output.vtu", MPI_COMM_WORLD);
 * @endcode
 *
 */
template <int dim, int spacedim = dim>
class DataOutPoints : public dealii::DataOutInterface<0, spacedim>
{
public:
  /**
   * Build the patches for a given set of points and optionally data in each
   * point.
   *
   * @param [in] locations The point locations.
   * @param [in] offset An offset to be applied to each point index to form the ID.
   * @param [in] data A vector of data values for each point.
   * @param [in] data_component_names An optional vector of strings that
   * describe the properties of each datum.
   * @param [in] data_component_interpretations An optional vector that
   * controls if the properties are interpreted as scalars, vectors,
   * or tensors. Has to be of the same length as @p data_component_names.
   */
  void
  build_patches(
    const std::vector<Point<spacedim>>     &locations,
    const std::size_t                       offset               = 0,
    const std::vector<std::vector<double>> &data                 = {},
    const std::vector<std::string>         &data_component_names = {},
    const std::vector<DataComponentInterpretation::DataComponentInterpretation>
      &data_component_interpretations = {});


protected:
  /**
   * Returns the patches previously built by the build_patches() function.
   */
  virtual const std::vector<DataOutBase::Patch<0, spacedim>> &
  get_patches() const override;

  /**
   * Virtual function through which the names of data sets are obtained from
   * this class
   */
  virtual std::vector<std::string>
  get_dataset_names() const override;

  /**
   * Overload of the respective DataOutInterface::get_nonscalar_data_ranges()
   * function. See there for a more extensive documentation.
   * This function is a reimplementation of the function
   * DataOut_DoFData::get_nonscalar_data_ranges().
   */
  virtual std::vector<
    std::tuple<unsigned int,
               unsigned int,
               std::string,
               DataComponentInterpretation::DataComponentInterpretation>>
  get_nonscalar_data_ranges() const override;

private:
  /**
   * This is a vector of patches that is created each time build_patches() is
   * called. These patches are used in the output routines of the base
   * classes.
   */
  std::vector<DataOutBase::Patch<0, spacedim>> patches;

  /**
   * A vector of field names for all data components stored in patches.
   */
  std::vector<std::string> dataset_names;

  /**
   * A vector that for each of the data components of the
   * current data set indicates whether they are scalar fields, parts of a
   * vector-field, or any of the other supported kinds of data.
   */
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretations;
};

DEAL_II_NAMESPACE_CLOSE

#endif
