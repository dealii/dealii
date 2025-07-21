
// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2000 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_export_lines_h
#define dealii_export_lines_h

#include <deal.II/base/point.h>

#include <deal.II/numerics/data_out.h>


DEAL_II_NAMESPACE_OPEN

/**
 * @brief Exports a collection of line segments to VTK files for visualization.
 *
 * This function takes a vector of line segments, each defined by a pair of
 * start and end points, and writes them to one or more `.vtu` files. In an MPI
 * parallel setting, each process writes its own file (e.g.,
 * `filename.proc0000.vtu`, `filename.proc0001.vtu`, etc.). Additionally,
 * process 0 writes a `.pvtu` file that aggregates all the individual `.vtu`
 * files, allowing visualization software to load the entire distributed data
 * set as a single entity.
 *
 * The function is intended for visualizing geometric entities like shift
 * vectors for Shifted Boundary Method or other line-based data. No data is
 * associated with the lines; only their geometry is exported.
 *
 * @param filename_without_extension The base name for the output files. The
 *        function will append suffixes like `.procXXXX.vtu` and `.pvtu`.
 * @param point_pairs A vector where each element is a `std::pair` of `Point<dim>`,
 *        representing the start and end points of a line segment to be
 *        exported.
 */
template <int dim>
void
export_line_segments(
  const std::string &filename_without_extension,
  const std::vector<std::pair<Point<dim>, Point<dim>>> &point_pairs,
  const MPI_Comm &mpi_communicator = MPI_COMM_WORLD)
{
  const unsigned int n_mpi_process =
    Utilities::MPI::n_mpi_processes(mpi_communicator);
  const unsigned int this_mpi_process =
    Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  using PointOutData = DataOutBase::Patch<1, dim>;
  std::vector<PointOutData> patches_output;
  patches_output.reserve(point_pairs.size());



  Triangulation<1, dim> tria_dummy;
  Triangulation<1, dim> tria_dummy;
  GridGenerator::hyper_cube(tria_dummy, -1., 1.);
  const ReferenceCell reference_cell =
    tria_dummy.begin_active()->reference_cell();

  for (unsigned int segment_index = 0; segment_index < point_pairs.size();
       ++segment_index)
    {
      PointOutData link_out;
      link_out.patch_index    = segment_index;
      link_out.reference_cell = reference_cell;

      link_out.vertices[0] = point_pairs[segment_index].first;
      link_out.vertices[1] = point_pairs[segment_index].second;
      patches_output.emplace_back(std::move(link_out));
    }
  std::vector<std::string> piece_names(n_mpi_process);
  for (unsigned int i = 0; i < n_mpi_process; ++i)
    piece_names[i] = filename_without_extension + ".proc" +
                     Utilities::int_to_string(i, 4) + ".vtu";
  const std::string new_file = piece_names[this_mpi_process];

  const std::string out_pvtu = filename_without_extension + ".pvtu";


  std::ofstream out(new_file);
  std::vector<
    std::tuple<unsigned int,
               unsigned int,
               std::string,
               DataComponentInterpretation::DataComponentInterpretation>>
    vector_data_ranges;

  DataOutBase::write_vtu(
    patches_output, {}, vector_data_ranges, vtu_flags, out);

  if (this_mpi_process == 0)
    {
      std::ofstream pvtu_output(out_pvtu);
      std::ostream &pvtu_out_stream = pvtu_output;
      DataOutBase::write_pvtu_record(pvtu_out_stream,
                                     piece_names,
                                     data_names,
                                     vector_data_ranges,
                                     vtu_flags);
    }
}

DEAL_II_NAMESPACE_CLOSE

#endif
