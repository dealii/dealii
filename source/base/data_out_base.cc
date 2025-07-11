// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1999 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/data_out_base.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/mpi_large_count.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/utilities.h>

#include <deal.II/numerics/data_component_interpretation.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <ctime>
#include <fstream>
#include <future>
#include <iomanip>
#include <limits>
#include <memory>
#include <numeric>
#include <set>
#include <sstream>
#include <vector>

#ifdef DEAL_II_WITH_ZLIB
#  include <zlib.h>
#endif

#ifdef DEAL_II_WITH_HDF5
#  include <hdf5.h>
#endif

#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#ifdef DEAL_II_WITH_ZLIB
#  include <boost/iostreams/filter/zlib.hpp>
#endif



DEAL_II_NAMESPACE_OPEN

#ifndef DOXYGEN
// we need the following exception from a global function, so can't declare it
// in the usual way inside a class
namespace
{
  DeclException2(ExcUnexpectedInput,
                 std::string,
                 std::string,
                 << "Unexpected input: expected line\n  <" << arg1
                 << ">\nbut got\n  <" << arg2 << ">");

#  ifdef DEAL_II_WITH_ZLIB
  constexpr bool deal_ii_with_zlib = true;
#  else
  constexpr bool deal_ii_with_zlib = false;
#  endif


#  ifdef DEAL_II_WITH_ZLIB
  /**
   * Convert between the CompressionLevel enum (used inside VtkFlags
   * for example) and the preprocessor constant defined by zlib.
   */
  int
  get_zlib_compression_level(const DataOutBase::CompressionLevel level)
  {
    switch (level)
      {
        case (DataOutBase::CompressionLevel::no_compression):
          return Z_NO_COMPRESSION;
        case (DataOutBase::CompressionLevel::best_speed):
          return Z_BEST_SPEED;
        case (DataOutBase::CompressionLevel::best_compression):
          return Z_BEST_COMPRESSION;
        case (DataOutBase::CompressionLevel::default_compression):
          return Z_DEFAULT_COMPRESSION;
        default:
          DEAL_II_NOT_IMPLEMENTED();
          return Z_NO_COMPRESSION;
      }
  }

#    ifdef DEAL_II_WITH_MPI
  /**
   * Convert between the CompressionLevel enum and the preprocessor
   * constant defined by boost::iostreams::zlib.
   */
  int
  get_boost_zlib_compression_level(const DataOutBase::CompressionLevel level)
  {
    switch (level)
      {
        case (DataOutBase::CompressionLevel::no_compression):
          return boost::iostreams::zlib::no_compression;
        case (DataOutBase::CompressionLevel::best_speed):
          return boost::iostreams::zlib::best_speed;
        case (DataOutBase::CompressionLevel::best_compression):
          return boost::iostreams::zlib::best_compression;
        case (DataOutBase::CompressionLevel::default_compression):
          return boost::iostreams::zlib::default_compression;
        default:
          DEAL_II_NOT_IMPLEMENTED();
          return boost::iostreams::zlib::no_compression;
      }
  }
#    endif
#  endif

  /**
   * Do a zlib compression followed by a base64 encoding of the given data. The
   * result is then returned as a string object.
   */
  template <typename T>
  std::string
  compress_array(const std::vector<T>               &data,
                 const DataOutBase::CompressionLevel compression_level)
  {
#  ifdef DEAL_II_WITH_ZLIB
    if (data.size() != 0)
      {
        const std::size_t uncompressed_size = (data.size() * sizeof(T));

        // While zlib's compress2 uses unsigned long (which is 64bits
        // on Linux), the vtu compression header stores the block size
        // as an std::uint32_t (see below). While we could implement
        // writing several smaller blocks, we haven't done that. Let's
        // trigger an error for the user instead:
        AssertThrow(uncompressed_size <=
                      std::numeric_limits<std::uint32_t>::max(),
                    ExcNotImplemented());

        // allocate a buffer for compressing data and do so
        auto compressed_data_length = compressBound(uncompressed_size);
        AssertThrow(compressed_data_length <=
                      std::numeric_limits<std::uint32_t>::max(),
                    ExcNotImplemented());

        std::vector<unsigned char> compressed_data(compressed_data_length);

        int err = compress2(&compressed_data[0],
                            &compressed_data_length,
                            reinterpret_cast<const Bytef *>(data.data()),
                            uncompressed_size,
                            get_zlib_compression_level(compression_level));
        (void)err;
        Assert(err == Z_OK, ExcInternalError());

        // Discard the unnecessary bytes
        compressed_data.resize(compressed_data_length);

        // now encode the compression header
        const std::uint32_t compression_header[4] = {
          1,                                             /* number of blocks */
          static_cast<std::uint32_t>(uncompressed_size), /* size of block */
          static_cast<std::uint32_t>(
            uncompressed_size), /* size of last block */
          static_cast<std::uint32_t>(
            compressed_data_length)}; /* list of compressed sizes of blocks */

        const auto *const header_start =
          reinterpret_cast<const unsigned char *>(&compression_header[0]);

        return (Utilities::encode_base64(
                  {header_start, header_start + 4 * sizeof(std::uint32_t)}) +
                Utilities::encode_base64(compressed_data));
      }
    else
      return {};
#  else
    (void)data;
    (void)compression_level;
    Assert(false,
           ExcMessage("This function can only be called if cmake found "
                      "a working libz installation."));
    return {};
#  endif
  }



  /**
   * Convert an array of data objects into a string that will form part of
   * what we then output as data into VTU objects.
   *
   * If libz was found during configuration, this function compresses and
   * encodes the entire data block. Otherwise, it simply writes it element by
   * element.
   */
  template <typename T>
  std::string
  vtu_stringize_array(const std::vector<T>               &data,
                      const DataOutBase::CompressionLevel compression_level,
                      const int                           precision)
  {
    if (deal_ii_with_zlib &&
        (compression_level != DataOutBase::CompressionLevel::plain_text))
      {
        // compress the data we have in memory
        return compress_array(data, compression_level);
      }
    else
      {
        std::ostringstream stream;
        stream.precision(precision);
        for (const T &el : data)
          stream << el << ' ';
        return stream.str();
      }
  }


  /**
   * The header in binary format that the parallel intermediate files
   * start with.
   *
   * @note We are using std::uint64_t for all variables for simplicity below,
   * so that we don't have to worry about packing/data member alignment
   * by the compiler.
   */
  struct ParallelIntermediateHeader
  {
    std::uint64_t magic;
    std::uint64_t version;
    std::uint64_t compression;
    std::uint64_t dimension;
    std::uint64_t space_dimension;
    std::uint64_t n_ranks;
    std::uint64_t n_patches;
  };
} // namespace
#endif


// some declarations of functions and locally used classes
namespace DataOutBase
{
#ifndef DOXYGEN
  namespace
  {
    /**
     * Class holding the data of one cell of a patch in two space dimensions for
     * output. It is the projection of a cell in three-dimensional space (two
     * coordinates, one height value) to the direction of sight.
     */
    class SvgCell
    {
    public:
      // Center of the cell (three-dimensional)
      Point<3> center;

      /**
       * Vector of vertices of this cell (three-dimensional)
       */
      Point<3> vertices[4];

      /**
       * Depth into the picture, which is defined as the distance from an
       * observer at an the origin in direction of the line of sight.
       */
      float depth;

      /**
       * Vector of vertices of this cell (projected, two-dimensional).
       */
      Point<2> projected_vertices[4];

      // Center of the cell (projected, two-dimensional)
      Point<2> projected_center;

      /**
       * Comparison operator for sorting.
       */
      bool
      operator<(const SvgCell &) const;
    };

    bool
    SvgCell::operator<(const SvgCell &e) const
    {
      // note the "wrong" order in which we sort the elements
      return depth > e.depth;
    }



    /**
     * Class holding the data of one cell of a patch in two space dimensions for
     * output. It is the projection of a cell in three-dimensional space (two
     * coordinates, one height value) to the direction of sight.
     */
    class EpsCell2d
    {
    public:
      /**
       * Vector of vertices of this cell.
       */
      Point<2> vertices[4];

      /**
       * Data value from which the actual colors will be computed by the
       * colorization function stated in the <tt>EpsFlags</tt> class.
       */
      float color_value;

      /**
       * Depth into the picture, which is defined as the distance from an
       * observer at an the origin in direction of the line of sight.
       */
      float depth;

      /**
       * Comparison operator for sorting.
       */
      bool
      operator<(const EpsCell2d &) const;
    };

    bool
    EpsCell2d::operator<(const EpsCell2d &e) const
    {
      // note the "wrong" order in which we sort the elements
      return depth > e.depth;
    }



    /**
     * This is a helper function that converts all of the data stored
     * in the `patches` array into one global data table. That data
     * table has as many rows as there are data sets in the patches,
     * and as many columns as there data points at which to output
     * data. In the end, each data set is then stored in one row of
     * this table, rather than scattered throughout the various patches.
     *
     * This function is used by all those output formats that write
     * data one data set at a time, rather than one cell at a time.
     */
    template <int dim, int spacedim, typename Number = double>
    std::unique_ptr<Table<2, Number>>
    create_global_data_table(const std::vector<Patch<dim, spacedim>> &patches)
    {
      // If there is nothing to write, just return
      if (patches.empty())
        return std::make_unique<Table<2, Number>>();

      // unlike in the main function, we don't have here the data_names field,
      // so we initialize it with the number of data sets in the first patch.
      // the equivalence of these two definitions is checked in the main
      // function.

      // we have to take care, however, whether the points are appended to the
      // end of the patch.data table
      const unsigned int n_data_sets = patches[0].points_are_available ?
                                         (patches[0].data.n_rows() - spacedim) :
                                         patches[0].data.n_rows();
      const unsigned int n_data_points =
        std::accumulate(patches.begin(),
                        patches.end(),
                        0U,
                        [](const unsigned int          count,
                           const Patch<dim, spacedim> &patch) {
                          return count + patch.data.n_cols();
                        });

      std::unique_ptr<Table<2, Number>> global_data_table =
        std::make_unique<Table<2, Number>>(n_data_sets, n_data_points);

      // loop over all patches
      unsigned int next_value = 0;
      for (const auto &patch : patches)
        {
          const unsigned int n_subdivisions = patch.n_subdivisions;
          (void)n_subdivisions;

          Assert((patch.data.n_rows() == n_data_sets &&
                  !patch.points_are_available) ||
                   (patch.data.n_rows() == n_data_sets + spacedim &&
                    patch.points_are_available),
                 ExcDimensionMismatch(patch.points_are_available ?
                                        (n_data_sets + spacedim) :
                                        n_data_sets,
                                      patch.data.n_rows()));
          Assert(patch.reference_cell != ReferenceCells::get_hypercube<dim>() ||
                   (n_data_sets == 0) ||
                   (patch.data.n_cols() ==
                    Utilities::fixed_power<dim>(n_subdivisions + 1)),
                 ExcInvalidDatasetSize(patch.data.n_cols(),
                                       n_subdivisions + 1));

          for (unsigned int i = 0; i < patch.data.n_cols(); ++i, ++next_value)
            for (unsigned int data_set = 0; data_set < n_data_sets; ++data_set)
              (*global_data_table)[data_set][next_value] =
                patch.data(data_set, i);
        }
      Assert(next_value == n_data_points, ExcInternalError());

      return global_data_table;
    }
  } // namespace


#endif


  DataOutFilter::DataOutFilter()
    : flags(false, true)
    , node_dim(numbers::invalid_unsigned_int)
    , num_cells(0)
  {}



  DataOutFilter::DataOutFilter(const DataOutBase::DataOutFilterFlags &flags)
    : flags(flags)
    , node_dim(numbers::invalid_unsigned_int)
    , num_cells(0)
  {}



  template <int dim>
  void
  DataOutFilter::write_point(const unsigned int index, const Point<dim> &p)
  {
    node_dim = dim;

    Point<3> int_pt;
    for (unsigned int d = 0; d < dim; ++d)
      int_pt[d] = p[d];

    const Map3DPoint::const_iterator it = existing_points.find(int_pt);
    unsigned int                     internal_ind;

    // If the point isn't in the set, or we're not filtering duplicate points,
    // add it
    if (it == existing_points.end() || !flags.filter_duplicate_vertices)
      {
        internal_ind = existing_points.size();
        existing_points.insert(std::make_pair(int_pt, internal_ind));
      }
    else
      {
        internal_ind = it->second;
      }
    // Now add the index to the list of filtered points
    filtered_points[index] = internal_ind;
  }



  void
  DataOutFilter::internal_add_cell(const unsigned int cell_index,
                                   const unsigned int pt_index)
  {
    filtered_cells[cell_index] = filtered_points[pt_index];

    // (Re)-initialize counter at any first call to this method.
    if (cell_index == 0)
      num_cells = 1;
  }



  void
  DataOutFilter::fill_node_data(std::vector<double> &node_data) const
  {
    node_data.resize(existing_points.size() * node_dim);

    for (const auto &existing_point : existing_points)
      {
        for (unsigned int d = 0; d < node_dim; ++d)
          node_data[node_dim * existing_point.second + d] =
            existing_point.first[d];
      }
  }



  void
  DataOutFilter::fill_cell_data(const unsigned int         local_node_offset,
                                std::vector<unsigned int> &cell_data) const
  {
    cell_data.resize(filtered_cells.size());

    for (const auto &filtered_cell : filtered_cells)
      {
        cell_data[filtered_cell.first] =
          filtered_cell.second + local_node_offset;
      }
  }



  std::string
  DataOutFilter::get_data_set_name(const unsigned int set_num) const
  {
    return data_set_names.at(set_num);
  }



  unsigned int
  DataOutFilter::get_data_set_dim(const unsigned int set_num) const
  {
    return data_set_dims.at(set_num);
  }



  const double *
  DataOutFilter::get_data_set(const unsigned int set_num) const
  {
    return data_sets[set_num].data();
  }



  unsigned int
  DataOutFilter::n_nodes() const
  {
    return existing_points.size();
  }



  unsigned int
  DataOutFilter::n_cells() const
  {
    return num_cells;
  }



  unsigned int
  DataOutFilter::n_data_sets() const
  {
    return data_set_names.size();
  }



  void
  DataOutFilter::flush_points()
  {}



  void
  DataOutFilter::flush_cells()
  {}



  template <int dim>
  void
  DataOutFilter::write_cell(const unsigned int                   index,
                            const unsigned int                   start,
                            const std::array<unsigned int, dim> &offsets)
  {
    ++num_cells;

    const unsigned int base_entry =
      index * GeometryInfo<dim>::vertices_per_cell;

    switch (dim)
      {
        case 0:
          {
            internal_add_cell(base_entry + 0, start);
            break;
          }

        case 1:
          {
            const unsigned int d1 = offsets[0];

            internal_add_cell(base_entry + 0, start);
            internal_add_cell(base_entry + 1, start + d1);
            break;
          }

        case 2:
          {
            const unsigned int d1 = offsets[0];
            const unsigned int d2 = offsets[1];

            internal_add_cell(base_entry + 0, start);
            internal_add_cell(base_entry + 1, start + d1);
            internal_add_cell(base_entry + 2, start + d2 + d1);
            internal_add_cell(base_entry + 3, start + d2);
            break;
          }

        case 3:
          {
            const unsigned int d1 = offsets[0];
            const unsigned int d2 = offsets[1];
            const unsigned int d3 = offsets[2];

            internal_add_cell(base_entry + 0, start);
            internal_add_cell(base_entry + 1, start + d1);
            internal_add_cell(base_entry + 2, start + d2 + d1);
            internal_add_cell(base_entry + 3, start + d2);
            internal_add_cell(base_entry + 4, start + d3);
            internal_add_cell(base_entry + 5, start + d3 + d1);
            internal_add_cell(base_entry + 6, start + d3 + d2 + d1);
            internal_add_cell(base_entry + 7, start + d3 + d2);
            break;
          }

        default:
          DEAL_II_NOT_IMPLEMENTED();
      }
  }



  void
  DataOutFilter::write_cell_single(const unsigned int   index,
                                   const unsigned int   start,
                                   const unsigned int   n_points,
                                   const ReferenceCell &reference_cell)
  {
    ++num_cells;

    const unsigned int base_entry = index * n_points;

    static const std::array<unsigned int, 5> table = {{0, 1, 3, 2, 4}};

    for (unsigned int i = 0; i < n_points; ++i)
      internal_add_cell(base_entry + i,
                        start + (reference_cell == ReferenceCells::Pyramid ?
                                   table[i] :
                                   i));
  }



  void
  DataOutFilter::write_data_set(const std::string      &name,
                                const unsigned int      dimension,
                                const unsigned int      set_num,
                                const Table<2, double> &data_vectors)
  {
    unsigned int new_dim;

    // HDF5/XDMF output only supports 1d or 3d output, so force rearrangement if
    // needed
    if (flags.xdmf_hdf5_output && dimension != 1)
      new_dim = 3;
    else
      new_dim = dimension;

    // Record the data set name, dimension, and allocate space for it
    data_set_names.push_back(name);
    data_set_dims.push_back(new_dim);
    data_sets.emplace_back(new_dim * existing_points.size());

    // TODO: averaging, min/max, etc for merged vertices
    for (unsigned int i = 0; i < filtered_points.size(); ++i)
      {
        const unsigned int r = filtered_points[i];

        for (unsigned int d = 0; d < new_dim; ++d)
          {
            if (d < dimension)
              data_sets.back()[r * new_dim + d] = data_vectors(set_num + d, i);
            else
              data_sets.back()[r * new_dim + d] = 0;
          }
      }
  }
} // namespace DataOutBase



//----------------------------------------------------------------------//
// Auxiliary data
//----------------------------------------------------------------------//

namespace
{
  const char *gmv_cell_type[4] = {"", "line 2", "quad 4", "hex 8"};

  const char *ucd_cell_type[4] = {"pt", "line", "quad", "hex"};

  const char *tecplot_cell_type[4] = {"", "lineseg", "quadrilateral", "brick"};

  /**
   * Return the tuple (vtk cell type, number of cells, number of nodes)
   * for a patch.
   *
   * The logic used here is as follows:
   * - If a cell is not subdivided or we don't use higher order cells,
   *   then we use linear cells
   * - For hypercubes, we support subdividing cells into sub-cells,
   *   which are then treated as each being linear
   * - For triangles and tetrahedra, we special-case the situation of
   *   n_subdivisions==2, in which case we treat the cell as a single
   *   quadratic cell (i.e., higher order)
   */
  template <int dim, int spacedim>
  std::array<unsigned int, 3>
  extract_vtk_patch_info(const DataOutBase::Patch<dim, spacedim> &patch,
                         const bool write_higher_order_cells)
  {
    std::array<unsigned int, 3> vtk_cell_id = {
      {/* cell type, tbd: */ numbers::invalid_unsigned_int,
       /* # of cells, default: just one cell */ 1,
       /* # of nodes, default: as many nodes as vertices */
       patch.reference_cell.n_vertices()}};

    if (write_higher_order_cells)
      {
        vtk_cell_id[0] = patch.reference_cell.vtk_lagrange_type();
        vtk_cell_id[2] = patch.data.n_cols();
      }
    else if (patch.data.n_cols() == patch.reference_cell.n_vertices())
      // One data set per vertex -> a linear cell
      vtk_cell_id[0] = patch.reference_cell.vtk_linear_type();
    else if (patch.reference_cell == ReferenceCells::Triangle &&
             patch.data.n_cols() == 6)
      {
        Assert(patch.n_subdivisions == 2, ExcInternalError());
        vtk_cell_id[0] = patch.reference_cell.vtk_quadratic_type();
        vtk_cell_id[2] = patch.data.n_cols();
      }
    else if (patch.reference_cell == ReferenceCells::Tetrahedron &&
             patch.data.n_cols() == 10)
      {
        Assert(patch.n_subdivisions == 2, ExcInternalError());
        vtk_cell_id[0] = patch.reference_cell.vtk_quadratic_type();
        vtk_cell_id[2] = patch.data.n_cols();
      }
    else if (patch.reference_cell.is_hyper_cube())
      {
        // For hypercubes, we support sub-divided linear cells
        vtk_cell_id[0] = patch.reference_cell.vtk_linear_type();
        vtk_cell_id[1] = Utilities::pow(patch.n_subdivisions, dim);
      }
    else if (patch.reference_cell.is_simplex())
      {
        vtk_cell_id[0] = patch.reference_cell.vtk_lagrange_type();
        vtk_cell_id[2] = patch.data.n_cols();
      }
    else
      {
        DEAL_II_NOT_IMPLEMENTED();
      }

    return vtk_cell_id;
  }

  //----------------------------------------------------------------------//
  // Auxiliary functions
  //----------------------------------------------------------------------//

  // For a given patch that corresponds to a hypercube cell, compute the
  // location of a node interpolating the corner nodes linearly
  // at the point lattice_location/n_subdivisions where lattice_location
  // is a dim-dimensional integer vector. If the points are
  // saved in the patch.data member, return the saved point instead.
  template <int dim, int spacedim>
  inline Point<spacedim>
  get_equispaced_location(
    const DataOutBase::Patch<dim, spacedim>   &patch,
    const std::initializer_list<unsigned int> &lattice_location,
    const unsigned int                         n_subdivisions)
  {
    // This function only makes sense when called on hypercube cells
    Assert(patch.reference_cell.is_hyper_cube(), ExcInternalError());

    Assert(lattice_location.size() == dim, ExcInternalError());

    const unsigned int xstep = (dim > 0 ? *(lattice_location.begin() + 0) : 0);
    const unsigned int ystep = (dim > 1 ? *(lattice_location.begin() + 1) : 0);
    const unsigned int zstep = (dim > 2 ? *(lattice_location.begin() + 2) : 0);

    // If the patch stores the locations of nodes (rather than of only the
    // vertices), then obtain the location by direct lookup.
    if (patch.points_are_available)
      {
        Assert(n_subdivisions == patch.n_subdivisions, ExcNotImplemented());

        unsigned int point_no = 0;
        switch (dim)
          {
            case 3:
              AssertIndexRange(zstep, n_subdivisions + 1);
              point_no += (n_subdivisions + 1) * (n_subdivisions + 1) * zstep;
              DEAL_II_FALLTHROUGH;
            case 2:
              AssertIndexRange(ystep, n_subdivisions + 1);
              point_no += (n_subdivisions + 1) * ystep;
              DEAL_II_FALLTHROUGH;
            case 1:
              AssertIndexRange(xstep, n_subdivisions + 1);
              point_no += xstep;
              DEAL_II_FALLTHROUGH;
            case 0:
              // break here for dim<=3
              break;

            default:
              DEAL_II_NOT_IMPLEMENTED();
          }
        Point<spacedim> node;
        for (unsigned int d = 0; d < spacedim; ++d)
          node[d] = patch.data(patch.data.size(0) - spacedim + d, point_no);
        return node;
      }
    else
      // The patch does not store node locations, so we have to interpolate
      // between its vertices:
      {
        if constexpr (dim == 0)
          return patch.vertices[0];
        else
          {
            // perform a dim-linear interpolation
            const double stepsize = 1. / n_subdivisions;
            const double xfrac    = xstep * stepsize;

            Point<spacedim> node =
              (patch.vertices[1] * xfrac) + (patch.vertices[0] * (1 - xfrac));
            if (dim > 1)
              {
                const double yfrac = ystep * stepsize;
                node *= 1 - yfrac;
                node += ((patch.vertices[3] * xfrac) +
                         (patch.vertices[2] * (1 - xfrac))) *
                        yfrac;
                if (dim > 2)
                  {
                    const double zfrac = zstep * stepsize;
                    node *= (1 - zfrac);
                    node += (((patch.vertices[5] * xfrac) +
                              (patch.vertices[4] * (1 - xfrac))) *
                               (1 - yfrac) +
                             ((patch.vertices[7] * xfrac) +
                              (patch.vertices[6] * (1 - xfrac))) *
                               yfrac) *
                            zfrac;
                  }
              }
            return node;
          }
      }
  }

  // For a given patch, compute the nodes for arbitrary (non-hypercube) cells.
  // If the points are saved in the patch.data member, return the saved point
  // instead.
  template <int dim, int spacedim>
  inline Point<spacedim>
  get_node_location(const DataOutBase::Patch<dim, spacedim> &patch,
                    const unsigned int                       node_index)
  {
    // Due to a historical accident, we are using a different indexing
    // for pyramids in this file than we do where we create patches.
    // So translate if necessary.
    unsigned int point_no_actual = node_index;
    if (patch.reference_cell == ReferenceCells::Pyramid)
      {
        AssertDimension(patch.n_subdivisions, 1);

        static const std::array<unsigned int, 5> table = {{0, 1, 3, 2, 4}};
        point_no_actual                                = table[node_index];
      }

    // If the patch stores the locations of nodes (rather than of only the
    // vertices), then obtain the location by direct lookup.
    if (patch.points_are_available)
      {
        Point<spacedim> node;
        for (unsigned int d = 0; d < spacedim; ++d)
          node[d] =
            patch.data(patch.data.size(0) - spacedim + d, point_no_actual);
        return node;
      }
    else
      // The patch does not store node locations, so we have to interpolate
      // between its vertices. This isn't currently implemented for anything
      // other than one subdivision, but would go here.
      //
      // For n_subdivisions==1, the locations are simply those of vertices, so
      // get the information from there.
      {
        AssertDimension(patch.n_subdivisions, 1);

        return patch.vertices[point_no_actual];
      }
  }



  /**
   * Count the number of nodes and cells referenced by the given
   * argument, and return these numbers (in order nodes, then cells)
   * as a tuple.
   */
  template <int dim, int spacedim>
  std::tuple<unsigned int, unsigned int>
  count_nodes_and_cells(
    const std::vector<DataOutBase::Patch<dim, spacedim>> &patches)
  {
    unsigned int n_nodes = 0;
    unsigned int n_cells = 0;
    for (const auto &patch : patches)
      {
        Assert(patch.reference_cell != ReferenceCells::Invalid,
               ExcMessage(
                 "The reference cell for this patch is set to 'Invalid', "
                 "but that is clearly not a valid choice. Did you forget "
                 "to set the reference cell for the patch?"));

        if (patch.reference_cell.is_hyper_cube())
          {
            n_nodes += Utilities::fixed_power<dim>(patch.n_subdivisions + 1);
            n_cells += Utilities::fixed_power<dim>(patch.n_subdivisions);
          }
        else
          {
            Assert(patch.n_subdivisions == 1, ExcNotImplemented());
            n_nodes += patch.reference_cell.n_vertices();
            n_cells += 1;
          }
      }

    return std::make_tuple(n_nodes, n_cells);
  }



  /**
   * Count the number of nodes and cells referenced by the given
   * argument, and return these numbers (in order nodes, then cells, then cells
   * plus points) as a tuple.
   */
  template <int dim, int spacedim>
  std::tuple<unsigned int, unsigned int, unsigned int>
  count_nodes_and_cells_and_points(
    const std::vector<DataOutBase::Patch<dim, spacedim>> &patches,
    const bool write_higher_order_cells)
  {
    unsigned int n_nodes              = 0;
    unsigned int n_cells              = 0;
    unsigned int n_points_and_n_cells = 0;

    for (const auto &patch : patches)
      {
        if (patch.reference_cell.is_hyper_cube())
          {
            n_nodes += Utilities::fixed_power<dim>(patch.n_subdivisions + 1);

            if (write_higher_order_cells)
              {
                // Write all of these nodes as a single higher-order cell. So
                // add one to the number of cells, and update the number of
                // points appropriately.
                n_cells += 1;
                n_points_and_n_cells +=
                  1 + Utilities::fixed_power<dim>(patch.n_subdivisions + 1);
              }
            else
              {
                // Write all of these nodes as a collection of d-linear
                // cells. Add the number of sub-cells to the total number of
                // cells, and then add one for each cell plus the number of
                // vertices per cell for each subcell to the number of points.
                const unsigned int n_subcells =
                  Utilities::fixed_power<dim>(patch.n_subdivisions);
                n_cells += n_subcells;
                n_points_and_n_cells +=
                  n_subcells * (1 + GeometryInfo<dim>::vertices_per_cell);
              }
          }
        else
          {
            n_nodes += patch.data.n_cols();
            n_cells += 1;
            n_points_and_n_cells += patch.data.n_cols() + 1;
          }
      }

    return std::make_tuple(n_nodes, n_cells, n_points_and_n_cells);
  }

  /**
   * Class describing common functionality between different output streams.
   *
   * @ingroup output
   */
  template <typename FlagsType>
  class StreamBase
  {
  public:
    /*
     * Constructor. Stores a reference to the output stream for immediate use.
     */
    StreamBase(std::ostream &stream, const FlagsType &flags)
      : selected_component(numbers::invalid_unsigned_int)
      , stream(stream)
      , flags(flags)
    {}

    /**
     * Output operator for points. All inheriting classes should implement this
     * function.
     */
    template <int dim>
    void
    write_point(const unsigned int, const Point<dim> &)
    {
      Assert(false,
             ExcMessage("The derived class you are using needs to "
                        "reimplement this function if you want to call "
                        "it."));
    }

    /**
     * Do whatever is necessary to terminate the list of points. The default
     * implementation does nothing; derived classes that do not require any
     * action do not need to reimplement this.
     */
    void
    flush_points()
    {}

    /**
     * Write dim-dimensional cell with first vertex at number start and further
     * vertices offset by the specified values. Values not needed are ignored.
     * All inheriting classes should implement this function.
     */
    template <int dim>
    void
    write_cell(const unsigned int /*index*/,
               const unsigned int /*start*/,
               std::array<unsigned int, dim> & /*offsets*/)
    {
      Assert(false,
             ExcMessage("The derived class you are using needs to "
                        "reimplement this function if you want to call "
                        "it."));
    }

    /**
     * Write dim-dimensional @p index cell with @p n_point vertices and first
     * vertex at number @p start.
     *
     * @note All inheriting classes should implement this function.
     */
    void
    write_cell_single(const unsigned int   index,
                      const unsigned int   start,
                      const unsigned int   n_points,
                      const ReferenceCell &reference_cell)
    {
      (void)index;
      (void)start;
      (void)n_points;
      (void)reference_cell;

      Assert(false,
             ExcMessage("The derived class you are using needs to "
                        "reimplement this function if you want to call "
                        "it."));
    }

    /**
     * Do whatever is necessary to terminate the list of cells. This function is
     * usually only reimplemented if deal.II is compiled with zlib. The default
     * implementation does nothing; derived classes that do not require any
     * action do not need to reimplement this.
     */
    void
    flush_cells()
    {}

    /**
     * Forwarding of an output stream. This function is usually only
     * reimplemented if inheriting classes use zlib.
     */
    template <typename T>
    std::ostream &
    operator<<(const T &t)
    {
      stream << t;
      return stream;
    }

    /**
     * Since the GMV and Tecplot formats read the x, y and z coordinates in
     * separate fields, we enable write() to output only a single selected
     * component at once and do this dim times for the whole set of nodes. This
     * integer can be used to select the component written.
     */
    unsigned int selected_component;

  protected:
    /**
     * The ostream to use. Since the life span of these objects is small, we use
     * a very simple storage technique.
     */
    std::ostream &stream;

    /**
     * The flags controlling the output.
     */
    const FlagsType flags;
  };

  /**
   * Class for writing basic entities in @ref SoftwareOpenDX format, depending on the flags.
   */
  class DXStream : public StreamBase<DataOutBase::DXFlags>
  {
  public:
    DXStream(std::ostream &stream, const DataOutBase::DXFlags &flags);

    template <int dim>
    void
    write_point(const unsigned int index, const Point<dim> &);

    /**
     * The order of vertices for these cells in different dimensions is
     * <ol>
     * <li> [0,1]
     * <li> [0,2,1,3]
     * <li> [0,4,2,6,1,5,3,7]
     * </ol>
     */
    template <int dim>
    void
    write_cell(const unsigned int                   index,
               const unsigned int                   start,
               const std::array<unsigned int, dim> &offsets);

    /**
     * Write a complete set of data for a single node.
     *
     * The index given as first argument indicates the number of a data set, as
     * some output formats require this number to be printed.
     */
    template <typename data>
    void
    write_dataset(const unsigned int index, const std::vector<data> &values);
  };

  /**
   * Class for writing basic entities in @ref SoftwareGMV format, depending on the flags.
   */
  class GmvStream : public StreamBase<DataOutBase::GmvFlags>
  {
  public:
    GmvStream(std::ostream &stream, const DataOutBase::GmvFlags &flags);

    template <int dim>
    void
    write_point(const unsigned int index, const Point<dim> &);

    /**
     * The order of vertices for these cells in different dimensions is
     * <ol>
     * <li> [0,1]
     * <li> [0,1,3,2]
     * <li> [0,1,3,2,4,5,7,6]
     * </ol>
     */
    template <int dim>
    void
    write_cell(const unsigned int                   index,
               const unsigned int                   start,
               const std::array<unsigned int, dim> &offsets);
  };

  /**
   * Class for writing basic entities in @ref SoftwareTecplot format, depending on the flags.
   */
  class TecplotStream : public StreamBase<DataOutBase::TecplotFlags>
  {
  public:
    TecplotStream(std::ostream &stream, const DataOutBase::TecplotFlags &flags);

    template <int dim>
    void
    write_point(const unsigned int index, const Point<dim> &);

    /**
     * The order of vertices for these cells in different dimensions is
     * <ol>
     * <li> [0,1]
     * <li> [0,1,3,2]
     * <li> [0,1,3,2,4,5,7,6]
     * </ol>
     */
    template <int dim>
    void
    write_cell(const unsigned int                   index,
               const unsigned int                   start,
               const std::array<unsigned int, dim> &offsets);
  };

  /**
   * Class for writing basic entities in UCD format for @ref SoftwareAVS, depending on the flags.
   */
  class UcdStream : public StreamBase<DataOutBase::UcdFlags>
  {
  public:
    UcdStream(std::ostream &stream, const DataOutBase::UcdFlags &flags);

    template <int dim>
    void
    write_point(const unsigned int index, const Point<dim> &);

    /**
     * The additional offset 1 is added inside this function.
     *
     * The order of vertices for these cells in different dimensions is
     * <ol>
     * <li> [0,1]
     * <li> [0,1,3,2]
     * <li> [0,1,5,4,2,3,7,6]
     * </ol>
     */
    template <int dim>
    void
    write_cell(const unsigned int                   index,
               const unsigned int                   start,
               const std::array<unsigned int, dim> &offsets);

    /**
     * Write a complete set of data for a single node.
     *
     * The index given as first argument indicates the number of a data set, as
     * some output formats require this number to be printed.
     */
    template <typename data>
    void
    write_dataset(const unsigned int index, const std::vector<data> &values);
  };

  /**
   * Class for writing basic entities in @ref SoftwareVTK format, depending on the flags.
   */
  class VtkStream : public StreamBase<DataOutBase::VtkFlags>
  {
  public:
    VtkStream(std::ostream &stream, const DataOutBase::VtkFlags &flags);

    template <int dim>
    void
    write_point(const unsigned int index, const Point<dim> &);

    /**
     * The order of vertices for these cells in different dimensions is
     * <ol>
     * <li> [0,1]
     * <li> []
     * <li> []
     * </ol>
     */
    template <int dim>
    void
    write_cell(const unsigned int                   index,
               const unsigned int                   start,
               const std::array<unsigned int, dim> &offsets);

    /**
     * Print vertices [start, start+n_points[
     */
    void
    write_cell_single(const unsigned int   index,
                      const unsigned int   start,
                      const unsigned int   n_points,
                      const ReferenceCell &reference_cell);

    /**
     * Write a high-order cell type, i.e., a Lagrange cell
     * in the VTK terminology.
     * The connectivity order of the points is given in the
     * @p connectivity array, which are offset
     * by the global index @p start.
     */
    template <int dim>
    void
    write_high_order_cell(const unsigned int           start,
                          const std::vector<unsigned> &connectivity);
  };


  //----------------------------------------------------------------------//

  DXStream::DXStream(std::ostream &out, const DataOutBase::DXFlags &f)
    : StreamBase<DataOutBase::DXFlags>(out, f)
  {}


  template <int dim>
  void
  DXStream::write_point(const unsigned int, const Point<dim> &p)
  {
    if (flags.coordinates_binary)
      {
        float data[dim];
        for (unsigned int d = 0; d < dim; ++d)
          data[d] = p[d];
        stream.write(reinterpret_cast<const char *>(data), dim * sizeof(*data));
      }
    else
      {
        for (unsigned int d = 0; d < dim; ++d)
          stream << p[d] << '\t';
        stream << '\n';
      }
  }



  // Separate these out to avoid an internal compiler error with intel 17
  namespace DataOutBaseImplementation
  {
    /**
     * Set up the node numbers for a given cell being written to an output
     * stream.
     */
    std::array<unsigned int, GeometryInfo<0>::vertices_per_cell>
    set_node_numbers(const unsigned int /*start*/,
                     const std::array<unsigned int, 0> & /*d1*/)
    {
      DEAL_II_ASSERT_UNREACHABLE();
      return {};
    }



    std::array<unsigned int, GeometryInfo<1>::vertices_per_cell>
    set_node_numbers(const unsigned int                 start,
                     const std::array<unsigned int, 1> &offsets)
    {
      std::array<unsigned int, GeometryInfo<1>::vertices_per_cell> nodes;
      nodes[0] = start;
      nodes[1] = start + offsets[0];
      return nodes;
    }



    std::array<unsigned int, GeometryInfo<2>::vertices_per_cell>
    set_node_numbers(const unsigned int                 start,
                     const std::array<unsigned int, 2> &offsets)

    {
      const unsigned int d1 = offsets[0];
      const unsigned int d2 = offsets[1];

      std::array<unsigned int, GeometryInfo<2>::vertices_per_cell> nodes;
      nodes[0] = start;
      nodes[1] = start + d1;
      nodes[2] = start + d2;
      nodes[3] = start + d2 + d1;
      return nodes;
    }



    std::array<unsigned int, GeometryInfo<3>::vertices_per_cell>
    set_node_numbers(const unsigned int                 start,
                     const std::array<unsigned int, 3> &offsets)
    {
      const unsigned int d1 = offsets[0];
      const unsigned int d2 = offsets[1];
      const unsigned int d3 = offsets[2];

      std::array<unsigned int, GeometryInfo<3>::vertices_per_cell> nodes;
      nodes[0] = start;
      nodes[1] = start + d1;
      nodes[2] = start + d2;
      nodes[3] = start + d2 + d1;
      nodes[4] = start + d3;
      nodes[5] = start + d3 + d1;
      nodes[6] = start + d3 + d2;
      nodes[7] = start + d3 + d2 + d1;
      return nodes;
    }
  } // namespace DataOutBaseImplementation



  template <int dim>
  void
  DXStream::write_cell(const unsigned int,
                       const unsigned int                   start,
                       const std::array<unsigned int, dim> &offsets)
  {
    const auto nodes =
      DataOutBaseImplementation::set_node_numbers(start, offsets);

    if (flags.int_binary)
      {
        std::array<unsigned int, GeometryInfo<dim>::vertices_per_cell> temp;
        for (unsigned int i = 0; i < nodes.size(); ++i)
          temp[i] = nodes[GeometryInfo<dim>::dx_to_deal[i]];
        stream.write(reinterpret_cast<const char *>(temp.data()),
                     temp.size() * sizeof(temp[0]));
      }
    else
      {
        for (unsigned int i = 0; i < nodes.size() - 1; ++i)
          stream << nodes[GeometryInfo<dim>::dx_to_deal[i]] << '\t';
        stream << nodes[GeometryInfo<dim>::dx_to_deal[nodes.size() - 1]]
               << '\n';
      }
  }



  template <typename data>
  inline void
  DXStream::write_dataset(const unsigned int, const std::vector<data> &values)
  {
    if (flags.data_binary)
      {
        stream.write(reinterpret_cast<const char *>(values.data()),
                     values.size() * sizeof(data));
      }
    else
      {
        for (unsigned int i = 0; i < values.size(); ++i)
          stream << '\t' << values[i];
        stream << '\n';
      }
  }



  //----------------------------------------------------------------------//

  GmvStream::GmvStream(std::ostream &out, const DataOutBase::GmvFlags &f)
    : StreamBase<DataOutBase::GmvFlags>(out, f)
  {}


  template <int dim>
  void
  GmvStream::write_point(const unsigned int, const Point<dim> &p)
  {
    Assert(selected_component != numbers::invalid_unsigned_int,
           ExcNotInitialized());
    stream << p[selected_component] << ' ';
  }



  template <int dim>
  void
  GmvStream::write_cell(const unsigned int,
                        const unsigned int                   s,
                        const std::array<unsigned int, dim> &offsets)
  {
    // Vertices are numbered starting with one.
    const unsigned int start = s + 1;
    stream << gmv_cell_type[dim] << '\n';

    switch (dim)
      {
        case 0:
          {
            stream << start;
            break;
          }

        case 1:
          {
            const unsigned int d1 = offsets[0];
            stream << start;
            stream << '\t' << start + d1;
            break;
          }

        case 2:
          {
            const unsigned int d1 = offsets[0];
            const unsigned int d2 = offsets[1];
            stream << start;
            stream << '\t' << start + d1;
            stream << '\t' << start + d2 + d1 << '\t' << start + d2;
            break;
          }

        case 3:
          {
            const unsigned int d1 = offsets[0];
            const unsigned int d2 = offsets[1];
            const unsigned int d3 = offsets[2];
            stream << start;
            stream << '\t' << start + d1;
            stream << '\t' << start + d2 + d1 << '\t' << start + d2;
            stream << '\t' << start + d3 << '\t' << start + d3 + d1 << '\t'
                   << start + d3 + d2 + d1 << '\t' << start + d3 + d2;
            break;
          }

        default:
          DEAL_II_NOT_IMPLEMENTED();
      }
    stream << '\n';
  }



  TecplotStream::TecplotStream(std::ostream                    &out,
                               const DataOutBase::TecplotFlags &f)
    : StreamBase<DataOutBase::TecplotFlags>(out, f)
  {}


  template <int dim>
  void
  TecplotStream::write_point(const unsigned int, const Point<dim> &p)
  {
    Assert(selected_component != numbers::invalid_unsigned_int,
           ExcNotInitialized());
    stream << p[selected_component] << '\n';
  }



  template <int dim>
  void
  TecplotStream::write_cell(const unsigned int,
                            const unsigned int                   s,
                            const std::array<unsigned int, dim> &offsets)
  {
    const unsigned int start = s + 1;

    switch (dim)
      {
        case 0:
          {
            stream << start;
            break;
          }

        case 1:
          {
            const unsigned int d1 = offsets[0];
            stream << start;
            stream << '\t' << start + d1;
            break;
          }

        case 2:
          {
            const unsigned int d1 = offsets[0];
            const unsigned int d2 = offsets[1];
            stream << start;
            stream << '\t' << start + d1;
            stream << '\t' << start + d2 + d1 << '\t' << start + d2;
            break;
          }

        case 3:
          {
            const unsigned int d1 = offsets[0];
            const unsigned int d2 = offsets[1];
            const unsigned int d3 = offsets[2];
            stream << start;
            stream << '\t' << start + d1;
            stream << '\t' << start + d2 + d1 << '\t' << start + d2;
            stream << '\t' << start + d3 << '\t' << start + d3 + d1 << '\t'
                   << start + d3 + d2 + d1 << '\t' << start + d3 + d2;
            break;
          }

        default:
          DEAL_II_NOT_IMPLEMENTED();
      }
    stream << '\n';
  }



  UcdStream::UcdStream(std::ostream &out, const DataOutBase::UcdFlags &f)
    : StreamBase<DataOutBase::UcdFlags>(out, f)
  {}


  template <int dim>
  void
  UcdStream::write_point(const unsigned int index, const Point<dim> &p)
  {
    stream << index + 1 << "   ";
    // write out coordinates
    for (unsigned int i = 0; i < dim; ++i)
      stream << p[i] << ' ';
    // fill with zeroes
    for (unsigned int i = dim; i < 3; ++i)
      stream << "0 ";
    stream << '\n';
  }



  template <int dim>
  void
  UcdStream::write_cell(const unsigned int                   index,
                        const unsigned int                   start,
                        const std::array<unsigned int, dim> &offsets)
  {
    const auto nodes =
      DataOutBaseImplementation::set_node_numbers(start, offsets);

    // Write out all cells and remember that all indices must be shifted by one.
    stream << index + 1 << "\t0 " << ucd_cell_type[dim];
    for (unsigned int i = 0; i < nodes.size(); ++i)
      stream << '\t' << nodes[GeometryInfo<dim>::ucd_to_deal[i]] + 1;
    stream << '\n';
  }



  template <typename data>
  inline void
  UcdStream::write_dataset(const unsigned int       index,
                           const std::vector<data> &values)
  {
    stream << index + 1;
    for (unsigned int i = 0; i < values.size(); ++i)
      stream << '\t' << values[i];
    stream << '\n';
  }



  //----------------------------------------------------------------------//

  VtkStream::VtkStream(std::ostream &out, const DataOutBase::VtkFlags &f)
    : StreamBase<DataOutBase::VtkFlags>(out, f)
  {}


  template <int dim>
  void
  VtkStream::write_point(const unsigned int, const Point<dim> &p)
  {
    // write out coordinates
    stream << p;
    // fill with zeroes
    for (unsigned int i = dim; i < 3; ++i)
      stream << " 0";
    stream << '\n';
  }



  template <int dim>
  void
  VtkStream::write_cell(const unsigned int,
                        const unsigned int                   start,
                        const std::array<unsigned int, dim> &offsets)
  {
    stream << GeometryInfo<dim>::vertices_per_cell << '\t';

    switch (dim)
      {
        case 0:
          {
            stream << start;
            break;
          }

        case 1:
          {
            const unsigned int d1 = offsets[0];
            stream << start;
            stream << '\t' << start + d1;
            break;
          }

        case 2:
          {
            const unsigned int d1 = offsets[0];
            const unsigned int d2 = offsets[1];
            stream << start;
            stream << '\t' << start + d1;
            stream << '\t' << start + d2 + d1 << '\t' << start + d2;
            break;
          }

        case 3:
          {
            const unsigned int d1 = offsets[0];
            const unsigned int d2 = offsets[1];
            const unsigned int d3 = offsets[2];
            stream << start;
            stream << '\t' << start + d1;
            stream << '\t' << start + d2 + d1 << '\t' << start + d2;
            stream << '\t' << start + d3 << '\t' << start + d3 + d1 << '\t'
                   << start + d3 + d2 + d1 << '\t' << start + d3 + d2;
            break;
          }

        default:
          DEAL_II_NOT_IMPLEMENTED();
      }
    stream << '\n';
  }



  void
  VtkStream::write_cell_single(const unsigned int   index,
                               const unsigned int   start,
                               const unsigned int   n_points,
                               const ReferenceCell &reference_cell)
  {
    (void)index;

    static const std::array<unsigned int, 5> table = {{0, 1, 3, 2, 4}};

    stream << '\t' << n_points;
    for (unsigned int i = 0; i < n_points; ++i)
      stream << '\t'
             << start +
                  (reference_cell == ReferenceCells::Pyramid ? table[i] : i);
    stream << '\n';
  }

  template <int dim>
  void
  VtkStream::write_high_order_cell(const unsigned int           start,
                                   const std::vector<unsigned> &connectivity)
  {
    stream << connectivity.size();
    for (const auto &c : connectivity)
      stream << '\t' << start + c;
    stream << '\n';
  }
} // namespace



namespace DataOutBase
{
  const unsigned int Deal_II_IntermediateFlags::format_version = 4;


  template <int dim, int spacedim>
  const unsigned int Patch<dim, spacedim>::space_dim;


  template <int dim, int spacedim>
  const unsigned int Patch<dim, spacedim>::no_neighbor;


  template <int dim, int spacedim>
  Patch<dim, spacedim>::Patch()
    : patch_index(no_neighbor)
    , n_subdivisions(1)
    , points_are_available(false)
    , reference_cell(ReferenceCells::Invalid)
  // all the other data has a constructor of its own, except for the "neighbors"
  // field, which we set to invalid values.
  {
    for (const unsigned int i : GeometryInfo<dim>::face_indices())
      neighbors[i] = no_neighbor;

    AssertIndexRange(dim, spacedim + 1);
    Assert(spacedim <= 3, ExcNotImplemented());
  }



  template <int dim, int spacedim>
  bool
  Patch<dim, spacedim>::operator==(const Patch &patch) const
  {
    if (reference_cell != patch.reference_cell)
      return false;

    // TODO: make tolerance relative
    const double epsilon = 3e-16;
    for (const unsigned int i : GeometryInfo<dim>::vertex_indices())
      if (vertices[i].distance(patch.vertices[i]) > epsilon)
        return false;

    for (const unsigned int i : GeometryInfo<dim>::face_indices())
      if (neighbors[i] != patch.neighbors[i])
        return false;

    if (patch_index != patch.patch_index)
      return false;

    if (n_subdivisions != patch.n_subdivisions)
      return false;

    if (points_are_available != patch.points_are_available)
      return false;

    if (data.n_rows() != patch.data.n_rows())
      return false;

    if (data.n_cols() != patch.data.n_cols())
      return false;

    for (unsigned int i = 0; i < data.n_rows(); ++i)
      for (unsigned int j = 0; j < data.n_cols(); ++j)
        if (data[i][j] != patch.data[i][j])
          return false;

    return true;
  }



  template <int dim, int spacedim>
  std::size_t
  Patch<dim, spacedim>::memory_consumption() const
  {
    return (sizeof(vertices) / sizeof(vertices[0]) *
              MemoryConsumption::memory_consumption(vertices[0]) +
            sizeof(neighbors) / sizeof(neighbors[0]) *
              MemoryConsumption::memory_consumption(neighbors[0]) +
            MemoryConsumption::memory_consumption(patch_index) +
            MemoryConsumption::memory_consumption(n_subdivisions) +
            MemoryConsumption::memory_consumption(data) +
            MemoryConsumption::memory_consumption(points_are_available) +
            sizeof(reference_cell));
  }



  template <int dim, int spacedim>
  void
  Patch<dim, spacedim>::swap(Patch<dim, spacedim> &other_patch) noexcept
  {
    std::swap(vertices, other_patch.vertices);
    std::swap(neighbors, other_patch.neighbors);
    std::swap(patch_index, other_patch.patch_index);
    std::swap(n_subdivisions, other_patch.n_subdivisions);
    data.swap(other_patch.data);
    std::swap(points_are_available, other_patch.points_are_available);
    std::swap(reference_cell, other_patch.reference_cell);
  }



  template <int spacedim>
  const unsigned int Patch<0, spacedim>::space_dim;


  template <int spacedim>
  const unsigned int Patch<0, spacedim>::no_neighbor;


  template <int spacedim>
  unsigned int Patch<0, spacedim>::neighbors[1] = {
    Patch<0, spacedim>::no_neighbor};

  template <int spacedim>
  const unsigned int Patch<0, spacedim>::n_subdivisions = 1;

  template <int spacedim>
  const ReferenceCell Patch<0, spacedim>::reference_cell =
    ReferenceCells::Vertex;

  template <int spacedim>
  Patch<0, spacedim>::Patch()
    : patch_index(no_neighbor)
    , points_are_available(false)
  {
    Assert(spacedim <= 3, ExcNotImplemented());
  }



  template <int spacedim>
  bool
  Patch<0, spacedim>::operator==(const Patch &patch) const
  {
    const unsigned int dim = 0;

    // TODO: make tolerance relative
    const double epsilon = 3e-16;
    for (const unsigned int i : GeometryInfo<dim>::vertex_indices())
      if (vertices[i].distance(patch.vertices[i]) > epsilon)
        return false;

    if (patch_index != patch.patch_index)
      return false;

    if (points_are_available != patch.points_are_available)
      return false;

    if (data.n_rows() != patch.data.n_rows())
      return false;

    if (data.n_cols() != patch.data.n_cols())
      return false;

    for (unsigned int i = 0; i < data.n_rows(); ++i)
      for (unsigned int j = 0; j < data.n_cols(); ++j)
        if (data[i][j] != patch.data[i][j])
          return false;

    return true;
  }



  template <int spacedim>
  std::size_t
  Patch<0, spacedim>::memory_consumption() const
  {
    return (sizeof(vertices) / sizeof(vertices[0]) *
              MemoryConsumption::memory_consumption(vertices[0]) +
            MemoryConsumption::memory_consumption(data) +
            MemoryConsumption::memory_consumption(points_are_available));
  }



  template <int spacedim>
  void
  Patch<0, spacedim>::swap(Patch<0, spacedim> &other_patch) noexcept
  {
    std::swap(vertices, other_patch.vertices);
    std::swap(patch_index, other_patch.patch_index);
    data.swap(other_patch.data);
    std::swap(points_are_available, other_patch.points_are_available);
  }



  UcdFlags::UcdFlags(const bool write_preamble)
    : write_preamble(write_preamble)
  {}



  GnuplotFlags::GnuplotFlags()
  {
    space_dimension_labels.emplace_back("x");
    space_dimension_labels.emplace_back("y");
    space_dimension_labels.emplace_back("z");
  }



  GnuplotFlags::GnuplotFlags(const std::vector<std::string> &labels)
    : space_dimension_labels(labels)
  {}



  std::size_t
  GnuplotFlags::memory_consumption() const
  {
    return MemoryConsumption::memory_consumption(space_dimension_labels);
  }



  PovrayFlags::PovrayFlags(const bool smooth,
                           const bool bicubic_patch,
                           const bool external_data)
    : smooth(smooth)
    , bicubic_patch(bicubic_patch)
    , external_data(external_data)
  {}


  DataOutFilterFlags::DataOutFilterFlags(const bool filter_duplicate_vertices,
                                         const bool xdmf_hdf5_output)
    : filter_duplicate_vertices(filter_duplicate_vertices)
    , xdmf_hdf5_output(xdmf_hdf5_output)
  {}


  void
  DataOutFilterFlags::declare_parameters(ParameterHandler &prm)
  {
    prm.declare_entry(
      "Filter duplicate vertices",
      "false",
      Patterns::Bool(),
      "Whether to remove duplicate vertex values. deal.II duplicates "
      "vertices once for each adjacent cell so that it can output "
      "discontinuous quantities for which there may be more than one "
      "value for each vertex position. Setting this flag to "
      "'true' will merge all of these values by selecting a "
      "random one and outputting this as 'the' value for the vertex. "
      "As long as the data to be output corresponds to continuous "
      "fields, merging vertices has no effect. On the other hand, "
      "if the data to be output corresponds to discontinuous fields "
      "(either because you are using a discontinuous finite element, "
      "or because you are using a DataPostprocessor that yields "
      "discontinuous data, or because the data to be output has been "
      "produced by entirely different means), then the data in the "
      "output file no longer faithfully represents the underlying data "
      "because the discontinuous field has been replaced by a "
      "continuous one. Note also that the filtering can not occur "
      "on processor boundaries. Thus, a filtered discontinuous field "
      "looks like a continuous field inside of a subdomain, "
      "but like a discontinuous field at the subdomain boundary."
      "\n\n"
      "In any case, filtering results in drastically smaller output "
      "files (smaller by about a factor of 2^dim).");
    prm.declare_entry(
      "XDMF HDF5 output",
      "false",
      Patterns::Bool(),
      "Whether the data will be used in an XDMF/HDF5 combination.");
  }



  void
  DataOutFilterFlags::parse_parameters(const ParameterHandler &prm)
  {
    filter_duplicate_vertices = prm.get_bool("Filter duplicate vertices");
    xdmf_hdf5_output          = prm.get_bool("XDMF HDF5 output");
  }



  DXFlags::DXFlags(const bool write_neighbors,
                   const bool int_binary,
                   const bool coordinates_binary,
                   const bool data_binary)
    : write_neighbors(write_neighbors)
    , int_binary(int_binary)
    , coordinates_binary(coordinates_binary)
    , data_binary(data_binary)
    , data_double(false)
  {}


  void
  DXFlags::declare_parameters(ParameterHandler &prm)
  {
    prm.declare_entry("Write neighbors",
                      "true",
                      Patterns::Bool(),
                      "A boolean field indicating whether neighborship "
                      "information between cells is to be written to the "
                      "OpenDX output file");
    prm.declare_entry("Integer format",
                      "ascii",
                      Patterns::Selection("ascii|32|64"),
                      "Output format of integer numbers, which is "
                      "either a text representation (ascii) or binary integer "
                      "values of 32 or 64 bits length");
    prm.declare_entry("Coordinates format",
                      "ascii",
                      Patterns::Selection("ascii|32|64"),
                      "Output format of vertex coordinates, which is "
                      "either a text representation (ascii) or binary "
                      "floating point values of 32 or 64 bits length");
    prm.declare_entry("Data format",
                      "ascii",
                      Patterns::Selection("ascii|32|64"),
                      "Output format of data values, which is "
                      "either a text representation (ascii) or binary "
                      "floating point values of 32 or 64 bits length");
  }



  void
  DXFlags::parse_parameters(const ParameterHandler &prm)
  {
    write_neighbors = prm.get_bool("Write neighbors");
    // TODO:[GK] Read the new  parameters
  }



  void
  UcdFlags::declare_parameters(ParameterHandler &prm)
  {
    prm.declare_entry("Write preamble",
                      "true",
                      Patterns::Bool(),
                      "A flag indicating whether a comment should be "
                      "written to the beginning of the output file "
                      "indicating date and time of creation as well "
                      "as the creating program");
  }



  void
  UcdFlags::parse_parameters(const ParameterHandler &prm)
  {
    write_preamble = prm.get_bool("Write preamble");
  }



  SvgFlags::SvgFlags(const unsigned int height_vector,
                     const int          azimuth_angle,
                     const int          polar_angle,
                     const unsigned int line_thickness,
                     const bool         margin,
                     const bool         draw_colorbar)
    : height(4000)
    , width(0)
    , height_vector(height_vector)
    , azimuth_angle(azimuth_angle)
    , polar_angle(polar_angle)
    , line_thickness(line_thickness)
    , margin(margin)
    , draw_colorbar(draw_colorbar)
  {}



  void
  PovrayFlags::declare_parameters(ParameterHandler &prm)
  {
    prm.declare_entry("Use smooth triangles",
                      "false",
                      Patterns::Bool(),
                      "A flag indicating whether POVRAY should use smoothed "
                      "triangles instead of the usual ones");
    prm.declare_entry("Use bicubic patches",
                      "false",
                      Patterns::Bool(),
                      "Whether POVRAY should use bicubic patches");
    prm.declare_entry("Include external file",
                      "true",
                      Patterns::Bool(),
                      "Whether camera and lighting information should "
                      "be put into an external file \"data.inc\" or into "
                      "the POVRAY input file");
  }



  void
  PovrayFlags::parse_parameters(const ParameterHandler &prm)
  {
    smooth        = prm.get_bool("Use smooth triangles");
    bicubic_patch = prm.get_bool("Use bicubic patches");
    external_data = prm.get_bool("Include external file");
  }



  EpsFlags::EpsFlags(const unsigned int  height_vector,
                     const unsigned int  color_vector,
                     const SizeType      size_type,
                     const unsigned int  size,
                     const double        line_width,
                     const double        azimut_angle,
                     const double        turn_angle,
                     const double        z_scaling,
                     const bool          draw_mesh,
                     const bool          draw_cells,
                     const bool          shade_cells,
                     const ColorFunction color_function)
    : height_vector(height_vector)
    , color_vector(color_vector)
    , size_type(size_type)
    , size(size)
    , line_width(line_width)
    , azimut_angle(azimut_angle)
    , turn_angle(turn_angle)
    , z_scaling(z_scaling)
    , draw_mesh(draw_mesh)
    , draw_cells(draw_cells)
    , shade_cells(shade_cells)
    , color_function(color_function)
  {}



  EpsFlags::RgbValues
  EpsFlags::default_color_function(const double x,
                                   const double xmin,
                                   const double xmax)
  {
    RgbValues rgb_values = {0, 0, 0};

    // A difficult color scale:
    //     xmin          = black  [1]
    // 3/4*xmin+1/4*xmax = blue   [2]
    // 1/2*xmin+1/2*xmax = green  (3)
    // 1/4*xmin+3/4*xmax = red    (4)
    //              xmax = white  (5)
    // Makes the following color functions:
    //
    // red      green    blue
    //       __
    //      /      /\  /  /\    /
    // ____/    __/  \/  /  \__/

    //     { 0                                [1] - (3)
    // r = { ( 4*x-2*xmin+2*xmax)/(xmax-xmin) (3) - (4)
    //     { 1                                (4) - (5)
    //
    //     { 0                                [1] - [2]
    // g = { ( 4*x-3*xmin-  xmax)/(xmax-xmin) [2] - (3)
    //     { (-4*x+  xmin+3*xmax)/(xmax-xmin) (3) - (4)
    //     { ( 4*x-  xmin-3*xmax)/(xmax-xmin) (4) - (5)
    //
    //     { ( 4*x-4*xmin       )/(xmax-xmin) [1] - [2]
    // b = { (-4*x+2*xmin+2*xmax)/(xmax-xmin) [2] - (3)
    //     { 0                                (3) - (4)
    //     { ( 4*x-  xmin-3*xmax)/(xmax-xmin) (4) - (5)

    double sum    = xmax + xmin;
    double sum13  = xmin + 3 * xmax;
    double sum22  = 2 * xmin + 2 * xmax;
    double sum31  = 3 * xmin + xmax;
    double dif    = xmax - xmin;
    double rezdif = 1.0 / dif;

    int where;

    if (x < (sum31) / 4)
      where = 0;
    else if (x < (sum22) / 4)
      where = 1;
    else if (x < (sum13) / 4)
      where = 2;
    else
      where = 3;

    if (dif != 0)
      {
        switch (where)
          {
            case 0:
              rgb_values.red   = 0;
              rgb_values.green = 0;
              rgb_values.blue  = (x - xmin) * 4. * rezdif;
              break;
            case 1:
              rgb_values.red   = 0;
              rgb_values.green = (4 * x - 3 * xmin - xmax) * rezdif;
              rgb_values.blue  = (sum22 - 4. * x) * rezdif;
              break;
            case 2:
              rgb_values.red   = (4 * x - 2 * sum) * rezdif;
              rgb_values.green = (xmin + 3 * xmax - 4 * x) * rezdif;
              rgb_values.blue  = 0;
              break;
            case 3:
              rgb_values.red   = 1;
              rgb_values.green = (4 * x - xmin - 3 * xmax) * rezdif;
              rgb_values.blue  = (4. * x - sum13) * rezdif;
              break;
            default:
              break;
          }
      }
    else // White
      rgb_values.red = rgb_values.green = rgb_values.blue = 1;

    return rgb_values;
  }



  EpsFlags::RgbValues
  EpsFlags::grey_scale_color_function(const double x,
                                      const double xmin,
                                      const double xmax)
  {
    EpsFlags::RgbValues rgb_values;
    rgb_values.red = rgb_values.blue = rgb_values.green =
      (x - xmin) / (xmax - xmin);
    return rgb_values;
  }



  EpsFlags::RgbValues
  EpsFlags::reverse_grey_scale_color_function(const double x,
                                              const double xmin,
                                              const double xmax)
  {
    EpsFlags::RgbValues rgb_values;
    rgb_values.red = rgb_values.blue = rgb_values.green =
      1 - (x - xmin) / (xmax - xmin);
    return rgb_values;
  }



  void
  EpsFlags::declare_parameters(ParameterHandler &prm)
  {
    prm.declare_entry("Index of vector for height",
                      "0",
                      Patterns::Integer(),
                      "Number of the input vector that is to be used to "
                      "generate height information");
    prm.declare_entry("Index of vector for color",
                      "0",
                      Patterns::Integer(),
                      "Number of the input vector that is to be used to "
                      "generate color information");
    prm.declare_entry("Scale to width or height",
                      "width",
                      Patterns::Selection("width|height"),
                      "Whether width or height should be scaled to match "
                      "the given size");
    prm.declare_entry("Size (width or height) in eps units",
                      "300",
                      Patterns::Integer(),
                      "The size (width or height) to which the eps output "
                      "file is to be scaled");
    prm.declare_entry("Line widths in eps units",
                      "0.5",
                      Patterns::Double(),
                      "The width in which the postscript renderer is to "
                      "plot lines");
    prm.declare_entry("Azimut angle",
                      "60",
                      Patterns::Double(0, 180),
                      "Angle of the viewing position against the vertical "
                      "axis");
    prm.declare_entry("Turn angle",
                      "30",
                      Patterns::Double(0, 360),
                      "Angle of the viewing direction against the y-axis");
    prm.declare_entry("Scaling for z-axis",
                      "1",
                      Patterns::Double(),
                      "Scaling for the z-direction relative to the scaling "
                      "used in x- and y-directions");
    prm.declare_entry("Draw mesh lines",
                      "true",
                      Patterns::Bool(),
                      "Whether the mesh lines, or only the surface should be "
                      "drawn");
    prm.declare_entry("Fill interior of cells",
                      "true",
                      Patterns::Bool(),
                      "Whether only the mesh lines, or also the interior of "
                      "cells should be plotted. If this flag is false, then "
                      "one can see through the mesh");
    prm.declare_entry("Color shading of interior of cells",
                      "true",
                      Patterns::Bool(),
                      "Whether the interior of cells shall be shaded");
    prm.declare_entry("Color function",
                      "default",
                      Patterns::Selection(
                        "default|grey scale|reverse grey scale"),
                      "Name of a color function used to colorize mesh lines "
                      "and/or cell interiors");
  }



  void
  EpsFlags::parse_parameters(const ParameterHandler &prm)
  {
    height_vector = prm.get_integer("Index of vector for height");
    color_vector  = prm.get_integer("Index of vector for color");
    if (prm.get("Scale to width or height") == "width")
      size_type = width;
    else
      size_type = height;
    size         = prm.get_integer("Size (width or height) in eps units");
    line_width   = prm.get_double("Line widths in eps units");
    azimut_angle = prm.get_double("Azimut angle");
    turn_angle   = prm.get_double("Turn angle");
    z_scaling    = prm.get_double("Scaling for z-axis");
    draw_mesh    = prm.get_bool("Draw mesh lines");
    draw_cells   = prm.get_bool("Fill interior of cells");
    shade_cells  = prm.get_bool("Color shading of interior of cells");
    if (prm.get("Color function") == "default")
      color_function = &default_color_function;
    else if (prm.get("Color function") == "grey scale")
      color_function = &grey_scale_color_function;
    else if (prm.get("Color function") == "reverse grey scale")
      color_function = &reverse_grey_scale_color_function;
    else
      // we shouldn't get here, since the parameter object should already have
      // checked that the given value is valid
      DEAL_II_ASSERT_UNREACHABLE();
  }


  Hdf5Flags::Hdf5Flags(const CompressionLevel compression_level)
    : compression_level(compression_level)
  {}


  TecplotFlags::TecplotFlags(const char *zone_name, const double solution_time)
    : zone_name(zone_name)
    , solution_time(solution_time)
  {}



  std::size_t
  TecplotFlags::memory_consumption() const
  {
    return sizeof(*this) + MemoryConsumption::memory_consumption(zone_name);
  }



  VtkFlags::VtkFlags(const double           time,
                     const unsigned int     cycle,
                     const bool             print_date_and_time,
                     const CompressionLevel compression_level,
                     const bool             write_higher_order_cells,
                     const std::map<std::string, std::string> &physical_units)
    : time(time)
    , cycle(cycle)
    , print_date_and_time(print_date_and_time)
    , compression_level(compression_level)
    , write_higher_order_cells(write_higher_order_cells)
    , physical_units(physical_units)
  {}



  OutputFormat
  parse_output_format(const std::string &format_name)
  {
    if (format_name == "none")
      return none;

    if (format_name == "dx")
      return dx;

    if (format_name == "ucd")
      return ucd;

    if (format_name == "gnuplot")
      return gnuplot;

    if (format_name == "povray")
      return povray;

    if (format_name == "eps")
      return eps;

    if (format_name == "gmv")
      return gmv;

    if (format_name == "tecplot")
      return tecplot;

    if (format_name == "vtk")
      return vtk;

    if (format_name == "vtu")
      return vtu;

    if (format_name == "deal.II intermediate")
      return deal_II_intermediate;

    if (format_name == "hdf5")
      return hdf5;

    AssertThrow(false,
                ExcMessage("The given file format name is not recognized: <" +
                           format_name + ">"));

    // return something invalid
    return OutputFormat(-1);
  }



  std::string
  get_output_format_names()
  {
    return "none|dx|ucd|gnuplot|povray|eps|gmv|tecplot|vtk|vtu|hdf5|svg|deal.II intermediate";
  }



  std::string
  default_suffix(const OutputFormat output_format)
  {
    switch (output_format)
      {
        case none:
          return "";
        case dx:
          return ".dx";
        case ucd:
          return ".inp";
        case gnuplot:
          return ".gnuplot";
        case povray:
          return ".pov";
        case eps:
          return ".eps";
        case gmv:
          return ".gmv";
        case tecplot:
          return ".dat";
        case vtk:
          return ".vtk";
        case vtu:
          return ".vtu";
        case deal_II_intermediate:
          return ".d2";
        case hdf5:
          return ".h5";
        case svg:
          return ".svg";
        default:
          DEAL_II_NOT_IMPLEMENTED();
          return "";
      }
  }


  //----------------------------------------------------------------------//


  /**
   * Obtain the positions of all nodes referenced in the patches given as
   * argument.
   */
  template <int dim, int spacedim>
  std::vector<Point<spacedim>>
  get_node_positions(const std::vector<Patch<dim, spacedim>> &patches)
  {
    Assert(dim <= 3, ExcNotImplemented());
    static const std::array<unsigned int, 5> table = {{0, 1, 3, 2, 4}};

    std::vector<Point<spacedim>> node_positions;
    for (const auto &patch : patches)
      {
        // special treatment of non-hypercube cells
        if (patch.reference_cell != ReferenceCells::get_hypercube<dim>())
          {
            for (unsigned int point_no = 0; point_no < patch.data.n_cols();
                 ++point_no)
              node_positions.emplace_back(get_node_location(
                patch,
                (patch.reference_cell == ReferenceCells::Pyramid ?
                   table[point_no] :
                   point_no)));
          }
        else
          {
            const unsigned int n_subdivisions = patch.n_subdivisions;
            const unsigned int n              = n_subdivisions + 1;

            switch (dim)
              {
                case 0:
                  node_positions.emplace_back(
                    get_equispaced_location(patch, {}, n_subdivisions));
                  break;
                case 1:
                  for (unsigned int i1 = 0; i1 < n; ++i1)
                    node_positions.emplace_back(
                      get_equispaced_location(patch, {i1}, n_subdivisions));
                  break;
                case 2:
                  for (unsigned int i2 = 0; i2 < n; ++i2)
                    for (unsigned int i1 = 0; i1 < n; ++i1)
                      node_positions.emplace_back(get_equispaced_location(
                        patch, {i1, i2}, n_subdivisions));
                  break;
                case 3:
                  for (unsigned int i3 = 0; i3 < n; ++i3)
                    for (unsigned int i2 = 0; i2 < n; ++i2)
                      for (unsigned int i1 = 0; i1 < n; ++i1)
                        node_positions.emplace_back(get_equispaced_location(
                          patch, {i1, i2, i3}, n_subdivisions));
                  break;

                default:
                  DEAL_II_ASSERT_UNREACHABLE();
              }
          }
      }

    return node_positions;
  }


  template <int dim, int spacedim, typename StreamType>
  void
  write_nodes(const std::vector<Patch<dim, spacedim>> &patches, StreamType &out)
  {
    // Obtain the node locations, and then output them via the given stream
    // object
    const std::vector<Point<spacedim>> node_positions =
      get_node_positions(patches);

    int count = 0;
    for (const auto &node : node_positions)
      out.write_point(count++, node);
    out.flush_points();
  }



  template <int dim, int spacedim, typename StreamType>
  void
  write_cells(const std::vector<Patch<dim, spacedim>> &patches, StreamType &out)
  {
    Assert(dim <= 3, ExcNotImplemented());
    unsigned int count                 = 0;
    unsigned int first_vertex_of_patch = 0;
    for (const auto &patch : patches)
      {
        // special treatment of simplices since they are not subdivided
        if (patch.reference_cell != ReferenceCells::get_hypercube<dim>())
          {
            out.write_cell_single(count++,
                                  first_vertex_of_patch,
                                  patch.data.n_cols(),
                                  patch.reference_cell);
            first_vertex_of_patch += patch.data.n_cols();
          }
        else // hypercube cell
          {
            const unsigned int n_subdivisions = patch.n_subdivisions;
            const unsigned int n              = n_subdivisions + 1;

            switch (dim)
              {
                case 0:
                  {
                    const unsigned int offset = first_vertex_of_patch;
                    out.template write_cell<0>(count++, offset, {});
                    break;
                  }

                case 1:
                  {
                    constexpr unsigned int d1 = 1;

                    for (unsigned int i1 = 0; i1 < n_subdivisions; ++i1)
                      {
                        const unsigned int offset =
                          first_vertex_of_patch + i1 * d1;
                        out.template write_cell<1>(count++, offset, {{d1}});
                      }

                    break;
                  }

                case 2:
                  {
                    constexpr unsigned int d1 = 1;
                    const unsigned int     d2 = n;

                    for (unsigned int i2 = 0; i2 < n_subdivisions; ++i2)
                      for (unsigned int i1 = 0; i1 < n_subdivisions; ++i1)
                        {
                          const unsigned int offset =
                            first_vertex_of_patch + i2 * d2 + i1 * d1;
                          out.template write_cell<2>(count++,
                                                     offset,
                                                     {{d1, d2}});
                        }

                    break;
                  }

                case 3:
                  {
                    constexpr unsigned int d1 = 1;
                    const unsigned int     d2 = n;
                    const unsigned int     d3 = n * n;

                    for (unsigned int i3 = 0; i3 < n_subdivisions; ++i3)
                      for (unsigned int i2 = 0; i2 < n_subdivisions; ++i2)
                        for (unsigned int i1 = 0; i1 < n_subdivisions; ++i1)
                          {
                            const unsigned int offset = first_vertex_of_patch +
                                                        i3 * d3 + i2 * d2 +
                                                        i1 * d1;
                            out.template write_cell<3>(count++,
                                                       offset,
                                                       {{d1, d2, d3}});
                          }

                    break;
                  }
                default:
                  DEAL_II_NOT_IMPLEMENTED();
              }

            // Update the number of the first vertex of this patch
            first_vertex_of_patch +=
              Utilities::fixed_power<dim>(n_subdivisions + 1);
          }
      }

    out.flush_cells();
  }



  template <int dim, int spacedim, typename StreamType>
  void
  write_high_order_cells(const std::vector<Patch<dim, spacedim>> &patches,
                         StreamType                              &out,
                         const bool                               legacy_format)
  {
    unsigned int first_vertex_of_patch = 0;
    // Array to hold all the node numbers of a cell
    std::vector<unsigned> connectivity;

    for (const auto &patch : patches)
      {
        if (patch.reference_cell != ReferenceCells::get_hypercube<dim>())
          {
            connectivity.resize(patch.data.n_cols());

            for (unsigned int i = 0; i < patch.data.n_cols(); ++i)
              connectivity[i] = i;

            out.template write_high_order_cell<dim>(first_vertex_of_patch,
                                                    connectivity);

            first_vertex_of_patch += patch.data.n_cols();
          }
        else
          {
            const unsigned int n_subdivisions = patch.n_subdivisions;
            const unsigned int n              = n_subdivisions + 1;

            connectivity.resize(Utilities::fixed_power<dim>(n));

            switch (dim)
              {
                case 0:
                  {
                    Assert(false,
                           ExcMessage("Point-like cells should not be possible "
                                      "when writing higher-order cells."));
                    break;
                  }
                case 1:
                  {
                    for (unsigned int i1 = 0; i1 < n_subdivisions + 1; ++i1)
                      {
                        const unsigned int local_index = i1;
                        const unsigned int connectivity_index =
                          patch.reference_cell
                            .template vtk_lexicographic_to_node_index<1>(
                              {{i1}}, {{n_subdivisions}}, legacy_format);
                        connectivity[connectivity_index] = local_index;
                      }

                    break;
                  }
                case 2:
                  {
                    for (unsigned int i2 = 0; i2 < n_subdivisions + 1; ++i2)
                      for (unsigned int i1 = 0; i1 < n_subdivisions + 1; ++i1)
                        {
                          const unsigned int local_index = i2 * n + i1;
                          const unsigned int connectivity_index =
                            patch.reference_cell
                              .template vtk_lexicographic_to_node_index<2>(
                                {{i1, i2}},
                                {{n_subdivisions, n_subdivisions}},
                                legacy_format);
                          connectivity[connectivity_index] = local_index;
                        }

                    break;
                  }
                case 3:
                  {
                    for (unsigned int i3 = 0; i3 < n_subdivisions + 1; ++i3)
                      for (unsigned int i2 = 0; i2 < n_subdivisions + 1; ++i2)
                        for (unsigned int i1 = 0; i1 < n_subdivisions + 1; ++i1)
                          {
                            const unsigned int local_index =
                              i3 * n * n + i2 * n + i1;
                            const unsigned int connectivity_index =
                              patch.reference_cell
                                .template vtk_lexicographic_to_node_index<3>(
                                  {{i1, i2, i3}},
                                  {{n_subdivisions,
                                    n_subdivisions,
                                    n_subdivisions}},
                                  legacy_format);
                            connectivity[connectivity_index] = local_index;
                          }

                    break;
                  }
                default:
                  DEAL_II_NOT_IMPLEMENTED();
              }

            // Having so set up the 'connectivity' data structure,
            // output it:
            out.template write_high_order_cell<dim>(first_vertex_of_patch,
                                                    connectivity);

            // Finally update the number of the first vertex of this patch
            first_vertex_of_patch += Utilities::fixed_power<dim>(n);
          }
      }

    out.flush_cells();
  }


  template <int dim, int spacedim, typename StreamType>
  void
  write_data(const std::vector<Patch<dim, spacedim>> &patches,
             unsigned int                             n_data_sets,
             const bool                               double_precision,
             StreamType                              &out)
  {
    Assert(dim <= 3, ExcNotImplemented());
    unsigned int count = 0;

    for (const auto &patch : patches)
      {
        const unsigned int n_subdivisions = patch.n_subdivisions;
        const unsigned int n              = n_subdivisions + 1;
        // Length of loops in all dimensions
        Assert((patch.data.n_rows() == n_data_sets &&
                !patch.points_are_available) ||
                 (patch.data.n_rows() == n_data_sets + spacedim &&
                  patch.points_are_available),
               ExcDimensionMismatch(patch.points_are_available ?
                                      (n_data_sets + spacedim) :
                                      n_data_sets,
                                    patch.data.n_rows()));
        Assert(patch.data.n_cols() == Utilities::fixed_power<dim>(n),
               ExcInvalidDatasetSize(patch.data.n_cols(), n));

        std::vector<float>  floats(n_data_sets);
        std::vector<double> doubles(n_data_sets);

        // Data is already in lexicographic ordering
        for (unsigned int i = 0; i < Utilities::fixed_power<dim>(n);
             ++i, ++count)
          if (double_precision)
            {
              for (unsigned int data_set = 0; data_set < n_data_sets;
                   ++data_set)
                doubles[data_set] = patch.data(data_set, i);
              out.write_dataset(count, doubles);
            }
          else
            {
              for (unsigned int data_set = 0; data_set < n_data_sets;
                   ++data_set)
                floats[data_set] = patch.data(data_set, i);
              out.write_dataset(count, floats);
            }
      }
  }



  namespace
  {
    /**
     * This function projects a three-dimensional point (Point<3> point)
     * onto a two-dimensional image plane, specified by the position of
     * the camera viewing system (Point<3> camera_position), camera
     * direction (Point<3> camera_position), camera horizontal (Point<3>
     * camera_horizontal, necessary for the correct alignment of the
     * later images), and the focus of the camera (float camera_focus).
     */
    Point<2>
    svg_project_point(Point<3> point,
                      Point<3> camera_position,
                      Point<3> camera_direction,
                      Point<3> camera_horizontal,
                      float    camera_focus)
    {
      Point<3> camera_vertical;
      camera_vertical[0] = camera_horizontal[1] * camera_direction[2] -
                           camera_horizontal[2] * camera_direction[1];
      camera_vertical[1] = camera_horizontal[2] * camera_direction[0] -
                           camera_horizontal[0] * camera_direction[2];
      camera_vertical[2] = camera_horizontal[0] * camera_direction[1] -
                           camera_horizontal[1] * camera_direction[0];

      float phi;
      phi = camera_focus;
      phi /= (point[0] - camera_position[0]) * camera_direction[0] +
             (point[1] - camera_position[1]) * camera_direction[1] +
             (point[2] - camera_position[2]) * camera_direction[2];

      Point<3> projection;
      projection[0] =
        camera_position[0] + phi * (point[0] - camera_position[0]);
      projection[1] =
        camera_position[1] + phi * (point[1] - camera_position[1]);
      projection[2] =
        camera_position[2] + phi * (point[2] - camera_position[2]);

      Point<2> projection_decomposition;
      projection_decomposition[0] = (projection[0] - camera_position[0] -
                                     camera_focus * camera_direction[0]) *
                                    camera_horizontal[0];
      projection_decomposition[0] += (projection[1] - camera_position[1] -
                                      camera_focus * camera_direction[1]) *
                                     camera_horizontal[1];
      projection_decomposition[0] += (projection[2] - camera_position[2] -
                                      camera_focus * camera_direction[2]) *
                                     camera_horizontal[2];

      projection_decomposition[1] = (projection[0] - camera_position[0] -
                                     camera_focus * camera_direction[0]) *
                                    camera_vertical[0];
      projection_decomposition[1] += (projection[1] - camera_position[1] -
                                      camera_focus * camera_direction[1]) *
                                     camera_vertical[1];
      projection_decomposition[1] += (projection[2] - camera_position[2] -
                                      camera_focus * camera_direction[2]) *
                                     camera_vertical[2];

      return projection_decomposition;
    }


    /**
     * Function to compute the gradient parameters for a triangle with given
     * values for the vertices.
     */
    Point<6>
    svg_get_gradient_parameters(Point<3> points[])
    {
      Point<3> v_min, v_max, v_inter;

      // Use the Bubblesort algorithm to sort the points with respect to the
      // third coordinate
      for (int i = 0; i < 2; ++i)
        {
          for (int j = 0; j < 2 - i; ++j)
            {
              if (points[j][2] > points[j + 1][2])
                {
                  Point<3> temp = points[j];
                  points[j]     = points[j + 1];
                  points[j + 1] = temp;
                }
            }
        }

      // save the related three-dimensional vectors v_min, v_inter, and v_max
      v_min   = points[0];
      v_inter = points[1];
      v_max   = points[2];

      Point<2> A[2];
      Point<2> b, gradient;

      // determine the plane offset c
      A[0][0] = v_max[0] - v_min[0];
      A[0][1] = v_inter[0] - v_min[0];
      A[1][0] = v_max[1] - v_min[1];
      A[1][1] = v_inter[1] - v_min[1];

      b[0] = -v_min[0];
      b[1] = -v_min[1];

      double x, sum;
      bool   col_change = false;

      if (A[0][0] == 0)
        {
          col_change = true;

          A[0][0] = A[0][1];
          A[0][1] = 0;

          double temp = A[1][0];
          A[1][0]     = A[1][1];
          A[1][1]     = temp;
        }

      for (unsigned int k = 0; k < 1; ++k)
        {
          for (unsigned int i = k + 1; i < 2; ++i)
            {
              x = A[i][k] / A[k][k];

              for (unsigned int j = k + 1; j < 2; ++j)
                A[i][j] = A[i][j] - A[k][j] * x;

              b[i] = b[i] - b[k] * x;
            }
        }

      b[1] = b[1] / A[1][1];

      for (int i = 0; i >= 0; i--)
        {
          sum = b[i];

          for (unsigned int j = i + 1; j < 2; ++j)
            sum = sum - A[i][j] * b[j];

          b[i] = sum / A[i][i];
        }

      if (col_change)
        {
          double temp = b[0];
          b[0]        = b[1];
          b[1]        = temp;
        }

      double c = b[0] * (v_max[2] - v_min[2]) + b[1] * (v_inter[2] - v_min[2]) +
                 v_min[2];

      // Determine the first entry of the gradient (phi, cf. documentation)
      A[0][0] = v_max[0] - v_min[0];
      A[0][1] = v_inter[0] - v_min[0];
      A[1][0] = v_max[1] - v_min[1];
      A[1][1] = v_inter[1] - v_min[1];

      b[0] = 1.0 - v_min[0];
      b[1] = -v_min[1];

      col_change = false;

      if (A[0][0] == 0)
        {
          col_change = true;

          A[0][0] = A[0][1];
          A[0][1] = 0;

          double temp = A[1][0];
          A[1][0]     = A[1][1];
          A[1][1]     = temp;
        }

      for (unsigned int k = 0; k < 1; ++k)
        {
          for (unsigned int i = k + 1; i < 2; ++i)
            {
              x = A[i][k] / A[k][k];

              for (unsigned int j = k + 1; j < 2; ++j)
                A[i][j] = A[i][j] - A[k][j] * x;

              b[i] = b[i] - b[k] * x;
            }
        }

      b[1] = b[1] / A[1][1];

      for (int i = 0; i >= 0; i--)
        {
          sum = b[i];

          for (unsigned int j = i + 1; j < 2; ++j)
            sum = sum - A[i][j] * b[j];

          b[i] = sum / A[i][i];
        }

      if (col_change)
        {
          double temp = b[0];
          b[0]        = b[1];
          b[1]        = temp;
        }

      gradient[0] = b[0] * (v_max[2] - v_min[2]) +
                    b[1] * (v_inter[2] - v_min[2]) - c + v_min[2];

      // determine the second entry of the gradient
      A[0][0] = v_max[0] - v_min[0];
      A[0][1] = v_inter[0] - v_min[0];
      A[1][0] = v_max[1] - v_min[1];
      A[1][1] = v_inter[1] - v_min[1];

      b[0] = -v_min[0];
      b[1] = 1.0 - v_min[1];

      col_change = false;

      if (A[0][0] == 0)
        {
          col_change = true;

          A[0][0] = A[0][1];
          A[0][1] = 0;

          double temp = A[1][0];
          A[1][0]     = A[1][1];
          A[1][1]     = temp;
        }

      for (unsigned int k = 0; k < 1; ++k)
        {
          for (unsigned int i = k + 1; i < 2; ++i)
            {
              x = A[i][k] / A[k][k];

              for (unsigned int j = k + 1; j < 2; ++j)
                A[i][j] = A[i][j] - A[k][j] * x;

              b[i] = b[i] - b[k] * x;
            }
        }

      b[1] = b[1] / A[1][1];

      for (int i = 0; i >= 0; i--)
        {
          sum = b[i];

          for (unsigned int j = i + 1; j < 2; ++j)
            sum = sum - A[i][j] * b[j];

          b[i] = sum / A[i][i];
        }

      if (col_change)
        {
          double temp = b[0];
          b[0]        = b[1];
          b[1]        = temp;
        }

      gradient[1] = b[0] * (v_max[2] - v_min[2]) +
                    b[1] * (v_inter[2] - v_min[2]) - c + v_min[2];

      // normalize the gradient
      gradient /= gradient.norm();

      const double lambda = -gradient[0] * (v_min[0] - v_max[0]) -
                            gradient[1] * (v_min[1] - v_max[1]);

      Point<6> gradient_parameters;

      gradient_parameters[0] = v_min[0];
      gradient_parameters[1] = v_min[1];

      gradient_parameters[2] = v_min[0] + lambda * gradient[0];
      gradient_parameters[3] = v_min[1] + lambda * gradient[1];

      gradient_parameters[4] = v_min[2];
      gradient_parameters[5] = v_max[2];

      return gradient_parameters;
    }
  } // namespace



  template <int dim, int spacedim>
  void
  write_ucd(
    const std::vector<Patch<dim, spacedim>> &patches,
    const std::vector<std::string>          &data_names,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>> &,
    const UcdFlags &flags,
    std::ostream   &out)
  {
    // Note that while in theory dim==0 should be implemented, this is not
    // tested, therefore currently not allowed.
    AssertThrow(dim > 0, ExcNotImplemented());

    AssertThrow(out.fail() == false, ExcIO());

#ifndef DEAL_II_WITH_MPI
    // verify that there are indeed patches to be written out. most of the
    // times, people just forget to call build_patches when there are no
    // patches, so a warning is in order. that said, the assertion is disabled
    // if we support MPI since then it can happen that on the coarsest mesh, a
    // processor simply has no cells it actually owns, and in that case it is
    // legit if there are no patches
    Assert(patches.size() > 0, ExcNoPatches());
#else
    if (patches.empty())
      return;
#endif

    const unsigned int n_data_sets = data_names.size();

    UcdStream ucd_out(out, flags);

    // first count the number of cells and cells for later use
    unsigned int n_nodes;
    unsigned int n_cells;
    std::tie(n_nodes, n_cells) = count_nodes_and_cells(patches);
    //---------------------
    // preamble
    if (flags.write_preamble)
      {
        out
          << "# This file was generated by the deal.II library." << '\n'
          << "# Date =  " << Utilities::System::get_date() << '\n'
          << "# Time =  " << Utilities::System::get_time() << '\n'
          << "#" << '\n'
          << "# For a description of the UCD format see the AVS Developer's guide."
          << '\n'
          << "#" << '\n';
      }

    // start with ucd data
    out << n_nodes << ' ' << n_cells << ' ' << n_data_sets << ' ' << 0
        << ' ' // no cell data at present
        << 0   // no model data
        << '\n';

    write_nodes(patches, ucd_out);
    out << '\n';

    write_cells(patches, ucd_out);
    out << '\n';

    //---------------------------
    // now write data
    if (n_data_sets != 0)
      {
        out << n_data_sets << "    "; // number of vectors
        for (unsigned int i = 0; i < n_data_sets; ++i)
          out << 1 << ' '; // number of components;
        // only 1 supported presently
        out << '\n';

        for (unsigned int data_set = 0; data_set < n_data_sets; ++data_set)
          out << data_names[data_set]
              << ",dimensionless" // no units supported at present
              << '\n';

        write_data(patches, n_data_sets, true, ucd_out);
      }
    // make sure everything now gets to disk
    out.flush();

    // assert the stream is still ok
    AssertThrow(out.fail() == false, ExcIO());
  }


  template <int dim, int spacedim>
  void
  write_dx(
    const std::vector<Patch<dim, spacedim>> &patches,
    const std::vector<std::string>          &data_names,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>> &,
    const DXFlags &flags,
    std::ostream  &out)
  {
    // Point output is currently not implemented.
    AssertThrow(dim > 0, ExcNotImplemented());

    AssertThrow(out.fail() == false, ExcIO());

#ifndef DEAL_II_WITH_MPI
    // verify that there are indeed patches to be written out. most of the
    // times, people just forget to call build_patches when there are no
    // patches, so a warning is in order. that said, the assertion is disabled
    // if we support MPI since then it can happen that on the coarsest mesh, a
    // processor simply has no cells it actually owns, and in that case it is
    // legit if there are no patches
    Assert(patches.size() > 0, ExcNoPatches());
#else
    if (patches.empty())
      return;
#endif
    // Stream with special features for dx output
    DXStream dx_out(out, flags);

    // Variable counting the offset of binary data.
    unsigned int offset = 0;

    const unsigned int n_data_sets = data_names.size();

    // first count the number of cells and cells for later use
    unsigned int n_nodes;
    unsigned int n_cells;
    std::tie(n_nodes, n_cells) = count_nodes_and_cells(patches);

    // start with vertices order is lexicographical, x varying fastest
    out << "object \"vertices\" class array type float rank 1 shape "
        << spacedim << " items " << n_nodes;

    if (flags.coordinates_binary)
      {
        out << " lsb ieee data 0" << '\n';
        offset += n_nodes * spacedim * sizeof(float);
      }
    else
      {
        out << " data follows" << '\n';
        write_nodes(patches, dx_out);
      }

    //-----------------------------
    // first write the coordinates of all vertices

    //---------------------------------------
    // write cells
    out << "object \"cells\" class array type int rank 1 shape "
        << GeometryInfo<dim>::vertices_per_cell << " items " << n_cells;

    if (flags.int_binary)
      {
        out << " lsb binary data " << offset << '\n';
        offset += n_cells * sizeof(int);
      }
    else
      {
        out << " data follows" << '\n';
        write_cells(patches, dx_out);
        out << '\n';
      }


    out << "attribute \"element type\" string \"";
    if constexpr (dim == 1)
      out << "lines";
    else if constexpr (dim == 2)
      out << "quads";
    else if constexpr (dim == 3)
      out << "cubes";
    out << "\"" << '\n' << "attribute \"ref\" string \"positions\"" << '\n';

    // TODO:[GK] Patches must be of same size!
    //---------------------------
    // write neighbor information
    if (flags.write_neighbors)
      {
        out << "object \"neighbors\" class array type int rank 1 shape "
            << GeometryInfo<dim>::faces_per_cell << " items " << n_cells
            << " data follows";

        for (const auto &patch : patches)
          {
            const unsigned int n               = patch.n_subdivisions;
            const unsigned int n1              = (dim > 0) ? n : 1;
            const unsigned int n2              = (dim > 1) ? n : 1;
            const unsigned int n3              = (dim > 2) ? n : 1;
            const unsigned int x_minus         = (dim > 0) ? 0 : 0;
            const unsigned int x_plus          = (dim > 0) ? 1 : 0;
            const unsigned int y_minus         = (dim > 1) ? 2 : 0;
            const unsigned int y_plus          = (dim > 1) ? 3 : 0;
            const unsigned int z_minus         = (dim > 2) ? 4 : 0;
            const unsigned int z_plus          = (dim > 2) ? 5 : 0;
            unsigned int       cells_per_patch = Utilities::fixed_power<dim>(n);
            unsigned int       dx              = 1;
            unsigned int       dy              = n;
            unsigned int       dz              = n * n;

            const unsigned int patch_start =
              patch.patch_index * cells_per_patch;

            for (unsigned int i3 = 0; i3 < n3; ++i3)
              for (unsigned int i2 = 0; i2 < n2; ++i2)
                for (unsigned int i1 = 0; i1 < n1; ++i1)
                  {
                    const unsigned int nx = i1 * dx;
                    const unsigned int ny = i2 * dy;
                    const unsigned int nz = i3 * dz;

                    // There are no neighbors for dim==0. Note that this case is
                    // caught by the AssertThrow at the beginning of this
                    // function anyway. This condition avoids compiler warnings.
                    if (dim < 1)
                      continue;

                    out << '\n';
                    // Direction -x Last cell in row of other patch
                    if (i1 == 0)
                      {
                        const unsigned int nn = patch.neighbors[x_minus];
                        out << '\t';
                        if (nn != patch.no_neighbor)
                          out
                            << (nn * cells_per_patch + ny + nz + dx * (n - 1));
                        else
                          out << "-1";
                      }
                    else
                      {
                        out << '\t' << patch_start + nx - dx + ny + nz;
                      }
                    // Direction +x First cell in row of other patch
                    if (i1 == n - 1)
                      {
                        const unsigned int nn = patch.neighbors[x_plus];
                        out << '\t';
                        if (nn != patch.no_neighbor)
                          out << (nn * cells_per_patch + ny + nz);
                        else
                          out << "-1";
                      }
                    else
                      {
                        out << '\t' << patch_start + nx + dx + ny + nz;
                      }
                    if (dim < 2)
                      continue;
                    // Direction -y
                    if (i2 == 0)
                      {
                        const unsigned int nn = patch.neighbors[y_minus];
                        out << '\t';
                        if (nn != patch.no_neighbor)
                          out
                            << (nn * cells_per_patch + nx + nz + dy * (n - 1));
                        else
                          out << "-1";
                      }
                    else
                      {
                        out << '\t' << patch_start + nx + ny - dy + nz;
                      }
                    // Direction +y
                    if (i2 == n - 1)
                      {
                        const unsigned int nn = patch.neighbors[y_plus];
                        out << '\t';
                        if (nn != patch.no_neighbor)
                          out << (nn * cells_per_patch + nx + nz);
                        else
                          out << "-1";
                      }
                    else
                      {
                        out << '\t' << patch_start + nx + ny + dy + nz;
                      }
                    if (dim < 3)
                      continue;

                    // Direction -z
                    if (i3 == 0)
                      {
                        const unsigned int nn = patch.neighbors[z_minus];
                        out << '\t';
                        if (nn != patch.no_neighbor)
                          out
                            << (nn * cells_per_patch + nx + ny + dz * (n - 1));
                        else
                          out << "-1";
                      }
                    else
                      {
                        out << '\t' << patch_start + nx + ny + nz - dz;
                      }
                    // Direction +z
                    if (i3 == n - 1)
                      {
                        const unsigned int nn = patch.neighbors[z_plus];
                        out << '\t';
                        if (nn != patch.no_neighbor)
                          out << (nn * cells_per_patch + nx + ny);
                        else
                          out << "-1";
                      }
                    else
                      {
                        out << '\t' << patch_start + nx + ny + nz + dz;
                      }
                  }
            out << '\n';
          }
      }
    //---------------------------
    // now write data
    if (n_data_sets != 0)
      {
        out << "object \"data\" class array type float rank 1 shape "
            << n_data_sets << " items " << n_nodes;

        if (flags.data_binary)
          {
            out << " lsb ieee data " << offset << '\n';
            offset += n_data_sets * n_nodes *
                      ((flags.data_double) ? sizeof(double) : sizeof(float));
          }
        else
          {
            out << " data follows" << '\n';
            write_data(patches, n_data_sets, flags.data_double, dx_out);
          }

        // loop over all patches
        out << "attribute \"dep\" string \"positions\"" << '\n';
      }
    else
      {
        out << "object \"data\" class constantarray type float rank 0 items "
            << n_nodes << " data follows" << '\n'
            << '0' << '\n';
      }

    // no model data

    out << "object \"deal data\" class field" << '\n'
        << "component \"positions\" value \"vertices\"" << '\n'
        << "component \"connections\" value \"cells\"" << '\n'
        << "component \"data\" value \"data\"" << '\n';

    if (flags.write_neighbors)
      out << "component \"neighbors\" value \"neighbors\"" << '\n';

    {
      out << "attribute \"created\" string \"" << Utilities::System::get_date()
          << ' ' << Utilities::System::get_time() << '"' << '\n';
    }

    out << "end" << '\n';
    // Write all binary data now
    if (flags.coordinates_binary)
      write_nodes(patches, dx_out);
    if (flags.int_binary)
      write_cells(patches, dx_out);
    if (flags.data_binary)
      write_data(patches, n_data_sets, flags.data_double, dx_out);

    // make sure everything now gets to disk
    out.flush();

    // assert the stream is still ok
    AssertThrow(out.fail() == false, ExcIO());
  }



  template <int dim, int spacedim>
  void
  write_gnuplot(
    const std::vector<Patch<dim, spacedim>> &patches,
    const std::vector<std::string>          &data_names,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>> &,
    const GnuplotFlags &flags,
    std::ostream       &out)
  {
    AssertThrow(out.fail() == false, ExcIO());

#ifndef DEAL_II_WITH_MPI
    // verify that there are indeed patches to be written out. most
    // of the times, people just forget to call build_patches when there
    // are no patches, so a warning is in order. that said, the
    // assertion is disabled if we support MPI since then it can
    // happen that on the coarsest mesh, a processor simply has no
    // cells it actually owns, and in that case it is legit if there
    // are no patches
    Assert(patches.size() > 0, ExcNoPatches());
#else
    if (patches.empty())
      return;
#endif

    const unsigned int n_data_sets = data_names.size();

    // write preamble
    {
      out << "# This file was generated by the deal.II library." << '\n'
          << "# Date =  " << Utilities::System::get_date() << '\n'
          << "# Time =  " << Utilities::System::get_time() << '\n'
          << "#" << '\n'
          << "# For a description of the GNUPLOT format see the GNUPLOT manual."
          << '\n'
          << "#" << '\n'
          << "# ";

      AssertThrow(spacedim <= flags.space_dimension_labels.size(),
                  GnuplotFlags::ExcNotEnoughSpaceDimensionLabels());
      for (unsigned int spacedim_n = 0; spacedim_n < spacedim; ++spacedim_n)
        {
          out << '<' << flags.space_dimension_labels.at(spacedim_n) << "> ";
        }

      for (const auto &data_name : data_names)
        out << '<' << data_name << "> ";
      out << '\n';
    }


    // loop over all patches
    for (const auto &patch : patches)
      {
        const unsigned int n_subdivisions         = patch.n_subdivisions;
        const unsigned int n_points_per_direction = n_subdivisions + 1;

        Assert((patch.data.n_rows() == n_data_sets &&
                !patch.points_are_available) ||
                 (patch.data.n_rows() == n_data_sets + spacedim &&
                  patch.points_are_available),
               ExcDimensionMismatch(patch.points_are_available ?
                                      (n_data_sets + spacedim) :
                                      n_data_sets,
                                    patch.data.n_rows()));

        auto output_point_data =
          [&out, &patch, n_data_sets](const unsigned int point_index) mutable {
            for (unsigned int data_set = 0; data_set < n_data_sets; ++data_set)
              out << patch.data(data_set, point_index) << ' ';
          };

        switch (dim)
          {
            case 0:
              {
                Assert(patch.reference_cell == ReferenceCells::Vertex,
                       ExcInternalError());
                Assert(patch.data.n_cols() == 1,
                       ExcInvalidDatasetSize(patch.data.n_cols(),
                                             n_subdivisions + 1));


                // compute coordinates for this patch point
                out << get_equispaced_location(patch, {}, n_subdivisions)
                    << ' ';
                output_point_data(0);
                out << '\n';
                out << '\n';
                break;
              }

            case 1:
              {
                Assert(patch.reference_cell == ReferenceCells::Line,
                       ExcInternalError());
                Assert(patch.data.n_cols() ==
                         Utilities::fixed_power<dim>(n_points_per_direction),
                       ExcInvalidDatasetSize(patch.data.n_cols(),
                                             n_subdivisions + 1));

                for (unsigned int i1 = 0; i1 < n_points_per_direction; ++i1)
                  {
                    // compute coordinates for this patch point
                    out << get_equispaced_location(patch, {i1}, n_subdivisions)
                        << ' ';

                    output_point_data(i1);
                    out << '\n';
                  }
                // end of patch
                out << '\n';
                out << '\n';
                break;
              }

            case 2:
              {
                if (patch.reference_cell == ReferenceCells::Quadrilateral)
                  {
                    Assert(patch.data.n_cols() == Utilities::fixed_power<dim>(
                                                    n_points_per_direction),
                           ExcInvalidDatasetSize(patch.data.n_cols(),
                                                 n_subdivisions + 1));

                    for (unsigned int i2 = 0; i2 < n_points_per_direction; ++i2)
                      {
                        for (unsigned int i1 = 0; i1 < n_points_per_direction;
                             ++i1)
                          {
                            // compute coordinates for this patch point
                            out << get_equispaced_location(patch,
                                                           {i1, i2},
                                                           n_subdivisions)
                                << ' ';

                            output_point_data(i1 + i2 * n_points_per_direction);
                            out << '\n';
                          }
                        // end of row in patch
                        out << '\n';
                      }
                  }
                else if (patch.reference_cell == ReferenceCells::Triangle)
                  {
                    Assert(n_subdivisions == 1, ExcNotImplemented());

                    Assert(patch.data.n_cols() == 3, ExcInternalError());

                    // Gnuplot can only plot surfaces if each facet of the
                    // surface is a bilinear patch, or a subdivided bilinear
                    // patch with equally many points along each row of the
                    // subdivision. This is what the code above for
                    // quadrilaterals does. We emulate this by repeating the
                    // third point of a triangle twice so that there are two
                    // points for that row as well -- i.e., we write a 2x2
                    // bilinear patch where two of the points are collapsed onto
                    // one vertex.
                    //
                    // This also matches the example here:
                    // https://stackoverflow.com/questions/42784369/drawing-triangular-mesh-using-gnuplot
                    out << get_node_location(patch, 0) << ' ';
                    output_point_data(0);
                    out << '\n';

                    out << get_node_location(patch, 1) << ' ';
                    output_point_data(1);
                    out << '\n';
                    out << '\n'; // end of one row of points

                    out << get_node_location(patch, 2) << ' ';
                    output_point_data(2);
                    out << '\n';

                    out << get_node_location(patch, 2) << ' ';
                    output_point_data(2);
                    out << '\n';
                    out << '\n'; // end of the second row of points
                    out << '\n'; // end of the entire patch
                  }
                else
                  // There aren't any other reference cells in 2d than the
                  // quadrilateral and the triangle. So whatever we got here
                  // can't be any good
                  DEAL_II_ASSERT_UNREACHABLE();
                // end of patch
                out << '\n';

                break;
              }

            case 3:
              {
                if (patch.reference_cell == ReferenceCells::Hexahedron)
                  {
                    Assert(patch.data.n_cols() == Utilities::fixed_power<dim>(
                                                    n_points_per_direction),
                           ExcInvalidDatasetSize(patch.data.n_cols(),
                                                 n_subdivisions + 1));

                    // for all grid points: draw lines into all positive
                    // coordinate directions if there is another grid point
                    // there
                    for (unsigned int i3 = 0; i3 < n_points_per_direction; ++i3)
                      for (unsigned int i2 = 0; i2 < n_points_per_direction;
                           ++i2)
                        for (unsigned int i1 = 0; i1 < n_points_per_direction;
                             ++i1)
                          {
                            // compute coordinates for this patch point
                            const Point<spacedim> this_point =
                              get_equispaced_location(patch,
                                                      {i1, i2, i3},
                                                      n_subdivisions);
                            // line into positive x-direction if possible
                            if (i1 < n_subdivisions)
                              {
                                // write point here and its data
                                out << this_point << ' ';
                                output_point_data(i1 +
                                                  i2 * n_points_per_direction +
                                                  i3 * n_points_per_direction *
                                                    n_points_per_direction);
                                out << '\n';

                                // write point there and its data
                                out << get_equispaced_location(patch,
                                                               {i1 + 1, i2, i3},
                                                               n_subdivisions)
                                    << ' ';

                                output_point_data((i1 + 1) +
                                                  i2 * n_points_per_direction +
                                                  i3 * n_points_per_direction *
                                                    n_points_per_direction);
                                out << '\n';

                                // end of line
                                out << '\n' << '\n';
                              }

                            // line into positive y-direction if possible
                            if (i2 < n_subdivisions)
                              {
                                // write point here and its data
                                out << this_point << ' ';
                                output_point_data(i1 +
                                                  i2 * n_points_per_direction +
                                                  i3 * n_points_per_direction *
                                                    n_points_per_direction);
                                out << '\n';

                                // write point there and its data
                                out << get_equispaced_location(patch,
                                                               {i1, i2 + 1, i3},
                                                               n_subdivisions)
                                    << ' ';

                                output_point_data(
                                  i1 + (i2 + 1) * n_points_per_direction +
                                  i3 * n_points_per_direction *
                                    n_points_per_direction);
                                out << '\n';

                                // end of line
                                out << '\n' << '\n';
                              }

                            // line into positive z-direction if possible
                            if (i3 < n_subdivisions)
                              {
                                // write point here and its data
                                out << this_point << ' ';
                                output_point_data(i1 +
                                                  i2 * n_points_per_direction +
                                                  i3 * n_points_per_direction *
                                                    n_points_per_direction);
                                out << '\n';

                                // write point there and its data
                                out << get_equispaced_location(patch,
                                                               {i1, i2, i3 + 1},
                                                               n_subdivisions)
                                    << ' ';

                                output_point_data(
                                  i1 + i2 * n_points_per_direction +
                                  (i3 + 1) * n_points_per_direction *
                                    n_points_per_direction);
                                out << '\n';
                                // end of line
                                out << '\n' << '\n';
                              }
                          }
                  }
                else if (patch.reference_cell == ReferenceCells::Tetrahedron)
                  {
                    Assert(n_subdivisions == 1, ExcNotImplemented());

                    // Draw the tetrahedron as a collection of two lines.
                    for (const unsigned int v : {0, 1, 2, 0, 3, 2})
                      {
                        out << get_node_location(patch, v) << ' ';
                        output_point_data(v);
                        out << '\n';
                      }
                    out << '\n'; // end of first line

                    for (const unsigned int v : {3, 1})
                      {
                        out << get_node_location(patch, v) << ' ';
                        output_point_data(v);
                        out << '\n';
                      }
                    out << '\n'; // end of second line
                  }
                else if (patch.reference_cell == ReferenceCells::Pyramid)
                  {
                    Assert(n_subdivisions == 1, ExcNotImplemented());

                    // Draw the pyramid as a collection of two lines.
                    for (const unsigned int v : {0, 1, 3, 2, 0, 4, 1})
                      {
                        out << get_node_location(patch, v) << ' ';
                        output_point_data(v);
                        out << '\n';
                      }
                    out << '\n'; // end of first line

                    for (const unsigned int v : {2, 4, 3})
                      {
                        out << get_node_location(patch, v) << ' ';
                        output_point_data(v);
                        out << '\n';
                      }
                    out << '\n'; // end of second line
                  }
                else if (patch.reference_cell == ReferenceCells::Wedge)
                  {
                    Assert(n_subdivisions == 1, ExcNotImplemented());

                    // Draw the wedge as a collection of three
                    // lines. The first one wraps around the base,
                    // goes up to the top, and wraps around that. The
                    // second and third are just individual lines
                    // going from base to top.
                    for (const unsigned int v : {0, 1, 2, 0, 3, 4, 5, 3})
                      {
                        out << get_node_location(patch, v) << ' ';
                        output_point_data(v);
                        out << '\n';
                      }
                    out << '\n'; // end of first line

                    for (const unsigned int v : {1, 4})
                      {
                        out << get_node_location(patch, v) << ' ';
                        output_point_data(v);
                        out << '\n';
                      }
                    out << '\n'; // end of second line

                    for (const unsigned int v : {2, 5})
                      {
                        out << get_node_location(patch, v) << ' ';
                        output_point_data(v);
                        out << '\n';
                      }
                    out << '\n'; // end of second line
                  }
                else
                  // No other reference cells are currently implemented
                  DEAL_II_NOT_IMPLEMENTED();

                break;
              }

            default:
              DEAL_II_NOT_IMPLEMENTED();
          }
      }
    // make sure everything now gets to disk
    out.flush();

    AssertThrow(out.fail() == false, ExcIO());
  }


  namespace
  {
    template <int dim, int spacedim>
    void
    do_write_povray(const std::vector<Patch<dim, spacedim>> &,
                    const std::vector<std::string> &,
                    const PovrayFlags &,
                    std::ostream &)
    {
      Assert(false,
             ExcMessage("Writing files in POVRAY format is only supported "
                        "for two-dimensional meshes."));
    }



    void
    do_write_povray(const std::vector<Patch<2, 2>> &patches,
                    const std::vector<std::string> &data_names,
                    const PovrayFlags              &flags,
                    std::ostream                   &out)
    {
      AssertThrow(out.fail() == false, ExcIO());

#ifndef DEAL_II_WITH_MPI
      // verify that there are indeed patches to be written out. most
      // of the times, people just forget to call build_patches when there
      // are no patches, so a warning is in order. that said, the
      // assertion is disabled if we support MPI since then it can
      // happen that on the coarsest mesh, a processor simply has no cells it
      // actually owns, and in that case it is legit if there are no patches
      Assert(patches.size() > 0, ExcNoPatches());
#else
      if (patches.empty())
        return;
#endif
      constexpr int dim = 2;
      (void)dim;
      constexpr int spacedim = 2;

      const unsigned int n_data_sets = data_names.size();
      (void)n_data_sets;

      // write preamble
      {
        out
          << "/* This file was generated by the deal.II library." << '\n'
          << "   Date =  " << Utilities::System::get_date() << '\n'
          << "   Time =  " << Utilities::System::get_time() << '\n'
          << '\n'
          << "   For a description of the POVRAY format see the POVRAY manual."
          << '\n'
          << "*/ " << '\n';

        // include files
        out << "#include \"colors.inc\" " << '\n'
            << "#include \"textures.inc\" " << '\n';


        // use external include file for textures, camera and light
        if (flags.external_data)
          out << "#include \"data.inc\" " << '\n';
        else // all definitions in data file
          {
            // camera
            out << '\n'
                << '\n'
                << "camera {" << '\n'
                << "  location <1,4,-7>" << '\n'
                << "  look_at <0,0,0>" << '\n'
                << "  angle 30" << '\n'
                << "}" << '\n';

            // light
            out << '\n'
                << "light_source {" << '\n'
                << "  <1,4,-7>" << '\n'
                << "  color Grey" << '\n'
                << "}" << '\n';
            out << '\n'
                << "light_source {" << '\n'
                << "  <0,20,0>" << '\n'
                << "  color White" << '\n'
                << "}" << '\n';
          }
      }

      // max. and min. height of solution
      Assert(patches.size() > 0, ExcNoPatches());
      double hmin = patches[0].data(0, 0);
      double hmax = patches[0].data(0, 0);

      for (const auto &patch : patches)
        {
          const unsigned int n_subdivisions = patch.n_subdivisions;

          Assert((patch.data.n_rows() == n_data_sets &&
                  !patch.points_are_available) ||
                   (patch.data.n_rows() == n_data_sets + spacedim &&
                    patch.points_are_available),
                 ExcDimensionMismatch(patch.points_are_available ?
                                        (n_data_sets + spacedim) :
                                        n_data_sets,
                                      patch.data.n_rows()));
          Assert(patch.data.n_cols() ==
                   Utilities::fixed_power<dim>(n_subdivisions + 1),
                 ExcInvalidDatasetSize(patch.data.n_cols(),
                                       n_subdivisions + 1));

          for (unsigned int i = 0; i < n_subdivisions + 1; ++i)
            for (unsigned int j = 0; j < n_subdivisions + 1; ++j)
              {
                const int dl = i * (n_subdivisions + 1) + j;
                if (patch.data(0, dl) < hmin)
                  hmin = patch.data(0, dl);
                if (patch.data(0, dl) > hmax)
                  hmax = patch.data(0, dl);
              }
        }

      out << "#declare HMIN=" << hmin << ";" << '\n'
          << "#declare HMAX=" << hmax << ";" << '\n'
          << '\n';

      if (!flags.external_data)
        {
          // texture with scaled niveau lines 10 lines in the surface
          out << "#declare Tex=texture{" << '\n'
              << "  pigment {" << '\n'
              << "    gradient y" << '\n'
              << "    scale y*(HMAX-HMIN)*" << 0.1 << '\n'
              << "    color_map {" << '\n'
              << "      [0.00 color Light_Purple] " << '\n'
              << "      [0.95 color Light_Purple] " << '\n'
              << "      [1.00 color White]    " << '\n'
              << "} } }" << '\n'
              << '\n';
        }

      if (!flags.bicubic_patch)
        {
          // start of mesh header
          out << '\n' << "mesh {" << '\n';
        }

      // loop over all patches
      for (const auto &patch : patches)
        {
          const unsigned int n_subdivisions = patch.n_subdivisions;
          const unsigned int n              = n_subdivisions + 1;
          const unsigned int d1             = 1;
          const unsigned int d2             = n;

          Assert((patch.data.n_rows() == n_data_sets &&
                  !patch.points_are_available) ||
                   (patch.data.n_rows() == n_data_sets + spacedim &&
                    patch.points_are_available),
                 ExcDimensionMismatch(patch.points_are_available ?
                                        (n_data_sets + spacedim) :
                                        n_data_sets,
                                      patch.data.n_rows()));
          Assert(patch.data.n_cols() == Utilities::fixed_power<dim>(n),
                 ExcInvalidDatasetSize(patch.data.n_cols(),
                                       n_subdivisions + 1));


          std::vector<Point<spacedim>> ver(n * n);

          for (unsigned int i2 = 0; i2 < n; ++i2)
            for (unsigned int i1 = 0; i1 < n; ++i1)
              {
                // compute coordinates for this patch point, storing in ver
                ver[i1 * d1 + i2 * d2] =
                  get_equispaced_location(patch, {i1, i2}, n_subdivisions);
              }


          if (!flags.bicubic_patch)
            {
              // approximate normal vectors in patch
              std::vector<Point<3>> nrml;
              // only if smooth triangles are used
              if (flags.smooth)
                {
                  nrml.resize(n * n);
                  // These are difference quotients of the surface
                  // mapping. We take them symmetric inside the
                  // patch and one-sided at the edges
                  Point<3> h1, h2;
                  // Now compute normals in every point
                  for (unsigned int i = 0; i < n; ++i)
                    for (unsigned int j = 0; j < n; ++j)
                      {
                        const unsigned int il = (i == 0) ? i : (i - 1);
                        const unsigned int ir =
                          (i == n_subdivisions) ? i : (i + 1);
                        const unsigned int jl = (j == 0) ? j : (j - 1);
                        const unsigned int jr =
                          (j == n_subdivisions) ? j : (j + 1);

                        h1[0] =
                          ver[ir * d1 + j * d2][0] - ver[il * d1 + j * d2][0];
                        h1[1] = patch.data(0, ir * d1 + j * d2) -
                                patch.data(0, il * d1 + j * d2);
                        h1[2] =
                          ver[ir * d1 + j * d2][1] - ver[il * d1 + j * d2][1];

                        h2[0] =
                          ver[i * d1 + jr * d2][0] - ver[i * d1 + jl * d2][0];
                        h2[1] = patch.data(0, i * d1 + jr * d2) -
                                patch.data(0, i * d1 + jl * d2);
                        h2[2] =
                          ver[i * d1 + jr * d2][1] - ver[i * d1 + jl * d2][1];

                        nrml[i * d1 + j * d2][0] =
                          h1[1] * h2[2] - h1[2] * h2[1];
                        nrml[i * d1 + j * d2][1] =
                          h1[2] * h2[0] - h1[0] * h2[2];
                        nrml[i * d1 + j * d2][2] =
                          h1[0] * h2[1] - h1[1] * h2[0];

                        // normalize Vector
                        double norm = std::hypot(nrml[i * d1 + j * d2][0],
                                                 nrml[i * d1 + j * d2][1],
                                                 nrml[i * d1 + j * d2][2]);

                        if (nrml[i * d1 + j * d2][1] < 0)
                          norm *= -1.;

                        for (unsigned int k = 0; k < 3; ++k)
                          nrml[i * d1 + j * d2][k] /= norm;
                      }
                }

              // setting up triangles
              for (unsigned int i = 0; i < n_subdivisions; ++i)
                for (unsigned int j = 0; j < n_subdivisions; ++j)
                  {
                    // down/left vertex of triangle
                    const int dl = i * d1 + j * d2;
                    if (flags.smooth)
                      {
                        // writing smooth_triangles

                        // down/right triangle
                        out << "smooth_triangle {" << '\n'
                            << "\t<" << ver[dl][0] << "," << patch.data(0, dl)
                            << "," << ver[dl][1] << ">, <" << nrml[dl][0]
                            << ", " << nrml[dl][1] << ", " << nrml[dl][2]
                            << ">," << '\n';
                        out << " \t<" << ver[dl + d1][0] << ","
                            << patch.data(0, dl + d1) << "," << ver[dl + d1][1]
                            << ">, <" << nrml[dl + d1][0] << ", "
                            << nrml[dl + d1][1] << ", " << nrml[dl + d1][2]
                            << ">," << '\n';
                        out << "\t<" << ver[dl + d1 + d2][0] << ","
                            << patch.data(0, dl + d1 + d2) << ","
                            << ver[dl + d1 + d2][1] << ">, <"
                            << nrml[dl + d1 + d2][0] << ", "
                            << nrml[dl + d1 + d2][1] << ", "
                            << nrml[dl + d1 + d2][2] << ">}" << '\n';

                        // upper/left triangle
                        out << "smooth_triangle {" << '\n'
                            << "\t<" << ver[dl][0] << "," << patch.data(0, dl)
                            << "," << ver[dl][1] << ">, <" << nrml[dl][0]
                            << ", " << nrml[dl][1] << ", " << nrml[dl][2]
                            << ">," << '\n';
                        out << "\t<" << ver[dl + d1 + d2][0] << ","
                            << patch.data(0, dl + d1 + d2) << ","
                            << ver[dl + d1 + d2][1] << ">, <"
                            << nrml[dl + d1 + d2][0] << ", "
                            << nrml[dl + d1 + d2][1] << ", "
                            << nrml[dl + d1 + d2][2] << ">," << '\n';
                        out << "\t<" << ver[dl + d2][0] << ","
                            << patch.data(0, dl + d2) << "," << ver[dl + d2][1]
                            << ">, <" << nrml[dl + d2][0] << ", "
                            << nrml[dl + d2][1] << ", " << nrml[dl + d2][2]
                            << ">}" << '\n';
                      }
                    else
                      {
                        // writing standard triangles down/right triangle
                        out << "triangle {" << '\n'
                            << "\t<" << ver[dl][0] << "," << patch.data(0, dl)
                            << "," << ver[dl][1] << ">," << '\n';
                        out << "\t<" << ver[dl + d1][0] << ","
                            << patch.data(0, dl + d1) << "," << ver[dl + d1][1]
                            << ">," << '\n';
                        out << "\t<" << ver[dl + d1 + d2][0] << ","
                            << patch.data(0, dl + d1 + d2) << ","
                            << ver[dl + d1 + d2][1] << ">}" << '\n';

                        // upper/left triangle
                        out << "triangle {" << '\n'
                            << "\t<" << ver[dl][0] << "," << patch.data(0, dl)
                            << "," << ver[dl][1] << ">," << '\n';
                        out << "\t<" << ver[dl + d1 + d2][0] << ","
                            << patch.data(0, dl + d1 + d2) << ","
                            << ver[dl + d1 + d2][1] << ">," << '\n';
                        out << "\t<" << ver[dl + d2][0] << ","
                            << patch.data(0, dl + d2) << "," << ver[dl + d2][1]
                            << ">}" << '\n';
                      }
                  }
            }
          else
            {
              // writing bicubic_patch
              Assert(n_subdivisions == 3,
                     ExcDimensionMismatch(n_subdivisions, 3));
              out << '\n'
                  << "bicubic_patch {" << '\n'
                  << "  type 0" << '\n'
                  << "  flatness 0" << '\n'
                  << "  u_steps 0" << '\n'
                  << "  v_steps 0" << '\n';
              for (int i = 0; i < 16; ++i)
                {
                  out << "\t<" << ver[i][0] << "," << patch.data(0, i) << ","
                      << ver[i][1] << ">";
                  if (i != 15)
                    out << ",";
                  out << '\n';
                }
              out << "  texture {Tex}" << '\n' << "}" << '\n';
            }
        }

      if (!flags.bicubic_patch)
        {
          // the end of the mesh
          out << "  texture {Tex}" << '\n' << "}" << '\n' << '\n';
        }

      // make sure everything now gets to disk
      out.flush();

      AssertThrow(out.fail() == false, ExcIO());
    }
  } // namespace



  template <int dim, int spacedim>
  void
  write_povray(
    const std::vector<Patch<dim, spacedim>> &patches,
    const std::vector<std::string>          &data_names,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>> &,
    const PovrayFlags &flags,
    std::ostream      &out)
  {
    do_write_povray(patches, data_names, flags, out);
  }



  template <int dim, int spacedim>
  void
  write_eps(
    const std::vector<Patch<dim, spacedim>> & /*patches*/,
    const std::vector<std::string> & /*data_names*/,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>> &,
    const EpsFlags & /*flags*/,
    std::ostream & /*out*/)
  {
    // not implemented, see the documentation of the function
    AssertThrow(dim == 2, ExcNotImplemented());
  }


  template <int spacedim>
  void
  write_eps(
    const std::vector<Patch<2, spacedim>> &patches,
    const std::vector<std::string> & /*data_names*/,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>> &,
    const EpsFlags &flags,
    std::ostream   &out)
  {
    AssertThrow(out.fail() == false, ExcIO());

#ifndef DEAL_II_WITH_MPI
    // verify that there are indeed patches to be written out. most of the
    // times, people just forget to call build_patches when there are no
    // patches, so a warning is in order. that said, the assertion is disabled
    // if we support MPI since then it can happen that on the coarsest mesh, a
    // processor simply has no cells it actually owns, and in that case it is
    // legit if there are no patches
    Assert(patches.size() > 0, ExcNoPatches());
#else
    if (patches.empty())
      return;
#endif

    // set up an array of cells to be written later. this array holds the cells
    // of all the patches as projected to the plane perpendicular to the line of
    // sight.
    //
    // note that they are kept sorted by the set, where we chose the value of
    // the center point of the cell along the line of sight as value for sorting
    std::multiset<EpsCell2d> cells;

    // two variables in which we will store the minimum and maximum values of
    // the field to be used for colorization
    float min_color_value = std::numeric_limits<float>::max();
    float max_color_value = std::numeric_limits<float>::min();

    // Array for z-coordinates of points. The elevation determined by a function
    // if spacedim=2 or the z-coordinate of the grid point if spacedim=3
    double heights[4] = {0, 0, 0, 0};

    // compute the cells for output and enter them into the set above note that
    // since dim==2, we have exactly four vertices per patch and per cell
    for (const auto &patch : patches)
      {
        const unsigned int n_subdivisions = patch.n_subdivisions;
        const unsigned int n              = n_subdivisions + 1;
        const unsigned int d1             = 1;
        const unsigned int d2             = n;

        for (unsigned int i2 = 0; i2 < n_subdivisions; ++i2)
          for (unsigned int i1 = 0; i1 < n_subdivisions; ++i1)
            {
              Point<spacedim> points[4];
              points[0] =
                get_equispaced_location(patch, {i1, i2}, n_subdivisions);
              points[1] =
                get_equispaced_location(patch, {i1 + 1, i2}, n_subdivisions);
              points[2] =
                get_equispaced_location(patch, {i1, i2 + 1}, n_subdivisions);
              points[3] = get_equispaced_location(patch,
                                                  {i1 + 1, i2 + 1},
                                                  n_subdivisions);

              switch (spacedim)
                {
                  case 2:
                    Assert((flags.height_vector < patch.data.n_rows()) ||
                             patch.data.n_rows() == 0,
                           ExcIndexRange(flags.height_vector,
                                         0,
                                         patch.data.n_rows()));
                    heights[0] =
                      patch.data.n_rows() != 0 ?
                        patch.data(flags.height_vector, i1 * d1 + i2 * d2) *
                          flags.z_scaling :
                        0;
                    heights[1] = patch.data.n_rows() != 0 ?
                                   patch.data(flags.height_vector,
                                              (i1 + 1) * d1 + i2 * d2) *
                                     flags.z_scaling :
                                   0;
                    heights[2] = patch.data.n_rows() != 0 ?
                                   patch.data(flags.height_vector,
                                              i1 * d1 + (i2 + 1) * d2) *
                                     flags.z_scaling :
                                   0;
                    heights[3] = patch.data.n_rows() != 0 ?
                                   patch.data(flags.height_vector,
                                              (i1 + 1) * d1 + (i2 + 1) * d2) *
                                     flags.z_scaling :
                                   0;

                    break;
                  case 3:
                    // Copy z-coordinates into the height vector
                    for (unsigned int i = 0; i < 4; ++i)
                      heights[i] = points[i][2];
                    break;
                  default:
                    DEAL_II_NOT_IMPLEMENTED();
                }


              // now compute the projection of the bilinear cell given by the
              // four vertices and their heights and write them to a proper cell
              // object. note that we only need the first two components of the
              // projected position for output, but we need the value along the
              // line of sight for sorting the cells for back-to- front-output
              //
              // this computation was first written by Stefan Nauber. please
              // no-one ask me why it works that way (or may be not), especially
              // not about the angles and the sign of the height field, I don't
              // know it.
              EpsCell2d    eps_cell;
              const double pi = numbers::PI;
              const double cx =
                             -std::cos(pi - flags.azimut_angle * 2 * pi / 360.),
                           cz = -std::cos(flags.turn_angle * 2 * pi / 360.),
                           sx =
                             std::sin(pi - flags.azimut_angle * 2 * pi / 360.),
                           sz = std::sin(flags.turn_angle * 2 * pi / 360.);
              for (unsigned int vertex = 0; vertex < 4; ++vertex)
                {
                  const double x = points[vertex][0], y = points[vertex][1],
                               z = -heights[vertex];

                  eps_cell.vertices[vertex][0] = -cz * x + sz * y;
                  eps_cell.vertices[vertex][1] =
                    -cx * sz * x - cx * cz * y - sx * z;

                  //      ( 1 0    0 )
                  // D1 = ( 0 cx -sx )
                  //      ( 0 sx  cx )

                  //      ( cy 0 sy )
                  // Dy = (  0 1  0 )
                  //      (-sy 0 cy )

                  //      ( cz -sz 0 )
                  // Dz = ( sz  cz 0 )
                  //      (  0   0 1 )

                  //       ( cz -sz 0 )( 1 0    0 )(x)   (
                  //       cz*x-sz*(cx*y-sx*z)+0*(sx*y+cx*z) )
                  // Dxz = ( sz  cz 0 )( 0 cx -sx )(y) = (
                  // sz*x+cz*(cx*y-sx*z)+0*(sx*y+cx*z) )
                  //       (  0   0 1 )( 0 sx  cx )(z)   (  0*x+
                  //       *(cx*y-sx*z)+1*(sx*y+cx*z) )
                }

              // compute coordinates of center of cell
              const Point<spacedim> center_point =
                (points[0] + points[1] + points[2] + points[3]) / 4;
              const double center_height =
                -(heights[0] + heights[1] + heights[2] + heights[3]) / 4;

              // compute the depth into the picture
              eps_cell.depth = -sx * sz * center_point[0] -
                               sx * cz * center_point[1] + cx * center_height;

              if (flags.draw_cells && flags.shade_cells)
                {
                  Assert((flags.color_vector < patch.data.n_rows()) ||
                           patch.data.n_rows() == 0,
                         ExcIndexRange(flags.color_vector,
                                       0,
                                       patch.data.n_rows()));
                  const double color_values[4] = {
                    patch.data.n_rows() != 0 ?
                      patch.data(flags.color_vector, i1 * d1 + i2 * d2) :
                      1,

                    patch.data.n_rows() != 0 ?
                      patch.data(flags.color_vector, (i1 + 1) * d1 + i2 * d2) :
                      1,

                    patch.data.n_rows() != 0 ?
                      patch.data(flags.color_vector, i1 * d1 + (i2 + 1) * d2) :
                      1,

                    patch.data.n_rows() != 0 ?
                      patch.data(flags.color_vector,
                                 (i1 + 1) * d1 + (i2 + 1) * d2) :
                      1};

                  // set color value to average of the value at the vertices
                  eps_cell.color_value = (color_values[0] + color_values[1] +
                                          color_values[3] + color_values[2]) /
                                         4;

                  // update bounds of color field
                  min_color_value =
                    std::min(min_color_value, eps_cell.color_value);
                  max_color_value =
                    std::max(max_color_value, eps_cell.color_value);
                }

              // finally add this cell
              cells.insert(eps_cell);
            }
      }

    // find out minimum and maximum x and y coordinates to compute offsets and
    // scaling factors
    double x_min = cells.begin()->vertices[0][0];
    double x_max = x_min;
    double y_min = cells.begin()->vertices[0][1];
    double y_max = y_min;

    for (const auto &cell : cells)
      for (const auto &vertex : cell.vertices)
        {
          x_min = std::min(x_min, vertex[0]);
          x_max = std::max(x_max, vertex[0]);
          y_min = std::min(y_min, vertex[1]);
          y_max = std::max(y_max, vertex[1]);
        }

    // scale in x-direction such that in the output 0 <= x <= 300. don't scale
    // in y-direction to preserve the shape of the triangulation
    const double scale =
      (flags.size /
       (flags.size_type == EpsFlags::width ? x_max - x_min : y_min - y_max));

    const Point<2> offset(x_min, y_min);


    // now write preamble
    {
      out << "%!PS-Adobe-2.0 EPSF-1.2" << '\n'
          << "%%Title: deal.II Output" << '\n'
          << "%%Creator: the deal.II library" << '\n'
          << "%%Creation Date: " << Utilities::System::get_date() << " - "
          << Utilities::System::get_time() << '\n'
          << "%%BoundingBox: "
          // lower left corner
          << "0 0 "
          // upper right corner
          << static_cast<unsigned int>((x_max - x_min) * scale + 0.5) << ' '
          << static_cast<unsigned int>((y_max - y_min) * scale + 0.5) << '\n';

      // define some abbreviations to keep the output small:
      // m=move turtle to
      // l=define a line
      // s=set rgb color
      // sg=set gray value
      // lx=close the line and plot the line
      // lf=close the line and fill the interior
      out << "/m {moveto} bind def" << '\n'
          << "/l {lineto} bind def" << '\n'
          << "/s {setrgbcolor} bind def" << '\n'
          << "/sg {setgray} bind def" << '\n'
          << "/lx {lineto closepath stroke} bind def" << '\n'
          << "/lf {lineto closepath fill} bind def" << '\n';

      out << "%%EndProlog" << '\n' << '\n';
      // set fine lines
      out << flags.line_width << " setlinewidth" << '\n';
    }

    // check if min and max values for the color are actually different. If
    // that is not the case (such things happen, for example, in the very first
    // time step of a time dependent problem, if the initial values are zero),
    // all values are equal, and then we can draw everything in an arbitrary
    // color. Thus, change one of the two values arbitrarily
    if (max_color_value == min_color_value)
      max_color_value = min_color_value + 1;

    // now we've got all the information we need. write the cells. note: due to
    // the ordering, we traverse the list of cells back-to-front
    for (const auto &cell : cells)
      {
        if (flags.draw_cells)
          {
            if (flags.shade_cells)
              {
                const EpsFlags::RgbValues rgb_values =
                  (*flags.color_function)(cell.color_value,
                                          min_color_value,
                                          max_color_value);

                // write out color
                if (rgb_values.is_grey())
                  out << rgb_values.red << " sg ";
                else
                  out << rgb_values.red << ' ' << rgb_values.green << ' '
                      << rgb_values.blue << " s ";
              }
            else
              out << "1 sg ";

            out << (cell.vertices[0] - offset) * scale << " m "
                << (cell.vertices[1] - offset) * scale << " l "
                << (cell.vertices[3] - offset) * scale << " l "
                << (cell.vertices[2] - offset) * scale << " lf" << '\n';
          }

        if (flags.draw_mesh)
          out << "0 sg " // draw lines in black
              << (cell.vertices[0] - offset) * scale << " m "
              << (cell.vertices[1] - offset) * scale << " l "
              << (cell.vertices[3] - offset) * scale << " l "
              << (cell.vertices[2] - offset) * scale << " lx" << '\n';
      }
    out << "showpage" << '\n';

    out.flush();

    AssertThrow(out.fail() == false, ExcIO());
  }



  template <int dim, int spacedim>
  void
  write_gmv(
    const std::vector<Patch<dim, spacedim>> &patches,
    const std::vector<std::string>          &data_names,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>> &,
    const GmvFlags &flags,
    std::ostream   &out)
  {
    // The gmv format does not support cells that only consist of a single
    // point. It does support the output of point data using the keyword
    // 'tracers' instead of 'nodes' and 'cells', but this output format is
    // currently not implemented.
    AssertThrow(dim > 0, ExcNotImplemented());

    Assert(dim <= 3, ExcNotImplemented());
    AssertThrow(out.fail() == false, ExcIO());

#ifndef DEAL_II_WITH_MPI
    // verify that there are indeed patches to be written out. most of the
    // times, people just forget to call build_patches when there are no
    // patches, so a warning is in order. that said, the assertion is disabled
    // if we support MPI since then it can happen that on the coarsest mesh, a
    // processor simply has no cells it actually owns, and in that case it is
    // legit if there are no patches
    Assert(patches.size() > 0, ExcNoPatches());
#else
    if (patches.empty())
      return;
#endif

    GmvStream          gmv_out(out, flags);
    const unsigned int n_data_sets = data_names.size();
    // check against # of data sets in first patch. checks against all other
    // patches are made in write_gmv_reorder_data_vectors
    Assert((patches[0].data.n_rows() == n_data_sets &&
            !patches[0].points_are_available) ||
             (patches[0].data.n_rows() == n_data_sets + spacedim &&
              patches[0].points_are_available),
           ExcDimensionMismatch(patches[0].points_are_available ?
                                  (n_data_sets + spacedim) :
                                  n_data_sets,
                                patches[0].data.n_rows()));

    //---------------------
    // preamble
    out << "gmvinput ascii" << '\n' << '\n';

    // first count the number of cells and cells for later use
    unsigned int n_nodes;
    unsigned int n_cells;
    std::tie(n_nodes, n_cells) = count_nodes_and_cells(patches);

    // For the format we write here, we need to write all node values relating
    // to one variable at a time. We could in principle do this by looping
    // over all patches and extracting the values corresponding to the one
    // variable we're dealing with right now, and then start the process over
    // for the next variable with another loop over all patches.
    //
    // An easier way is to create a global table that for each variable
    // lists all values. This copying of data vectors can be done in the
    // background while we're already working on vertices and cells,
    // so do this on a separate task and when wanting to write out the
    // data, we wait for that task to finish.
    Threads::Task<std::unique_ptr<Table<2, double>>>
      create_global_data_table_task = Threads::new_task(
        [&patches]() { return create_global_data_table(patches); });

    //-----------------------------
    // first make up a list of used vertices along with their coordinates
    //
    // note that we have to print 3 dimensions
    out << "nodes " << n_nodes << '\n';
    for (unsigned int d = 0; d < spacedim; ++d)
      {
        gmv_out.selected_component = d;
        write_nodes(patches, gmv_out);
        out << '\n';
      }
    gmv_out.selected_component = numbers::invalid_unsigned_int;

    for (unsigned int d = spacedim; d < 3; ++d)
      {
        for (unsigned int i = 0; i < n_nodes; ++i)
          out << "0 ";
        out << '\n';
      }

    //-------------------------------
    // now for the cells. note that vertices are counted from 1 onwards
    out << "cells " << n_cells << '\n';
    write_cells(patches, gmv_out);

    //-------------------------------------
    // data output.
    out << "variable" << '\n';

    // Wait for the reordering to be done and retrieve the reordered data:
    const Table<2, double> data_vectors =
      std::move(*create_global_data_table_task.return_value());

    // then write data. the '1' means: node data (as opposed to cell data, which
    // we do not support explicitly here)
    for (unsigned int data_set = 0; data_set < n_data_sets; ++data_set)
      {
        out << data_names[data_set] << " 1" << '\n';
        std::copy(data_vectors[data_set].begin(),
                  data_vectors[data_set].end(),
                  std::ostream_iterator<double>(out, " "));
        out << '\n' << '\n';
      }



    // end of variable section
    out << "endvars" << '\n';

    // end of output
    out << "endgmv" << '\n';

    // make sure everything now gets to disk
    out.flush();

    // assert the stream is still ok
    AssertThrow(out.fail() == false, ExcIO());
  }



  template <int dim, int spacedim>
  void
  write_tecplot(
    const std::vector<Patch<dim, spacedim>> &patches,
    const std::vector<std::string>          &data_names,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>> &,
    const TecplotFlags &flags,
    std::ostream       &out)
  {
    AssertThrow(out.fail() == false, ExcIO());

    // The FEBLOCK or FEPOINT formats of tecplot only allows full elements (e.g.
    // triangles), not single points. Other tecplot format allow point output,
    // but they are currently not implemented.
    AssertThrow(dim > 0, ExcNotImplemented());

#ifndef DEAL_II_WITH_MPI
    // verify that there are indeed patches to be written out. most of the
    // times, people just forget to call build_patches when there are no
    // patches, so a warning is in order. that said, the assertion is disabled
    // if we support MPI since then it can happen that on the coarsest mesh, a
    // processor simply has no cells it actually owns, and in that case it is
    // legit if there are no patches
    Assert(patches.size() > 0, ExcNoPatches());
#else
    if (patches.empty())
      return;
#endif

    TecplotStream tecplot_out(out, flags);

    const unsigned int n_data_sets = data_names.size();
    // check against # of data sets in first patch. checks against all other
    // patches are made in write_gmv_reorder_data_vectors
    Assert((patches[0].data.n_rows() == n_data_sets &&
            !patches[0].points_are_available) ||
             (patches[0].data.n_rows() == n_data_sets + spacedim &&
              patches[0].points_are_available),
           ExcDimensionMismatch(patches[0].points_are_available ?
                                  (n_data_sets + spacedim) :
                                  n_data_sets,
                                patches[0].data.n_rows()));

    // first count the number of cells and cells for later use
    unsigned int n_nodes;
    unsigned int n_cells;
    std::tie(n_nodes, n_cells) = count_nodes_and_cells(patches);

    //---------
    // preamble
    {
      out
        << "# This file was generated by the deal.II library." << '\n'
        << "# Date =  " << Utilities::System::get_date() << '\n'
        << "# Time =  " << Utilities::System::get_time() << '\n'
        << "#" << '\n'
        << "# For a description of the Tecplot format see the Tecplot documentation."
        << '\n'
        << "#" << '\n';


      out << "Variables=";

      switch (spacedim)
        {
          case 1:
            out << "\"x\"";
            break;
          case 2:
            out << "\"x\", \"y\"";
            break;
          case 3:
            out << "\"x\", \"y\", \"z\"";
            break;
          default:
            DEAL_II_NOT_IMPLEMENTED();
        }

      for (unsigned int data_set = 0; data_set < n_data_sets; ++data_set)
        out << ", \"" << data_names[data_set] << "\"";

      out << '\n';

      out << "zone ";
      if (flags.zone_name)
        out << "t=\"" << flags.zone_name << "\" ";

      if (flags.solution_time >= 0.0)
        out << "strandid=1, solutiontime=" << flags.solution_time << ", ";

      out << "f=feblock, n=" << n_nodes << ", e=" << n_cells
          << ", et=" << tecplot_cell_type[dim] << '\n';
    }


    // For the format we write here, we need to write all node values relating
    // to one variable at a time. We could in principle do this by looping
    // over all patches and extracting the values corresponding to the one
    // variable we're dealing with right now, and then start the process over
    // for the next variable with another loop over all patches.
    //
    // An easier way is to create a global table that for each variable
    // lists all values. This copying of data vectors can be done in the
    // background while we're already working on vertices and cells,
    // so do this on a separate task and when wanting to write out the
    // data, we wait for that task to finish.
    Threads::Task<std::unique_ptr<Table<2, double>>>
      create_global_data_table_task = Threads::new_task(
        [&patches]() { return create_global_data_table(patches); });

    //-----------------------------
    // first make up a list of used vertices along with their coordinates


    for (unsigned int d = 0; d < spacedim; ++d)
      {
        tecplot_out.selected_component = d;
        write_nodes(patches, tecplot_out);
        out << '\n';
      }


    //-------------------------------------
    // data output.
    //
    // Wait for the reordering to be done and retrieve the reordered data:
    const Table<2, double> data_vectors =
      std::move(*create_global_data_table_task.return_value());

    // then write data.
    for (unsigned int data_set = 0; data_set < n_data_sets; ++data_set)
      {
        std::copy(data_vectors[data_set].begin(),
                  data_vectors[data_set].end(),
                  std::ostream_iterator<double>(out, "\n"));
        out << '\n';
      }

    write_cells(patches, tecplot_out);

    // make sure everything now gets to disk
    out.flush();

    // assert the stream is still ok
    AssertThrow(out.fail() == false, ExcIO());
  }



  template <int dim, int spacedim>
  void
  write_vtk(
    const std::vector<Patch<dim, spacedim>> &patches,
    const std::vector<std::string>          &data_names,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>
                   &nonscalar_data_ranges,
    const VtkFlags &flags,
    std::ostream   &out)
  {
    AssertThrow(out.fail() == false, ExcIO());

#ifndef DEAL_II_WITH_MPI
    // verify that there are indeed patches to be written out. most of the
    // times, people just forget to call build_patches when there are no
    // patches, so a warning is in order. that said, the assertion is disabled
    // if we support MPI since then it can happen that on the coarsest mesh, a
    // processor simply has no cells it actually owns, and in that case it is
    // legit if there are no patches
    Assert(patches.size() > 0, ExcNoPatches());
#else
    if (patches.empty())
      return;
#endif

    VtkStream vtk_out(out, flags);

    const unsigned int n_data_sets = data_names.size();
    // check against # of data sets in first patch.
    if (patches[0].points_are_available)
      {
        AssertDimension(n_data_sets + spacedim, patches[0].data.n_rows());
      }
    else
      {
        AssertDimension(n_data_sets, patches[0].data.n_rows());
      }

    //---------------------
    // preamble
    {
      out << "# vtk DataFile Version 3.0" << '\n'
          << "#This file was generated by the deal.II library";
      if (flags.print_date_and_time)
        {
          out << " on " << Utilities::System::get_date() << " at "
              << Utilities::System::get_time();
        }
      else
        out << '.';
      out << '\n' << "ASCII" << '\n';
      // now output the data header
      out << "DATASET UNSTRUCTURED_GRID\n" << '\n';
    }

    // if desired, output time and cycle of the simulation, following the
    // instructions at
    // http://www.visitusers.org/index.php?title=Time_and_Cycle_in_VTK_files
    {
      const unsigned int n_metadata =
        ((flags.cycle != std::numeric_limits<unsigned int>::min() ? 1 : 0) +
         (flags.time != std::numeric_limits<double>::min() ? 1 : 0));
      if (n_metadata > 0)
        {
          out << "FIELD FieldData " << n_metadata << '\n';

          if (flags.cycle != std::numeric_limits<unsigned int>::min())
            {
              out << "CYCLE 1 1 int\n" << flags.cycle << '\n';
            }
          if (flags.time != std::numeric_limits<double>::min())
            {
              out << "TIME 1 1 double\n" << flags.time << '\n';
            }
        }
    }

    // first count the number of cells and cells for later use
    unsigned int n_nodes;
    unsigned int n_cells;
    unsigned int n_points_and_n_cells;
    std::tie(n_nodes, n_cells, n_points_and_n_cells) =
      count_nodes_and_cells_and_points(patches, flags.write_higher_order_cells);

    // For the format we write here, we need to write all node values relating
    // to one variable at a time. We could in principle do this by looping
    // over all patches and extracting the values corresponding to the one
    // variable we're dealing with right now, and then start the process over
    // for the next variable with another loop over all patches.
    //
    // An easier way is to create a global table that for each variable
    // lists all values. This copying of data vectors can be done in the
    // background while we're already working on vertices and cells,
    // so do this on a separate task and when wanting to write out the
    // data, we wait for that task to finish.
    Threads::Task<std::unique_ptr<Table<2, double>>>
      create_global_data_table_task = Threads::new_task(
        [&patches]() { return create_global_data_table(patches); });

    //-----------------------------
    // first make up a list of used vertices along with their coordinates
    //
    // note that we have to print d=1..3 dimensions
    out << "POINTS " << n_nodes << " double" << '\n';
    write_nodes(patches, vtk_out);
    out << '\n';
    //-------------------------------
    // now for the cells
    out << "CELLS " << n_cells << ' ' << n_points_and_n_cells << '\n';
    if (flags.write_higher_order_cells)
      write_high_order_cells(patches, vtk_out, /* legacy_format = */ true);
    else
      write_cells(patches, vtk_out);
    out << '\n';
    // next output the types of the cells. since all cells are the same, this is
    // simple
    out << "CELL_TYPES " << n_cells << '\n';

    // need to distinguish between linear cells, simplex cells (linear or
    // quadratic), and  high order cells
    for (const auto &patch : patches)
      {
        const auto vtk_cell_id =
          extract_vtk_patch_info(patch, flags.write_higher_order_cells);

        for (unsigned int i = 0; i < vtk_cell_id[1]; ++i)
          out << ' ' << vtk_cell_id[0];
      }

    out << '\n';
    //-------------------------------------
    // data output.

    // Wait for the reordering to be done and retrieve the reordered data:
    const Table<2, double> data_vectors =
      std::move(*create_global_data_table_task.return_value());

    // then write data.  the 'POINT_DATA' means: node data (as opposed to cell
    // data, which we do not support explicitly here). all following data sets
    // are point data
    out << "POINT_DATA " << n_nodes << '\n';

    // when writing, first write out all vector data, then handle the scalar
    // data sets that have been left over
    std::vector<bool> data_set_written(n_data_sets, false);
    for (const auto &nonscalar_data_range : nonscalar_data_ranges)
      {
        AssertThrow(std::get<3>(nonscalar_data_range) !=
                      DataComponentInterpretation::component_is_part_of_tensor,
                    ExcMessage(
                      "The VTK writer does not currently support outputting "
                      "tensor data. Use the VTU writer instead."));

        AssertThrow(std::get<1>(nonscalar_data_range) >=
                      std::get<0>(nonscalar_data_range),
                    ExcLowerRange(std::get<1>(nonscalar_data_range),
                                  std::get<0>(nonscalar_data_range)));
        AssertThrow(std::get<1>(nonscalar_data_range) < n_data_sets,
                    ExcIndexRange(std::get<1>(nonscalar_data_range),
                                  0,
                                  n_data_sets));
        AssertThrow(std::get<1>(nonscalar_data_range) + 1 -
                        std::get<0>(nonscalar_data_range) <=
                      3,
                    ExcMessage(
                      "Can't declare a vector with more than 3 components "
                      "in VTK"));

        // mark these components as already written:
        for (unsigned int i = std::get<0>(nonscalar_data_range);
             i <= std::get<1>(nonscalar_data_range);
             ++i)
          data_set_written[i] = true;

        // write the header. concatenate all the component names with double
        // underscores unless a vector name has been specified
        out << "VECTORS ";

        if (!std::get<2>(nonscalar_data_range).empty())
          out << std::get<2>(nonscalar_data_range);
        else
          {
            for (unsigned int i = std::get<0>(nonscalar_data_range);
                 i < std::get<1>(nonscalar_data_range);
                 ++i)
              out << data_names[i] << "__";
            out << data_names[std::get<1>(nonscalar_data_range)];
          }

        out << " double" << '\n';

        // now write data. pad all vectors to have three components
        for (unsigned int n = 0; n < n_nodes; ++n)
          {
            switch (std::get<1>(nonscalar_data_range) -
                    std::get<0>(nonscalar_data_range))
              {
                case 0:
                  out << data_vectors(std::get<0>(nonscalar_data_range), n)
                      << " 0 0" << '\n';
                  break;

                case 1:
                  out << data_vectors(std::get<0>(nonscalar_data_range), n)
                      << ' '
                      << data_vectors(std::get<0>(nonscalar_data_range) + 1, n)
                      << " 0" << '\n';
                  break;
                case 2:
                  out << data_vectors(std::get<0>(nonscalar_data_range), n)
                      << ' '
                      << data_vectors(std::get<0>(nonscalar_data_range) + 1, n)
                      << ' '
                      << data_vectors(std::get<0>(nonscalar_data_range) + 2, n)
                      << '\n';
                  break;

                default:
                  // VTK doesn't support anything else than vectors with 1, 2,
                  // or 3 components
                  DEAL_II_ASSERT_UNREACHABLE();
              }
          }
      }

    // now do the left over scalar data sets
    for (unsigned int data_set = 0; data_set < n_data_sets; ++data_set)
      if (data_set_written[data_set] == false)
        {
          out << "SCALARS " << data_names[data_set] << " double 1" << '\n'
              << "LOOKUP_TABLE default" << '\n';
          std::copy(data_vectors[data_set].begin(),
                    data_vectors[data_set].end(),
                    std::ostream_iterator<double>(out, " "));
          out << '\n';
        }

    // make sure everything now gets to disk
    out.flush();

    // assert the stream is still ok
    AssertThrow(out.fail() == false, ExcIO());
  }


  void
  write_vtu_header(std::ostream &out, const VtkFlags &flags)
  {
    AssertThrow(out.fail() == false, ExcIO());
    out << "<?xml version=\"1.0\" ?> \n";
    out << "<!-- \n";
    out << "# vtk DataFile Version 3.0" << '\n'
        << "#This file was generated by the deal.II library";
    if (flags.print_date_and_time)
      {
        out << " on " << Utilities::System::get_time() << " at "
            << Utilities::System::get_date();
      }
    else
      out << '.';
    out << "\n-->\n";

    if (flags.write_higher_order_cells)
      out << "<VTKFile type=\"UnstructuredGrid\" version=\"2.2\"";
    else
      out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"";
    if (deal_ii_with_zlib &&
        (flags.compression_level != CompressionLevel::plain_text))
      out << " compressor=\"vtkZLibDataCompressor\"";
#ifdef DEAL_II_WORDS_BIGENDIAN
    out << " byte_order=\"BigEndian\"";
#else
    out << " byte_order=\"LittleEndian\"";
#endif
    out << ">";
    out << '\n';
    out << "<UnstructuredGrid>";
    out << '\n';
  }



  void
  write_vtu_footer(std::ostream &out)
  {
    AssertThrow(out.fail() == false, ExcIO());
    out << " </UnstructuredGrid>\n";
    out << "</VTKFile>\n";
  }



  template <int dim, int spacedim>
  void
  write_vtu(
    const std::vector<Patch<dim, spacedim>> &patches,
    const std::vector<std::string>          &data_names,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>
                   &nonscalar_data_ranges,
    const VtkFlags &flags,
    std::ostream   &out)
  {
    write_vtu_header(out, flags);
    write_vtu_main(patches, data_names, nonscalar_data_ranges, flags, out);
    write_vtu_footer(out);

    out << std::flush;
  }


  template <int dim, int spacedim>
  void
  write_vtu_main(
    const std::vector<Patch<dim, spacedim>> &patches,
    const std::vector<std::string>          &data_names,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>
                   &nonscalar_data_ranges,
    const VtkFlags &flags,
    std::ostream   &out)
  {
    AssertThrow(out.fail() == false, ExcIO());

    // If the user provided physical units, make sure that they don't contain
    // quote characters as this would make the VTU file invalid XML and
    // probably lead to all sorts of difficult error messages. Other than that,
    // trust the user that whatever they provide makes sense somehow.
    for (const auto &unit : flags.physical_units)
      {
        (void)unit;
        Assert(
          unit.second.find('\"') == std::string::npos,
          ExcMessage(
            "A physical unit you provided, <" + unit.second +
            ">, contained a quotation mark character. This is not allowed."));
      }

#ifndef DEAL_II_WITH_MPI
    // verify that there are indeed patches to be written out. most of the
    // times, people just forget to call build_patches when there are no
    // patches, so a warning is in order. that said, the assertion is disabled
    // if we support MPI since then it can happen that on the coarsest mesh, a
    // processor simply has no cells it actually owns, and in that case it is
    // legit if there are no patches
    Assert(patches.size() > 0, ExcNoPatches());
#else
    if (patches.empty())
      {
        // we still need to output a valid vtu file, because other CPUs might
        // output data. This is the minimal file that is accepted by paraview
        // and visit. if we remove the field definitions, visit is complaining.
        out << "<Piece NumberOfPoints=\"0\" NumberOfCells=\"0\" >\n"
            << "<Cells>\n"
            << "<DataArray type=\"UInt8\" Name=\"types\"></DataArray>\n"
            << "</Cells>\n"
            << "  <PointData Scalars=\"scalars\">\n";
        std::vector<bool> data_set_written(data_names.size(), false);
        for (const auto &nonscalar_data_range : nonscalar_data_ranges)
          {
            // mark these components as already written:
            for (unsigned int i = std::get<0>(nonscalar_data_range);
                 i <= std::get<1>(nonscalar_data_range);
                 ++i)
              data_set_written[i] = true;

            // write the header. concatenate all the component names with double
            // underscores unless a vector name has been specified
            out << "    <DataArray type=\"Float32\" Name=\"";

            if (!std::get<2>(nonscalar_data_range).empty())
              out << std::get<2>(nonscalar_data_range);
            else
              {
                for (unsigned int i = std::get<0>(nonscalar_data_range);
                     i < std::get<1>(nonscalar_data_range);
                     ++i)
                  out << data_names[i] << "__";
                out << data_names[std::get<1>(nonscalar_data_range)];
              }

            out << "\" NumberOfComponents=\"3\"></DataArray>\n";
          }

        for (unsigned int data_set = 0; data_set < data_names.size();
             ++data_set)
          if (data_set_written[data_set] == false)
            {
              out << "    <DataArray type=\"Float32\" Name=\""
                  << data_names[data_set] << "\"></DataArray>\n";
            }

        out << "  </PointData>\n";
        out << "</Piece>\n";

        out << std::flush;

        return;
      }
#endif

    // first up: metadata
    //
    // if desired, output time and cycle of the simulation, following the
    // instructions at
    // http://www.visitusers.org/index.php?title=Time_and_Cycle_in_VTK_files
    {
      const unsigned int n_metadata =
        ((flags.cycle != std::numeric_limits<unsigned int>::min() ? 1 : 0) +
         (flags.time != std::numeric_limits<double>::min() ? 1 : 0));
      if (n_metadata > 0)
        out << "<FieldData>\n";

      if (flags.cycle != std::numeric_limits<unsigned int>::min())
        {
          out
            << "<DataArray type=\"Float32\" Name=\"CYCLE\" NumberOfTuples=\"1\" format=\"ascii\">"
            << flags.cycle << "</DataArray>\n";
        }
      if (flags.time != std::numeric_limits<double>::min())
        {
          out
            << "<DataArray type=\"Float32\" Name=\"TIME\" NumberOfTuples=\"1\" format=\"ascii\">"
            << flags.time << "</DataArray>\n";
        }

      if (n_metadata > 0)
        out << "</FieldData>\n";
    }


    const unsigned int n_data_sets = data_names.size();
    // check against # of data sets in first patch. checks against all other
    // patches are made in write_gmv_reorder_data_vectors
    if (patches[0].points_are_available)
      {
        AssertDimension(n_data_sets + spacedim, patches[0].data.n_rows());
      }
    else
      {
        AssertDimension(n_data_sets, patches[0].data.n_rows());
      }

    const char *ascii_or_binary =
      (deal_ii_with_zlib &&
       (flags.compression_level != CompressionLevel::plain_text)) ?
        "binary" :
        "ascii";


    // first count the number of cells and cells for later use
    unsigned int n_nodes;
    unsigned int n_cells;
    std::tie(n_nodes, n_cells, std::ignore) =
      count_nodes_and_cells_and_points(patches, flags.write_higher_order_cells);

    // -----------------
    // In the following, let us first set up a number of lambda functions that
    // will be used in building the different parts of the VTU file. We will
    // later call them in turn on different tasks.
    // first make up a list of used vertices along with their coordinates
    const auto stringize_vertex_information = [&patches,
                                               &flags,
                                               output_precision =
                                                 out.precision(),
                                               ascii_or_binary]() {
      std::ostringstream o;
      o << "  <Points>\n";
      o << "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\""
        << ascii_or_binary << "\">\n";
      const std::vector<Point<spacedim>> node_positions =
        get_node_positions(patches);

      // VTK/VTU always wants to see three coordinates, even if we are
      // in 1d or 2d. So pad node positions with zeros as appropriate.
      std::vector<float> node_coordinates_3d;
      node_coordinates_3d.reserve(node_positions.size() * 3);
      for (const auto &node_position : node_positions)
        {
          for (unsigned int d = 0; d < 3; ++d)
            if (d < spacedim)
              node_coordinates_3d.emplace_back(node_position[d]);
            else
              node_coordinates_3d.emplace_back(0.0f);
        }
      o << vtu_stringize_array(node_coordinates_3d,
                               flags.compression_level,
                               output_precision)
        << '\n';
      o << "    </DataArray>\n";
      o << "  </Points>\n\n";

      return o.str();
    };


    //-------------------------------
    // Now for the cells. The first part of this is how vertices
    // build cells.
    const auto stringize_cell_to_vertex_information = [&patches,
                                                       &flags,
                                                       ascii_or_binary,
                                                       output_precision =
                                                         out.precision()]() {
      std::ostringstream o;

      o << "  <Cells>\n";
      o << "    <DataArray type=\"Int32\" Name=\"connectivity\" format=\""
        << ascii_or_binary << "\">\n";

      std::vector<std::int32_t> cells;
      Assert(dim <= 3, ExcNotImplemented());

      unsigned int first_vertex_of_patch = 0;

      for (const auto &patch : patches)
        {
          // First treat a slight oddball case: For triangles and tetrahedra,
          // the case with n_subdivisions==2 is treated as if the cell was
          // output as a single, quadratic, cell rather than as one would
          // expect as 4 sub-cells (for triangles; and the corresponding
          // number of sub-cells for tetrahedra). This is courtesy of some
          // special-casing in the function extract_vtk_patch_info().
          if ((dim >= 2) &&
              (patch.reference_cell == ReferenceCells::get_simplex<dim>()) &&
              (patch.n_subdivisions == 2))
            {
              const unsigned int n_points = patch.data.n_cols();
              Assert((dim == 2 && n_points == 6) ||
                       (dim == 3 && n_points == 10),
                     ExcInternalError());

              if (deal_ii_with_zlib &&
                  (flags.compression_level !=
                   DataOutBase::CompressionLevel::plain_text))
                {
                  for (unsigned int i = 0; i < n_points; ++i)
                    cells.push_back(first_vertex_of_patch + i);
                }
              else
                {
                  for (unsigned int i = 0; i < n_points; ++i)
                    o << '\t' << first_vertex_of_patch + i;
                  o << '\n';
                }

              first_vertex_of_patch += n_points;
            }
          // Then treat all of the other non-hypercube cases since they can
          // currently not be subdivided (into sub-cells, or into higher-order
          // cells):
          else if (patch.reference_cell != ReferenceCells::get_hypercube<dim>())
            {
              Assert(patch.n_subdivisions == 1, ExcNotImplemented());

              const unsigned int n_points = patch.data.n_cols();

              if (deal_ii_with_zlib &&
                  (flags.compression_level !=
                   DataOutBase::CompressionLevel::plain_text))
                {
                  for (unsigned int i = 0; i < n_points; ++i)
                    cells.push_back(
                      first_vertex_of_patch +
                      patch.reference_cell.vtk_vertex_to_deal_vertex(i));
                }
              else
                {
                  for (unsigned int i = 0; i < n_points; ++i)
                    o << '\t'
                      << (first_vertex_of_patch +
                          patch.reference_cell.vtk_vertex_to_deal_vertex(i));
                  o << '\n';
                }

              first_vertex_of_patch += n_points;
            }
          else // a hypercube cell
            {
              const unsigned int n_subdivisions         = patch.n_subdivisions;
              const unsigned int n_points_per_direction = n_subdivisions + 1;

              std::vector<unsigned> local_vertex_order;

              // Output the current state of the local_vertex_order array,
              // then clear it:
              const auto flush_current_cell = [&flags,
                                               &o,
                                               &cells,
                                               first_vertex_of_patch,
                                               &local_vertex_order]() {
                if (deal_ii_with_zlib &&
                    (flags.compression_level !=
                     DataOutBase::CompressionLevel::plain_text))
                  {
                    for (const auto &c : local_vertex_order)
                      cells.push_back(first_vertex_of_patch + c);
                  }
                else
                  {
                    for (const auto &c : local_vertex_order)
                      o << '\t' << first_vertex_of_patch + c;
                    o << '\n';
                  }

                local_vertex_order.clear();
              };

              if (flags.write_higher_order_cells == false)
                {
                  local_vertex_order.reserve(Utilities::fixed_power<dim>(2));

                  switch (dim)
                    {
                      case 0:
                        {
                          local_vertex_order.emplace_back(0);
                          flush_current_cell();
                          break;
                        }

                      case 1:
                        {
                          for (unsigned int i1 = 0; i1 < n_subdivisions; ++i1)
                            {
                              const unsigned int starting_offset = i1;
                              local_vertex_order.emplace_back(starting_offset);
                              local_vertex_order.emplace_back(starting_offset +
                                                              1);
                              flush_current_cell();
                            }
                          break;
                        }

                      case 2:
                        {
                          for (unsigned int i2 = 0; i2 < n_subdivisions; ++i2)
                            for (unsigned int i1 = 0; i1 < n_subdivisions; ++i1)
                              {
                                const unsigned int starting_offset =
                                  i2 * n_points_per_direction + i1;
                                local_vertex_order.emplace_back(
                                  starting_offset);
                                local_vertex_order.emplace_back(
                                  starting_offset + 1);
                                local_vertex_order.emplace_back(
                                  starting_offset + n_points_per_direction + 1);
                                local_vertex_order.emplace_back(
                                  starting_offset + n_points_per_direction);
                                flush_current_cell();
                              }
                          break;
                        }

                      case 3:
                        {
                          for (unsigned int i3 = 0; i3 < n_subdivisions; ++i3)
                            for (unsigned int i2 = 0; i2 < n_subdivisions; ++i2)
                              for (unsigned int i1 = 0; i1 < n_subdivisions;
                                   ++i1)
                                {
                                  const unsigned int starting_offset =
                                    i3 * n_points_per_direction *
                                      n_points_per_direction +
                                    i2 * n_points_per_direction + i1;
                                  local_vertex_order.emplace_back(
                                    starting_offset);
                                  local_vertex_order.emplace_back(
                                    starting_offset + 1);
                                  local_vertex_order.emplace_back(
                                    starting_offset + n_points_per_direction +
                                    1);
                                  local_vertex_order.emplace_back(
                                    starting_offset + n_points_per_direction);
                                  local_vertex_order.emplace_back(
                                    starting_offset + n_points_per_direction *
                                                        n_points_per_direction);
                                  local_vertex_order.emplace_back(
                                    starting_offset +
                                    n_points_per_direction *
                                      n_points_per_direction +
                                    1);
                                  local_vertex_order.emplace_back(
                                    starting_offset +
                                    n_points_per_direction *
                                      n_points_per_direction +
                                    n_points_per_direction + 1);
                                  local_vertex_order.emplace_back(
                                    starting_offset +
                                    n_points_per_direction *
                                      n_points_per_direction +
                                    n_points_per_direction);
                                  flush_current_cell();
                                }
                          break;
                        }

                      default:
                        DEAL_II_NOT_IMPLEMENTED();
                    }
                }
              else // use higher-order output
                {
                  local_vertex_order.resize(
                    Utilities::fixed_power<dim>(n_points_per_direction));

                  switch (dim)
                    {
                      case 0:
                        {
                          Assert(false,
                                 ExcMessage(
                                   "Point-like cells should not be possible "
                                   "when writing higher-order cells."));
                          break;
                        }
                      case 1:
                        {
                          for (unsigned int i1 = 0; i1 < n_subdivisions + 1;
                               ++i1)
                            {
                              const unsigned int local_index = i1;
                              const unsigned int connectivity_index =
                                patch.reference_cell
                                  .template vtk_lexicographic_to_node_index<1>(
                                    {{i1}},
                                    {{n_subdivisions}},
                                    /* use VTU, not VTK: */ false);
                              local_vertex_order[connectivity_index] =
                                local_index;
                            }
                          flush_current_cell();

                          break;
                        }
                      case 2:
                        {
                          for (unsigned int i2 = 0; i2 < n_subdivisions + 1;
                               ++i2)
                            for (unsigned int i1 = 0; i1 < n_subdivisions + 1;
                                 ++i1)
                              {
                                const unsigned int local_index =
                                  i2 * n_points_per_direction + i1;
                                const unsigned int connectivity_index =
                                  patch.reference_cell
                                    .template vtk_lexicographic_to_node_index<
                                      2>({{i1, i2}},
                                         {{n_subdivisions, n_subdivisions}},
                                         /* use VTU, not VTK: */ false);
                                local_vertex_order[connectivity_index] =
                                  local_index;
                              }
                          flush_current_cell();

                          break;
                        }
                      case 3:
                        {
                          for (unsigned int i3 = 0; i3 < n_subdivisions + 1;
                               ++i3)
                            for (unsigned int i2 = 0; i2 < n_subdivisions + 1;
                                 ++i2)
                              for (unsigned int i1 = 0; i1 < n_subdivisions + 1;
                                   ++i1)
                                {
                                  const unsigned int local_index =
                                    i3 * n_points_per_direction *
                                      n_points_per_direction +
                                    i2 * n_points_per_direction + i1;
                                  const unsigned int connectivity_index =
                                    patch.reference_cell
                                      .template vtk_lexicographic_to_node_index<
                                        3>({{i1, i2, i3}},
                                           {{n_subdivisions,
                                             n_subdivisions,
                                             n_subdivisions}},
                                           /* use VTU, not VTK: */ false);
                                  local_vertex_order[connectivity_index] =
                                    local_index;
                                }

                          flush_current_cell();
                          break;
                        }
                      default:
                        DEAL_II_NOT_IMPLEMENTED();
                    }
                }

              // Finally update the number of the first vertex of this
              // patch
              first_vertex_of_patch +=
                Utilities::fixed_power<dim>(patch.n_subdivisions + 1);
            }
        }

      // Flush the 'cells' object we created herein.
      if (deal_ii_with_zlib && (flags.compression_level !=
                                DataOutBase::CompressionLevel::plain_text))
        {
          o << vtu_stringize_array(cells,
                                   flags.compression_level,
                                   output_precision)
            << '\n';
        }
      o << "    </DataArray>\n";

      return o.str();
    };


    //-------------------------------
    // The second part of cell information is the offsets in
    // the array built by the previous lambda function that indicate
    // individual cells.
    //
    // Note that this separates XML VTU format from the VTK format; the latter
    // puts the number of nodes per cell in front of the connectivity list for
    // each cell, whereas the VTU format uses one large list of vertex indices
    // and a separate array of offsets.
    //
    // The third piece to cell information is that we need to
    // output the types of the cells.
    //
    // The following function does both of these pieces.
    const auto stringize_cell_offset_and_type_information =
      [&patches,
       &flags,
       ascii_or_binary,
       n_cells,
       output_precision = out.precision()]() {
        std::ostringstream o;

        o << "    <DataArray type=\"Int32\" Name=\"offsets\" format=\""
          << ascii_or_binary << "\">\n";

        std::vector<std::int32_t> offsets;
        offsets.reserve(n_cells);

        // std::uint8_t might be an alias to unsigned char which is then not
        // printed as ascii integers
        std::vector<unsigned int> cell_types;
        cell_types.reserve(n_cells);

        unsigned int first_vertex_of_patch = 0;

        for (const auto &patch : patches)
          {
            const auto vtk_cell_id =
              extract_vtk_patch_info(patch, flags.write_higher_order_cells);

            for (unsigned int i = 0; i < vtk_cell_id[1]; ++i)
              {
                cell_types.push_back(vtk_cell_id[0]);
                first_vertex_of_patch += vtk_cell_id[2];
                offsets.push_back(first_vertex_of_patch);
              }
          }

        o << vtu_stringize_array(offsets,
                                 flags.compression_level,
                                 output_precision);
        o << '\n';
        o << "    </DataArray>\n";

        o << "    <DataArray type=\"UInt8\" Name=\"types\" format=\""
          << ascii_or_binary << "\">\n";

        if (deal_ii_with_zlib &&
            (flags.compression_level != CompressionLevel::plain_text))
          {
            std::vector<std::uint8_t> cell_types_uint8_t(cell_types.size());
            for (unsigned int i = 0; i < cell_types.size(); ++i)
              cell_types_uint8_t[i] = static_cast<std::uint8_t>(cell_types[i]);

            o << vtu_stringize_array(cell_types_uint8_t,
                                     flags.compression_level,
                                     output_precision);
          }
        else
          {
            o << vtu_stringize_array(cell_types,
                                     flags.compression_level,
                                     output_precision);
          }

        o << '\n';
        o << "    </DataArray>\n";
        o << "  </Cells>\n";

        return o.str();
      };


    //-------------------------------------
    // data output.

    const auto stringize_nonscalar_data_range =
      [&flags,
       &data_names,
       ascii_or_binary,
       n_data_sets,
       n_nodes,
       output_precision = out.precision()](const Table<2, float> &data_vectors,
                                           const auto            &range) {
        std::ostringstream o;

        const auto  first_component = std::get<0>(range);
        const auto  last_component  = std::get<1>(range);
        const auto &name            = std::get<2>(range);
        const bool  is_tensor =
          (std::get<3>(range) ==
           DataComponentInterpretation::component_is_part_of_tensor);
        const unsigned int n_components = (is_tensor ? 9 : 3);
        AssertThrow(last_component >= first_component,
                    ExcLowerRange(last_component, first_component));
        AssertThrow(last_component < n_data_sets,
                    ExcIndexRange(last_component, 0, n_data_sets));
        if (is_tensor)
          {
            AssertThrow((last_component + 1 - first_component <= 9),
                        ExcMessage(
                          "Can't declare a tensor with more than 9 components "
                          "in VTK/VTU format."));
          }
        else
          {
            AssertThrow((last_component + 1 - first_component <= 3),
                        ExcMessage(
                          "Can't declare a vector with more than 3 components "
                          "in VTK/VTU format."));
          }

        // write the header. concatenate all the component names with double
        // underscores unless a vector name has been specified
        o << "    <DataArray type=\"Float32\" Name=\"";

        if (!name.empty())
          o << name;
        else
          {
            for (unsigned int i = first_component; i < last_component; ++i)
              o << data_names[i] << "__";
            o << data_names[last_component];
          }

        o << "\" NumberOfComponents=\"" << n_components << "\" format=\""
          << ascii_or_binary << "\"";
        // If present, also list the physical units for this quantity. Look
        // this up for either the name of the whole vector/tensor, or if that
        // isn't listed, via its first component.
        if (!name.empty())
          {
            if (flags.physical_units.find(name) != flags.physical_units.end())
              o << " units=\"" << flags.physical_units.at(name) << "\"";
          }
        else
          {
            if (flags.physical_units.find(data_names[first_component]) !=
                flags.physical_units.end())
              o << " units=\""
                << flags.physical_units.at(data_names[first_component]) << "\"";
          }
        o << ">\n";

        // now write data. pad all vectors to have three components
        std::vector<float> data;
        data.reserve(n_nodes * n_components);

        for (unsigned int n = 0; n < n_nodes; ++n)
          {
            if (!is_tensor)
              {
                switch (last_component - first_component)
                  {
                    case 0:
                      data.push_back(data_vectors(first_component, n));
                      data.push_back(0);
                      data.push_back(0);
                      break;

                    case 1:
                      data.push_back(data_vectors(first_component, n));
                      data.push_back(data_vectors(first_component + 1, n));
                      data.push_back(0);
                      break;

                    case 2:
                      data.push_back(data_vectors(first_component, n));
                      data.push_back(data_vectors(first_component + 1, n));
                      data.push_back(data_vectors(first_component + 2, n));
                      break;

                    default:
                      // Anything else is not yet implemented
                      DEAL_II_ASSERT_UNREACHABLE();
                  }
              }
            else
              {
                Tensor<2, 3> vtk_data;
                vtk_data = 0.;

                const unsigned int size = last_component - first_component + 1;
                if (size == 1)
                  // 1d, 1 element
                  {
                    vtk_data[0][0] = data_vectors(first_component, n);
                  }
                else if (size == 4)
                  // 2d, 4 elements
                  {
                    for (unsigned int c = 0; c < size; ++c)
                      {
                        const auto ind =
                          Tensor<2, 2>::unrolled_to_component_indices(c);
                        vtk_data[ind[0]][ind[1]] =
                          data_vectors(first_component + c, n);
                      }
                  }
                else if (size == 9)
                  // 3d 9 elements
                  {
                    for (unsigned int c = 0; c < size; ++c)
                      {
                        const auto ind =
                          Tensor<2, 3>::unrolled_to_component_indices(c);
                        vtk_data[ind[0]][ind[1]] =
                          data_vectors(first_component + c, n);
                      }
                  }
                else
                  {
                    DEAL_II_ASSERT_UNREACHABLE();
                  }

                // now put the tensor into data
                // note we pad with zeros because VTK format always wants to
                // see a 3x3 tensor, regardless of dimension
                for (unsigned int i = 0; i < 3; ++i)
                  for (unsigned int j = 0; j < 3; ++j)
                    data.push_back(vtk_data[i][j]);
              }
          } // loop over nodes

        o << vtu_stringize_array(data,
                                 flags.compression_level,
                                 output_precision);
        o << '\n';
        o << "    </DataArray>\n";

        return o.str();
      };

    const auto stringize_scalar_data_set =
      [&flags,
       &data_names,
       ascii_or_binary,
       output_precision = out.precision()](const Table<2, float> &data_vectors,
                                           const unsigned int     data_set) {
        std::ostringstream o;

        o << "    <DataArray type=\"Float32\" Name=\"" << data_names[data_set]
          << "\" format=\"" << ascii_or_binary << "\"";
        // If present, also list the physical units for this quantity.
        if (flags.physical_units.find(data_names[data_set]) !=
            flags.physical_units.end())
          o << " units=\"" << flags.physical_units.at(data_names[data_set])
            << "\"";

        o << ">\n";

        const std::vector<float> data(data_vectors[data_set].begin(),
                                      data_vectors[data_set].end());
        o << vtu_stringize_array(data,
                                 flags.compression_level,
                                 output_precision);
        o << '\n';
        o << "    </DataArray>\n";

        return o.str();
      };


    // For the format we write here, we need to write all node values relating
    // to one variable at a time. We could in principle do this by looping
    // over all patches and extracting the values corresponding to the one
    // variable we're dealing with right now, and then start the process over
    // for the next variable with another loop over all patches.
    //
    // An easier way is to create a global table that for each variable
    // lists all values. This copying of data vectors can be done in the
    // background while we're already working on vertices and cells,
    // so do this on a separate task and when wanting to write out the
    // data, we wait for that task to finish.
    Threads::Task<std::unique_ptr<Table<2, float>>>
      create_global_data_table_task = Threads::new_task([&patches]() {
        return create_global_data_table<dim, spacedim, float>(patches);
      });

    // -----------------------------
    // Now finally get around to actually doing anything. Let's start with
    // running the first three tasks generating the vertex and cell information:
    Threads::TaskGroup<std::string> mesh_tasks;
    mesh_tasks += Threads::new_task(stringize_vertex_information);
    mesh_tasks += Threads::new_task(stringize_cell_to_vertex_information);
    mesh_tasks += Threads::new_task(stringize_cell_offset_and_type_information);

    // For what follows, we have to have the reordered data available. So wait
    // for that task to conclude and get the resulting data table:
    const Table<2, float> data_vectors =
      std::move(*create_global_data_table_task.return_value());

    // Then create the strings for the actual values of the solution vectors,
    // again on separate tasks:
    Threads::TaskGroup<std::string> data_tasks;
    // When writing, first write out all vector and tensor data
    std::vector<bool> data_set_handled(n_data_sets, false);
    for (const auto &range : nonscalar_data_ranges)
      {
        // Mark these components as already handled:
        const auto first_component = std::get<0>(range);
        const auto last_component  = std::get<1>(range);
        for (unsigned int i = first_component; i <= last_component; ++i)
          data_set_handled[i] = true;

        data_tasks += Threads::new_task([&, range]() {
          return stringize_nonscalar_data_range(data_vectors, range);
        });
      }

    // Now do the left over scalar data sets
    for (unsigned int data_set = 0; data_set < n_data_sets; ++data_set)
      if (data_set_handled[data_set] == false)
        {
          data_tasks += Threads::new_task([&, data_set]() {
            return stringize_scalar_data_set(data_vectors, data_set);
          });
        }

    // Alright, all tasks are now running. Wait for their conclusion and output
    // all of the data they have produced:
    out << "<Piece NumberOfPoints=\"" << n_nodes << "\" NumberOfCells=\""
        << n_cells << "\" >\n";
    for (const auto &s : mesh_tasks.return_values())
      out << s;
    out << "  <PointData Scalars=\"scalars\">\n";
    for (const auto &s : data_tasks.return_values())
      out << s;
    out << "  </PointData>\n";
    out << " </Piece>\n";

    // make sure everything now gets to disk
    out.flush();

    // assert the stream is still ok
    AssertThrow(out.fail() == false, ExcIO());
  }



  void
  write_pvtu_record(
    std::ostream                   &out,
    const std::vector<std::string> &piece_names,
    const std::vector<std::string> &data_names,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>
                   &nonscalar_data_ranges,
    const VtkFlags &flags)
  {
    AssertThrow(out.fail() == false, ExcIO());

    // If the user provided physical units, make sure that they don't contain
    // quote characters as this would make the VTU file invalid XML and
    // probably lead to all sorts of difficult error messages. Other than that,
    // trust the user that whatever they provide makes sense somehow.
    for (const auto &unit : flags.physical_units)
      {
        (void)unit;
        Assert(
          unit.second.find('\"') == std::string::npos,
          ExcMessage(
            "A physical unit you provided, <" + unit.second +
            ">, contained a quotation mark character. This is not allowed."));
      }

    const unsigned int n_data_sets = data_names.size();

    out << "<?xml version=\"1.0\"?>\n";

    out << "<!--\n";
    out << "#This file was generated by the deal.II library"
        << " on " << Utilities::System::get_date() << " at "
        << Utilities::System::get_time() << "\n-->\n";

    out
      << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out << "  <PUnstructuredGrid GhostLevel=\"0\">\n";
    out << "    <PPointData Scalars=\"scalars\">\n";

    // We need to output in the same order as the write_vtu function does:
    std::vector<bool> data_set_written(n_data_sets, false);
    for (const auto &nonscalar_data_range : nonscalar_data_ranges)
      {
        const auto first_component = std::get<0>(nonscalar_data_range);
        const auto last_component  = std::get<1>(nonscalar_data_range);
        const bool is_tensor =
          (std::get<3>(nonscalar_data_range) ==
           DataComponentInterpretation::component_is_part_of_tensor);
        const unsigned int n_components = (is_tensor ? 9 : 3);
        AssertThrow(last_component >= first_component,
                    ExcLowerRange(last_component, first_component));
        AssertThrow(last_component < n_data_sets,
                    ExcIndexRange(last_component, 0, n_data_sets));
        if (is_tensor)
          {
            AssertThrow((last_component + 1 - first_component <= 9),
                        ExcMessage(
                          "Can't declare a tensor with more than 9 components "
                          "in VTK"));
          }
        else
          {
            Assert((last_component + 1 - first_component <= 3),
                   ExcMessage(
                     "Can't declare a vector with more than 3 components "
                     "in VTK"));
          }

        // mark these components as already written:
        for (unsigned int i = std::get<0>(nonscalar_data_range);
             i <= std::get<1>(nonscalar_data_range);
             ++i)
          data_set_written[i] = true;

        // write the header. concatenate all the component names with double
        // underscores unless a vector name has been specified
        out << "    <PDataArray type=\"Float32\" Name=\"";

        const std::string &name = std::get<2>(nonscalar_data_range);
        if (!name.empty())
          out << name;
        else
          {
            for (unsigned int i = std::get<0>(nonscalar_data_range);
                 i < std::get<1>(nonscalar_data_range);
                 ++i)
              out << data_names[i] << "__";
            out << data_names[std::get<1>(nonscalar_data_range)];
          }

        out << "\" NumberOfComponents=\"" << n_components
            << "\" format=\"ascii\"";
        // If present, also list the physical units for this quantity. Look this
        // up for either the name of the whole vector/tensor, or if that isn't
        // listed, via its first component.
        if (!name.empty())
          {
            if (flags.physical_units.find(name) != flags.physical_units.end())
              out << " units=\"" << flags.physical_units.at(name) << "\"";
          }
        else
          {
            if (flags.physical_units.find(
                  data_names[std::get<1>(nonscalar_data_range)]) !=
                flags.physical_units.end())
              out << " units=\""
                  << flags.physical_units.at(
                       data_names[std::get<1>(nonscalar_data_range)])
                  << "\"";
          }

        out << "/>\n";
      }

    // Now for the scalar fields
    for (unsigned int data_set = 0; data_set < n_data_sets; ++data_set)
      if (data_set_written[data_set] == false)
        {
          out << "    <PDataArray type=\"Float32\" Name=\""
              << data_names[data_set] << "\" format=\"ascii\"";

          if (flags.physical_units.find(data_names[data_set]) !=
              flags.physical_units.end())
            out << " units=\"" << flags.physical_units.at(data_names[data_set])
                << "\"";

          out << "/>\n";
        }

    out << "    </PPointData>\n";

    out << "    <PPoints>\n";
    out << "      <PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>\n";
    out << "    </PPoints>\n";

    for (const auto &piece_name : piece_names)
      out << "    <Piece Source=\"" << piece_name << "\"/>\n";

    out << "  </PUnstructuredGrid>\n";
    out << "</VTKFile>\n";

    out.flush();

    // assert the stream is still ok
    AssertThrow(out.fail() == false, ExcIO());
  }



  void
  write_pvd_record(
    std::ostream                                      &out,
    const std::vector<std::pair<double, std::string>> &times_and_names)
  {
    AssertThrow(out.fail() == false, ExcIO());

    out << "<?xml version=\"1.0\"?>\n";

    out << "<!--\n";
    out << "#This file was generated by the deal.II library"
        << " on " << Utilities::System::get_date() << " at "
        << Utilities::System::get_time() << "\n-->\n";

    out
      << "<VTKFile type=\"Collection\" version=\"0.1\" ByteOrder=\"LittleEndian\">\n";
    out << "  <Collection>\n";

    std::streamsize ss = out.precision();
    out.precision(12);

    for (const auto &time_and_name : times_and_names)
      out << "    <DataSet timestep=\"" << time_and_name.first
          << "\" group=\"\" part=\"0\" file=\"" << time_and_name.second
          << "\"/>\n";

    out << "  </Collection>\n";
    out << "</VTKFile>\n";

    out.flush();
    out.precision(ss);

    AssertThrow(out.fail() == false, ExcIO());
  }



  void
  write_visit_record(std::ostream                   &out,
                     const std::vector<std::string> &piece_names)
  {
    out << "!NBLOCKS " << piece_names.size() << '\n';
    for (const auto &piece_name : piece_names)
      out << piece_name << '\n';

    out << std::flush;
  }



  void
  write_visit_record(std::ostream                                &out,
                     const std::vector<std::vector<std::string>> &piece_names)
  {
    AssertThrow(out.fail() == false, ExcIO());

    if (piece_names.empty())
      return;

    const double nblocks = piece_names[0].size();
    Assert(nblocks > 0,
           ExcMessage("piece_names should be a vector of nonempty vectors."));

    out << "!NBLOCKS " << nblocks << '\n';
    for (const auto &domain : piece_names)
      {
        Assert(domain.size() == nblocks,
               ExcMessage(
                 "piece_names should be a vector of equal sized vectors."));
        for (const auto &subdomain : domain)
          out << subdomain << '\n';
      }

    out << std::flush;
  }



  void
  write_visit_record(
    std::ostream &out,
    const std::vector<std::pair<double, std::vector<std::string>>>
      &times_and_piece_names)
  {
    AssertThrow(out.fail() == false, ExcIO());

    if (times_and_piece_names.empty())
      return;

    const double nblocks = times_and_piece_names[0].second.size();
    Assert(
      nblocks > 0,
      ExcMessage(
        "time_and_piece_names should contain nonempty vectors of filenames for every timestep."));

    for (const auto &domain : times_and_piece_names)
      out << "!TIME " << domain.first << '\n';

    out << "!NBLOCKS " << nblocks << '\n';
    for (const auto &domain : times_and_piece_names)
      {
        Assert(domain.second.size() == nblocks,
               ExcMessage(
                 "piece_names should be a vector of equal sized vectors."));
        for (const auto &subdomain : domain.second)
          out << subdomain << '\n';
      }

    out << std::flush;
  }



  template <int dim, int spacedim>
  void
  write_svg(
    const std::vector<Patch<dim, spacedim>> &,
    const std::vector<std::string> &,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>> &,
    const SvgFlags &,
    std::ostream &)
  {
    DEAL_II_NOT_IMPLEMENTED();
  }

  template <int spacedim>
  void
  write_svg(
    const std::vector<Patch<2, spacedim>> &patches,
    const std::vector<std::string> & /*data_names*/,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>
      & /*nonscalar_data_ranges*/,
    const SvgFlags &flags,
    std::ostream   &out)
  {
    const unsigned int height = flags.height;
    unsigned int       width  = flags.width;

    // margin around the plotted area
    unsigned int margin_in_percent = 0;
    if (flags.margin)
      margin_in_percent = 5;


    // determine the bounding box in the model space
    double x_dimension, y_dimension, z_dimension;

    const auto &first_patch = patches[0];

    unsigned int       n_subdivisions = first_patch.n_subdivisions;
    unsigned int       n              = n_subdivisions + 1;
    const unsigned int d1             = 1;
    const unsigned int d2             = n;

    Point<spacedim>                projected_point;
    std::array<Point<spacedim>, 4> projected_points;

    Point<2>                projection_decomposition;
    std::array<Point<2>, 4> projection_decompositions;

    projected_point =
      get_equispaced_location(first_patch, {0, 0}, n_subdivisions);

    if (first_patch.data.n_rows() != 0)
      {
        AssertIndexRange(flags.height_vector, first_patch.data.n_rows());
      }

    double x_min = projected_point[0];
    double x_max = x_min;
    double y_min = projected_point[1];
    double y_max = y_min;
    double z_min = first_patch.data.n_rows() != 0 ?
                     first_patch.data(flags.height_vector, 0) :
                     0;
    double z_max = z_min;

    // iterate over the patches
    for (const auto &patch : patches)
      {
        n_subdivisions = patch.n_subdivisions;
        n              = n_subdivisions + 1;

        for (unsigned int i2 = 0; i2 < n_subdivisions; ++i2)
          {
            for (unsigned int i1 = 0; i1 < n_subdivisions; ++i1)
              {
                projected_points[0] =
                  get_equispaced_location(patch, {i1, i2}, n_subdivisions);
                projected_points[1] =
                  get_equispaced_location(patch, {i1 + 1, i2}, n_subdivisions);
                projected_points[2] =
                  get_equispaced_location(patch, {i1, i2 + 1}, n_subdivisions);
                projected_points[3] = get_equispaced_location(patch,
                                                              {i1 + 1, i2 + 1},
                                                              n_subdivisions);

                x_min = std::min(x_min, projected_points[0][0]);
                x_min = std::min(x_min, projected_points[1][0]);
                x_min = std::min(x_min, projected_points[2][0]);
                x_min = std::min(x_min, projected_points[3][0]);

                x_max = std::max(x_max, projected_points[0][0]);
                x_max = std::max(x_max, projected_points[1][0]);
                x_max = std::max(x_max, projected_points[2][0]);
                x_max = std::max(x_max, projected_points[3][0]);

                y_min = std::min(y_min, projected_points[0][1]);
                y_min = std::min(y_min, projected_points[1][1]);
                y_min = std::min(y_min, projected_points[2][1]);
                y_min = std::min(y_min, projected_points[3][1]);

                y_max = std::max(y_max, projected_points[0][1]);
                y_max = std::max(y_max, projected_points[1][1]);
                y_max = std::max(y_max, projected_points[2][1]);
                y_max = std::max(y_max, projected_points[3][1]);

                Assert((flags.height_vector < patch.data.n_rows()) ||
                         patch.data.n_rows() == 0,
                       ExcIndexRange(flags.height_vector,
                                     0,
                                     patch.data.n_rows()));

                z_min = std::min<double>(z_min,
                                         patch.data(flags.height_vector,
                                                    i1 * d1 + i2 * d2));
                z_min = std::min<double>(z_min,
                                         patch.data(flags.height_vector,
                                                    (i1 + 1) * d1 + i2 * d2));
                z_min = std::min<double>(z_min,
                                         patch.data(flags.height_vector,
                                                    i1 * d1 + (i2 + 1) * d2));
                z_min =
                  std::min<double>(z_min,
                                   patch.data(flags.height_vector,
                                              (i1 + 1) * d1 + (i2 + 1) * d2));

                z_max = std::max<double>(z_max,
                                         patch.data(flags.height_vector,
                                                    i1 * d1 + i2 * d2));
                z_max = std::max<double>(z_max,
                                         patch.data(flags.height_vector,
                                                    (i1 + 1) * d1 + i2 * d2));
                z_max = std::max<double>(z_max,
                                         patch.data(flags.height_vector,
                                                    i1 * d1 + (i2 + 1) * d2));
                z_max =
                  std::max<double>(z_max,
                                   patch.data(flags.height_vector,
                                              (i1 + 1) * d1 + (i2 + 1) * d2));
              }
          }
      }

    x_dimension = x_max - x_min;
    y_dimension = y_max - y_min;
    z_dimension = z_max - z_min;


    // set initial camera position
    Point<3> camera_position;
    Point<3> camera_direction;
    Point<3> camera_horizontal;
    float    camera_focus = 0;

    // translate camera from the origin to the initial position
    camera_position[0] = 0.;
    camera_position[1] = 0.;
    camera_position[2] = z_min + 2. * z_dimension;

    camera_direction[0] = 0.;
    camera_direction[1] = 0.;
    camera_direction[2] = -1.;

    camera_horizontal[0] = 1.;
    camera_horizontal[1] = 0.;
    camera_horizontal[2] = 0.;

    camera_focus = .5 * z_dimension;

    Point<3> camera_position_temp;
    Point<3> camera_direction_temp;
    Point<3> camera_horizontal_temp;

    const float angle_factor = 3.14159265f / 180.f;

    // (I) rotate the camera to the chosen polar angle
    camera_position_temp[1] =
      std::cos(angle_factor * flags.polar_angle) * camera_position[1] -
      std::sin(angle_factor * flags.polar_angle) * camera_position[2];
    camera_position_temp[2] =
      std::sin(angle_factor * flags.polar_angle) * camera_position[1] +
      std::cos(angle_factor * flags.polar_angle) * camera_position[2];

    camera_direction_temp[1] =
      std::cos(angle_factor * flags.polar_angle) * camera_direction[1] -
      std::sin(angle_factor * flags.polar_angle) * camera_direction[2];
    camera_direction_temp[2] =
      std::sin(angle_factor * flags.polar_angle) * camera_direction[1] +
      std::cos(angle_factor * flags.polar_angle) * camera_direction[2];

    camera_horizontal_temp[1] =
      std::cos(angle_factor * flags.polar_angle) * camera_horizontal[1] -
      std::sin(angle_factor * flags.polar_angle) * camera_horizontal[2];
    camera_horizontal_temp[2] =
      std::sin(angle_factor * flags.polar_angle) * camera_horizontal[1] +
      std::cos(angle_factor * flags.polar_angle) * camera_horizontal[2];

    camera_position[1] = camera_position_temp[1];
    camera_position[2] = camera_position_temp[2];

    camera_direction[1] = camera_direction_temp[1];
    camera_direction[2] = camera_direction_temp[2];

    camera_horizontal[1] = camera_horizontal_temp[1];
    camera_horizontal[2] = camera_horizontal_temp[2];

    // (II) rotate the camera to the chosen azimuth angle
    camera_position_temp[0] =
      std::cos(angle_factor * flags.azimuth_angle) * camera_position[0] -
      std::sin(angle_factor * flags.azimuth_angle) * camera_position[1];
    camera_position_temp[1] =
      std::sin(angle_factor * flags.azimuth_angle) * camera_position[0] +
      std::cos(angle_factor * flags.azimuth_angle) * camera_position[1];

    camera_direction_temp[0] =
      std::cos(angle_factor * flags.azimuth_angle) * camera_direction[0] -
      std::sin(angle_factor * flags.azimuth_angle) * camera_direction[1];
    camera_direction_temp[1] =
      std::sin(angle_factor * flags.azimuth_angle) * camera_direction[0] +
      std::cos(angle_factor * flags.azimuth_angle) * camera_direction[1];

    camera_horizontal_temp[0] =
      std::cos(angle_factor * flags.azimuth_angle) * camera_horizontal[0] -
      std::sin(angle_factor * flags.azimuth_angle) * camera_horizontal[1];
    camera_horizontal_temp[1] =
      std::sin(angle_factor * flags.azimuth_angle) * camera_horizontal[0] +
      std::cos(angle_factor * flags.azimuth_angle) * camera_horizontal[1];

    camera_position[0] = camera_position_temp[0];
    camera_position[1] = camera_position_temp[1];

    camera_direction[0] = camera_direction_temp[0];
    camera_direction[1] = camera_direction_temp[1];

    camera_horizontal[0] = camera_horizontal_temp[0];
    camera_horizontal[1] = camera_horizontal_temp[1];

    // (III) translate the camera
    camera_position[0] = x_min + .5 * x_dimension;
    camera_position[1] = y_min + .5 * y_dimension;

    camera_position[0] += (z_min + 2. * z_dimension) *
                          std::sin(angle_factor * flags.polar_angle) *
                          std::sin(angle_factor * flags.azimuth_angle);
    camera_position[1] -= (z_min + 2. * z_dimension) *
                          std::sin(angle_factor * flags.polar_angle) *
                          std::cos(angle_factor * flags.azimuth_angle);


    // determine the bounding box on the projection plane
    double x_min_perspective, y_min_perspective;
    double x_max_perspective, y_max_perspective;
    double x_dimension_perspective, y_dimension_perspective;

    n_subdivisions = first_patch.n_subdivisions;
    n              = n_subdivisions + 1;

    Point<3> point;

    projected_point =
      get_equispaced_location(first_patch, {0, 0}, n_subdivisions);

    if (first_patch.data.n_rows() != 0)
      {
        AssertIndexRange(flags.height_vector, first_patch.data.n_rows());
      }

    point[0] = projected_point[0];
    point[1] = projected_point[1];
    point[2] = first_patch.data.n_rows() != 0 ?
                 first_patch.data(flags.height_vector, 0) :
                 0;

    projection_decomposition = svg_project_point(point,
                                                 camera_position,
                                                 camera_direction,
                                                 camera_horizontal,
                                                 camera_focus);

    x_min_perspective = projection_decomposition[0];
    x_max_perspective = projection_decomposition[0];
    y_min_perspective = projection_decomposition[1];
    y_max_perspective = projection_decomposition[1];

    // iterate over the patches
    for (const auto &patch : patches)
      {
        n_subdivisions = patch.n_subdivisions;
        for (unsigned int i2 = 0; i2 < n_subdivisions; ++i2)
          {
            for (unsigned int i1 = 0; i1 < n_subdivisions; ++i1)
              {
                const std::array<Point<spacedim>, 4> projected_vertices{
                  {get_equispaced_location(patch, {i1, i2}, n_subdivisions),
                   get_equispaced_location(patch, {i1 + 1, i2}, n_subdivisions),
                   get_equispaced_location(patch, {i1, i2 + 1}, n_subdivisions),
                   get_equispaced_location(patch,
                                           {i1 + 1, i2 + 1},
                                           n_subdivisions)}};

                Assert((flags.height_vector < patch.data.n_rows()) ||
                         patch.data.n_rows() == 0,
                       ExcIndexRange(flags.height_vector,
                                     0,
                                     patch.data.n_rows()));

                const std::array<Point<3>, 4> vertices = {
                  {Point<3>{projected_vertices[0][0],
                            projected_vertices[0][1],
                            patch.data.n_rows() != 0 ?
                              patch.data(0, i1 * d1 + i2 * d2) :
                              0},
                   Point<3>{projected_vertices[1][0],
                            projected_vertices[1][1],
                            patch.data.n_rows() != 0 ?
                              patch.data(0, (i1 + 1) * d1 + i2 * d2) :
                              0},
                   Point<3>{projected_vertices[2][0],
                            projected_vertices[2][1],
                            patch.data.n_rows() != 0 ?
                              patch.data(0, i1 * d1 + (i2 + 1) * d2) :
                              0},
                   Point<3>{projected_vertices[3][0],
                            projected_vertices[3][1],
                            patch.data.n_rows() != 0 ?
                              patch.data(0, (i1 + 1) * d1 + (i2 + 1) * d2) :
                              0}}};

                projection_decompositions = {
                  {svg_project_point(vertices[0],
                                     camera_position,
                                     camera_direction,
                                     camera_horizontal,
                                     camera_focus),
                   svg_project_point(vertices[1],
                                     camera_position,
                                     camera_direction,
                                     camera_horizontal,
                                     camera_focus),
                   svg_project_point(vertices[2],
                                     camera_position,
                                     camera_direction,
                                     camera_horizontal,
                                     camera_focus),
                   svg_project_point(vertices[3],
                                     camera_position,
                                     camera_direction,
                                     camera_horizontal,
                                     camera_focus)}};

                x_min_perspective =
                  std::min(x_min_perspective,
                           static_cast<double>(
                             projection_decompositions[0][0]));
                x_min_perspective =
                  std::min(x_min_perspective,
                           static_cast<double>(
                             projection_decompositions[1][0]));
                x_min_perspective =
                  std::min(x_min_perspective,
                           static_cast<double>(
                             projection_decompositions[2][0]));
                x_min_perspective =
                  std::min(x_min_perspective,
                           static_cast<double>(
                             projection_decompositions[3][0]));

                x_max_perspective =
                  std::max(x_max_perspective,
                           static_cast<double>(
                             projection_decompositions[0][0]));
                x_max_perspective =
                  std::max(x_max_perspective,
                           static_cast<double>(
                             projection_decompositions[1][0]));
                x_max_perspective =
                  std::max(x_max_perspective,
                           static_cast<double>(
                             projection_decompositions[2][0]));
                x_max_perspective =
                  std::max(x_max_perspective,
                           static_cast<double>(
                             projection_decompositions[3][0]));

                y_min_perspective =
                  std::min(y_min_perspective,
                           static_cast<double>(
                             projection_decompositions[0][1]));
                y_min_perspective =
                  std::min(y_min_perspective,
                           static_cast<double>(
                             projection_decompositions[1][1]));
                y_min_perspective =
                  std::min(y_min_perspective,
                           static_cast<double>(
                             projection_decompositions[2][1]));
                y_min_perspective =
                  std::min(y_min_perspective,
                           static_cast<double>(
                             projection_decompositions[3][1]));

                y_max_perspective =
                  std::max(y_max_perspective,
                           static_cast<double>(
                             projection_decompositions[0][1]));
                y_max_perspective =
                  std::max(y_max_perspective,
                           static_cast<double>(
                             projection_decompositions[1][1]));
                y_max_perspective =
                  std::max(y_max_perspective,
                           static_cast<double>(
                             projection_decompositions[2][1]));
                y_max_perspective =
                  std::max(y_max_perspective,
                           static_cast<double>(
                             projection_decompositions[3][1]));
              }
          }
      }

    x_dimension_perspective = x_max_perspective - x_min_perspective;
    y_dimension_perspective = y_max_perspective - y_min_perspective;

    std::multiset<SvgCell> cells;

    // iterate over the patches
    for (const auto &patch : patches)
      {
        n_subdivisions = patch.n_subdivisions;

        for (unsigned int i2 = 0; i2 < n_subdivisions; ++i2)
          {
            for (unsigned int i1 = 0; i1 < n_subdivisions; ++i1)
              {
                const std::array<Point<spacedim>, 4> projected_vertices = {
                  {get_equispaced_location(patch, {i1, i2}, n_subdivisions),
                   get_equispaced_location(patch, {i1 + 1, i2}, n_subdivisions),
                   get_equispaced_location(patch, {i1, i2 + 1}, n_subdivisions),
                   get_equispaced_location(patch,
                                           {i1 + 1, i2 + 1},
                                           n_subdivisions)}};

                Assert((flags.height_vector < patch.data.n_rows()) ||
                         patch.data.n_rows() == 0,
                       ExcIndexRange(flags.height_vector,
                                     0,
                                     patch.data.n_rows()));

                SvgCell cell;

                cell.vertices[0][0] = projected_vertices[0][0];
                cell.vertices[0][1] = projected_vertices[0][1];
                cell.vertices[0][2] = patch.data.n_rows() != 0 ?
                                        patch.data(0, i1 * d1 + i2 * d2) :
                                        0;

                cell.vertices[1][0] = projected_vertices[1][0];
                cell.vertices[1][1] = projected_vertices[1][1];
                cell.vertices[1][2] = patch.data.n_rows() != 0 ?
                                        patch.data(0, (i1 + 1) * d1 + i2 * d2) :
                                        0;

                cell.vertices[2][0] = projected_vertices[2][0];
                cell.vertices[2][1] = projected_vertices[2][1];
                cell.vertices[2][2] = patch.data.n_rows() != 0 ?
                                        patch.data(0, i1 * d1 + (i2 + 1) * d2) :
                                        0;

                cell.vertices[3][0] = projected_vertices[3][0];
                cell.vertices[3][1] = projected_vertices[3][1];
                cell.vertices[3][2] =
                  patch.data.n_rows() != 0 ?
                    patch.data(0, (i1 + 1) * d1 + (i2 + 1) * d2) :
                    0;

                cell.projected_vertices[0] =
                  svg_project_point(cell.vertices[0],
                                    camera_position,
                                    camera_direction,
                                    camera_horizontal,
                                    camera_focus);
                cell.projected_vertices[1] =
                  svg_project_point(cell.vertices[1],
                                    camera_position,
                                    camera_direction,
                                    camera_horizontal,
                                    camera_focus);
                cell.projected_vertices[2] =
                  svg_project_point(cell.vertices[2],
                                    camera_position,
                                    camera_direction,
                                    camera_horizontal,
                                    camera_focus);
                cell.projected_vertices[3] =
                  svg_project_point(cell.vertices[3],
                                    camera_position,
                                    camera_direction,
                                    camera_horizontal,
                                    camera_focus);

                cell.center = .25 * (cell.vertices[0] + cell.vertices[1] +
                                     cell.vertices[2] + cell.vertices[3]);
                cell.projected_center = svg_project_point(cell.center,
                                                          camera_position,
                                                          camera_direction,
                                                          camera_horizontal,
                                                          camera_focus);

                cell.depth = cell.center.distance(camera_position);

                cells.insert(cell);
              }
          }
      }


    // write the svg file
    if (width == 0)
      width = static_cast<unsigned int>(
        .5 + height * (x_dimension_perspective / y_dimension_perspective));
    unsigned int additional_width = 0;

    if (flags.draw_colorbar)
      additional_width = static_cast<unsigned int>(
        .5 + height * .3); // additional width for colorbar

    // basic svg header and background rectangle
    out << "<svg width=\"" << width + additional_width << "\" height=\""
        << height << "\" xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">"
        << '\n'
        << " <rect width=\"" << width + additional_width << "\" height=\""
        << height << "\" style=\"fill:white\"/>" << '\n'
        << '\n';

    unsigned int triangle_counter = 0;

    // write the cells in the correct order
    for (const auto &cell : cells)
      {
        Point<3> points3d_triangle[3];

        for (unsigned int triangle_index = 0; triangle_index < 4;
             triangle_index++)
          {
            switch (triangle_index)
              {
                case 0:
                  points3d_triangle[0] = cell.vertices[0],
                  points3d_triangle[1] = cell.vertices[1],
                  points3d_triangle[2] = cell.center;
                  break;
                case 1:
                  points3d_triangle[0] = cell.vertices[1],
                  points3d_triangle[1] = cell.vertices[3],
                  points3d_triangle[2] = cell.center;
                  break;
                case 2:
                  points3d_triangle[0] = cell.vertices[3],
                  points3d_triangle[1] = cell.vertices[2],
                  points3d_triangle[2] = cell.center;
                  break;
                case 3:
                  points3d_triangle[0] = cell.vertices[2],
                  points3d_triangle[1] = cell.vertices[0],
                  points3d_triangle[2] = cell.center;
                  break;
                default:
                  break;
              }

            Point<6> gradient_param =
              svg_get_gradient_parameters(points3d_triangle);

            double start_h =
              .667 - ((gradient_param[4] - z_min) / z_dimension) * .667;
            double stop_h =
              .667 - ((gradient_param[5] - z_min) / z_dimension) * .667;

            unsigned int start_r = 0;
            unsigned int start_g = 0;
            unsigned int start_b = 0;

            unsigned int stop_r = 0;
            unsigned int stop_g = 0;
            unsigned int stop_b = 0;

            unsigned int start_i = static_cast<unsigned int>(start_h * 6.);
            unsigned int stop_i  = static_cast<unsigned int>(stop_h * 6.);

            double start_f = start_h * 6. - start_i;
            double start_q = 1. - start_f;

            double stop_f = stop_h * 6. - stop_i;
            double stop_q = 1. - stop_f;

            switch (start_i % 6)
              {
                case 0:
                  start_r = 255,
                  start_g = static_cast<unsigned int>(.5 + 255. * start_f);
                  break;
                case 1:
                  start_r = static_cast<unsigned int>(.5 + 255. * start_q),
                  start_g = 255;
                  break;
                case 2:
                  start_g = 255,
                  start_b = static_cast<unsigned int>(.5 + 255. * start_f);
                  break;
                case 3:
                  start_g = static_cast<unsigned int>(.5 + 255. * start_q),
                  start_b = 255;
                  break;
                case 4:
                  start_r = static_cast<unsigned int>(.5 + 255. * start_f),
                  start_b = 255;
                  break;
                case 5:
                  start_r = 255,
                  start_b = static_cast<unsigned int>(.5 + 255. * start_q);
                  break;
                default:
                  break;
              }

            switch (stop_i % 6)
              {
                case 0:
                  stop_r = 255,
                  stop_g = static_cast<unsigned int>(.5 + 255. * stop_f);
                  break;
                case 1:
                  stop_r = static_cast<unsigned int>(.5 + 255. * stop_q),
                  stop_g = 255;
                  break;
                case 2:
                  stop_g = 255,
                  stop_b = static_cast<unsigned int>(.5 + 255. * stop_f);
                  break;
                case 3:
                  stop_g = static_cast<unsigned int>(.5 + 255. * stop_q),
                  stop_b = 255;
                  break;
                case 4:
                  stop_r = static_cast<unsigned int>(.5 + 255. * stop_f),
                  stop_b = 255;
                  break;
                case 5:
                  stop_r = 255,
                  stop_b = static_cast<unsigned int>(.5 + 255. * stop_q);
                  break;
                default:
                  break;
              }

            Point<3> gradient_start_point_3d, gradient_stop_point_3d;

            gradient_start_point_3d[0] = gradient_param[0];
            gradient_start_point_3d[1] = gradient_param[1];
            gradient_start_point_3d[2] = gradient_param[4];

            gradient_stop_point_3d[0] = gradient_param[2];
            gradient_stop_point_3d[1] = gradient_param[3];
            gradient_stop_point_3d[2] = gradient_param[5];

            Point<2> gradient_start_point =
              svg_project_point(gradient_start_point_3d,
                                camera_position,
                                camera_direction,
                                camera_horizontal,
                                camera_focus);
            Point<2> gradient_stop_point =
              svg_project_point(gradient_stop_point_3d,
                                camera_position,
                                camera_direction,
                                camera_horizontal,
                                camera_focus);

            // define linear gradient
            out << "  <linearGradient id=\"" << triangle_counter
                << "\" gradientUnits=\"userSpaceOnUse\" "
                << "x1=\""
                << static_cast<unsigned int>(
                     .5 +
                     ((gradient_start_point[0] - x_min_perspective) /
                      x_dimension_perspective) *
                       (width - (width / 100.) * 2. * margin_in_percent) +
                     ((width / 100.) * margin_in_percent))
                << "\" "
                << "y1=\""
                << static_cast<unsigned int>(
                     .5 + height - (height / 100.) * margin_in_percent -
                     ((gradient_start_point[1] - y_min_perspective) /
                      y_dimension_perspective) *
                       (height - (height / 100.) * 2. * margin_in_percent))
                << "\" "
                << "x2=\""
                << static_cast<unsigned int>(
                     .5 +
                     ((gradient_stop_point[0] - x_min_perspective) /
                      x_dimension_perspective) *
                       (width - (width / 100.) * 2. * margin_in_percent) +
                     ((width / 100.) * margin_in_percent))
                << "\" "
                << "y2=\""
                << static_cast<unsigned int>(
                     .5 + height - (height / 100.) * margin_in_percent -
                     ((gradient_stop_point[1] - y_min_perspective) /
                      y_dimension_perspective) *
                       (height - (height / 100.) * 2. * margin_in_percent))
                << "\""
                << ">" << '\n'
                << "   <stop offset=\"0\" style=\"stop-color:rgb(" << start_r
                << "," << start_g << "," << start_b << ")\"/>" << '\n'
                << "   <stop offset=\"1\" style=\"stop-color:rgb(" << stop_r
                << "," << stop_g << "," << stop_b << ")\"/>" << '\n'
                << "  </linearGradient>" << '\n';

            // draw current triangle
            double x1 = 0, y1 = 0, x2 = 0, y2 = 0;
            double x3 = cell.projected_center[0];
            double y3 = cell.projected_center[1];

            switch (triangle_index)
              {
                case 0:
                  x1 = cell.projected_vertices[0][0],
                  y1 = cell.projected_vertices[0][1],
                  x2 = cell.projected_vertices[1][0],
                  y2 = cell.projected_vertices[1][1];
                  break;
                case 1:
                  x1 = cell.projected_vertices[1][0],
                  y1 = cell.projected_vertices[1][1],
                  x2 = cell.projected_vertices[3][0],
                  y2 = cell.projected_vertices[3][1];
                  break;
                case 2:
                  x1 = cell.projected_vertices[3][0],
                  y1 = cell.projected_vertices[3][1],
                  x2 = cell.projected_vertices[2][0],
                  y2 = cell.projected_vertices[2][1];
                  break;
                case 3:
                  x1 = cell.projected_vertices[2][0],
                  y1 = cell.projected_vertices[2][1],
                  x2 = cell.projected_vertices[0][0],
                  y2 = cell.projected_vertices[0][1];
                  break;
                default:
                  break;
              }

            out << "  <path d=\"M "
                << static_cast<unsigned int>(
                     .5 +
                     ((x1 - x_min_perspective) / x_dimension_perspective) *
                       (width - (width / 100.) * 2. * margin_in_percent) +
                     ((width / 100.) * margin_in_percent))
                << ' '
                << static_cast<unsigned int>(
                     .5 + height - (height / 100.) * margin_in_percent -
                     ((y1 - y_min_perspective) / y_dimension_perspective) *
                       (height - (height / 100.) * 2. * margin_in_percent))
                << " L "
                << static_cast<unsigned int>(
                     .5 +
                     ((x2 - x_min_perspective) / x_dimension_perspective) *
                       (width - (width / 100.) * 2. * margin_in_percent) +
                     ((width / 100.) * margin_in_percent))
                << ' '
                << static_cast<unsigned int>(
                     .5 + height - (height / 100.) * margin_in_percent -
                     ((y2 - y_min_perspective) / y_dimension_perspective) *
                       (height - (height / 100.) * 2. * margin_in_percent))
                << " L "
                << static_cast<unsigned int>(
                     .5 +
                     ((x3 - x_min_perspective) / x_dimension_perspective) *
                       (width - (width / 100.) * 2. * margin_in_percent) +
                     ((width / 100.) * margin_in_percent))
                << ' '
                << static_cast<unsigned int>(
                     .5 + height - (height / 100.) * margin_in_percent -
                     ((y3 - y_min_perspective) / y_dimension_perspective) *
                       (height - (height / 100.) * 2. * margin_in_percent))
                << " L "
                << static_cast<unsigned int>(
                     .5 +
                     ((x1 - x_min_perspective) / x_dimension_perspective) *
                       (width - (width / 100.) * 2. * margin_in_percent) +
                     ((width / 100.) * margin_in_percent))
                << ' '
                << static_cast<unsigned int>(
                     .5 + height - (height / 100.) * margin_in_percent -
                     ((y1 - y_min_perspective) / y_dimension_perspective) *
                       (height - (height / 100.) * 2. * margin_in_percent))
                << "\" style=\"stroke:black; fill:url(#" << triangle_counter
                << "); stroke-width:" << flags.line_thickness << "\"/>" << '\n';

            ++triangle_counter;
          }
      }


    // draw the colorbar
    if (flags.draw_colorbar)
      {
        out << '\n' << " <!-- colorbar -->" << '\n';

        unsigned int element_height = static_cast<unsigned int>(
          ((height / 100.) * (71. - 2. * margin_in_percent)) / 4);
        unsigned int element_width =
          static_cast<unsigned int>(.5 + (height / 100.) * 2.5);

        additional_width = 0;
        if (!flags.margin)
          additional_width =
            static_cast<unsigned int>(.5 + (height / 100.) * 2.5);

        for (unsigned int index = 0; index < 4; ++index)
          {
            double start_h = .667 - ((index + 1) / 4.) * .667;
            double stop_h  = .667 - (index / 4.) * .667;

            unsigned int start_r = 0;
            unsigned int start_g = 0;
            unsigned int start_b = 0;

            unsigned int stop_r = 0;
            unsigned int stop_g = 0;
            unsigned int stop_b = 0;

            unsigned int start_i = static_cast<unsigned int>(start_h * 6.);
            unsigned int stop_i  = static_cast<unsigned int>(stop_h * 6.);

            double start_f = start_h * 6. - start_i;
            double start_q = 1. - start_f;

            double stop_f = stop_h * 6. - stop_i;
            double stop_q = 1. - stop_f;

            switch (start_i % 6)
              {
                case 0:
                  start_r = 255,
                  start_g = static_cast<unsigned int>(.5 + 255. * start_f);
                  break;
                case 1:
                  start_r = static_cast<unsigned int>(.5 + 255. * start_q),
                  start_g = 255;
                  break;
                case 2:
                  start_g = 255,
                  start_b = static_cast<unsigned int>(.5 + 255. * start_f);
                  break;
                case 3:
                  start_g = static_cast<unsigned int>(.5 + 255. * start_q),
                  start_b = 255;
                  break;
                case 4:
                  start_r = static_cast<unsigned int>(.5 + 255. * start_f),
                  start_b = 255;
                  break;
                case 5:
                  start_r = 255,
                  start_b = static_cast<unsigned int>(.5 + 255. * start_q);
                  break;
                default:
                  break;
              }

            switch (stop_i % 6)
              {
                case 0:
                  stop_r = 255,
                  stop_g = static_cast<unsigned int>(.5 + 255. * stop_f);
                  break;
                case 1:
                  stop_r = static_cast<unsigned int>(.5 + 255. * stop_q),
                  stop_g = 255;
                  break;
                case 2:
                  stop_g = 255,
                  stop_b = static_cast<unsigned int>(.5 + 255. * stop_f);
                  break;
                case 3:
                  stop_g = static_cast<unsigned int>(.5 + 255. * stop_q),
                  stop_b = 255;
                  break;
                case 4:
                  stop_r = static_cast<unsigned int>(.5 + 255. * stop_f),
                  stop_b = 255;
                  break;
                case 5:
                  stop_r = 255,
                  stop_b = static_cast<unsigned int>(.5 + 255. * stop_q);
                  break;
                default:
                  break;
              }

            // define gradient
            out << "  <linearGradient id=\"colorbar_" << index
                << "\" gradientUnits=\"userSpaceOnUse\" "
                << "x1=\"" << width + additional_width << "\" "
                << "y1=\""
                << static_cast<unsigned int>(.5 + (height / 100.) *
                                                    (margin_in_percent + 29)) +
                     (3 - index) * element_height
                << "\" "
                << "x2=\"" << width + additional_width << "\" "
                << "y2=\""
                << static_cast<unsigned int>(.5 + (height / 100.) *
                                                    (margin_in_percent + 29)) +
                     (4 - index) * element_height
                << "\""
                << ">" << '\n'
                << "   <stop offset=\"0\" style=\"stop-color:rgb(" << start_r
                << "," << start_g << "," << start_b << ")\"/>" << '\n'
                << "   <stop offset=\"1\" style=\"stop-color:rgb(" << stop_r
                << "," << stop_g << "," << stop_b << ")\"/>" << '\n'
                << "  </linearGradient>" << '\n';

            // draw box corresponding to the gradient above
            out
              << "  <rect"
              << " x=\"" << width + additional_width << "\" y=\""
              << static_cast<unsigned int>(.5 + (height / 100.) *
                                                  (margin_in_percent + 29)) +
                   (3 - index) * element_height
              << "\" width=\"" << element_width << "\" height=\""
              << element_height
              << "\" style=\"stroke:black; stroke-width:2; fill:url(#colorbar_"
              << index << ")\"/>" << '\n';
          }

        for (unsigned int index = 0; index < 5; ++index)
          {
            out
              << "  <text x=\""
              << width + additional_width +
                   static_cast<unsigned int>(1.5 * element_width)
              << "\" y=\""
              << static_cast<unsigned int>(
                   .5 + (height / 100.) * (margin_in_percent + 29) +
                   (4. - index) * element_height + 30.)
              << "\""
              << " style=\"text-anchor:start; font-size:80; font-family:Helvetica";

            if (index == 0 || index == 4)
              out << "; font-weight:bold";

            out << "\">"
                << static_cast<float>(
                     (static_cast<int>((z_min + index * (z_dimension / 4.)) *
                                       10000)) /
                     10000.);

            if (index == 4)
              out << " max";
            if (index == 0)
              out << " min";

            out << "</text>" << '\n';
          }
      }

    // finalize the svg file
    out << '\n' << "</svg>";
    out.flush();
  }



  template <int dim, int spacedim>
  void
  write_deal_II_intermediate(
    const std::vector<Patch<dim, spacedim>> &patches,
    const std::vector<std::string>          &data_names,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>
      &nonscalar_data_ranges,
    const Deal_II_IntermediateFlags & /*flags*/,
    std::ostream &out)
  {
    AssertThrow(out.fail() == false, ExcIO());

    // first write tokens indicating the template parameters. we need this in
    // here because we may want to read in data again even if we don't know in
    // advance the template parameters:
    out << dim << ' ' << spacedim << '\n';

    // then write a header
    out << "[deal.II intermediate format graphics data]" << '\n'
        << "[written by " << DEAL_II_PACKAGE_NAME << " "
        << DEAL_II_PACKAGE_VERSION << "]" << '\n'
        << "[Version: " << Deal_II_IntermediateFlags::format_version << "]"
        << '\n';

    out << data_names.size() << '\n';
    for (const auto &data_name : data_names)
      out << data_name << '\n';

    out << patches.size() << '\n';
    for (unsigned int i = 0; i < patches.size(); ++i)
      out << patches[i] << '\n';

    out << nonscalar_data_ranges.size() << '\n';
    for (const auto &nonscalar_data_range : nonscalar_data_ranges)
      out << std::get<0>(nonscalar_data_range) << ' '
          << std::get<1>(nonscalar_data_range) << '\n'
          << std::get<2>(nonscalar_data_range) << '\n';

    out << '\n';
    // make sure everything now gets to disk
    out.flush();
  }


  template <int dim, int spacedim>
  void
  write_deal_II_intermediate_in_parallel(
    const std::vector<Patch<dim, spacedim>> &patches,
    const std::vector<std::string>          &data_names,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>
                                    &nonscalar_data_ranges,
    const Deal_II_IntermediateFlags &flags,
    const std::string               &filename,
    const MPI_Comm                   comm,
    const CompressionLevel           compression)
  {
#ifndef DEAL_II_WITH_MPI
    (void)patches;
    (void)data_names;
    (void)nonscalar_data_ranges;
    (void)flags;
    (void)filename;
    (void)comm;
    (void)compression;

    AssertThrow(false,
                ExcMessage("This functionality requires MPI to be enabled."));

#else

    // We write a simple format based on the text format of
    // write_deal_II_intermediate() on each MPI rank. The text format
    // is quite verbose and we should probably change this to a more
    // efficient binary representation at some point. The file layout
    // is as follows:
    //
    // 1. A binary header with layout struct
    //    ParallelIntermediateHeaderType.
    // 2. A list of uint64_t with one value per rank denoting the
    //    compressed size of the chunks of the next step.
    // 3. The (potentially compressed) chunks as generated by
    //    write_deal_II_intermediate() on each MPI rank.

    // First generate my data by writing (optionally compressed) data into
    // my_buffer:
    std::vector<char> my_buffer;
    {
      boost::iostreams::filtering_ostream f;

      AssertThrow(compression != CompressionLevel::plain_text,
                  ExcNotImplemented());

      if (compression != CompressionLevel::no_compression)
#  ifdef DEAL_II_WITH_ZLIB
        f.push(boost::iostreams::zlib_compressor(
          get_boost_zlib_compression_level(compression)));
#  else
        AssertThrow(
          false,
          ExcMessage(
            "Compression requires deal.II to be configured with ZLIB support."));
#  endif

      boost::iostreams::back_insert_device<std::vector<char>> inserter(
        my_buffer);
      f.push(inserter);

      write_deal_II_intermediate<dim, spacedim>(
        patches, data_names, nonscalar_data_ranges, flags, f);
    }
    const std::uint64_t my_size = my_buffer.size();

    const unsigned int  my_rank   = Utilities::MPI::this_mpi_process(comm);
    const std::uint64_t n_ranks   = Utilities::MPI::n_mpi_processes(comm);
    const std::uint64_t n_patches = Utilities::MPI::sum(patches.size(), comm);

    const ParallelIntermediateHeader header{
      0x00dea111,
      Deal_II_IntermediateFlags::format_version,
      static_cast<std::uint64_t>(compression),
      dim,
      spacedim,
      n_ranks,
      n_patches};

    // Rank 0 also collects and writes the size of the data from each
    // rank in bytes. The static_cast for the destination buffer looks
    // useless, but without it clang-tidy will complain about a wrong
    // MPI type.
    std::vector<std::uint64_t> chunk_sizes(n_ranks);
    int                        ierr = MPI_Gather(&my_size,
                          1,
                          Utilities::MPI::mpi_type_id_for_type<std::uint64_t>,
                          static_cast<std::uint64_t *>(chunk_sizes.data()),
                          1,
                          Utilities::MPI::mpi_type_id_for_type<std::uint64_t>,
                          0,
                          comm);
    AssertThrowMPI(ierr);

    MPI_Info info;
    ierr = MPI_Info_create(&info);
    AssertThrowMPI(ierr);
    MPI_File fh;
    ierr = MPI_File_open(
      comm, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &fh);
    AssertThrow(ierr == MPI_SUCCESS, ExcFileNotOpen(filename));
    ierr = MPI_Info_free(&info);
    AssertThrowMPI(ierr);

    // Delete the file contents:
    ierr = MPI_File_set_size(fh, 0);
    AssertThrowMPI(ierr);
    // This barrier is necessary, because otherwise others might already write
    // while one core is still setting the size to zero.
    ierr = MPI_Barrier(comm);
    AssertThrowMPI(ierr);

    // Write the two parts of the header on rank 0:
    if (my_rank == 0)
      {
        ierr = Utilities::MPI::LargeCount::File_write_at_c(
          fh, 0, &header, sizeof(header), MPI_CHAR, MPI_STATUS_IGNORE);
        AssertThrowMPI(ierr);

        ierr = Utilities::MPI::LargeCount::File_write_at_c(
          fh,
          /* offset = */ sizeof(header),
          chunk_sizes.data(),
          chunk_sizes.size(),
          Utilities::MPI::mpi_type_id_for_type<std::uint64_t>,
          MPI_STATUS_IGNORE);
        AssertThrowMPI(ierr);
      }

    // Write the main part on each rank:
    {
      std::uint64_t prefix_sum = 0;
      ierr                     = MPI_Exscan(&my_size,
                        &prefix_sum,
                        1,
                        Utilities::MPI::mpi_type_id_for_type<std::uint64_t>,
                        MPI_SUM,
                        comm);
      AssertThrowMPI(ierr);

      // Locate specific offset for each processor.
      const MPI_Offset offset = static_cast<MPI_Offset>(sizeof(header)) +
                                n_ranks * sizeof(std::uint64_t) + prefix_sum;

      ierr = Utilities::MPI::LargeCount::File_write_at_all_c(
        fh, offset, my_buffer.data(), my_size, MPI_CHAR, MPI_STATUS_IGNORE);
      AssertThrowMPI(ierr);
    }

    // Make sure we sync to disk. As written in the standard,
    // MPI_File_close() actually already implies a sync but there seems
    // to be a bug on at least one configuration (running with multiple
    // nodes using OpenMPI 4.1) that requires it. Without this call, the
    // footer is sometimes missing.
    ierr = MPI_File_sync(fh);
    AssertThrowMPI(ierr);

    ierr = MPI_File_close(&fh);
    AssertThrowMPI(ierr);
#endif
  }



  std::pair<unsigned int, unsigned int>
  determine_intermediate_format_dimensions(std::istream &input)
  {
    AssertThrow(input.fail() == false, ExcIO());

    unsigned int dim, spacedim;
    input >> dim >> spacedim;

    return std::make_pair(dim, spacedim);
  }
} // namespace DataOutBase



/* --------------------------- class DataOutInterface ---------------------- */


template <int dim, int spacedim>
DataOutInterface<dim, spacedim>::DataOutInterface()
  : default_subdivisions(1)
  , default_fmt(DataOutBase::default_format)
{}



template <int dim, int spacedim>
void
DataOutInterface<dim, spacedim>::write_dx(std::ostream &out) const
{
  DataOutBase::write_dx(get_patches(),
                        get_dataset_names(),
                        get_nonscalar_data_ranges(),
                        dx_flags,
                        out);
}



template <int dim, int spacedim>
void
DataOutInterface<dim, spacedim>::write_ucd(std::ostream &out) const
{
  DataOutBase::write_ucd(get_patches(),
                         get_dataset_names(),
                         get_nonscalar_data_ranges(),
                         ucd_flags,
                         out);
}



template <int dim, int spacedim>
void
DataOutInterface<dim, spacedim>::write_gnuplot(std::ostream &out) const
{
  DataOutBase::write_gnuplot(get_patches(),
                             get_dataset_names(),
                             get_nonscalar_data_ranges(),
                             gnuplot_flags,
                             out);
}



template <int dim, int spacedim>
void
DataOutInterface<dim, spacedim>::write_povray(std::ostream &out) const
{
  DataOutBase::write_povray(get_patches(),
                            get_dataset_names(),
                            get_nonscalar_data_ranges(),
                            povray_flags,
                            out);
}



template <int dim, int spacedim>
void
DataOutInterface<dim, spacedim>::write_eps(std::ostream &out) const
{
  DataOutBase::write_eps(get_patches(),
                         get_dataset_names(),
                         get_nonscalar_data_ranges(),
                         eps_flags,
                         out);
}



template <int dim, int spacedim>
void
DataOutInterface<dim, spacedim>::write_gmv(std::ostream &out) const
{
  DataOutBase::write_gmv(get_patches(),
                         get_dataset_names(),
                         get_nonscalar_data_ranges(),
                         gmv_flags,
                         out);
}



template <int dim, int spacedim>
void
DataOutInterface<dim, spacedim>::write_tecplot(std::ostream &out) const
{
  DataOutBase::write_tecplot(get_patches(),
                             get_dataset_names(),
                             get_nonscalar_data_ranges(),
                             tecplot_flags,
                             out);
}



template <int dim, int spacedim>
void
DataOutInterface<dim, spacedim>::write_vtk(std::ostream &out) const
{
  DataOutBase::write_vtk(get_patches(),
                         get_dataset_names(),
                         get_nonscalar_data_ranges(),
                         vtk_flags,
                         out);
}

template <int dim, int spacedim>
void
DataOutInterface<dim, spacedim>::write_vtu(std::ostream &out) const
{
  DataOutBase::write_vtu(get_patches(),
                         get_dataset_names(),
                         get_nonscalar_data_ranges(),
                         vtk_flags,
                         out);
}

template <int dim, int spacedim>
void
DataOutInterface<dim, spacedim>::write_svg(std::ostream &out) const
{
  DataOutBase::write_svg(get_patches(),
                         get_dataset_names(),
                         get_nonscalar_data_ranges(),
                         svg_flags,
                         out);
}


template <int dim, int spacedim>
void
DataOutInterface<dim, spacedim>::write_vtu_in_parallel(
  const std::string &filename,
  const MPI_Comm     comm) const
{
#ifndef DEAL_II_WITH_MPI
  // without MPI fall back to the normal way to write a vtu file:
  (void)comm;

  std::ofstream f(filename);
  AssertThrow(f, ExcFileNotOpen(filename));
  write_vtu(f);
#else

  const unsigned int myrank  = Utilities::MPI::this_mpi_process(comm);
  const unsigned int n_ranks = Utilities::MPI::n_mpi_processes(comm);
  MPI_Info           info;
  int                ierr = MPI_Info_create(&info);
  AssertThrowMPI(ierr);
  MPI_File fh;
  ierr = MPI_File_open(
    comm, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &fh);
  AssertThrow(ierr == MPI_SUCCESS, ExcFileNotOpen(filename));

  ierr = MPI_File_set_size(fh, 0); // delete the file contents
  AssertThrowMPI(ierr);
  // this barrier is necessary, because otherwise others might already write
  // while one core is still setting the size to zero.
  ierr = MPI_Barrier(comm);
  AssertThrowMPI(ierr);
  ierr = MPI_Info_free(&info);
  AssertThrowMPI(ierr);

  // Define header size so we can broadcast later.
  unsigned int  header_size;
  std::uint64_t footer_offset;

  // write header
  if (myrank == 0)
    {
      std::stringstream ss;
      DataOutBase::write_vtu_header(ss, vtk_flags);
      header_size = ss.str().size();
      // Write the header on rank 0 at the start of a file, i.e., offset 0.
      ierr = Utilities::MPI::LargeCount::File_write_at_c(
        fh, 0, ss.str().c_str(), header_size, MPI_CHAR, MPI_STATUS_IGNORE);
      AssertThrowMPI(ierr);
    }

  ierr = MPI_Bcast(&header_size, 1, MPI_UNSIGNED, 0, comm);
  AssertThrowMPI(ierr);

  {
    const auto                   &patches      = get_patches();
    const types::global_dof_index my_n_patches = patches.size();
    const types::global_dof_index global_n_patches =
      Utilities::MPI::sum(my_n_patches, comm);

    // Do not write pieces with 0 cells as this will crash paraview if this is
    // the first piece written. But if nobody has any pieces to write (file is
    // empty), let processor 0 write their empty data, otherwise the vtk file is
    // invalid.
    std::stringstream ss;
    if (my_n_patches > 0 || (global_n_patches == 0 && myrank == 0))
      DataOutBase::write_vtu_main(patches,
                                  get_dataset_names(),
                                  get_nonscalar_data_ranges(),
                                  vtk_flags,
                                  ss);

    // Use prefix sum to find specific offset to write at.
    const std::uint64_t size_on_proc = ss.str().size();
    std::uint64_t       prefix_sum   = 0;
    ierr                             = MPI_Exscan(&size_on_proc,
                      &prefix_sum,
                      1,
                      Utilities::MPI::mpi_type_id_for_type<std::uint64_t>,
                      MPI_SUM,
                      comm);
    AssertThrowMPI(ierr);

    // Locate specific offset for each processor.
    const MPI_Offset offset = static_cast<MPI_Offset>(header_size) + prefix_sum;

    ierr = Utilities::MPI::LargeCount::File_write_at_all_c(fh,
                                                           offset,
                                                           ss.str().c_str(),
                                                           ss.str().size(),
                                                           MPI_CHAR,
                                                           MPI_STATUS_IGNORE);
    AssertThrowMPI(ierr);

    if (myrank == n_ranks - 1)
      {
        // Locating Footer with offset on last rank.
        footer_offset = size_on_proc + offset;

        std::stringstream ss;
        DataOutBase::write_vtu_footer(ss);
        const unsigned int footer_size = ss.str().size();

        // Writing footer:
        ierr = Utilities::MPI::LargeCount::File_write_at_c(fh,
                                                           footer_offset,
                                                           ss.str().c_str(),
                                                           footer_size,
                                                           MPI_CHAR,
                                                           MPI_STATUS_IGNORE);
        AssertThrowMPI(ierr);
      }
  }

  // Make sure we sync to disk. As written in the standard,
  // MPI_File_close() actually already implies a sync but there seems
  // to be a bug on at least one configuration (running with multiple
  // nodes using OpenMPI 4.1) that requires it. Without this call, the
  // footer is sometimes missing.
  ierr = MPI_File_sync(fh);
  AssertThrowMPI(ierr);

  ierr = MPI_File_close(&fh);
  AssertThrowMPI(ierr);
#endif
}



template <int dim, int spacedim>
void
DataOutInterface<dim, spacedim>::write_pvtu_record(
  std::ostream                   &out,
  const std::vector<std::string> &piece_names) const
{
  DataOutBase::write_pvtu_record(out,
                                 piece_names,
                                 get_dataset_names(),
                                 get_nonscalar_data_ranges(),
                                 vtk_flags);
}



template <int dim, int spacedim>
std::string
DataOutInterface<dim, spacedim>::write_vtu_with_pvtu_record(
  const std::string &directory,
  const std::string &filename_without_extension,
  const unsigned int counter,
  const MPI_Comm     mpi_communicator,
  const unsigned int n_digits_for_counter,
  const unsigned int n_groups) const
{
  const unsigned int rank = Utilities::MPI::this_mpi_process(mpi_communicator);
  const unsigned int n_ranks =
    Utilities::MPI::n_mpi_processes(mpi_communicator);
  const unsigned int n_files_written =
    (n_groups == 0 || n_groups > n_ranks) ? n_ranks : n_groups;

  Assert(n_files_written >= 1, ExcInternalError());
  // the "-1" is needed since we use C++ style counting starting with 0, so
  // writing 10 files means the filename runs from 0 to 9
  const unsigned int n_digits =
    Utilities::needed_digits(std::max(0, int(n_files_written) - 1));

  const unsigned int color = rank % n_files_written;
  const std::string  filename =
    directory + filename_without_extension + "_" +
    Utilities::int_to_string(counter, n_digits_for_counter) + "." +
    Utilities::int_to_string(color, n_digits) + ".vtu";

  if (n_groups == 0 || n_groups > n_ranks)
    {
      // every processor writes one file
      std::ofstream output(filename);
      AssertThrow(output, ExcFileNotOpen(filename));
      this->write_vtu(output);
    }
  else if (n_groups == 1)
    {
      // write only a single data file in parallel
      this->write_vtu_in_parallel(filename, mpi_communicator);
    }
  else
    {
#ifdef DEAL_II_WITH_MPI
      // write n_groups data files
      MPI_Comm comm_group;
      int ierr = MPI_Comm_split(mpi_communicator, color, rank, &comm_group);
      AssertThrowMPI(ierr);
      this->write_vtu_in_parallel(filename, comm_group);
      Utilities::MPI::free_communicator(comm_group);
#else
      AssertThrow(false, ExcMessage("Logical error. Should not arrive here."));
#endif
    }

  // write pvtu record
  const std::string pvtu_filename =
    filename_without_extension + "_" +
    Utilities::int_to_string(counter, n_digits_for_counter) + ".pvtu";

  if (rank == 0)
    {
      std::vector<std::string> filename_vector;
      for (unsigned int i = 0; i < n_files_written; ++i)
        {
          const std::string filename =
            filename_without_extension + "_" +
            Utilities::int_to_string(counter, n_digits_for_counter) + "." +
            Utilities::int_to_string(i, n_digits) + ".vtu";

          filename_vector.emplace_back(filename);
        }

      std::ofstream pvtu_output(directory + pvtu_filename);
      this->write_pvtu_record(pvtu_output, filename_vector);
    }

  return pvtu_filename;
}



template <int dim, int spacedim>
void
DataOutInterface<dim, spacedim>::write_deal_II_intermediate(
  std::ostream &out) const
{
  DataOutBase::write_deal_II_intermediate(get_patches(),
                                          get_dataset_names(),
                                          get_nonscalar_data_ranges(),
                                          deal_II_intermediate_flags,
                                          out);
}



template <int dim, int spacedim>
void
DataOutInterface<dim, spacedim>::write_deal_II_intermediate_in_parallel(
  const std::string                  &filename,
  const MPI_Comm                      comm,
  const DataOutBase::CompressionLevel compression) const
{
  DataOutBase::write_deal_II_intermediate_in_parallel(
    get_patches(),
    get_dataset_names(),
    get_nonscalar_data_ranges(),
    deal_II_intermediate_flags,
    filename,
    comm,
    compression);
}



template <int dim, int spacedim>
XDMFEntry
DataOutInterface<dim, spacedim>::create_xdmf_entry(
  const DataOutBase::DataOutFilter &data_filter,
  const std::string                &h5_filename,
  const double                      cur_time,
  const MPI_Comm                    comm) const
{
  return create_xdmf_entry(
    data_filter, h5_filename, h5_filename, cur_time, comm);
}



template <int dim, int spacedim>
XDMFEntry
DataOutInterface<dim, spacedim>::create_xdmf_entry(
  const DataOutBase::DataOutFilter &data_filter,
  const std::string                &h5_mesh_filename,
  const std::string                &h5_solution_filename,
  const double                      cur_time,
  const MPI_Comm                    comm) const
{
  AssertThrow(spacedim == 2 || spacedim == 3,
              ExcMessage("XDMF only supports 2 or 3 space dimensions."));

#ifndef DEAL_II_WITH_HDF5
  // throw an exception, but first make sure the compiler does not warn about
  // the now unused function arguments
  (void)data_filter;
  (void)h5_mesh_filename;
  (void)h5_solution_filename;
  (void)cur_time;
  (void)comm;
  AssertThrow(false, ExcMessage("XDMF support requires HDF5 to be turned on."));

  return {};

#else

  std::uint64_t local_node_cell_count[2], global_node_cell_count[2];

  local_node_cell_count[0] = data_filter.n_nodes();
  local_node_cell_count[1] = data_filter.n_cells();

  const int myrank = Utilities::MPI::this_mpi_process(comm);
  // And compute the global total
  int ierr = MPI_Allreduce(local_node_cell_count,
                           global_node_cell_count,
                           2,
                           Utilities::MPI::mpi_type_id_for_type<std::uint64_t>,
                           MPI_SUM,
                           comm);
  AssertThrowMPI(ierr);

  // The implementation is a bit complicated because we are supposed to return
  // the correct data on rank 0 and an empty object on all other ranks but all
  // information (for example the attributes) are only available on ranks that
  // have any cells.
  // We will identify the smallest rank that has data and then communicate
  // from this rank to rank 0 (if they are different ranks).

  const bool have_data = (data_filter.n_nodes() > 0);
  MPI_Comm   split_comm;
  {
    const int key   = myrank;
    const int color = (have_data ? 1 : 0);
    const int ierr  = MPI_Comm_split(comm, color, key, &split_comm);
    AssertThrowMPI(ierr);
  }

  const bool am_i_first_rank_with_data =
    have_data && (Utilities::MPI::this_mpi_process(split_comm) == 0);

  ierr = MPI_Comm_free(&split_comm);
  AssertThrowMPI(ierr);

  const int tag = 47381;

  // Output the XDMF file only on the root process of all ranks with data:
  if (am_i_first_rank_with_data)
    {
      const auto &patches = get_patches();
      Assert(patches.size() > 0, DataOutBase::ExcNoPatches());

      // We currently don't support writing mixed meshes:
      if constexpr (running_in_debug_mode())
        {
          for (const auto &patch : patches)
            Assert(patch.reference_cell == patches[0].reference_cell,
                   ExcNotImplemented());
        }

      XDMFEntry          entry(h5_mesh_filename,
                      h5_solution_filename,
                      cur_time,
                      global_node_cell_count[0],
                      global_node_cell_count[1],
                      dim,
                      spacedim,
                      patches[0].reference_cell);
      const unsigned int n_data_sets = data_filter.n_data_sets();

      // The vector names generated here must match those generated in
      // the HDF5 file
      for (unsigned int i = 0; i < n_data_sets; ++i)
        {
          entry.add_attribute(data_filter.get_data_set_name(i),
                              data_filter.get_data_set_dim(i));
        }

      if (myrank != 0)
        {
          // send to rank 0
          const std::vector<char> buffer = Utilities::pack(entry, false);
          ierr = MPI_Send(buffer.data(), buffer.size(), MPI_BYTE, 0, tag, comm);
          AssertThrowMPI(ierr);

          return {};
        }

      return entry;
    }

  if (myrank == 0 && !am_i_first_rank_with_data)
    {
      // receive the XDMF data on rank 0 if we don't have it...

      MPI_Status status;
      int        ierr = MPI_Probe(MPI_ANY_SOURCE, tag, comm, &status);
      AssertThrowMPI(ierr);

      int len;
      ierr = MPI_Get_count(&status, MPI_BYTE, &len);
      AssertThrowMPI(ierr);

      std::vector<char> buffer(len);
      ierr = MPI_Recv(buffer.data(),
                      len,
                      MPI_BYTE,
                      status.MPI_SOURCE,
                      tag,
                      comm,
                      MPI_STATUS_IGNORE);
      AssertThrowMPI(ierr);

      return Utilities::unpack<XDMFEntry>(buffer, false);
    }

  // default case for any other rank is to return an empty object
  return {};
#endif
}

template <int dim, int spacedim>
void
DataOutInterface<dim, spacedim>::write_xdmf_file(
  const std::vector<XDMFEntry> &entries,
  const std::string            &filename,
  const MPI_Comm                comm) const
{
#ifdef DEAL_II_WITH_MPI
  const int myrank = Utilities::MPI::this_mpi_process(comm);
#else
  (void)comm;
  const int myrank = 0;
#endif

  // Only rank 0 process writes the XDMF file
  if (myrank == 0)
    {
      std::ofstream xdmf_file(filename);

      xdmf_file << "<?xml version=\"1.0\" ?>\n";
      xdmf_file << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
      xdmf_file << "<Xdmf Version=\"2.0\">\n";
      xdmf_file << "  <Domain>\n";
      xdmf_file
        << "    <Grid Name=\"CellTime\" GridType=\"Collection\" CollectionType=\"Temporal\">\n";

      for (const auto &entry : entries)
        {
          xdmf_file << entry.get_xdmf_content(3);
        }

      xdmf_file << "    </Grid>\n";
      xdmf_file << "  </Domain>\n";
      xdmf_file << "</Xdmf>\n";

      xdmf_file.close();
    }
}



/*
 * Write the data in this DataOutInterface to a DataOutFilter object. Filtering
 * is performed based on the DataOutFilter flags.
 */
template <int dim, int spacedim>
void
DataOutInterface<dim, spacedim>::write_filtered_data(
  DataOutBase::DataOutFilter &filtered_data) const
{
  DataOutBase::write_filtered_data(get_patches(),
                                   get_dataset_names(),
                                   get_nonscalar_data_ranges(),
                                   filtered_data);
}


namespace
{
#ifdef DEAL_II_WITH_HDF5
  /**
   * Helper function to actually perform the HDF5 output.
   */
  template <int dim, int spacedim>
  void
  do_write_hdf5(const std::vector<DataOutBase::Patch<dim, spacedim>> &patches,
                const DataOutBase::DataOutFilter &data_filter,
                const DataOutBase::Hdf5Flags     &flags,
                const bool                        write_mesh_file,
                const std::string                &mesh_filename,
                const std::string                &solution_filename,
                const MPI_Comm                    comm)
  {
    hid_t h5_mesh_file_id = -1, h5_solution_file_id, file_plist_id, plist_id;
    hid_t node_dataspace, node_dataset, node_file_dataspace,
      node_memory_dataspace, node_dataset_id;
    hid_t cell_dataspace, cell_dataset, cell_file_dataspace,
      cell_memory_dataspace;
    hid_t pt_data_dataspace, pt_data_dataset, pt_data_file_dataspace,
      pt_data_memory_dataspace;
    herr_t              status;
    std::uint64_t       local_node_cell_count[2];
    hsize_t             count[2], offset[2], node_ds_dim[2], cell_ds_dim[2];
    std::vector<double> node_data_vec;
    std::vector<unsigned int> cell_data_vec;



    local_node_cell_count[0] = data_filter.n_nodes();
    local_node_cell_count[1] = data_filter.n_cells();

    // Create file access properties
    file_plist_id = H5Pcreate(H5P_FILE_ACCESS);
    AssertThrow(file_plist_id != -1, ExcIO());
    // If MPI is enabled *and* HDF5 is parallel, we can do parallel output
#  ifdef DEAL_II_WITH_MPI
#    ifdef H5_HAVE_PARALLEL
    // Set the access to use the specified MPI_Comm object
    status = H5Pset_fapl_mpio(file_plist_id, comm, MPI_INFO_NULL);
    AssertThrow(status >= 0, ExcIO());
#    endif
#  endif
    // if zlib support is disabled flags are unused
#  ifndef DEAL_II_WITH_ZLIB
    (void)flags;
#  endif

    // Compute the global total number of nodes/cells and determine the offset
    // of the data for this process

    std::uint64_t global_node_cell_count[2]   = {0, 0};
    std::uint64_t global_node_cell_offsets[2] = {0, 0};

#  ifdef DEAL_II_WITH_MPI
    int ierr =
      MPI_Allreduce(local_node_cell_count,
                    global_node_cell_count,
                    2,
                    Utilities::MPI::mpi_type_id_for_type<std::uint64_t>,
                    MPI_SUM,
                    comm);
    AssertThrowMPI(ierr);
    ierr = MPI_Exscan(local_node_cell_count,
                      global_node_cell_offsets,
                      2,
                      Utilities::MPI::mpi_type_id_for_type<std::uint64_t>,
                      MPI_SUM,
                      comm);
    AssertThrowMPI(ierr);
#  else
    global_node_cell_count[0]   = local_node_cell_count[0];
    global_node_cell_count[1]   = local_node_cell_count[1];
    global_node_cell_offsets[0] = global_node_cell_offsets[1] = 0;
#  endif

    // Create the property list for a collective write
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    AssertThrow(plist_id >= 0, ExcIO());
#  ifdef DEAL_II_WITH_MPI
#    ifdef H5_HAVE_PARALLEL
    status = H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    AssertThrow(status >= 0, ExcIO());
#    endif
#  endif

    if (write_mesh_file)
      {
        // Overwrite any existing files (change this to an option?)
        h5_mesh_file_id = H5Fcreate(mesh_filename.c_str(),
                                    H5F_ACC_TRUNC,
                                    H5P_DEFAULT,
                                    file_plist_id);
        AssertThrow(h5_mesh_file_id >= 0, ExcIO());

        // Create the dataspace for the nodes and cells. HDF5 only supports 2-
        // or 3-dimensional coordinates
        node_ds_dim[0] = global_node_cell_count[0];
        node_ds_dim[1] = (spacedim < 2) ? 2 : spacedim;
        node_dataspace = H5Screate_simple(2, node_ds_dim, nullptr);
        AssertThrow(node_dataspace >= 0, ExcIO());

        cell_ds_dim[0] = global_node_cell_count[1];
        cell_ds_dim[1] = patches[0].reference_cell.n_vertices();
        cell_dataspace = H5Screate_simple(2, cell_ds_dim, nullptr);
        AssertThrow(cell_dataspace >= 0, ExcIO());

        // Create the dataset for the nodes and cells
#  if H5Gcreate_vers == 1
        node_dataset = H5Dcreate(h5_mesh_file_id,
                                 "nodes",
                                 H5T_NATIVE_DOUBLE,
                                 node_dataspace,
                                 H5P_DEFAULT);
#  else
        node_dataset_id = H5Pcreate(H5P_DATASET_CREATE);
#    ifdef DEAL_II_WITH_ZLIB
        H5Pset_deflate(node_dataset_id,
                       get_zlib_compression_level(flags.compression_level));
        H5Pset_chunk(node_dataset_id, 2, node_ds_dim);
#    endif
        node_dataset = H5Dcreate(h5_mesh_file_id,
                                 "nodes",
                                 H5T_NATIVE_DOUBLE,
                                 node_dataspace,
                                 H5P_DEFAULT,
                                 node_dataset_id,
                                 H5P_DEFAULT);
        H5Pclose(node_dataset_id);
#  endif
        AssertThrow(node_dataset >= 0, ExcIO());
#  if H5Gcreate_vers == 1
        cell_dataset = H5Dcreate(h5_mesh_file_id,
                                 "cells",
                                 H5T_NATIVE_UINT,
                                 cell_dataspace,
                                 H5P_DEFAULT);
#  else
        node_dataset_id = H5Pcreate(H5P_DATASET_CREATE);
#    ifdef DEAL_II_WITH_ZLIB
        H5Pset_deflate(node_dataset_id,
                       get_zlib_compression_level(flags.compression_level));
        H5Pset_chunk(node_dataset_id, 2, cell_ds_dim);
#    endif
        cell_dataset = H5Dcreate(h5_mesh_file_id,
                                 "cells",
                                 H5T_NATIVE_UINT,
                                 cell_dataspace,
                                 H5P_DEFAULT,
                                 node_dataset_id,
                                 H5P_DEFAULT);
        H5Pclose(node_dataset_id);
#  endif
        AssertThrow(cell_dataset >= 0, ExcIO());

        // Close the node and cell dataspaces since we're done with them
        status = H5Sclose(node_dataspace);
        AssertThrow(status >= 0, ExcIO());
        status = H5Sclose(cell_dataspace);
        AssertThrow(status >= 0, ExcIO());

        // Create the data subset we'll use to read from memory. HDF5 only
        // supports 2- or 3-dimensional coordinates
        count[0] = local_node_cell_count[0];
        count[1] = (spacedim < 2) ? 2 : spacedim;

        offset[0] = global_node_cell_offsets[0];
        offset[1] = 0;

        node_memory_dataspace = H5Screate_simple(2, count, nullptr);
        AssertThrow(node_memory_dataspace >= 0, ExcIO());

        // Select the hyperslab in the file
        node_file_dataspace = H5Dget_space(node_dataset);
        AssertThrow(node_file_dataspace >= 0, ExcIO());
        status = H5Sselect_hyperslab(
          node_file_dataspace, H5S_SELECT_SET, offset, nullptr, count, nullptr);
        AssertThrow(status >= 0, ExcIO());

        // And repeat for cells
        count[0]              = local_node_cell_count[1];
        count[1]              = patches[0].reference_cell.n_vertices();
        offset[0]             = global_node_cell_offsets[1];
        offset[1]             = 0;
        cell_memory_dataspace = H5Screate_simple(2, count, nullptr);
        AssertThrow(cell_memory_dataspace >= 0, ExcIO());

        cell_file_dataspace = H5Dget_space(cell_dataset);
        AssertThrow(cell_file_dataspace >= 0, ExcIO());
        status = H5Sselect_hyperslab(
          cell_file_dataspace, H5S_SELECT_SET, offset, nullptr, count, nullptr);
        AssertThrow(status >= 0, ExcIO());

        // And finally, write the node data
        data_filter.fill_node_data(node_data_vec);
        status = H5Dwrite(node_dataset,
                          H5T_NATIVE_DOUBLE,
                          node_memory_dataspace,
                          node_file_dataspace,
                          plist_id,
                          node_data_vec.data());
        AssertThrow(status >= 0, ExcIO());
        node_data_vec.clear();

        // And the cell data
        data_filter.fill_cell_data(global_node_cell_offsets[0], cell_data_vec);
        status = H5Dwrite(cell_dataset,
                          H5T_NATIVE_UINT,
                          cell_memory_dataspace,
                          cell_file_dataspace,
                          plist_id,
                          cell_data_vec.data());
        AssertThrow(status >= 0, ExcIO());
        cell_data_vec.clear();

        // Close the file dataspaces
        status = H5Sclose(node_file_dataspace);
        AssertThrow(status >= 0, ExcIO());
        status = H5Sclose(cell_file_dataspace);
        AssertThrow(status >= 0, ExcIO());

        // Close the memory dataspaces
        status = H5Sclose(node_memory_dataspace);
        AssertThrow(status >= 0, ExcIO());
        status = H5Sclose(cell_memory_dataspace);
        AssertThrow(status >= 0, ExcIO());

        // Close the datasets
        status = H5Dclose(node_dataset);
        AssertThrow(status >= 0, ExcIO());
        status = H5Dclose(cell_dataset);
        AssertThrow(status >= 0, ExcIO());

        // If the filenames are different, we need to close the mesh file
        if (mesh_filename != solution_filename)
          {
            status = H5Fclose(h5_mesh_file_id);
            AssertThrow(status >= 0, ExcIO());
          }
      }

    // If the filenames are identical, continue with the same file
    if (mesh_filename == solution_filename && write_mesh_file)
      {
        h5_solution_file_id = h5_mesh_file_id;
      }
    else
      {
        // Otherwise we need to open a new file
        h5_solution_file_id = H5Fcreate(solution_filename.c_str(),
                                        H5F_ACC_TRUNC,
                                        H5P_DEFAULT,
                                        file_plist_id);
        AssertThrow(h5_solution_file_id >= 0, ExcIO());
      }

    // when writing, first write out all vector data, then handle the scalar
    // data sets that have been left over
    unsigned int i;
    std::string  vector_name;
    for (i = 0; i < data_filter.n_data_sets(); ++i)
      {
        // Allocate space for the point data
        // Must be either 1d or 3d
        const unsigned int pt_data_vector_dim = data_filter.get_data_set_dim(i);
        vector_name = data_filter.get_data_set_name(i);

        // Create the dataspace for the point data
        node_ds_dim[0]    = global_node_cell_count[0];
        node_ds_dim[1]    = pt_data_vector_dim;
        pt_data_dataspace = H5Screate_simple(2, node_ds_dim, nullptr);
        AssertThrow(pt_data_dataspace >= 0, ExcIO());

#  if H5Gcreate_vers == 1
        pt_data_dataset = H5Dcreate(h5_solution_file_id,
                                    vector_name.c_str(),
                                    H5T_NATIVE_DOUBLE,
                                    pt_data_dataspace,
                                    H5P_DEFAULT);
#  else
        node_dataset_id = H5Pcreate(H5P_DATASET_CREATE);
#    ifdef DEAL_II_WITH_ZLIB
        H5Pset_deflate(node_dataset_id,
                       get_zlib_compression_level(flags.compression_level));
        H5Pset_chunk(node_dataset_id, 2, node_ds_dim);
#    endif
        pt_data_dataset = H5Dcreate(h5_solution_file_id,
                                    vector_name.c_str(),
                                    H5T_NATIVE_DOUBLE,
                                    pt_data_dataspace,
                                    H5P_DEFAULT,
                                    node_dataset_id,
                                    H5P_DEFAULT);
        H5Pclose(node_dataset_id);
#  endif
        AssertThrow(pt_data_dataset >= 0, ExcIO());

        // Create the data subset we'll use to read from memory
        count[0]                 = local_node_cell_count[0];
        count[1]                 = pt_data_vector_dim;
        offset[0]                = global_node_cell_offsets[0];
        offset[1]                = 0;
        pt_data_memory_dataspace = H5Screate_simple(2, count, nullptr);
        AssertThrow(pt_data_memory_dataspace >= 0, ExcIO());

        // Select the hyperslab in the file
        pt_data_file_dataspace = H5Dget_space(pt_data_dataset);
        AssertThrow(pt_data_file_dataspace >= 0, ExcIO());
        status = H5Sselect_hyperslab(pt_data_file_dataspace,
                                     H5S_SELECT_SET,
                                     offset,
                                     nullptr,
                                     count,
                                     nullptr);
        AssertThrow(status >= 0, ExcIO());

        // And finally, write the data
        status = H5Dwrite(pt_data_dataset,
                          H5T_NATIVE_DOUBLE,
                          pt_data_memory_dataspace,
                          pt_data_file_dataspace,
                          plist_id,
                          data_filter.get_data_set(i));
        AssertThrow(status >= 0, ExcIO());

        // Close the dataspaces
        status = H5Sclose(pt_data_dataspace);
        AssertThrow(status >= 0, ExcIO());
        status = H5Sclose(pt_data_memory_dataspace);
        AssertThrow(status >= 0, ExcIO());
        status = H5Sclose(pt_data_file_dataspace);
        AssertThrow(status >= 0, ExcIO());
        // Close the dataset
        status = H5Dclose(pt_data_dataset);
        AssertThrow(status >= 0, ExcIO());
      }

    // Close the file property list
    status = H5Pclose(file_plist_id);
    AssertThrow(status >= 0, ExcIO());

    // Close the parallel access
    status = H5Pclose(plist_id);
    AssertThrow(status >= 0, ExcIO());

    // Close the file
    status = H5Fclose(h5_solution_file_id);
    AssertThrow(status >= 0, ExcIO());
  }
#endif
} // namespace



template <int dim, int spacedim>
void
DataOutBase::write_filtered_data(
  const std::vector<Patch<dim, spacedim>> &patches,
  const std::vector<std::string>          &data_names,
  const std::vector<
    std::tuple<unsigned int,
               unsigned int,
               std::string,
               DataComponentInterpretation::DataComponentInterpretation>>
                             &nonscalar_data_ranges,
  DataOutBase::DataOutFilter &filtered_data)
{
  const unsigned int n_data_sets = data_names.size();

#ifndef DEAL_II_WITH_MPI
  // verify that there are indeed patches to be written out. most of the times,
  // people just forget to call build_patches when there are no patches, so a
  // warning is in order. that said, the assertion is disabled if we support MPI
  // since then it can happen that on the coarsest mesh, a processor simply has
  // no cells it actually owns, and in that case it is legit if there are no
  // patches
  Assert(patches.size() > 0, ExcNoPatches());
#else
  if (patches.empty())
    return;
#endif

  unsigned int n_nodes;
  std::tie(n_nodes, std::ignore) = count_nodes_and_cells(patches);

  // For the format we write here, we need to write all node values relating
  // to one variable at a time. We could in principle do this by looping
  // over all patches and extracting the values corresponding to the one
  // variable we're dealing with right now, and then start the process over
  // for the next variable with another loop over all patches.
  //
  // An easier way is to create a global table that for each variable
  // lists all values. This copying of data vectors can be done in the
  // background while we're already working on vertices and cells,
  // so do this on a separate task and when wanting to write out the
  // data, we wait for that task to finish.
  Threads::Task<std::unique_ptr<Table<2, double>>>
    create_global_data_table_task = Threads::new_task(
      [&patches]() { return create_global_data_table(patches); });

  // Write the nodes/cells to the DataOutFilter object.
  write_nodes(patches, filtered_data);
  write_cells(patches, filtered_data);

  // Wait for the reordering to be done and retrieve the reordered data:
  const Table<2, double> data_vectors =
    std::move(*create_global_data_table_task.return_value());

  // when writing, first write out all vector data, then handle the scalar data
  // sets that have been left over
  unsigned int i, n_th_vector, data_set, pt_data_vector_dim;
  std::string  vector_name;
  for (n_th_vector = 0, data_set = 0; data_set < n_data_sets;)
    {
      // Advance n_th_vector to at least the current data set we are on
      while (n_th_vector < nonscalar_data_ranges.size() &&
             std::get<0>(nonscalar_data_ranges[n_th_vector]) < data_set)
        ++n_th_vector;

      // Determine the dimension of this data
      if (n_th_vector < nonscalar_data_ranges.size() &&
          std::get<0>(nonscalar_data_ranges[n_th_vector]) == data_set)
        {
          // Multiple dimensions
          pt_data_vector_dim = std::get<1>(nonscalar_data_ranges[n_th_vector]) -
                               std::get<0>(nonscalar_data_ranges[n_th_vector]) +
                               1;

          // Ensure the dimensionality of the data is correct
          AssertThrow(
            std::get<1>(nonscalar_data_ranges[n_th_vector]) >=
              std::get<0>(nonscalar_data_ranges[n_th_vector]),
            ExcLowerRange(std::get<1>(nonscalar_data_ranges[n_th_vector]),
                          std::get<0>(nonscalar_data_ranges[n_th_vector])));
          AssertThrow(
            std::get<1>(nonscalar_data_ranges[n_th_vector]) < n_data_sets,
            ExcIndexRange(std::get<1>(nonscalar_data_ranges[n_th_vector]),
                          0,
                          n_data_sets));

          // Determine the vector name. Concatenate all the component names with
          // double underscores unless a vector name has been specified
          if (!std::get<2>(nonscalar_data_ranges[n_th_vector]).empty())
            {
              vector_name = std::get<2>(nonscalar_data_ranges[n_th_vector]);
            }
          else
            {
              vector_name = "";
              for (i = std::get<0>(nonscalar_data_ranges[n_th_vector]);
                   i < std::get<1>(nonscalar_data_ranges[n_th_vector]);
                   ++i)
                vector_name += data_names[i] + "__";
              vector_name +=
                data_names[std::get<1>(nonscalar_data_ranges[n_th_vector])];
            }
        }
      else
        {
          // One dimension
          pt_data_vector_dim = 1;
          vector_name        = data_names[data_set];
        }

      // Write data to the filter object
      filtered_data.write_data_set(vector_name,
                                   pt_data_vector_dim,
                                   data_set,
                                   data_vectors);

      // Advance the current data set
      data_set += pt_data_vector_dim;
    }
}



template <int dim, int spacedim>
void
DataOutInterface<dim, spacedim>::write_hdf5_parallel(
  const DataOutBase::DataOutFilter &data_filter,
  const std::string                &filename,
  const MPI_Comm                    comm) const
{
  DataOutBase::write_hdf5_parallel(
    get_patches(), data_filter, hdf5_flags, filename, comm);
}



template <int dim, int spacedim>
void
DataOutInterface<dim, spacedim>::write_hdf5_parallel(
  const DataOutBase::DataOutFilter &data_filter,
  const bool                        write_mesh_file,
  const std::string                &mesh_filename,
  const std::string                &solution_filename,
  const MPI_Comm                    comm) const
{
  DataOutBase::write_hdf5_parallel(get_patches(),
                                   data_filter,
                                   hdf5_flags,
                                   write_mesh_file,
                                   mesh_filename,
                                   solution_filename,
                                   comm);
}



template <int dim, int spacedim>
void
DataOutBase::write_hdf5_parallel(
  const std::vector<Patch<dim, spacedim>> &patches,
  const DataOutBase::DataOutFilter        &data_filter,
  const DataOutBase::Hdf5Flags            &flags,
  const std::string                       &filename,
  const MPI_Comm                           comm)
{
  write_hdf5_parallel(
    patches, data_filter, flags, true, filename, filename, comm);
}



template <int dim, int spacedim>
void
DataOutBase::write_hdf5_parallel(
  const std::vector<Patch<dim, spacedim>> &patches,
  const DataOutBase::DataOutFilter        &data_filter,
  const DataOutBase::Hdf5Flags            &flags,
  const bool                               write_mesh_file,
  const std::string                       &mesh_filename,
  const std::string                       &solution_filename,
  const MPI_Comm                           comm)
{
  AssertThrow(
    spacedim >= 2,
    ExcMessage(
      "DataOutBase was asked to write HDF5 output for a space dimension of 1. "
      "HDF5 only supports datasets that live in 2 or 3 dimensions."));

#ifndef DEAL_II_WITH_HDF5
  // throw an exception, but first make sure the compiler does not warn about
  // the now unused function arguments
  (void)patches;
  (void)data_filter;
  (void)flags;
  (void)write_mesh_file;
  (void)mesh_filename;
  (void)solution_filename;
  (void)comm;
  AssertThrow(false, ExcNeedsHDF5());
#else

  const unsigned int n_ranks = Utilities::MPI::n_mpi_processes(comm);
  (void)n_ranks;

  // If HDF5 is not parallel and we're using multiple processes, abort:
#  ifndef H5_HAVE_PARALLEL
  AssertThrow(
    n_ranks <= 1,
    ExcMessage(
      "Serial HDF5 output on multiple processes is not yet supported."));
#  endif

  // Verify that there are indeed patches to be written out. most of
  // the times, people just forget to call build_patches when there
  // are no patches, so a warning is in order. That said, the
  // assertion is disabled if we run with more than one MPI rank,
  // since then it can happen that, on coarse meshes, a processor
  // simply has no cells it actually owns, and in that case it is
  // legit if there are no patches.
  Assert((patches.size() > 0) || (n_ranks > 1), ExcNoPatches());

  // The HDF5 routines perform a bunch of collective calls that expect all
  // ranks to participate. One ranks without any patches we are missing
  // critical information, so rather than broadcasting that information, just
  // create a new communicator that only contains ranks with cells and
  // use that to perform the write operations:
  const bool have_patches = (patches.size() > 0);
  MPI_Comm   split_comm;
  {
    const int key   = Utilities::MPI::this_mpi_process(comm);
    const int color = (have_patches ? 1 : 0);
    const int ierr  = MPI_Comm_split(comm, color, key, &split_comm);
    AssertThrowMPI(ierr);
  }

  if (have_patches)
    {
      do_write_hdf5<dim, spacedim>(patches,
                                   data_filter,
                                   flags,
                                   write_mesh_file,
                                   mesh_filename,
                                   solution_filename,
                                   split_comm);
    }

  const int ierr = MPI_Comm_free(&split_comm);
  AssertThrowMPI(ierr);

#endif
}



template <int dim, int spacedim>
void
DataOutInterface<dim, spacedim>::write(
  std::ostream                   &out,
  const DataOutBase::OutputFormat output_format_) const
{
  DataOutBase::OutputFormat output_format = output_format_;
  if (output_format == DataOutBase::default_format)
    output_format = default_fmt;

  switch (output_format)
    {
      case DataOutBase::none:
        break;

      case DataOutBase::dx:
        write_dx(out);
        break;

      case DataOutBase::ucd:
        write_ucd(out);
        break;

      case DataOutBase::gnuplot:
        write_gnuplot(out);
        break;

      case DataOutBase::povray:
        write_povray(out);
        break;

      case DataOutBase::eps:
        write_eps(out);
        break;

      case DataOutBase::gmv:
        write_gmv(out);
        break;

      case DataOutBase::tecplot:
        write_tecplot(out);
        break;

      case DataOutBase::vtk:
        write_vtk(out);
        break;

      case DataOutBase::vtu:
        write_vtu(out);
        break;

      case DataOutBase::svg:
        write_svg(out);
        break;

      case DataOutBase::deal_II_intermediate:
        write_deal_II_intermediate(out);
        break;

      default:
        DEAL_II_NOT_IMPLEMENTED();
    }
}



template <int dim, int spacedim>
void
DataOutInterface<dim, spacedim>::set_default_format(
  const DataOutBase::OutputFormat fmt)
{
  Assert(fmt != DataOutBase::default_format, ExcNotImplemented());
  default_fmt = fmt;
}

template <int dim, int spacedim>
template <typename FlagType>
void
DataOutInterface<dim, spacedim>::set_flags(const FlagType &flags)
{
  if constexpr (std::is_same_v<FlagType, DataOutBase::DXFlags>)
    dx_flags = flags;
  else if constexpr (std::is_same_v<FlagType, DataOutBase::UcdFlags>)
    ucd_flags = flags;
  else if constexpr (std::is_same_v<FlagType, DataOutBase::PovrayFlags>)
    povray_flags = flags;
  else if constexpr (std::is_same_v<FlagType, DataOutBase::EpsFlags>)
    eps_flags = flags;
  else if constexpr (std::is_same_v<FlagType, DataOutBase::GmvFlags>)
    gmv_flags = flags;
  else if constexpr (std::is_same_v<FlagType, DataOutBase::Hdf5Flags>)
    hdf5_flags = flags;
  else if constexpr (std::is_same_v<FlagType, DataOutBase::TecplotFlags>)
    tecplot_flags = flags;
  else if constexpr (std::is_same_v<FlagType, DataOutBase::VtkFlags>)
    vtk_flags = flags;
  else if constexpr (std::is_same_v<FlagType, DataOutBase::SvgFlags>)
    svg_flags = flags;
  else if constexpr (std::is_same_v<FlagType, DataOutBase::GnuplotFlags>)
    gnuplot_flags = flags;
  else if constexpr (std::is_same_v<FlagType,
                                    DataOutBase::Deal_II_IntermediateFlags>)
    deal_II_intermediate_flags = flags;
  else
    DEAL_II_NOT_IMPLEMENTED();
}



template <int dim, int spacedim>
std::string
DataOutInterface<dim, spacedim>::default_suffix(
  const DataOutBase::OutputFormat output_format) const
{
  if (output_format == DataOutBase::default_format)
    return DataOutBase::default_suffix(default_fmt);
  else
    return DataOutBase::default_suffix(output_format);
}



template <int dim, int spacedim>
void
DataOutInterface<dim, spacedim>::declare_parameters(ParameterHandler &prm)
{
  prm.declare_entry("Output format",
                    "gnuplot",
                    Patterns::Selection(DataOutBase::get_output_format_names()),
                    "A name for the output format to be used");
  prm.declare_entry("Subdivisions",
                    "1",
                    Patterns::Integer(),
                    "Number of subdivisions of each mesh cell");

  prm.enter_subsection("DX output parameters");
  DataOutBase::DXFlags::declare_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("UCD output parameters");
  DataOutBase::UcdFlags::declare_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("Gnuplot output parameters");
  DataOutBase::GnuplotFlags::declare_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("Povray output parameters");
  DataOutBase::PovrayFlags::declare_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("Eps output parameters");
  DataOutBase::EpsFlags::declare_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("Gmv output parameters");
  DataOutBase::GmvFlags::declare_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("HDF5 output parameters");
  DataOutBase::Hdf5Flags::declare_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("Tecplot output parameters");
  DataOutBase::TecplotFlags::declare_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("Vtk output parameters");
  DataOutBase::VtkFlags::declare_parameters(prm);
  prm.leave_subsection();


  prm.enter_subsection("deal.II intermediate output parameters");
  DataOutBase::Deal_II_IntermediateFlags::declare_parameters(prm);
  prm.leave_subsection();
}



template <int dim, int spacedim>
void
DataOutInterface<dim, spacedim>::parse_parameters(ParameterHandler &prm)
{
  const std::string &output_name = prm.get("Output format");
  default_fmt          = DataOutBase::parse_output_format(output_name);
  default_subdivisions = prm.get_integer("Subdivisions");

  prm.enter_subsection("DX output parameters");
  dx_flags.parse_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("UCD output parameters");
  ucd_flags.parse_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("Gnuplot output parameters");
  gnuplot_flags.parse_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("Povray output parameters");
  povray_flags.parse_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("Eps output parameters");
  eps_flags.parse_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("Gmv output parameters");
  gmv_flags.parse_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("HDF5 output parameters");
  hdf5_flags.parse_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("Tecplot output parameters");
  tecplot_flags.parse_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("Vtk output parameters");
  vtk_flags.parse_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("deal.II intermediate output parameters");
  deal_II_intermediate_flags.parse_parameters(prm);
  prm.leave_subsection();
}



template <int dim, int spacedim>
std::size_t
DataOutInterface<dim, spacedim>::memory_consumption() const
{
  return (sizeof(default_fmt) +
          MemoryConsumption::memory_consumption(dx_flags) +
          MemoryConsumption::memory_consumption(ucd_flags) +
          MemoryConsumption::memory_consumption(gnuplot_flags) +
          MemoryConsumption::memory_consumption(povray_flags) +
          MemoryConsumption::memory_consumption(eps_flags) +
          MemoryConsumption::memory_consumption(gmv_flags) +
          MemoryConsumption::memory_consumption(hdf5_flags) +
          MemoryConsumption::memory_consumption(tecplot_flags) +
          MemoryConsumption::memory_consumption(vtk_flags) +
          MemoryConsumption::memory_consumption(svg_flags) +
          MemoryConsumption::memory_consumption(deal_II_intermediate_flags));
}



template <int dim, int spacedim>
std::vector<
  std::tuple<unsigned int,
             unsigned int,
             std::string,
             DataComponentInterpretation::DataComponentInterpretation>>
DataOutInterface<dim, spacedim>::get_nonscalar_data_ranges() const
{
  return std::vector<
    std::tuple<unsigned int,
               unsigned int,
               std::string,
               DataComponentInterpretation::DataComponentInterpretation>>();
}


template <int dim, int spacedim>
void
DataOutInterface<dim, spacedim>::validate_dataset_names() const
{
  if constexpr (running_in_debug_mode())
    {
      {
        // Check that names for datasets are only used once. This is somewhat
        // complicated, because vector ranges might have a name or not.
        std::set<std::string> all_names;

        const std::vector<
          std::tuple<unsigned int,
                     unsigned int,
                     std::string,
                     DataComponentInterpretation::DataComponentInterpretation>>
          ranges = this->get_nonscalar_data_ranges();
        const std::vector<std::string> data_names  = this->get_dataset_names();
        const unsigned int             n_data_sets = data_names.size();
        std::vector<bool>              data_set_written(n_data_sets, false);

        for (const auto &range : ranges)
          {
            const std::string &name = std::get<2>(range);
            if (!name.empty())
              {
                Assert(all_names.find(name) == all_names.end(),
                       ExcMessage(
                         "Error: names of fields in DataOut need to be unique, "
                         "but '" +
                         name + "' is used more than once."));
                all_names.insert(name);
                for (unsigned int i = std::get<0>(range);
                     i <= std::get<1>(range);
                     ++i)
                  data_set_written[i] = true;
              }
          }

        for (unsigned int data_set = 0; data_set < n_data_sets; ++data_set)
          if (data_set_written[data_set] == false)
            {
              const std::string &name = data_names[data_set];
              Assert(all_names.find(name) == all_names.end(),
                     ExcMessage(
                       "Error: names of fields in DataOut need to be unique, "
                       "but '" +
                       name + "' is used more than once."));
              all_names.insert(name);
            }
      }
    }
}



// ---------------------------------------------- DataOutReader ----------

template <int dim, int spacedim>
void
DataOutReader<dim, spacedim>::read(std::istream &in)
{
  AssertThrow(in.fail() == false, ExcIO());

  // first empty previous content
  {
    std::vector<typename dealii::DataOutBase::Patch<dim, spacedim>> tmp;
    tmp.swap(patches);
  }
  {
    std::vector<std::string> tmp;
    tmp.swap(dataset_names);
  }
  {
    std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>
      tmp;
    tmp.swap(nonscalar_data_ranges);
  }

  // then check that we have the correct header of this file. both the first and
  // second real lines have to match, as well as the dimension information
  // written before that and the Version information written in the third line
  {
    std::pair<unsigned int, unsigned int> dimension_info =
      DataOutBase::determine_intermediate_format_dimensions(in);
    AssertThrow((dimension_info.first == dim) &&
                  (dimension_info.second == spacedim),
                ExcIncompatibleDimensions(
                  dimension_info.first, dim, dimension_info.second, spacedim));

    // read to the end of the line
    std::string tmp;
    getline(in, tmp);
  }

  {
    std::string header;
    getline(in, header);

    std::ostringstream s;
    s << "[deal.II intermediate format graphics data]";

    Assert(header == s.str(), ExcUnexpectedInput(s.str(), header));
  }
  {
    std::string header;
    getline(in, header);

    std::ostringstream s;
    s << "[written by " << DEAL_II_PACKAGE_NAME << " "
      << DEAL_II_PACKAGE_VERSION << "]";

    Assert(header == s.str(), ExcUnexpectedInput(s.str(), header));
  }
  {
    std::string header;
    getline(in, header);

    std::ostringstream s;
    s << "[Version: "
      << dealii::DataOutBase::Deal_II_IntermediateFlags::format_version << "]";

    Assert(header == s.str(),
           ExcMessage(
             "Invalid or incompatible file format. Intermediate format "
             "files can only be read by the same deal.II version as they "
             "are written by."));
  }

  // then read the rest of the data
  unsigned int n_datasets;
  in >> n_datasets;
  dataset_names.resize(n_datasets);
  for (unsigned int i = 0; i < n_datasets; ++i)
    in >> dataset_names[i];

  unsigned int n_patches;
  in >> n_patches;
  patches.resize(n_patches);
  for (unsigned int i = 0; i < n_patches; ++i)
    in >> patches[i];

  unsigned int n_nonscalar_data_ranges;
  in >> n_nonscalar_data_ranges;
  nonscalar_data_ranges.resize(n_nonscalar_data_ranges);
  for (unsigned int i = 0; i < n_nonscalar_data_ranges; ++i)
    {
      in >> std::get<0>(nonscalar_data_ranges[i]) >>
        std::get<1>(nonscalar_data_ranges[i]);

      // read in the name of that vector range. because it is on a separate
      // line, we first need to read to the end of the previous line (nothing
      // should be there any more after we've read the previous two integers)
      // and then read the entire next line for the name
      std::string name;
      getline(in, name);
      getline(in, name);
      std::get<2>(nonscalar_data_ranges[i]) = name;
    }

  AssertThrow(in.fail() == false, ExcIO());
}



template <int dim, int spacedim>
void
DataOutReader<dim, spacedim>::read_whole_parallel_file(std::istream &in)
{
  AssertThrow(in.fail() == false, ExcIO());

  ParallelIntermediateHeader header;
  in.read(reinterpret_cast<char *>(&header), sizeof(header));
  AssertThrow(
    header.magic == 0x00dea111,
    ExcMessage(
      "Invalid header of parallel deal.II intermediate format encountered."));
  AssertThrow(
    header.version == DataOutBase::Deal_II_IntermediateFlags::format_version,
    ExcMessage(
      "Incorrect header version of parallel deal.II intermediate format."));

  std::vector<std::uint64_t> chunk_sizes(header.n_ranks);
  in.read(reinterpret_cast<char *>(chunk_sizes.data()),
          header.n_ranks * sizeof(std::uint64_t));

  for (unsigned int n = 0; n < header.n_ranks; ++n)
    {
      // First read the compressed data into temp_buffer and then
      // decompress and put into datastream
      std::vector<char> temp_buffer(chunk_sizes[n]);
      in.read(temp_buffer.data(), chunk_sizes[n]);

      AssertThrow(static_cast<DataOutBase::CompressionLevel>(
                    header.compression) !=
                    DataOutBase::CompressionLevel::plain_text,
                  ExcNotImplemented());

      boost::iostreams::filtering_istreambuf f;
      if (static_cast<DataOutBase::CompressionLevel>(header.compression) !=
          DataOutBase::CompressionLevel::no_compression)
#ifdef DEAL_II_WITH_ZLIB
        f.push(boost::iostreams::zlib_decompressor());
#else
        AssertThrow(
          false,
          ExcMessage(
            "Decompression requires deal.II to be configured with ZLIB support."));
#endif

      boost::iostreams::basic_array_source<char> source(temp_buffer.data(),
                                                        temp_buffer.size());
      f.push(source);

      std::stringstream datastream;
      boost::iostreams::copy(f, datastream);

      // Now we can load the data and merge this chunk into *this
      if (n == 0)
        {
          read(datastream);
        }
      else
        {
          DataOutReader<dim, spacedim> temp_reader;
          temp_reader.read(datastream);
          merge(temp_reader);
        }
    }
}



template <int dim, int spacedim>
void
DataOutReader<dim, spacedim>::merge(const DataOutReader<dim, spacedim> &source)
{
  using Patch = typename dealii::DataOutBase::Patch<dim, spacedim>;


  const std::vector<Patch> &source_patches = source.get_patches();
  Assert(patches.size() != 0, DataOutBase::ExcNoPatches());
  Assert(source_patches.size() != 0, DataOutBase::ExcNoPatches());
  // check equality of component names
  Assert(get_dataset_names() == source.get_dataset_names(),
         ExcIncompatibleDatasetNames());

  // check equality of the vector data specifications
  Assert(get_nonscalar_data_ranges().size() ==
           source.get_nonscalar_data_ranges().size(),
         ExcMessage("Both sources need to declare the same components "
                    "as vectors."));
  for (unsigned int i = 0; i < get_nonscalar_data_ranges().size(); ++i)
    {
      Assert(std::get<0>(get_nonscalar_data_ranges()[i]) ==
               std::get<0>(source.get_nonscalar_data_ranges()[i]),
             ExcMessage("Both sources need to declare the same components "
                        "as vectors."));
      Assert(std::get<1>(get_nonscalar_data_ranges()[i]) ==
               std::get<1>(source.get_nonscalar_data_ranges()[i]),
             ExcMessage("Both sources need to declare the same components "
                        "as vectors."));
      Assert(std::get<2>(get_nonscalar_data_ranges()[i]) ==
               std::get<2>(source.get_nonscalar_data_ranges()[i]),
             ExcMessage("Both sources need to declare the same components "
                        "as vectors."));
    }

  // make sure patches are compatible
  Assert(patches[0].n_subdivisions == source_patches[0].n_subdivisions,
         ExcIncompatiblePatchLists());
  Assert(patches[0].data.n_rows() == source_patches[0].data.n_rows(),
         ExcIncompatiblePatchLists());
  Assert(patches[0].data.n_cols() == source_patches[0].data.n_cols(),
         ExcIncompatiblePatchLists());

  // merge patches. store old number of elements, since we need to adjust patch
  // numbers, etc afterwards
  const unsigned int old_n_patches = patches.size();
  patches.insert(patches.end(), source_patches.begin(), source_patches.end());

  // adjust patch numbers
  for (unsigned int i = old_n_patches; i < patches.size(); ++i)
    patches[i].patch_index += old_n_patches;

  // adjust patch neighbors
  for (unsigned int i = old_n_patches; i < patches.size(); ++i)
    for (const unsigned int n : GeometryInfo<dim>::face_indices())
      if (patches[i].neighbors[n] !=
          dealii::DataOutBase::Patch<dim, spacedim>::no_neighbor)
        patches[i].neighbors[n] += old_n_patches;
}



template <int dim, int spacedim>
const std::vector<typename dealii::DataOutBase::Patch<dim, spacedim>> &
DataOutReader<dim, spacedim>::get_patches() const
{
  return patches;
}



template <int dim, int spacedim>
std::vector<std::string>
DataOutReader<dim, spacedim>::get_dataset_names() const
{
  return dataset_names;
}



template <int dim, int spacedim>
std::vector<
  std::tuple<unsigned int,
             unsigned int,
             std::string,
             DataComponentInterpretation::DataComponentInterpretation>>
DataOutReader<dim, spacedim>::get_nonscalar_data_ranges() const
{
  return nonscalar_data_ranges;
}



// ---------------------------------------------- XDMFEntry ----------

XDMFEntry::XDMFEntry()
  : valid(false)
  , h5_sol_filename("")
  , h5_mesh_filename("")
  , entry_time(0.0)
  , num_nodes(numbers::invalid_unsigned_int)
  , num_cells(numbers::invalid_unsigned_int)
  , dimension(numbers::invalid_unsigned_int)
  , space_dimension(numbers::invalid_unsigned_int)
  , cell_type()
{}



XDMFEntry::XDMFEntry(const std::string   &filename,
                     const double         time,
                     const std::uint64_t  nodes,
                     const std::uint64_t  cells,
                     const unsigned int   dim,
                     const ReferenceCell &cell_type)
  : XDMFEntry(filename, filename, time, nodes, cells, dim, dim, cell_type)
{}



XDMFEntry::XDMFEntry(const std::string   &mesh_filename,
                     const std::string   &solution_filename,
                     const double         time,
                     const std::uint64_t  nodes,
                     const std::uint64_t  cells,
                     const unsigned int   dim,
                     const ReferenceCell &cell_type)
  : XDMFEntry(mesh_filename,
              solution_filename,
              time,
              nodes,
              cells,
              dim,
              dim,
              cell_type)
{}



namespace
{
  /**
   * Deprecated XDMFEntry constructors do not fill the cell_type, so we use this
   * little helper to convert it to the appropriate hex cell.
   */
  ReferenceCell
  cell_type_hex_if_invalid(const ReferenceCell &cell_type,
                           const unsigned int   dimension)
  {
    if (cell_type == ReferenceCells::Invalid)
      {
        switch (dimension)
          {
            case 0:
              return ReferenceCells::get_hypercube<0>();
            case 1:
              return ReferenceCells::get_hypercube<1>();
            case 2:
              return ReferenceCells::get_hypercube<2>();
            case 3:
              return ReferenceCells::get_hypercube<3>();
            default:
              AssertThrow(false, ExcMessage("Invalid dimension"));
          }
      }
    else
      return cell_type;
  }
} // namespace



XDMFEntry::XDMFEntry(const std::string   &mesh_filename,
                     const std::string   &solution_filename,
                     const double         time,
                     const std::uint64_t  nodes,
                     const std::uint64_t  cells,
                     const unsigned int   dim,
                     const unsigned int   spacedim,
                     const ReferenceCell &cell_type_)
  : valid(true)
  , h5_sol_filename(solution_filename)
  , h5_mesh_filename(mesh_filename)
  , entry_time(time)
  , num_nodes(nodes)
  , num_cells(cells)
  , dimension(dim)
  , space_dimension(spacedim)
  , cell_type(cell_type_hex_if_invalid(cell_type_, dim))
{}



void
XDMFEntry::add_attribute(const std::string &attr_name,
                         const unsigned int dimension)
{
  attribute_dims[attr_name] = dimension;
}



namespace
{
  /**
   * Small function to create indentation for XML file.
   */
  std::string
  indent(const unsigned int indent_level)
  {
    std::string res = "";
    for (unsigned int i = 0; i < indent_level; ++i)
      res += "  ";
    return res;
  }
} // namespace



std::string
XDMFEntry::get_xdmf_content(const unsigned int indent_level) const
{
  if (!valid)
    return "";

  std::stringstream ss;
  ss.precision(12);
  ss << indent(indent_level + 0)
     << "<Grid Name=\"mesh\" GridType=\"Uniform\">\n";
  ss << indent(indent_level + 1) << "<Time Value=\"" << entry_time << "\"/>\n";
  ss << indent(indent_level + 1) << "<Geometry GeometryType=\""
     << (space_dimension <= 2 ? "XY" : "XYZ") << "\">\n";
  ss << indent(indent_level + 2) << "<DataItem Dimensions=\"" << num_nodes
     << " " << (space_dimension <= 2 ? 2 : space_dimension)
     << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n";
  ss << indent(indent_level + 3) << h5_mesh_filename << ":/nodes\n";
  ss << indent(indent_level + 2) << "</DataItem>\n";
  ss << indent(indent_level + 1) << "</Geometry>\n";

  // If we have cells defined, use the topology corresponding to the dimension
  if (num_cells > 0)
    {
      ss << indent(indent_level + 1) << "<Topology TopologyType=\"";

      if (dimension == 0)
        {
          ss << "Polyvertex";
        }
      else if (dimension == 1)
        {
          ss << "Polyline";
        }
      else if (dimension == 2)
        {
          Assert(cell_type == ReferenceCells::Quadrilateral ||
                   cell_type == ReferenceCells::Triangle,
                 ExcNotImplemented());

          if (cell_type == ReferenceCells::Quadrilateral)
            {
              ss << "Quadrilateral";
            }
          else // if (cell_type == ReferenceCells::Triangle)
            {
              ss << "Triangle";
            }
        }
      else if (dimension == 3)
        {
          Assert(cell_type == ReferenceCells::Hexahedron ||
                   cell_type == ReferenceCells::Tetrahedron,
                 ExcNotImplemented());

          if (cell_type == ReferenceCells::Hexahedron)
            {
              ss << "Hexahedron";
            }
          else // if (reference_cell == ReferenceCells::Tetrahedron)
            {
              ss << "Tetrahedron";
            }
        }

      ss << "\" NumberOfElements=\"" << num_cells;
      if (dimension == 0)
        ss << "\" NodesPerElement=\"1\">\n";
      else if (dimension == 1)
        ss << "\" NodesPerElement=\"2\">\n";
      else
        // no "NodesPerElement" for dimension 2 and higher
        ss << "\">\n";

      ss << indent(indent_level + 2) << "<DataItem Dimensions=\"" << num_cells
         << " " << cell_type.n_vertices()
         << "\" NumberType=\"UInt\" Format=\"HDF\">\n";

      ss << indent(indent_level + 3) << h5_mesh_filename << ":/cells\n";
      ss << indent(indent_level + 2) << "</DataItem>\n";
      ss << indent(indent_level + 1) << "</Topology>\n";
    }
  // Otherwise, we assume the points are isolated in space and use a Polyvertex
  // topology
  else
    {
      ss << indent(indent_level + 1)
         << "<Topology TopologyType=\"Polyvertex\" NumberOfElements=\""
         << num_nodes << "\">\n";
      ss << indent(indent_level + 1) << "</Topology>\n";
    }

  for (const auto &attribute_dim : attribute_dims)
    {
      ss << indent(indent_level + 1) << "<Attribute Name=\""
         << attribute_dim.first << "\" AttributeType=\""
         << (attribute_dim.second > 1 ? "Vector" : "Scalar")
         << "\" Center=\"Node\">\n";
      // Vectors must have 3 elements even for 2d models
      ss << indent(indent_level + 2) << "<DataItem Dimensions=\"" << num_nodes
         << " " << (attribute_dim.second > 1 ? 3 : 1)
         << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n";
      ss << indent(indent_level + 3) << h5_sol_filename << ":/"
         << attribute_dim.first << '\n';
      ss << indent(indent_level + 2) << "</DataItem>\n";
      ss << indent(indent_level + 1) << "</Attribute>\n";
    }

  ss << indent(indent_level + 0) << "</Grid>\n";

  return ss.str();
}



namespace DataOutBase
{
  template <int dim, int spacedim>
  std::ostream &
  operator<<(std::ostream &out, const Patch<dim, spacedim> &patch)
  {
    // write a header line
    out << "[deal.II intermediate Patch<" << dim << ',' << spacedim << ">]"
        << '\n';

    // First export what kind of reference cell we are looking at:
    out << patch.reference_cell << '\n';

    // then write all the data that is in this patch
    for (const unsigned int i : patch.reference_cell.vertex_indices())
      out << patch.vertices[i] << ' ';
    out << '\n';

    for (const unsigned int i : patch.reference_cell.face_indices())
      out << patch.neighbors[i] << ' ';
    out << '\n';

    out << patch.patch_index << ' ' << patch.n_subdivisions << '\n';

    out << patch.points_are_available << '\n';

    out << patch.data.n_rows() << ' ' << patch.data.n_cols() << '\n';
    for (unsigned int i = 0; i < patch.data.n_rows(); ++i)
      for (unsigned int j = 0; j < patch.data.n_cols(); ++j)
        out << patch.data[i][j] << ' ';
    out << '\n';
    out << '\n';

    return out;
  }



  template <int dim, int spacedim>
  std::istream &
  operator>>(std::istream &in, Patch<dim, spacedim> &patch)
  {
    AssertThrow(in.fail() == false, ExcIO());

    // read a header line and compare it to what we usually write. skip all
    // lines that contain only blanks at the start
    {
      std::string header;
      do
        {
          getline(in, header);
          while ((header.size() != 0) && (header.back() == ' '))
            header.erase(header.size() - 1);
        }
      while ((header.empty()) && in);

      std::ostringstream s;
      s << "[deal.II intermediate Patch<" << dim << ',' << spacedim << ">]";

      Assert(header == s.str(), ExcUnexpectedInput(s.str(), header));
    }

    // First import what kind of reference cell we are looking at:
    if constexpr (dim > 0)
      in >> patch.reference_cell;

    // then read all the data that is in this patch
    for (const unsigned int i : patch.reference_cell.vertex_indices())
      in >> patch.vertices[i];

    for (const unsigned int i : patch.reference_cell.face_indices())
      in >> patch.neighbors[i];

    in >> patch.patch_index;

    // If dim>1, we also need to set the number of subdivisions, whereas
    // in dim==1, this is a const variable equal to one that can't be changed.
    unsigned int n_subdivisions;
    in >> n_subdivisions;
    if constexpr (dim > 1)
      patch.n_subdivisions = n_subdivisions;

    in >> patch.points_are_available;

    unsigned int n_rows, n_cols;
    in >> n_rows >> n_cols;
    patch.data.reinit(n_rows, n_cols);
    for (unsigned int i = 0; i < patch.data.n_rows(); ++i)
      for (unsigned int j = 0; j < patch.data.n_cols(); ++j)
        in >> patch.data[i][j];

    AssertThrow(in.fail() == false, ExcIO());

    return in;
  }
} // namespace DataOutBase



// explicit instantiations
#include "base/data_out_base.inst"

DEAL_II_NAMESPACE_CLOSE
