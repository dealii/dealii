// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/bounding_box.h>
#include <deal.II/base/floating_point_comparator.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/base/utilities.h>

#include <deal.II/grid/grid_tools_geometry.h>
#include <deal.II/grid/grid_tools_topology.h>

#include <boost/container/small_vector.hpp>

#include <algorithm>
#include <map>
#include <numeric>
#include <set>
#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace GridTools
{
  // Generic functions for appending face data in 2d or 3d. TODO: we can
  // remove these once we have 'if constexpr'.
  namespace internal
  {
    inline void
    append_face_data(const CellData<1> &face_data, SubCellData &subcell_data)
    {
      subcell_data.boundary_lines.push_back(face_data);
    }



    inline void
    append_face_data(const CellData<2> &face_data, SubCellData &subcell_data)
    {
      subcell_data.boundary_quads.push_back(face_data);
    }



    // Lexical comparison for sorting CellData objects.
    template <int structdim>
    struct CellDataComparator
    {
      bool
      operator()(const CellData<structdim> &a,
                 const CellData<structdim> &b) const
      {
        // Check vertices:
        if (std::lexicographical_compare(std::begin(a.vertices),
                                         std::end(a.vertices),
                                         std::begin(b.vertices),
                                         std::end(b.vertices)))
          return true;
        // it should never be necessary to check the material or manifold
        // ids as a 'tiebreaker' (since they must be equal if the vertex
        // indices are equal). Assert it anyway:
        if constexpr (running_in_debug_mode())
          {
            if (std::equal(std::begin(a.vertices),
                           std::end(a.vertices),
                           std::begin(b.vertices)))
              {
                Assert(a.material_id == b.material_id &&
                         a.manifold_id == b.manifold_id,
                       ExcMessage(
                         "Two CellData objects with equal vertices must "
                         "have the same material/boundary ids and manifold "
                         "ids."));
              }
          }
        return false;
      }
    };


    /**
     * get_coarse_mesh_description() needs to store face data for dim>1, but
     * we can not have this code in the function, as this requires either an
     * instantiation of CellData<dim-1>, or constexpr if. We use a class with
     * specialization instead for now.
     *
     * Data on faces is added with insert_face_data() and then retrieved with
     * get().
     */
    template <int dim>
    class FaceDataHelper
    {
    public:
      /**
       * Store the data about the given face @p face.
       */
      template <typename FaceIteratorType>
      void
      insert_face_data(const FaceIteratorType &face)
      {
        CellData<dim - 1> face_cell_data(face->n_vertices());
        for (unsigned int vertex_n = 0; vertex_n < face->n_vertices();
             ++vertex_n)
          face_cell_data.vertices[vertex_n] = face->vertex_index(vertex_n);
        face_cell_data.boundary_id = face->boundary_id();
        face_cell_data.manifold_id = face->manifold_id();

        face_data.insert(std::move(face_cell_data));
      }

      /**
       * Return  the @p subcell_data with the stored information.
       */
      SubCellData
      get()
      {
        SubCellData subcell_data;

        for (const CellData<dim - 1> &face_cell_data : face_data)
          internal::append_face_data(face_cell_data, subcell_data);
        return subcell_data;
      }


    private:
      std::set<CellData<dim - 1>, internal::CellDataComparator<dim - 1>>
        face_data;
    };


    // Do nothing for dim=1:
    template <>
    class FaceDataHelper<1>
    {
    public:
      template <typename FaceIteratorType>
      void
      insert_face_data(const FaceIteratorType &)
      {}

      SubCellData
      get()
      {
        return SubCellData();
      }
    };
  } // namespace internal



  template <int dim, int spacedim>
  std::
    tuple<std::vector<Point<spacedim>>, std::vector<CellData<dim>>, SubCellData>
    get_coarse_mesh_description(const Triangulation<dim, spacedim> &tria)
  {
    Assert(tria.n_levels() >= 1,
           ExcMessage("The input triangulation must be non-empty."));

    std::vector<Point<spacedim>> vertices = tria.get_vertices();
    std::vector<CellData<dim>>   cells;

    internal::FaceDataHelper<dim> face_data;
    std::set<CellData<1>, internal::CellDataComparator<1>>
      line_data; // only used in 3d

    for (const auto &cell : tria.cell_iterators_on_level(0))
      {
        // Save cell data
        CellData<dim> cell_data(cell->n_vertices());
        for (const unsigned int cell_vertex_n : cell->vertex_indices())
          {
            Assert(cell->vertex_index(cell_vertex_n) < vertices.size(),
                   ExcInternalError());
            cell_data.vertices[cell_vertex_n] =
              cell->vertex_index(cell_vertex_n);
          }
        cell_data.material_id = cell->material_id();
        cell_data.manifold_id = cell->manifold_id();
        cells.emplace_back(std::move(cell_data));

        // Save face data
        if (dim > 1)
          {
            for (const unsigned int face_n : cell->face_indices())
              // We don't need to insert anything if we have default values
              {
                const auto face = cell->face(face_n);
                if (face->boundary_id() != numbers::internal_face_boundary_id ||
                    face->manifold_id() != numbers::flat_manifold_id)
                  face_data.insert_face_data(face);
              }
          }
        // Save line data
        if (dim == 3)
          {
            for (unsigned int line_n = 0; line_n < cell->n_lines(); ++line_n)
              {
                const auto line = cell->line(line_n);
                // We don't need to insert anything if we have default values
                if (line->boundary_id() != numbers::internal_face_boundary_id ||
                    line->manifold_id() != numbers::flat_manifold_id)
                  {
                    CellData<1> line_cell_data(line->n_vertices());
                    for (const unsigned int vertex_n : line->vertex_indices())
                      line_cell_data.vertices[vertex_n] =
                        line->vertex_index(vertex_n);
                    line_cell_data.boundary_id = line->boundary_id();
                    line_cell_data.manifold_id = line->manifold_id();
                    line_data.insert(std::move(line_cell_data));
                  }
              }
          }
      }

    SubCellData subcell_data = face_data.get();

    if (dim == 3)
      for (const CellData<1> &face_line_data : line_data)
        subcell_data.boundary_lines.push_back(face_line_data);

    // We end up with a 'vertices' array that uses some of the entries,
    // but not all -- specifically, all vertices referenced by level-0
    // cells. We can compress the array:
    GridTools::delete_unused_vertices(vertices, cells, subcell_data);

    return std::tuple<std::vector<Point<spacedim>>,
                      std::vector<CellData<dim>>,
                      SubCellData>(std::move(vertices),
                                   std::move(cells),
                                   std::move(subcell_data));
  }



  template <int dim, int spacedim>
  void
  delete_unused_vertices(std::vector<Point<spacedim>> &vertices,
                         std::vector<CellData<dim>>   &cells,
                         SubCellData                  &subcelldata)
  {
    Assert(
      subcelldata.check_consistency(dim),
      ExcMessage(
        "Invalid SubCellData supplied according to ::check_consistency(). "
        "This is caused by data containing objects for the wrong dimension."));

    // first check which vertices are actually used
    std::vector<bool> vertex_used(vertices.size(), false);
    for (unsigned int c = 0; c < cells.size(); ++c)
      for (unsigned int v = 0; v < cells[c].vertices.size(); ++v)
        {
          Assert(cells[c].vertices[v] < vertices.size(),
                 ExcMessage("Invalid vertex index encountered! cells[" +
                            Utilities::int_to_string(c) + "].vertices[" +
                            Utilities::int_to_string(v) + "]=" +
                            Utilities::int_to_string(cells[c].vertices[v]) +
                            " is invalid, because only " +
                            Utilities::int_to_string(vertices.size()) +
                            " vertices were supplied."));
          vertex_used[cells[c].vertices[v]] = true;
        }


    // then renumber the vertices that are actually used in the same order as
    // they were beforehand
    const unsigned int        invalid_vertex = numbers::invalid_unsigned_int;
    std::vector<unsigned int> new_vertex_numbers(vertices.size(),
                                                 invalid_vertex);
    unsigned int              next_free_number = 0;
    for (unsigned int i = 0; i < vertices.size(); ++i)
      if (vertex_used[i] == true)
        {
          new_vertex_numbers[i] = next_free_number;
          ++next_free_number;
        }

    // next replace old vertex numbers by the new ones
    for (unsigned int c = 0; c < cells.size(); ++c)
      for (auto &v : cells[c].vertices)
        v = new_vertex_numbers[v];

    // same for boundary data
    for (unsigned int c = 0; c < subcelldata.boundary_lines.size(); // NOLINT
         ++c)
      for (unsigned int v = 0;
           v < subcelldata.boundary_lines[c].vertices.size();
           ++v)
        {
          Assert(subcelldata.boundary_lines[c].vertices[v] <
                   new_vertex_numbers.size(),
                 ExcMessage(
                   "Invalid vertex index in subcelldata.boundary_lines. "
                   "subcelldata.boundary_lines[" +
                   Utilities::int_to_string(c) + "].vertices[" +
                   Utilities::int_to_string(v) + "]=" +
                   Utilities::int_to_string(
                     subcelldata.boundary_lines[c].vertices[v]) +
                   " is invalid, because only " +
                   Utilities::int_to_string(vertices.size()) +
                   " vertices were supplied."));
          subcelldata.boundary_lines[c].vertices[v] =
            new_vertex_numbers[subcelldata.boundary_lines[c].vertices[v]];
        }

    for (unsigned int c = 0; c < subcelldata.boundary_quads.size(); // NOLINT
         ++c)
      for (unsigned int v = 0;
           v < subcelldata.boundary_quads[c].vertices.size();
           ++v)
        {
          Assert(subcelldata.boundary_quads[c].vertices[v] <
                   new_vertex_numbers.size(),
                 ExcMessage(
                   "Invalid vertex index in subcelldata.boundary_quads. "
                   "subcelldata.boundary_quads[" +
                   Utilities::int_to_string(c) + "].vertices[" +
                   Utilities::int_to_string(v) + "]=" +
                   Utilities::int_to_string(
                     subcelldata.boundary_quads[c].vertices[v]) +
                   " is invalid, because only " +
                   Utilities::int_to_string(vertices.size()) +
                   " vertices were supplied."));

          subcelldata.boundary_quads[c].vertices[v] =
            new_vertex_numbers[subcelldata.boundary_quads[c].vertices[v]];
        }

    // finally copy over the vertices which we really need to a new array and
    // replace the old one by the new one
    std::vector<Point<spacedim>> tmp;
    tmp.reserve(std::count(vertex_used.begin(), vertex_used.end(), true));
    for (unsigned int v = 0; v < vertices.size(); ++v)
      if (vertex_used[v] == true)
        tmp.push_back(vertices[v]);
    swap(vertices, tmp);
  }



  template <int dim, int spacedim>
  void
  delete_duplicated_vertices(std::vector<Point<spacedim>> &vertices,
                             std::vector<CellData<dim>>   &cells,
                             SubCellData                  &subcelldata,
                             std::vector<unsigned int>    &considered_vertices,
                             const double                  tol)
  {
    if (tol == 0.0)
      return; // nothing to do per definition

    AssertIndexRange(2, vertices.size());
    std::vector<unsigned int> new_vertex_numbers(vertices.size());
    std::iota(new_vertex_numbers.begin(), new_vertex_numbers.end(), 0);

    // if the considered_vertices vector is empty, consider all vertices
    if (considered_vertices.empty())
      considered_vertices = new_vertex_numbers;
    Assert(considered_vertices.size() <= vertices.size(), ExcInternalError());

    // The algorithm below improves upon the naive O(n^2) algorithm by first
    // sorting vertices by their value in one component and then only
    // comparing vertices for equality which are nearly equal in that
    // component. For example, if @p vertices form a cube, then we will only
    // compare points that have the same x coordinate when we try to find
    // duplicated vertices.

    // Start by finding the longest coordinate direction. This minimizes the
    // number of points that need to be compared against each-other in a
    // single set for typical geometries.
    const BoundingBox<spacedim> bbox(vertices);

    unsigned int longest_coordinate_direction = 0;
    double       longest_coordinate_length    = bbox.side_length(0);
    for (unsigned int d = 1; d < spacedim; ++d)
      {
        const double coordinate_length = bbox.side_length(d);
        if (longest_coordinate_length < coordinate_length)
          {
            longest_coordinate_length    = coordinate_length;
            longest_coordinate_direction = d;
          }
      }

    // Sort vertices (while preserving their vertex numbers) along that
    // coordinate direction:
    std::vector<std::pair<unsigned int, Point<spacedim>>> sorted_vertices;
    sorted_vertices.reserve(vertices.size());
    for (const unsigned int vertex_n : considered_vertices)
      {
        AssertIndexRange(vertex_n, vertices.size());
        sorted_vertices.emplace_back(vertex_n, vertices[vertex_n]);
      }
    std::sort(sorted_vertices.begin(),
              sorted_vertices.end(),
              [&](const std::pair<unsigned int, Point<spacedim>> &a,
                  const std::pair<unsigned int, Point<spacedim>> &b) {
                return a.second[longest_coordinate_direction] <
                       b.second[longest_coordinate_direction];
              });

    auto within_tolerance = [=](const Point<spacedim> &a,
                                const Point<spacedim> &b) {
      for (unsigned int d = 0; d < spacedim; ++d)
        if (std::abs(a[d] - b[d]) > tol)
          return false;
      return true;
    };

    // Find a range of numbers that have the same component in the longest
    // coordinate direction:
    auto range_start = sorted_vertices.begin();
    while (range_start != sorted_vertices.end())
      {
        auto range_end = range_start + 1;
        while (range_end != sorted_vertices.end() &&
               std::abs(range_end->second[longest_coordinate_direction] -
                        range_start->second[longest_coordinate_direction]) <
                 tol)
          ++range_end;

        // preserve behavior with older versions of this function by replacing
        // higher vertex numbers by lower vertex numbers
        std::sort(range_start,
                  range_end,
                  [](const std::pair<unsigned int, Point<spacedim>> &a,
                     const std::pair<unsigned int, Point<spacedim>> &b) {
                    return a.first < b.first;
                  });

        // Now de-duplicate [range_start, range_end)
        //
        // We have identified all points that are within a strip of width 'tol'
        // in one coordinate direction. Now we need to figure out which of these
        // are also close in other coordinate directions. If two are close, we
        // can mark the second one for deletion.
        for (auto reference = range_start; reference != range_end; ++reference)
          {
            if (reference->first != numbers::invalid_unsigned_int)
              for (auto it = reference + 1; it != range_end; ++it)
                {
                  if (within_tolerance(reference->second, it->second))
                    {
                      new_vertex_numbers[it->first] = reference->first;
                      // skip the replaced vertex in the future
                      it->first = numbers::invalid_unsigned_int;
                    }
                }
          }
        range_start = range_end;
      }

    // now we got a renumbering list. simply renumber all vertices
    // (non-duplicate vertices get renumbered to themselves, so nothing bad
    // happens). after that, the duplicate vertices will be unused, so call
    // delete_unused_vertices() to do that part of the job.
    for (auto &cell : cells)
      for (auto &vertex_index : cell.vertices)
        vertex_index = new_vertex_numbers[vertex_index];
    for (auto &quad : subcelldata.boundary_quads)
      for (auto &vertex_index : quad.vertices)
        vertex_index = new_vertex_numbers[vertex_index];
    for (auto &line : subcelldata.boundary_lines)
      for (auto &vertex_index : line.vertices)
        vertex_index = new_vertex_numbers[vertex_index];

    delete_unused_vertices(vertices, cells, subcelldata);
  }



  template <int dim>
  void
  delete_duplicated_vertices(std::vector<Point<dim>> &vertices,
                             const double             tol)
  {
    if (vertices.empty())
      return;

    // 1) map point to local vertex index
    std::map<Point<dim>, unsigned int, FloatingPointComparator<double>>
      map_point_to_local_vertex_index{FloatingPointComparator<double>(tol)};

    // 2) initialize map with existing points uniquely
    for (unsigned int i = 0; i < vertices.size(); ++i)
      map_point_to_local_vertex_index[vertices[i]] = i;

    // no duplicate points are found
    if (map_point_to_local_vertex_index.size() == vertices.size())
      return;

    // 3) remove duplicate entries from vertices
    vertices.resize(map_point_to_local_vertex_index.size());
    {
      unsigned int j = 0;
      for (const auto &p : map_point_to_local_vertex_index)
        vertices[j++] = p.first;
    }
  }



  template <int dim, int spacedim>
  std::size_t
  invert_cells_with_negative_measure(
    const std::vector<Point<spacedim>> &all_vertices,
    std::vector<CellData<dim>>         &cells)
  {
    // This function is presently only implemented for volumetric (codimension
    // 0) elements.

    if (dim == 1)
      return 0;
    if (dim == 2 && spacedim == 3)
      DEAL_II_NOT_IMPLEMENTED();

    std::size_t n_negative_cells = 0;
    std::size_t cell_no          = 0;
    for (auto &cell : cells)
      {
        const ArrayView<const unsigned int> vertices(cell.vertices);
        // Some pathologically twisted cells can have exactly zero measure but
        // we can still fix them
        if (GridTools::cell_measure(all_vertices, vertices) <= 0)
          {
            ++n_negative_cells;
            const auto reference_cell =
              ReferenceCell::n_vertices_to_type(dim, vertices.size());

            if (reference_cell.is_hyper_cube())
              {
                if (dim == 2)
                  {
                    // flip the cell across the y = x line in 2d
                    std::swap(cell.vertices[1], cell.vertices[2]);
                  }
                else if (dim == 3)
                  {
                    // swap the front and back faces in 3d
                    std::swap(cell.vertices[0], cell.vertices[2]);
                    std::swap(cell.vertices[1], cell.vertices[3]);
                    std::swap(cell.vertices[4], cell.vertices[6]);
                    std::swap(cell.vertices[5], cell.vertices[7]);
                  }
              }
            else if (reference_cell.is_simplex())
              {
                // By basic rules for computing determinants we can just swap
                // two vertices to fix a negative volume. Arbitrarily pick the
                // last two.
                std::swap(cell.vertices[cell.vertices.size() - 2],
                          cell.vertices[cell.vertices.size() - 1]);
              }
            else if (reference_cell == ReferenceCells::Wedge)
              {
                // swap the two triangular faces
                std::swap(cell.vertices[0], cell.vertices[3]);
                std::swap(cell.vertices[1], cell.vertices[4]);
                std::swap(cell.vertices[2], cell.vertices[5]);
              }
            else if (reference_cell == ReferenceCells::Pyramid)
              {
                // Try swapping two vertices in the base - perhaps things were
                // read in the UCD (counter-clockwise) order instead of lexical
                std::swap(cell.vertices[2], cell.vertices[3]);
              }
            else
              {
                AssertThrow(false, ExcNotImplemented());
              }
            // Check whether the resulting cell is now ok.
            // If not, then the grid is seriously broken and
            // we just give up.
            AssertThrow(GridTools::cell_measure(all_vertices, vertices) > 0,
                        ExcGridHasInvalidCell(cell_no));
          }
        ++cell_no;
      }
    return n_negative_cells;
  }


  template <int dim, int spacedim>
  void
  invert_all_negative_measure_cells(
    const std::vector<Point<spacedim>> &all_vertices,
    std::vector<CellData<dim>>         &cells)
  {
    const std::size_t n_negative_cells =
      invert_cells_with_negative_measure(all_vertices, cells);

    // We assume that all cells of a grid have
    // either positive or negative volumes but
    // not both mixed. Although above reordering
    // might work also on single cells, grids
    // with both kind of cells are very likely to
    // be broken. Check for this here.
    AssertThrow(n_negative_cells == 0 || n_negative_cells == cells.size(),
                ExcMessage(
                  std::string(
                    "This function assumes that either all cells have positive "
                    "volume, or that all cells have been specified in an "
                    "inverted vertex order so that their volume is negative. "
                    "(In the latter case, this class automatically inverts "
                    "every cell.) However, the mesh you have specified "
                    "appears to have both cells with positive and cells with "
                    "negative volume. You need to check your mesh which "
                    "cells these are and how they got there.\n"
                    "As a hint, of the total ") +
                  std::to_string(cells.size()) + " cells in the mesh, " +
                  std::to_string(n_negative_cells) +
                  " appear to have a negative volume."));
  }



  // Functions and classes for consistently_order_cells
  namespace
  {
    /**
     * A simple data structure denoting an edge, i.e., the ordered pair
     * of its vertex indices. This is only used in the is_consistent()
     * function.
     */
    struct CheapEdge
    {
      /**
       * Construct an edge from the global indices of its two vertices.
       */
      CheapEdge(const unsigned int v0, const unsigned int v1)
        : v0(v0)
        , v1(v1)
      {}

      /**
       * Comparison operator for edges. It compares based on the
       * lexicographic ordering of the two vertex indices.
       */
      bool
      operator<(const CheapEdge &e) const
      {
        return ((v0 < e.v0) || ((v0 == e.v0) && (v1 < e.v1)));
      }

    private:
      /**
       * The global indices of the vertices that define the edge.
       */
      const unsigned int v0, v1;
    };


    /**
     * A function that determines whether the edges in a mesh are
     * already consistently oriented. It does so by adding all edges
     * of all cells into a set (which automatically eliminates
     * duplicates) but before that checks whether the reverse edge is
     * already in the set -- which would imply that a neighboring cell
     * is inconsistently oriented.
     */
    template <int dim>
    bool
    is_consistent(const std::vector<CellData<dim>> &cells)
    {
      std::set<CheapEdge> edges;

      for (typename std::vector<CellData<dim>>::const_iterator c =
             cells.begin();
           c != cells.end();
           ++c)
        {
          // construct the edges in reverse order. for each of them,
          // ensure that the reverse edge is not yet in the list of
          // edges (return false if the reverse edge already *is* in
          // the list) and then add the actual edge to it; std::set
          // eliminates duplicates automatically
          for (unsigned int l = 0; l < GeometryInfo<dim>::lines_per_cell; ++l)
            {
              const CheapEdge reverse_edge(
                c->vertices[GeometryInfo<dim>::line_to_cell_vertices(l, 1)],
                c->vertices[GeometryInfo<dim>::line_to_cell_vertices(l, 0)]);
              if (edges.find(reverse_edge) != edges.end())
                return false;


              // ok, not. insert edge in correct order
              const CheapEdge correct_edge(
                c->vertices[GeometryInfo<dim>::line_to_cell_vertices(l, 0)],
                c->vertices[GeometryInfo<dim>::line_to_cell_vertices(l, 1)]);
              edges.insert(correct_edge);
            }
        }

      // no conflicts found, so return true
      return true;
    }


    /**
     * A structure that describes some properties of parallel edges
     * such as what starter edges are (i.e., representative elements
     * of the sets of parallel edges within a cell) and what the set
     * of parallel edges to each edge is.
     */
    template <int dim>
    struct ParallelEdges
    {
      /**
       * An array that contains the indices of dim edges that can
       * serve as (arbitrarily chosen) starting points for the
       * dim sets of parallel edges within each cell.
       */
      static const unsigned int starter_edges[dim];

      /**
       * Number and indices of all of those edges parallel to each of the
       * edges in a cell.
       */
      static const unsigned int n_other_parallel_edges = (1 << (dim - 1)) - 1;
      static const unsigned int
        parallel_edges[GeometryInfo<dim>::lines_per_cell]
                      [n_other_parallel_edges];
    };

    template <>
    const unsigned int ParallelEdges<2>::starter_edges[2] = {0, 2};

    template <>
    const unsigned int ParallelEdges<2>::parallel_edges[4][1] = {{1},
                                                                 {0},
                                                                 {3},
                                                                 {2}};

    template <>
    const unsigned int ParallelEdges<3>::starter_edges[3] = {0, 2, 8};

    template <>
    const unsigned int ParallelEdges<3>::parallel_edges[12][3] = {
      {1, 4, 5},   // line 0
      {0, 4, 5},   // line 1
      {3, 6, 7},   // line 2
      {2, 6, 7},   // line 3
      {0, 1, 5},   // line 4
      {0, 1, 4},   // line 5
      {2, 3, 7},   // line 6
      {2, 3, 6},   // line 7
      {9, 10, 11}, // line 8
      {8, 10, 11}, // line 9
      {8, 9, 11},  // line 10
      {8, 9, 10}   // line 11
    };


    /**
     * A structure that store the index of a cell and, crucially, how a
     * given edge relates to this cell.
     */
    struct AdjacentCell
    {
      /**
       * Default constructor. Initialize the fields with invalid values.
       */
      AdjacentCell()
        : cell_index(numbers::invalid_unsigned_int)
        , edge_within_cell(numbers::invalid_unsigned_int)
      {}

      /**
       * Constructor. Initialize the fields with the given values.
       */
      AdjacentCell(const unsigned int cell_index,
                   const unsigned int edge_within_cell)
        : cell_index(cell_index)
        , edge_within_cell(edge_within_cell)
      {}


      unsigned int cell_index;
      unsigned int edge_within_cell;
    };



    template <int dim>
    class AdjacentCells;

    /**
     * A class that represents all of the cells adjacent to a given edge.
     * This class corresponds to the 2d case where each edge has at most
     * two adjacent cells.
     */
    template <>
    class AdjacentCells<2>
    {
    public:
      /**
       * An iterator that allows iterating over all cells adjacent
       * to the edge represented by the current object.
       */
      using const_iterator = const AdjacentCell *;

      /**
       * Add the given cell to the collection of cells adjacent to
       * the edge this object corresponds to. Since we are covering
       * the 2d case, the set of adjacent cells currently
       * represented by this object must have either zero or
       * one element already, since we can not add more than two
       * adjacent cells for each edge.
       */
      void
      push_back(const AdjacentCell &adjacent_cell)
      {
        if (adjacent_cells[0].cell_index == numbers::invalid_unsigned_int)
          adjacent_cells[0] = adjacent_cell;
        else
          {
            Assert(adjacent_cells[1].cell_index ==
                     numbers::invalid_unsigned_int,
                   ExcInternalError());
            adjacent_cells[1] = adjacent_cell;
          }
      }


      /**
       * Return an iterator to the first valid cell stored as adjacent to the
       * edge represented by the current object.
       */
      const_iterator
      begin() const
      {
        return adjacent_cells;
      }


      /**
       * Return an iterator to the element past the last valid cell stored
       * as adjacent to the edge represented by the current object.
       * @return
       */
      const_iterator
      end() const
      {
        // check whether the current object stores zero, one, or two
        // adjacent cells, and use this to point to the element past the
        // last valid one
        if (adjacent_cells[0].cell_index == numbers::invalid_unsigned_int)
          return adjacent_cells;
        else if (adjacent_cells[1].cell_index == numbers::invalid_unsigned_int)
          return adjacent_cells + 1;
        else
          return adjacent_cells + 2;
      }

    private:
      /**
       * References to the (at most) two cells that are adjacent to
       * the edge this object corresponds to. Unused elements are
       * default-initialized and have invalid values; in particular,
       * their cell_index field equals numbers::invalid_unsigned_int.
       */
      AdjacentCell adjacent_cells[2];
    };



    /**
     * A class that represents all of the cells adjacent to a given edge.
     * This class corresponds to the 3d case where each edge can have an
     * arbitrary number of adjacent cells. We represent this as a
     * std::vector<AdjacentCell>, from which class the current one is
     * derived and from which it inherits all of its member functions.
     */
    template <>
    class AdjacentCells<3> : public std::vector<AdjacentCell>
    {};


    /**
     * A class that describes all of the relevant properties of an
     * edge. For the purpose of what we do here, that includes the
     * indices of the two vertices, and the indices of the adjacent
     * cells (together with a description *where* in each of the
     * adjacent cells the edge is located). It also includes the
     * (global) direction of the edge: either from the first vertex to
     * the second, the other way around, or so far undetermined.
     */
    template <int dim>
    class Edge
    {
    public:
      /**
       * Constructor. Create the edge based on the information given
       * in @p cell, and selecting the edge with number @p edge_number
       * within this cell. Initialize the edge as unoriented.
       */
      Edge(const CellData<dim> &cell, const unsigned int edge_number)
        : orientation_status(not_oriented)
      {
        Assert(edge_number < GeometryInfo<dim>::lines_per_cell,
               ExcInternalError());

        // copy vertices for this particular line
        vertex_indices[0] =
          cell
            .vertices[GeometryInfo<dim>::line_to_cell_vertices(edge_number, 0)];
        vertex_indices[1] =
          cell
            .vertices[GeometryInfo<dim>::line_to_cell_vertices(edge_number, 1)];

        // bring them into standard orientation
        if (vertex_indices[0] > vertex_indices[1])
          std::swap(vertex_indices[0], vertex_indices[1]);
      }

      /**
       * Comparison operator for edges. It compares based on the
       * lexicographic ordering of the two vertex indices.
       */
      bool
      operator<(const Edge<dim> &e) const
      {
        return ((vertex_indices[0] < e.vertex_indices[0]) ||
                ((vertex_indices[0] == e.vertex_indices[0]) &&
                 (vertex_indices[1] < e.vertex_indices[1])));
      }

      /**
       * Compare two edges for equality based on their vertex indices.
       */
      bool
      operator==(const Edge<dim> &e) const
      {
        return ((vertex_indices[0] == e.vertex_indices[0]) &&
                (vertex_indices[1] == e.vertex_indices[1]));
      }

      /**
       * The global indices of the two vertices that bound this edge. These
       * will be ordered so that the first index is less than the second.
       */
      unsigned int vertex_indices[2];

      /**
       * An enum that indicates the direction of this edge with
       * regard to the two vertices that bound it.
       */
      enum OrientationStatus
      {
        not_oriented,
        forward,
        backward
      };

      OrientationStatus orientation_status;

      /**
       * Store the set of cells adjacent to this edge (these cells then
       * also store *where* in the cell the edge is located).
       */
      AdjacentCells<dim> adjacent_cells;
    };



    /**
     * A data structure that represents a cell with all of its vertices
     * and edges.
     */
    template <int dim>
    struct Cell
    {
      /**
       * Construct a Cell object from a CellData object. Also take a
       * (sorted) list of edges and to point the edges of the current
       * object into this list of edges.
       */
      Cell(const CellData<dim> &c, const std::vector<Edge<dim>> &edge_list)
      {
        for (const unsigned int i : GeometryInfo<dim>::vertex_indices())
          vertex_indices[i] = c.vertices[i];

        // now for each of the edges of this cell, find the location inside the
        // given edge_list array and store than index
        for (unsigned int l = 0; l < GeometryInfo<dim>::lines_per_cell; ++l)
          {
            const Edge<dim> e(c, l);
            edge_indices[l] =
              (std::lower_bound(edge_list.begin(), edge_list.end(), e) -
               edge_list.begin());
            Assert(edge_indices[l] < edge_list.size(), ExcInternalError());
            Assert(edge_list[edge_indices[l]] == e, ExcInternalError());
          }
      }

      /**
       * A list of global indices for the vertices that bound this cell.
       */
      unsigned int vertex_indices[GeometryInfo<dim>::vertices_per_cell];

      /**
       * A list of indices into the 'edge_list' array passed to the constructor
       * for the edges of the current cell.
       */
      unsigned int edge_indices[GeometryInfo<dim>::lines_per_cell];
    };



    template <int dim>
    class EdgeDeltaSet;

    /**
     * A class that represents by how much the set of parallel edges
     * grows in each step. In the graph orientation paper, this set is
     * called $\Delta_k$, thus the name.
     *
     * In 2d, this set can only include zero, one, or two elements.
     * Consequently, the appropriate data structure is one in which
     * we store at most 2 elements in a fixed sized data structure.
     */
    template <>
    class EdgeDeltaSet<2>
    {
    public:
      /**
       * Iterator type for the elements of the set.
       */
      using const_iterator = const unsigned int *;

      /**
       * Default constructor. Initialize both slots as unused, corresponding
       * to an empty set.
       */
      EdgeDeltaSet()
      {
        edge_indices[0] = edge_indices[1] = numbers::invalid_unsigned_int;
      }


      /**
       * Delete the elements of the set by marking both slots as unused.
       */
      void
      clear()
      {
        edge_indices[0] = edge_indices[1] = numbers::invalid_unsigned_int;
      }

      /**
       * Insert one element into the set. This will fail if the set already
       * has two elements.
       */
      void
      insert(const unsigned int edge_index)
      {
        if (edge_indices[0] == numbers::invalid_unsigned_int)
          edge_indices[0] = edge_index;
        else
          {
            Assert(edge_indices[1] == numbers::invalid_unsigned_int,
                   ExcInternalError());
            edge_indices[1] = edge_index;
          }
      }


      /**
       * Return an iterator pointing to the first element of the set.
       */
      const_iterator
      begin() const
      {
        return edge_indices;
      }


      /**
       * Return an iterator pointing to the element past the last used one.
       */
      const_iterator
      end() const
      {
        // check whether the current object stores zero, one, or two
        // indices, and use this to point to the element past the
        // last valid one
        if (edge_indices[0] == numbers::invalid_unsigned_int)
          return edge_indices;
        else if (edge_indices[1] == numbers::invalid_unsigned_int)
          return edge_indices + 1;
        else
          return edge_indices + 2;
      }

    private:
      /**
       * Storage space to store the indices of at most two edges.
       */
      unsigned int edge_indices[2];
    };



    /**
     * A class that represents by how much the set of parallel edges
     * grows in each step. In the graph orientation paper, this set is
     * called $\Delta_k$, thus the name.
     *
     * In 3d, this set can have arbitrarily many elements, unlike the
     * 2d case specialized above. Consequently, we simply represent
     * the data structure with a std::set. Class derivation ensures
     * that we simply inherit all of the member functions of the
     * base class.
     */
    template <>
    class EdgeDeltaSet<3> : public std::set<unsigned int>
    {};



    /**
     * From a list of cells, build a sorted vector that contains all of the
     * edges that exist in the mesh.
     */
    template <int dim>
    std::vector<Edge<dim>>
    build_edges(const std::vector<CellData<dim>> &cells)
    {
      // build the edge list for all cells. because each cell has
      // GeometryInfo<dim>::lines_per_cell edges, the total number
      // of edges is this many times the number of cells. of course
      // some of them will be duplicates, and we throw them out below
      std::vector<Edge<dim>> edge_list;
      edge_list.reserve(cells.size() * GeometryInfo<dim>::lines_per_cell);
      for (unsigned int i = 0; i < cells.size(); ++i)
        for (unsigned int l = 0; l < GeometryInfo<dim>::lines_per_cell; ++l)
          edge_list.emplace_back(cells[i], l);

      // next sort the edge list and then remove duplicates
      std::sort(edge_list.begin(), edge_list.end());
      edge_list.erase(std::unique(edge_list.begin(), edge_list.end()),
                      edge_list.end());

      return edge_list;
    }



    /**
     * Build the cell list. Update the edge array to let edges know
     * which cells are adjacent to them.
     */
    template <int dim>
    std::vector<Cell<dim>>
    build_cells_and_connect_edges(const std::vector<CellData<dim>> &cells,
                                  std::vector<Edge<dim>>           &edges)
    {
      std::vector<Cell<dim>> cell_list;
      cell_list.reserve(cells.size());
      for (unsigned int i = 0; i < cells.size(); ++i)
        {
          // create our own data structure for the cells and let it
          // connect to the edges array
          cell_list.emplace_back(cells[i], edges);

          // then also inform the edges that they are adjacent
          // to the current cell, and where within this cell
          for (unsigned int l = 0; l < GeometryInfo<dim>::lines_per_cell; ++l)
            edges[cell_list.back().edge_indices[l]].adjacent_cells.push_back(
              AdjacentCell(i, l));
        }
      Assert(cell_list.size() == cells.size(), ExcInternalError());

      return cell_list;
    }



    /**
     * Return the index within 'cells' of the first cell that has at least one
     * edge that is not yet oriented.
     */
    template <int dim>
    unsigned int
    get_next_unoriented_cell(const std::vector<Cell<dim>> &cells,
                             const std::vector<Edge<dim>> &edges,
                             const unsigned int            current_cell)
    {
      for (unsigned int c = current_cell; c < cells.size(); ++c)
        for (unsigned int l = 0; l < GeometryInfo<dim>::lines_per_cell; ++l)
          if (edges[cells[c].edge_indices[l]].orientation_status ==
              Edge<dim>::not_oriented)
            return c;

      return numbers::invalid_unsigned_int;
    }



    /**
     * Given a set of cells and edges, orient all edges that are
     * (global) parallel to the one identified by the @p cell and
     * within it the one with index @p local_edge.
     */
    template <int dim>
    void
    orient_one_set_of_parallel_edges(const std::vector<Cell<dim>> &cells,
                                     std::vector<Edge<dim>>       &edges,
                                     const unsigned int            cell,
                                     const unsigned int            local_edge)
    {
      // choose the direction of the first edge. we have free choice
      // here and could simply choose "forward" if that's what pleases
      // us. however, for backward compatibility with the previous
      // implementation used till 2016, let us just choose the
      // direction so that it matches what we have in the given cell.
      //
      // in fact, in what can only be assumed to be a bug in the
      // original implementation, after orienting all edges, the code
      // that rotates the cells so that they match edge orientations
      // (see the rotate_cell() function below) rotated the cell two
      // more times by 90 degrees. this is ok -- it simply flips all
      // edge orientations, which leaves them valid. rather than do
      // the same in the current implementation, we can achieve the
      // same effect by modifying the rule above to choose the
      // direction of the starting edge of this parallel set
      // *opposite* to what it looks like in the current cell
      //
      // this bug only existed in the 2d implementation since there
      // were different implementations for 2d and 3d. consequently,
      // only replicate it for the 2d case and be "intuitive" in 3d.
      if (edges[cells[cell].edge_indices[local_edge]].vertex_indices[0] ==
          cells[cell].vertex_indices[GeometryInfo<dim>::line_to_cell_vertices(
            local_edge, 0)])
        // orient initial edge *opposite* to the way it is in the cell
        // (see above for the reason)
        edges[cells[cell].edge_indices[local_edge]].orientation_status =
          (dim == 2 ? Edge<dim>::backward : Edge<dim>::forward);
      else
        {
          Assert(
            edges[cells[cell].edge_indices[local_edge]].vertex_indices[0] ==
              cells[cell].vertex_indices
                [GeometryInfo<dim>::line_to_cell_vertices(local_edge, 1)],
            ExcInternalError());
          Assert(
            edges[cells[cell].edge_indices[local_edge]].vertex_indices[1] ==
              cells[cell].vertex_indices
                [GeometryInfo<dim>::line_to_cell_vertices(local_edge, 0)],
            ExcInternalError());

          // orient initial edge *opposite* to the way it is in the cell
          // (see above for the reason)
          edges[cells[cell].edge_indices[local_edge]].orientation_status =
            (dim == 2 ? Edge<dim>::forward : Edge<dim>::backward);
        }

      // walk outward from the given edge as described in
      // the algorithm in the paper that documents all of
      // this
      //
      // note that in 2d, each of the Deltas can at most
      // contain two elements, whereas in 3d it can be arbitrarily many
      EdgeDeltaSet<dim> Delta_k;
      EdgeDeltaSet<dim> Delta_k_minus_1;
      Delta_k_minus_1.insert(cells[cell].edge_indices[local_edge]);

      while (Delta_k_minus_1.begin() !=
             Delta_k_minus_1.end()) // while set is not empty
        {
          Delta_k.clear();

          for (typename EdgeDeltaSet<dim>::const_iterator delta =
                 Delta_k_minus_1.begin();
               delta != Delta_k_minus_1.end();
               ++delta)
            {
              Assert(edges[*delta].orientation_status !=
                       Edge<dim>::not_oriented,
                     ExcInternalError());

              // now go through the cells adjacent to this edge
              for (typename AdjacentCells<dim>::const_iterator adjacent_cell =
                     edges[*delta].adjacent_cells.begin();
                   adjacent_cell != edges[*delta].adjacent_cells.end();
                   ++adjacent_cell)
                {
                  const unsigned int K = adjacent_cell->cell_index;
                  const unsigned int delta_is_edge_in_K =
                    adjacent_cell->edge_within_cell;

                  // figure out the direction of delta with respect to the cell
                  // K (in the orientation in which the user has given it to us)
                  const unsigned int first_edge_vertex =
                    (edges[*delta].orientation_status == Edge<dim>::forward ?
                       edges[*delta].vertex_indices[0] :
                       edges[*delta].vertex_indices[1]);
                  const unsigned int first_edge_vertex_in_K =
                    cells[K]
                      .vertex_indices[GeometryInfo<dim>::line_to_cell_vertices(
                        delta_is_edge_in_K, 0)];
                  Assert(
                    first_edge_vertex == first_edge_vertex_in_K ||
                      first_edge_vertex ==
                        cells[K].vertex_indices[GeometryInfo<
                          dim>::line_to_cell_vertices(delta_is_edge_in_K, 1)],
                    ExcInternalError());

                  // now figure out which direction the each of the "opposite"
                  // edges needs to be oriented into.
                  for (unsigned int o_e = 0;
                       o_e < ParallelEdges<dim>::n_other_parallel_edges;
                       ++o_e)
                    {
                      // get the index of the opposite edge and select which its
                      // first vertex needs to be based on how the current edge
                      // is oriented in the current cell
                      const unsigned int opposite_edge =
                        cells[K].edge_indices[ParallelEdges<
                          dim>::parallel_edges[delta_is_edge_in_K][o_e]];
                      const unsigned int first_opposite_edge_vertex =
                        cells[K].vertex_indices
                          [GeometryInfo<dim>::line_to_cell_vertices(
                            ParallelEdges<
                              dim>::parallel_edges[delta_is_edge_in_K][o_e],
                            (first_edge_vertex == first_edge_vertex_in_K ? 0 :
                                                                           1))];

                      // then determine the orientation of the edge based on
                      // whether the vertex we want to be the edge's first
                      // vertex is already the first vertex of the edge, or
                      // whether it points in the opposite direction
                      const typename Edge<dim>::OrientationStatus
                        opposite_edge_orientation =
                          (edges[opposite_edge].vertex_indices[0] ==
                               first_opposite_edge_vertex ?
                             Edge<dim>::forward :
                             Edge<dim>::backward);

                      // see if the opposite edge (there is only one in 2d) has
                      // already been oriented.
                      if (edges[opposite_edge].orientation_status ==
                          Edge<dim>::not_oriented)
                        {
                          // the opposite edge is not yet oriented. do orient it
                          // and add it to Delta_k
                          edges[opposite_edge].orientation_status =
                            opposite_edge_orientation;
                          Delta_k.insert(opposite_edge);
                        }
                      else
                        {
                          // this opposite edge has already been oriented. it
                          // should be consistent with the current one in 2d,
                          // while in 3d it may in fact be mis-oriented, and in
                          // that case the mesh will not be orientable. indicate
                          // this by throwing an exception that we can catch
                          // further up; this has the advantage that we can
                          // propagate through a couple of functions without
                          // having to do error checking and without modifying
                          // the 'cells' array that the user gave us
                          if (dim == 2)
                            {
                              Assert(edges[opposite_edge].orientation_status ==
                                       opposite_edge_orientation,
                                     ExcMeshNotOrientable());
                            }
                          else if (dim == 3)
                            {
                              if (edges[opposite_edge].orientation_status !=
                                  opposite_edge_orientation)
                                throw ExcMeshNotOrientable();
                            }
                          else
                            DEAL_II_NOT_IMPLEMENTED();
                        }
                    }
                }
            }

          // finally copy the new set to the previous one
          // (corresponding to increasing 'k' by one in the
          // algorithm)
          Delta_k_minus_1 = Delta_k;
        }
    }


    /**
     * Given data structures @p cell_list and @p edge_list, where
     * all edges are already oriented, rotate the cell with
     * index @p cell_index in such a way that its local coordinate
     * system matches the ones of the adjacent edges. Store the
     * rotated order of vertices in <code>raw_cells[cell_index]</code>.
     */
    template <int dim>
    void
    rotate_cell(const std::vector<Cell<dim>> &cell_list,
                const std::vector<Edge<dim>> &edge_list,
                const unsigned int            cell_index,
                std::vector<CellData<dim>>   &raw_cells)
    {
      // find the first vertex of the cell. this is the vertex where dim edges
      // originate, so for each of the edges record which the starting vertex is
      unsigned int starting_vertex_of_edge[GeometryInfo<dim>::lines_per_cell];
      for (unsigned int e = 0; e < GeometryInfo<dim>::lines_per_cell; ++e)
        {
          Assert(edge_list[cell_list[cell_index].edge_indices[e]]
                     .orientation_status != Edge<dim>::not_oriented,
                 ExcInternalError());
          if (edge_list[cell_list[cell_index].edge_indices[e]]
                .orientation_status == Edge<dim>::forward)
            starting_vertex_of_edge[e] =
              edge_list[cell_list[cell_index].edge_indices[e]]
                .vertex_indices[0];
          else
            starting_vertex_of_edge[e] =
              edge_list[cell_list[cell_index].edge_indices[e]]
                .vertex_indices[1];
        }

      // find the vertex number that appears dim times. this will then be
      // the vertex at which we want to locate the origin of the cell's
      // coordinate system (i.e., vertex 0)
      unsigned int origin_vertex_of_cell = numbers::invalid_unsigned_int;
      switch (dim)
        {
          case 2:
            {
              // in 2d, we can simply enumerate the possibilities where the
              // origin may be located because edges zero and one don't share
              // any vertices, and the same for edges two and three
              if ((starting_vertex_of_edge[0] == starting_vertex_of_edge[2]) ||
                  (starting_vertex_of_edge[0] == starting_vertex_of_edge[3]))
                origin_vertex_of_cell = starting_vertex_of_edge[0];
              else if ((starting_vertex_of_edge[1] ==
                        starting_vertex_of_edge[2]) ||
                       (starting_vertex_of_edge[1] ==
                        starting_vertex_of_edge[3]))
                origin_vertex_of_cell = starting_vertex_of_edge[1];
              else
                DEAL_II_ASSERT_UNREACHABLE();

              break;
            }

          case 3:
            {
              // one could probably do something similar in 3d, but that seems
              // more complicated than one wants to write down. just go
              // through the list of possible starting vertices and check
              for (origin_vertex_of_cell = 0;
                   origin_vertex_of_cell < GeometryInfo<dim>::vertices_per_cell;
                   ++origin_vertex_of_cell)
                if (std::count(starting_vertex_of_edge,
                               starting_vertex_of_edge +
                                 GeometryInfo<dim>::lines_per_cell,
                               cell_list[cell_index]
                                 .vertex_indices[origin_vertex_of_cell]) == dim)
                  break;
              Assert(origin_vertex_of_cell <
                       GeometryInfo<dim>::vertices_per_cell,
                     ExcInternalError());

              break;
            }

          default:
            DEAL_II_NOT_IMPLEMENTED();
        }

      // now rotate raw_cells[cell_index] in such a way that its orientation
      // matches that of cell_list[cell_index]
      switch (dim)
        {
          case 2:
            {
              // in 2d, we can literally rotate the cell until its origin
              // matches the one that we have determined above should be
              // the origin vertex
              //
              // when doing a rotation, take into account the ordering of
              // vertices (not in clockwise or counter-clockwise sense)
              while (raw_cells[cell_index].vertices[0] != origin_vertex_of_cell)
                {
                  const unsigned int tmp = raw_cells[cell_index].vertices[0];
                  raw_cells[cell_index].vertices[0] =
                    raw_cells[cell_index].vertices[1];
                  raw_cells[cell_index].vertices[1] =
                    raw_cells[cell_index].vertices[3];
                  raw_cells[cell_index].vertices[3] =
                    raw_cells[cell_index].vertices[2];
                  raw_cells[cell_index].vertices[2] = tmp;
                }
              break;
            }

          case 3:
            {
              // in 3d, the situation is a bit more complicated. from above, we
              // now know which vertex is at the origin (because 3 edges
              // originate from it), but that still leaves 3 possible rotations
              // of the cube. the important realization is that we can choose
              // any of them: in all 3 rotations, all edges originate from the
              // one vertex, and that fixes the directions of all 12 edges in
              // the cube because these 3 cover all 3 equivalence classes!
              // consequently, we can select an arbitrary one among the
              // permutations -- for example the following ones:
              static const unsigned int cube_permutations[8][8] = {
                {0, 1, 2, 3, 4, 5, 6, 7},
                {1, 5, 3, 7, 0, 4, 2, 6},
                {2, 6, 0, 4, 3, 7, 1, 5},
                {3, 2, 1, 0, 7, 6, 5, 4},
                {4, 0, 6, 2, 5, 1, 7, 3},
                {5, 4, 7, 6, 1, 0, 3, 2},
                {6, 7, 4, 5, 2, 3, 0, 1},
                {7, 3, 5, 1, 6, 2, 4, 0}};

              unsigned int
                temp_vertex_indices[GeometryInfo<dim>::vertices_per_cell];
              for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
                temp_vertex_indices[v] =
                  raw_cells[cell_index]
                    .vertices[cube_permutations[origin_vertex_of_cell][v]];
              for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
                raw_cells[cell_index].vertices[v] = temp_vertex_indices[v];

              break;
            }

          default:
            {
              DEAL_II_NOT_IMPLEMENTED();
            }
        }
    }


    /**
     * Given a set of cells, find globally unique edge orientations
     * and then rotate cells so that the coordinate system of the cell
     * coincides with the coordinate systems of the adjacent edges.
     */
    template <int dim>
    void
    reorient(std::vector<CellData<dim>> &cells)
    {
      // first build the arrays that connect cells to edges and the other
      // way around
      std::vector<Edge<dim>> edge_list = build_edges(cells);
      std::vector<Cell<dim>> cell_list =
        build_cells_and_connect_edges(cells, edge_list);

      // then loop over all cells and start orienting parallel edge sets
      // of cells that still have non-oriented edges
      unsigned int next_cell_with_unoriented_edge = 0;
      while ((next_cell_with_unoriented_edge = get_next_unoriented_cell(
                cell_list, edge_list, next_cell_with_unoriented_edge)) !=
             numbers::invalid_unsigned_int)
        {
          // see which edge sets are still not oriented
          //
          // we do not need to look at each edge because if we orient edge
          // 0, we will end up with edge 1 also oriented (in 2d; in 3d, there
          // will be 3 other edges that are also oriented). there are only
          // dim independent sets of edges, so loop over these.
          //
          // we need to check whether each one of these starter edges may
          // already be oriented because the line (sheet) that connects
          // globally parallel edges may be self-intersecting in the
          // current cell
          for (unsigned int l = 0; l < dim; ++l)
            if (edge_list[cell_list[next_cell_with_unoriented_edge]
                            .edge_indices[ParallelEdges<dim>::starter_edges[l]]]
                  .orientation_status == Edge<dim>::not_oriented)
              orient_one_set_of_parallel_edges(
                cell_list,
                edge_list,
                next_cell_with_unoriented_edge,
                ParallelEdges<dim>::starter_edges[l]);

          // ensure that we have really oriented all edges now, not just
          // the starter edges
          for (unsigned int l = 0; l < GeometryInfo<dim>::lines_per_cell; ++l)
            Assert(edge_list[cell_list[next_cell_with_unoriented_edge]
                               .edge_indices[l]]
                       .orientation_status != Edge<dim>::not_oriented,
                   ExcInternalError());
        }

      // now that we have oriented all edges, we need to rotate cells
      // so that the edges point in the right direction with the now
      // rotated coordinate system
      for (unsigned int c = 0; c < cells.size(); ++c)
        rotate_cell(cell_list, edge_list, c, cells);
    }


    // overload of the function above for 1d -- there is nothing
    // to orient in that case
    void
    reorient(std::vector<CellData<1>> &)
    {}
  } // namespace



  template <int dim>
  void
  consistently_order_cells(std::vector<CellData<dim>> &cells)
  {
    Assert(cells.size() != 0,
           ExcMessage(
             "List of elements to orient must have at least one cell"));

    // there is nothing for us to do in 1d
    if (dim == 1)
      return;

    // check if grids are already consistent. if so, do
    // nothing. if not, then do the reordering
    if (!is_consistent(cells))
      try
        {
          reorient(cells);
        }
      catch (const ExcMeshNotOrientable &)
        {
          // the mesh is not orientable. this is acceptable if we are in 3d,
          // as class Triangulation knows how to handle this, but it is
          // not in 2d; in that case, re-throw the exception
          if (dim < 3)
            throw;
        }
  }



  template <int dim, int spacedim>
  std::map<unsigned int, Point<spacedim>>
  get_all_vertices_at_boundary(const Triangulation<dim, spacedim> &tria)
  {
    std::map<unsigned int, Point<spacedim>> vertex_map;
    typename Triangulation<dim, spacedim>::active_cell_iterator
      cell = tria.begin_active(),
      endc = tria.end();
    for (; cell != endc; ++cell)
      {
        for (const unsigned int i : cell->face_indices())
          {
            const typename Triangulation<dim, spacedim>::face_iterator &face =
              cell->face(i);
            if (face->at_boundary())
              {
                for (unsigned j = 0; j < face->n_vertices(); ++j)
                  {
                    const Point<spacedim> &vertex       = face->vertex(j);
                    const unsigned int     vertex_index = face->vertex_index(j);
                    vertex_map[vertex_index]            = vertex;
                  }
              }
          }
      }
    return vertex_map;
  }



  template <int dim, int spacedim>
  std::vector<std::vector<std::pair<unsigned int, Point<spacedim>>>>
  extract_ordered_boundary_vertices(const Triangulation<dim, spacedim> &tria,
                                    const Mapping<dim, spacedim>       &mapping)
  {
    Assert(dim == 2, ExcMessage("Only implemented for 2D triangulations"));
    // This map holds the two vertex indices of each face.
    // Counterclockwise first vertex index on first position,
    // counterclockwise second vertex index on second position.
    std::map<unsigned int, unsigned int>    face_vertex_indices;
    std::map<unsigned int, Point<spacedim>> vertex_to_point;

    // Iterate over all active cells at the boundary
    for (const auto &cell : tria.active_cell_iterators())
      {
        for (const unsigned int f : cell->face_indices())
          {
            if (cell->face(f)->at_boundary())
              {
                // get mapped vertices of the cell
                const auto         v_mapped = mapping.get_vertices(cell);
                const unsigned int v0       = cell->face(f)->vertex_index(0);
                const unsigned int v1       = cell->face(f)->vertex_index(1);

                if (cell->reference_cell() == ReferenceCells::Triangle)
                  {
                    // add indices and first mapped vertex of the face
                    vertex_to_point[v0]     = v_mapped[f];
                    face_vertex_indices[v0] = v1;
                  }
                else if (cell->reference_cell() ==
                         ReferenceCells::Quadrilateral)
                  {
                    // Ensure that vertex indices of the face are in
                    // counterclockwise order inserted in the map.
                    if (f == 0 || f == 3)
                      {
                        // add indices and first mapped vertex of the face
                        vertex_to_point[v1] =
                          v_mapped[GeometryInfo<2>::face_to_cell_vertices(f,
                                                                          1)];
                        face_vertex_indices[v1] = v0;
                      }
                    else
                      {
                        // add indices and first mapped vertex of the face
                        vertex_to_point[v0] =
                          v_mapped[GeometryInfo<2>::face_to_cell_vertices(f,
                                                                          0)];
                        face_vertex_indices[v0] = v1;
                      }
                  }
                else
                  {
                    DEAL_II_ASSERT_UNREACHABLE();
                  }
              }
          }
      }

    std::vector<std::vector<std::pair<unsigned int, Point<spacedim>>>>
                                                          boundaries;
    std::vector<std::pair<unsigned int, Point<spacedim>>> current_boundary;

    // Vertex to start counterclockwise insertion
    unsigned int start_index   = face_vertex_indices.begin()->first;
    unsigned int current_index = start_index;

    // As long as still entries in the map, use last vertex index to
    // find next vertex index
    while (face_vertex_indices.size() > 0)
      {
        const auto vertex_it = vertex_to_point.find(current_index);
        Assert(vertex_it != vertex_to_point.end(),
               ExcMessage("This should not occur, please report bug"));
        current_boundary.emplace_back(vertex_it->first, vertex_it->second);
        vertex_to_point.erase(vertex_it);

        const auto it = face_vertex_indices.find(current_index);
        // If the boundary is one closed loop, the next vertex index
        // must exist as key until the map is empty.
        Assert(it != face_vertex_indices.end(),
               ExcMessage("Triangulation might contain holes"));

        current_index = it->second;
        face_vertex_indices.erase(it);

        // traversed one closed boundary loop
        if (current_index == start_index)
          {
            boundaries.push_back(current_boundary);
            current_boundary.clear();

            if (face_vertex_indices.size() == 0)
              {
                break;
              }

            // Take arbitrary remaining vertex as new start
            // for next boundary loop
            start_index   = face_vertex_indices.begin()->first;
            current_index = start_index;
          }
      }
    return boundaries;
  }



  template <int dim, int spacedim>
  void
  remove_hanging_nodes(Triangulation<dim, spacedim> &tria,
                       const bool                    isotropic,
                       const unsigned int            max_iterations)
  {
    unsigned int iter                = 0;
    bool         continue_refinement = true;

    while (continue_refinement && (iter < max_iterations))
      {
        if (max_iterations != numbers::invalid_unsigned_int)
          ++iter;
        continue_refinement = false;

        for (const auto &cell : tria.active_cell_iterators())
          for (const unsigned int j : cell->face_indices())
            if (cell->at_boundary(j) == false &&
                cell->neighbor(j)->has_children())
              {
                if (isotropic)
                  {
                    cell->set_refine_flag();
                    continue_refinement = true;
                  }
                else
                  continue_refinement |= cell->flag_for_face_refinement(j);
              }

        tria.execute_coarsening_and_refinement();
      }
  }



  template <int dim, int spacedim>
  void
  remove_anisotropy(Triangulation<dim, spacedim> &tria,
                    const double                  max_ratio,
                    const unsigned int            max_iterations)
  {
    unsigned int iter                = 0;
    bool         continue_refinement = true;

    while (continue_refinement && (iter < max_iterations))
      {
        ++iter;
        continue_refinement = false;
        for (const auto &cell : tria.active_cell_iterators())
          {
            std::pair<unsigned int, double> info =
              GridTools::get_longest_direction<dim, spacedim>(cell);
            if (info.second > max_ratio)
              {
                cell->set_refine_flag(
                  RefinementCase<dim>::cut_axis(info.first));
                continue_refinement = true;
              }
          }
        tria.execute_coarsening_and_refinement();
      }
  }



  template <int dim, int spacedim>
  std::map<unsigned int, Point<spacedim>>
  extract_used_vertices(const Triangulation<dim, spacedim> &container,
                        const Mapping<dim, spacedim>       &mapping)
  {
    std::map<unsigned int, Point<spacedim>> result;
    for (const auto &cell : container.active_cell_iterators())
      {
        if (!cell->is_artificial())
          {
            const auto vs = mapping.get_vertices(cell);
            for (unsigned int i = 0; i < vs.size(); ++i)
              result[cell->vertex_index(i)] = vs[i];
          }
      }
    return result;
  }



  template <int dim, int spacedim>
  std::vector<
    std::set<typename Triangulation<dim, spacedim>::active_cell_iterator>>
  vertex_to_cell_map(const Triangulation<dim, spacedim> &triangulation)
  {
    std::vector<
      std::set<typename Triangulation<dim, spacedim>::active_cell_iterator>>
      vertex_to_cell_map(triangulation.n_vertices());
    typename Triangulation<dim, spacedim>::active_cell_iterator
      cell = triangulation.begin_active(),
      endc = triangulation.end();
    for (; cell != endc; ++cell)
      for (const unsigned int i : cell->vertex_indices())
        vertex_to_cell_map[cell->vertex_index(i)].insert(cell);

    // Check if mesh has hanging nodes. Do this only locally to
    // prevent communication and possible deadlock.
    if (triangulation.Triangulation<dim, spacedim>::has_hanging_nodes())
      {
        Assert(triangulation.all_reference_cells_are_hyper_cube(),
               ExcNotImplemented());

        // Take care of hanging nodes
        cell = triangulation.begin_active();
        for (; cell != endc; ++cell)
          {
            for (const unsigned int i : cell->face_indices())
              {
                if ((cell->at_boundary(i) == false) &&
                    (cell->neighbor(i)->is_active()))
                  {
                    typename Triangulation<dim, spacedim>::active_cell_iterator
                      adjacent_cell = cell->neighbor(i);
                    for (unsigned int j = 0; j < cell->face(i)->n_vertices();
                         ++j)
                      vertex_to_cell_map[cell->face(i)->vertex_index(j)].insert(
                        adjacent_cell);
                  }
              }

            // in 3d also loop over the edges
            if (dim == 3)
              {
                for (unsigned int i = 0; i < cell->n_lines(); ++i)
                  if (cell->line(i)->has_children())
                    // the only place where this vertex could have been
                    // hiding is on the mid-edge point of the edge we
                    // are looking at
                    vertex_to_cell_map[cell->line(i)->child(0)->vertex_index(1)]
                      .insert(cell);
              }
          }
      }

    return vertex_to_cell_map;
  }



  template <int dim, int spacedim>
  void
  get_face_connectivity_of_cells(
    const Triangulation<dim, spacedim> &triangulation,
    DynamicSparsityPattern             &cell_connectivity)
  {
    cell_connectivity.reinit(triangulation.n_active_cells(),
                             triangulation.n_active_cells());

    // loop over all cells and their neighbors to build the sparsity
    // pattern. note that it's a bit hard to enter all the connections when a
    // neighbor has children since we would need to find out which of its
    // children is adjacent to the current cell. this problem can be omitted
    // if we only do something if the neighbor has no children -- in that case
    // it is either on the same or a coarser level than we are. in return, we
    // have to add entries in both directions for both cells
    for (const auto &cell : triangulation.active_cell_iterators())
      {
        const unsigned int index = cell->active_cell_index();
        cell_connectivity.add(index, index);
        for (auto f : cell->face_indices())
          if ((cell->at_boundary(f) == false) &&
              (cell->neighbor(f)->has_children() == false))
            {
              const unsigned int other_index =
                cell->neighbor(f)->active_cell_index();
              cell_connectivity.add(index, other_index);
              cell_connectivity.add(other_index, index);
            }
      }
  }



  template <int dim, int spacedim>
  void
  get_vertex_connectivity_of_cells(
    const Triangulation<dim, spacedim> &triangulation,
    DynamicSparsityPattern             &cell_connectivity)
  {
    // The choice of 16 or fewer neighbors here is based on empirical
    // measurements.
    //
    // Vertices in a structured hexahedral mesh have 8 adjacent cells. In a
    // structured tetrahedral mesh, about 98% of vertices have 16 neighbors or
    // fewer. Similarly, in an unstructured tetrahedral mesh, if we count the
    // number of neighbors each vertex has we obtain the following distribution:
    //
    // 3, 1
    // 4, 728
    // 5, 4084
    // 6, 7614
    // 7, 17329
    // 8, 31145
    // 9, 46698
    // 10, 64193
    // 11, 68269
    // 12, 63574
    // 13, 57016
    // 14, 50476
    // 15, 41886
    // 16, 31820
    // 17, 21269
    // 18, 12217
    // 19, 6072
    // 20, 2527
    // 21, 825
    // 22, 262
    // 23, 61
    // 24, 12
    // 26, 1
    //
    // so about 86% of vertices have 16 neighbors or fewer. Hence, we picked 16
    // neighbors here to cover most cases without allocation.
    std::vector<boost::container::small_vector<unsigned int, 16>>
      vertex_to_cell(triangulation.n_vertices());
    for (const auto &cell : triangulation.active_cell_iterators())
      {
        for (const unsigned int v : cell->vertex_indices())
          vertex_to_cell[cell->vertex_index(v)].push_back(
            cell->active_cell_index());
      }

    cell_connectivity.reinit(triangulation.n_active_cells(),
                             triangulation.n_active_cells());
    std::vector<types::global_dof_index> neighbors;
    for (const auto &cell : triangulation.active_cell_iterators())
      {
        neighbors.clear();
        for (const unsigned int v : cell->vertex_indices())
          neighbors.insert(neighbors.end(),
                           vertex_to_cell[cell->vertex_index(v)].begin(),
                           vertex_to_cell[cell->vertex_index(v)].end());
        std::sort(neighbors.begin(), neighbors.end());
        cell_connectivity.add_entries(cell->active_cell_index(),
                                      neighbors.begin(),
                                      std::unique(neighbors.begin(),
                                                  neighbors.end()),
                                      true);
      }
  }


  template <int dim, int spacedim>
  void
  get_vertex_connectivity_of_cells_on_level(
    const Triangulation<dim, spacedim> &triangulation,
    const unsigned int                  level,
    DynamicSparsityPattern             &cell_connectivity)
  {
    std::vector<boost::container::small_vector<unsigned int, 16>>
      vertex_to_cell(triangulation.n_vertices());
    for (typename Triangulation<dim, spacedim>::cell_iterator cell =
           triangulation.begin(level);
         cell != triangulation.end(level);
         ++cell)
      {
        for (const unsigned int v : cell->vertex_indices())
          vertex_to_cell[cell->vertex_index(v)].push_back(cell->index());
      }

    cell_connectivity.reinit(triangulation.n_cells(level),
                             triangulation.n_cells(level));
    std::vector<types::global_dof_index> neighbors;
    for (const auto &cell : triangulation.cell_iterators_on_level(level))
      {
        neighbors.clear();
        for (const unsigned int v : cell->vertex_indices())
          neighbors.insert(neighbors.end(),
                           vertex_to_cell[cell->vertex_index(v)].begin(),
                           vertex_to_cell[cell->vertex_index(v)].end());
        std::sort(neighbors.begin(), neighbors.end());
        cell_connectivity.add_entries(cell->index(),
                                      neighbors.begin(),
                                      std::unique(neighbors.begin(),
                                                  neighbors.end()),
                                      true);
      }
  }
} /* namespace GridTools */

// explicit instantiations
#include "grid/grid_tools_topology.inst"

DEAL_II_NAMESPACE_CLOSE
