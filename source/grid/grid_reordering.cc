// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2016 by the deal.II authors
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


#include <deal.II/grid/grid_reordering.h>
#include <deal.II/grid/grid_reordering_internal.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/std_cxx11/bind.h>

#include <algorithm>
#include <set>
#include <iostream>
#include <fstream>
#include <functional>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace GridReordering2d
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
      CheapEdge (const unsigned int v0,
                 const unsigned int v1)
        :
        v0(v0), v1(v1)
      {}

      /**
       * Comparison operator for edges. It compares based on the
       * lexicographic ordering of the two vertex indices.
       */
      bool operator < (const CheapEdge &e) const
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
    is_consistent  (const std::vector<CellData<dim> > &cells)
    {
      std::set<CheapEdge> edges;

      for (typename std::vector<CellData<dim> >::const_iterator c = cells.begin();
           c != cells.end(); ++c)
        {
          // construct the edges in reverse order. for each of them,
          // ensure that the reverse edge is not yet in the list of
          // edges (return false if the reverse edge already *is* in
          // the list) and then add the actual edge to it; std::set
          // eliminates duplicates automatically
          for (unsigned int l=0; l<GeometryInfo<dim>::lines_per_cell; ++l)
            {
              const CheapEdge reverse_edge (c->vertices[GeometryInfo<dim>::line_to_cell_vertices(l, 1)],
                                            c->vertices[GeometryInfo<dim>::line_to_cell_vertices(l, 0)]);
              if (edges.find (reverse_edge) != edges.end())
                return false;


              // ok, not. insert edge in correct order
              const CheapEdge correct_edge (c->vertices[GeometryInfo<dim>::line_to_cell_vertices(l, 0)],
                                            c->vertices[GeometryInfo<dim>::line_to_cell_vertices(l, 1)]);
              edges.insert (correct_edge);
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
      static const unsigned int starter_edges[dim];
      static const unsigned int parallel_edges[GeometryInfo<dim>::lines_per_cell][(1<<(dim-1)) - 1];
    };

    template <>
    const unsigned int ParallelEdges<2>::starter_edges[2] = { 0, 2 };

    template <>
    const unsigned int ParallelEdges<2>::parallel_edges[4][1] = { {1}, {0}, {3}, {2} };

    template <>
    const unsigned int ParallelEdges<3>::starter_edges[3] = { 0, 2, 8 };

    template <>
    const unsigned int ParallelEdges<3>::parallel_edges[12][3] = { {1, 4, 5},    // line 0
      {0, 4, 5},    // line 1
      {3, 6, 7},    // line 2
      {2, 6, 7},    // line 3
      {0, 1, 5},    // line 4
      {0, 1, 4},    // line 5
      {2, 3, 7},    // line 6
      {2, 3, 6},    // line 7
      {9, 10, 11},  // line 8
      {8, 10, 11},  // line 9
      {8, 9, 11},   // line 10
      {8, 9, 10}   // line 11
    };


    template <int dim>
    struct Edge;

    /**
     * A class that describes all of the relevant properties of an
     * edge. For the purpose of what we do here, that includes the
     * indices of the two vertices, and the indices of the adjacent
     * cells (together with a description *where* in each of the
     * adjacent cells the edge is located). It also includes the
     * (global) direction of the edge: either from the first vertex to
     * the second, the other way around, or so far undetermined.
     */
    template <>
    struct Edge<2>
    {
      static const unsigned int dim = 2;

      /**
       * Default constructor. Creates an invalid edge.
       */
      Edge ()
        :
        orientation_status (not_oriented)
      {
        for (unsigned int i=0; i<2; ++i)
          vertex_indices[i] = numbers::invalid_unsigned_int;

        for (unsigned int i=0; i<2; ++i)
          {
            static const AdjacentCell invalid_cell = { numbers::invalid_unsigned_int,
                                                       numbers::invalid_unsigned_int
                                                     };
            adjacent_cells[i] = invalid_cell;
          }
      }

      /**
       * Constructor. Create the edge based on the information given
       * in @p cell, and selecting the edge with number @p edge_number
       * within this cell. Initialize the edge as unoriented.
       */
      Edge (const CellData<dim> &cell,
            const unsigned int   edge_number)
        :
        orientation_status (not_oriented)
      {
        Assert (edge_number < GeometryInfo<dim>::lines_per_cell, ExcInternalError());

        // copy vertices for this particular line
        vertex_indices[0] = cell.vertices[GeometryInfo<dim>::line_to_cell_vertices(edge_number, 0)];
        vertex_indices[1] = cell.vertices[GeometryInfo<dim>::line_to_cell_vertices(edge_number, 1)];

        // bring them into standard orientation
        if (vertex_indices[0] > vertex_indices[1])
          std::swap (vertex_indices[0], vertex_indices[1]);

        for (unsigned int i=0; i<2; ++i)
          {
            static const AdjacentCell invalid_cell = { numbers::invalid_unsigned_int,
                                                       numbers::invalid_unsigned_int
                                                     };
            adjacent_cells[i] = invalid_cell;
          }
      }

      /**
       * Comparison operator for edges. It compares based on the
       * lexicographic ordering of the two vertex indices.
       */
      bool operator< (const Edge<dim> &e) const
      {
        return ((vertex_indices[0] < e.vertex_indices[0])
                ||
                ((vertex_indices[0] == e.vertex_indices[0]) && (vertex_indices[1] < e.vertex_indices[1])));
      }

      /**
       * Compare two edges for equality based on their vertex indices.
       */
      bool operator== (const Edge<dim> &e) const
      {
        return ((vertex_indices[0] == e.vertex_indices[0])
                &&
                (vertex_indices[1] == e.vertex_indices[1]));
      }

      /**
       * Store the given cell index as one of the cells as adjacent to
       * the current edge. Also store that the current edge has index
       * @p edge_within_cell within the given cell.
       */
      void add_adjacent_cell (const unsigned int cell_index,
                              const unsigned int edge_within_cell)
      {
        const AdjacentCell adjacent_cell = { cell_index, edge_within_cell };
        if (adjacent_cells[0].cell_index == numbers::invalid_unsigned_int)
          adjacent_cells[0] = adjacent_cell;
        else
          {
            Assert (adjacent_cells[1].cell_index == numbers::invalid_unsigned_int,
                    ExcInternalError());
            adjacent_cells[1] = adjacent_cell;
          }
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
       * A structure that stores how this edge relates to the at most
       * two cells that are next to it.
       */
      struct AdjacentCell
      {
        unsigned int cell_index;
        unsigned int edge_within_cell;
      };

      AdjacentCell adjacent_cells[2];
    };



    /**
     * A data structure that represents a cell with all of its vertices
     * and edges.
     */
    template <int dim>
    struct Cell
    {
      /**
       * Default construct a cell.
       */
      Cell ()
      {
        for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
          vertex_indices[i] = numbers::invalid_unsigned_int;
        for (unsigned int i=0; i<GeometryInfo<dim>::lines_per_cell; ++i)
          edge_indices[i] = numbers::invalid_unsigned_int;
      }

      /**
       * Construct a Cell object from a CellData object. Also take a (sorted)
       * list of edges and to point into from the current object.
       * @param c
       * @param edge_list
       */
      Cell (const CellData<dim>           &c,
            const std::vector<Edge<dim> > &edge_list)
      {
        for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
          vertex_indices[i] = c.vertices[i];

        // now for each of the edges of this cell, find the location inside the
        // given edge_list array and store than index
        for (unsigned int l=0; l<GeometryInfo<dim>::lines_per_cell; ++l)
          {
            const Edge<dim> e (c, l);
            edge_indices[l] = (std::lower_bound (edge_list.begin(), edge_list.end(), e)
                               -
                               edge_list.begin());
            Assert (edge_indices[l] < edge_list.size(), ExcInternalError());
            Assert (edge_list[edge_indices[l]] == e, ExcInternalError())
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




    /**
     * From a list of cells, build a sorted vector that contains all of the edges
     * that exist in the mesh.
     */
    template <int dim>
    std::vector<Edge<dim> >
    build_edges (const std::vector<CellData<dim> > &cells)
    {
      // build the edge list for all cells. because each cell has
      // GeometryInfo<dim>::lines_per_cell edges, the total number
      // of edges is this many times the number of cells. of course
      // some of them will be duplicates, and we throw them out below
      std::vector<Edge<dim> > edge_list;
      edge_list.reserve(cells.size()*GeometryInfo<dim>::lines_per_cell);
      for (unsigned int i=0; i<cells.size(); ++i)
        for (unsigned int l=0; l<GeometryInfo<dim>::lines_per_cell; ++l)
          edge_list.push_back (Edge<dim>(cells[i], l));

      // next sort the edge list and then remove duplicates
      std::sort (edge_list.begin(), edge_list.end());
      edge_list.erase(std::unique(edge_list.begin(),edge_list.end()),
                      edge_list.end());

      return edge_list;
    }

    /**
     * Build the cell list. Update the edge array to let edges know
     * which cells are adjacent to them.
     */
    template <int dim>
    std::vector<Cell<dim> >
    build_cells_and_connect_edges (const std::vector<CellData<dim> > &cells,
                                   std::vector<Edge<dim> > &edges)
    {
      std::vector<Cell<dim> > cell_list;
      cell_list.reserve(cells.size());
      for (unsigned int i=0; i<cells.size(); ++i)
        {
          // create our own data structure for the cells and let it
          // connect to the edges array
          cell_list.push_back (Cell<dim>(cells[i], edges));

          // then also inform the edges that they are adjacent
          // to the current cell, and where within this cell
          for (unsigned int l=0; l<GeometryInfo<dim>::lines_per_cell; ++l)
            edges[cell_list.back().edge_indices[l]].add_adjacent_cell (i, l);
        }
      Assert (cell_list.size() == cells.size(), ExcInternalError());

      return cell_list;
    }


    /**
     * Return the index within 'cells' of the first cell that has at least one
     * edge that is not yet oriented.
     */
    template <int dim>
    unsigned int
    get_next_unoriented_quad(const std::vector<Cell<dim> > &cells,
                             const std::vector<Edge<dim> > &edges)
    {
      for (unsigned int c=0; c<cells.size(); ++c)
        for (unsigned int l=0; l<GeometryInfo<dim>::lines_per_cell; ++l)
          if (edges[cells[c].edge_indices[l]].orientation_status == Edge<dim>::not_oriented)
            return c;

      return numbers::invalid_unsigned_int;
    }


    /**
     * Given a set of cells and edges, orient all edges that are
     * (globall) parallel to the one identified by the @p cell and
     * within it the one with index @p local_edge.
     */
    void
    orient_one_set_of_parallel_edges (const std::vector<Cell<2> > &cells,
                                      std::vector<Edge<2> >       &edges,
                                      const unsigned int           cell,
                                      const unsigned int           local_edge)
    {
      const unsigned int dim = 2;

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
      if (edges[cells[cell].edge_indices[local_edge]].vertex_indices[0]
          ==
          cells[cell].vertex_indices[GeometryInfo<dim>::line_to_cell_vertices (local_edge, 0)])
        // orient initial edge *opposite* to the way it is in the cell
        // (see above for the reason)
        edges[cells[cell].edge_indices[local_edge]].orientation_status = Edge<dim>::backward;
      else
        {
          Assert (edges[cells[cell].edge_indices[local_edge]].vertex_indices[0]
                  ==
                  cells[cell].vertex_indices[GeometryInfo<dim>::line_to_cell_vertices (local_edge, 1)],
                  ExcInternalError());
          Assert (edges[cells[cell].edge_indices[local_edge]].vertex_indices[1]
                  ==
                  cells[cell].vertex_indices[GeometryInfo<dim>::line_to_cell_vertices (local_edge, 0)],
                  ExcInternalError());

          // orient initial edge *opposite* to the way it is in the cell
          // (see above for the reason)
          edges[cells[cell].edge_indices[local_edge]].orientation_status = Edge<dim>::forward;
        }

      // walk outward from the given edge as described in
      // the algorithm in the paper that documents all of
      // this
      //
      // note that in 2d, each of the Deltas can at most
      // contain two elements. we indicate non-used elements
      // of these sets by invalid unsigned ints; if the set has
      // only one element, then we use the first
      unsigned int Delta_k[2] = { numbers::invalid_unsigned_int,
                                  numbers::invalid_unsigned_int
                                };
      unsigned int Delta_k_minus_1[2] = { cells[cell].edge_indices[local_edge],
                                          numbers::invalid_unsigned_int
                                        };
      while (Delta_k_minus_1[0] != numbers::invalid_unsigned_int)   // while set is not empty
        {
          Delta_k[0] = Delta_k[1] = numbers::invalid_unsigned_int;

          for (unsigned int delta_element=0; delta_element<2; ++delta_element)
            if (Delta_k_minus_1[delta_element] != numbers::invalid_unsigned_int)
              {
                // get the edge we are currently looking at. it must already
                // have been oriented
                const unsigned int delta = Delta_k_minus_1[delta_element];
                Assert (edges[delta].orientation_status != Edge<dim>::not_oriented,
                        ExcInternalError());

                // now go through the cells adjacent to this edge
                for (unsigned int K_element=0; K_element<2; ++K_element)
                  if (edges[delta].adjacent_cells[K_element].cell_index != numbers::invalid_unsigned_int)
                    {
                      const unsigned int K = edges[delta].adjacent_cells[K_element].cell_index;
                      const unsigned int delta_is_edge_in_K = edges[delta].adjacent_cells[K_element].edge_within_cell;

                      // figure out the direction of delta with respect to the cell K
                      // (in the orientation in which the user has given it to us)
                      const unsigned int first_edge_vertex
                        = (edges[delta].orientation_status == Edge<dim>::forward
                           ?
                           edges[delta].vertex_indices[0]
                           :
                           edges[delta].vertex_indices[1]);
                      const unsigned int first_edge_vertex_in_K = cells[K].vertex_indices[GeometryInfo<dim>::face_to_cell_vertices(delta_is_edge_in_K, 0)];
                      Assert (first_edge_vertex == first_edge_vertex_in_K
                              ||
                              first_edge_vertex == cells[K].vertex_indices[GeometryInfo<dim>::face_to_cell_vertices(delta_is_edge_in_K, 1)],
                              ExcInternalError());

                      // now figure out which direction the opposite edge needs to be into.
                      const unsigned int opposite_edge
                        = cells[K].edge_indices[ParallelEdges<2>::parallel_edges[delta_is_edge_in_K][0]];
                      const unsigned int first_opposite_edge_vertex
                        =  cells[K].vertex_indices[GeometryInfo<dim>::face_to_cell_vertices(
                                                     ParallelEdges<dim>::parallel_edges[delta_is_edge_in_K][0],
                                                     (first_edge_vertex == first_edge_vertex_in_K
                                                      ?
                                                      0
                                                      :
                                                      1))];

                      const Edge<dim>::OrientationStatus opposite_edge_orientation
                        = (edges[opposite_edge].vertex_indices[0]
                           ==
                           first_opposite_edge_vertex
                           ?
                           Edge<dim>::forward
                           :
                           Edge<dim>::backward);

                      // see if the opposite edge (there is only one in 2d) has already been
                      // oriented.
                      if (edges[opposite_edge].orientation_status == Edge<dim>::not_oriented)
                        {
                          // the opposite edge is not yet oriented. do orient it and add it to
                          // Delta_k;
                          edges[opposite_edge].orientation_status = opposite_edge_orientation;
                          if (Delta_k[0] == numbers::invalid_unsigned_int)
                            Delta_k[0] = opposite_edge;
                          else
                            {
                              Assert (Delta_k[1] == numbers::invalid_unsigned_int, ExcInternalError());
                              Delta_k[1] = opposite_edge;
                            }
                        }
                      else
                        {
                          // the opposite edge has already been oriented. assert that it is
                          // consistent with the current one
                          Assert (edges[opposite_edge].orientation_status == opposite_edge_orientation,
                                  ExcInternalError());
                        }
                    }
              }

          // finally copy the new set to the previous one
          // (corresponding to increasing 'k' by one in the
          // algorithm)
          Delta_k_minus_1[0] = Delta_k[0];
          Delta_k_minus_1[1] = Delta_k[1];
        }
    }


    /**
     * Given data structures @p cell_list and @p edge_list, where
     * all edges are already oriented, rotate the cell with
     * index @p cell_index in such a way that its local coordinate
     * system matches the ones of the adjacent edges. Store the
     * rotated order of vertices in <code>raw_cells[cell_index]</code>.
     */
    void
    rotate_cell (const std::vector<Cell<2> > &cell_list,
                 const std::vector<Edge<2> > &edge_list,
                 const unsigned int           cell_index,
                 std::vector<CellData<2> >   &raw_cells)
    {
      // find the first vertex of the cell. this is the
      // vertex where two edges originate, so for
      // each of the four edges record which the
      // starting vertex is
      unsigned int starting_vertex_of_edge[4];
      for (unsigned int e=0; e<4; ++e)
        {
          Assert (edge_list[cell_list[cell_index].edge_indices[e]].orientation_status
                  != Edge<2>::not_oriented,
                  ExcInternalError());
          if (edge_list[cell_list[cell_index].edge_indices[e]].orientation_status == Edge<2>::forward)
            starting_vertex_of_edge[e] = edge_list[cell_list[cell_index].edge_indices[e]].vertex_indices[0];
          else
            starting_vertex_of_edge[e] = edge_list[cell_list[cell_index].edge_indices[e]].vertex_indices[1];
        }

      // find the vertex number that appears twice. this must either be
      // the first, second, or third vertex in the list. because edges
      // zero and one don't share any vertices, and the same for edges
      // two and three, the possibilities can easily be enumerated
      unsigned int starting_vertex_of_cell = numbers::invalid_unsigned_int;
      if ((starting_vertex_of_edge[0] == starting_vertex_of_edge[2])
          ||
          (starting_vertex_of_edge[0] == starting_vertex_of_edge[3]))
        starting_vertex_of_cell = starting_vertex_of_edge[0];
      else if ((starting_vertex_of_edge[1] == starting_vertex_of_edge[2])
               ||
               (starting_vertex_of_edge[1] == starting_vertex_of_edge[3]))
        starting_vertex_of_cell = starting_vertex_of_edge[1];
      else
        Assert (false, ExcInternalError());

      // now rotate raw_cells[cell_index] until the starting indices match.
      // take into account the ordering of vertices (not in clockwise
      // or counter-clockwise sense)
      while (raw_cells[cell_index].vertices[0] != starting_vertex_of_cell)
        {
          const unsigned int tmp = raw_cells[cell_index].vertices[0];
          raw_cells[cell_index].vertices[0] = raw_cells[cell_index].vertices[1];
          raw_cells[cell_index].vertices[1] = raw_cells[cell_index].vertices[3];
          raw_cells[cell_index].vertices[3] = raw_cells[cell_index].vertices[2];
          raw_cells[cell_index].vertices[2] = tmp;
        }
    }


    /**
     * Given a set of cells, find globally unique edge orientations
     * and then rotate cells so that the coordinate system of the cell
     * coincides with the coordinate systems of the adjacent edges.
     */
    template <int dim>
    void reorient (std::vector<CellData<dim> > &cells)
    {
      // first build the arrays that connect cells to edges and the other
      // way around
      std::vector<Edge<dim> > edge_list = build_edges(cells);
      std::vector<Cell<dim> > cell_list = build_cells_and_connect_edges(cells, edge_list);

      // then loop over all cells and start orienting parallel edge sets
      // of cells that still have non-oriented edges
      unsigned int next_cell_with_unoriented_edge;
      while ((next_cell_with_unoriented_edge = get_next_unoriented_quad(cell_list, edge_list)) !=
             numbers::invalid_unsigned_int)
        {
          // see which edge sets are still not oriented
          //
          // we do not need to look at each edge because if we orient edge
          // 0, we will end up with edge 1 also oriented. there are only
          // dim independent sets of edges
          for (unsigned int l=0; l<dim; ++l)
            if (edge_list[cell_list[next_cell_with_unoriented_edge].edge_indices[ParallelEdges<dim>::starter_edges[l]]].orientation_status
                == Edge<dim>::not_oriented)
              orient_one_set_of_parallel_edges (cell_list,
                                                edge_list,
                                                next_cell_with_unoriented_edge,
                                                ParallelEdges<dim>::starter_edges[l]);

          // ensure that we have really oriented all edges now, not just
          // the starter edges
          for (unsigned int l=0; l<GeometryInfo<dim>::lines_per_cell; ++l)
            Assert (edge_list[cell_list[next_cell_with_unoriented_edge].edge_indices[l]].orientation_status
                    != Edge<dim>::not_oriented,
                    ExcInternalError());
        }

      // now that we have oriented all edges, we need to rotate cells
      // so that the edges point in the right direction with the now
      // rotated coordinate system
      for (unsigned int c=0; c<cells.size(); ++c)
        rotate_cell (cell_list, edge_list, c, cells);
    }

  }
}


template<>
void
GridReordering<1>::reorder_cells (std::vector<CellData<1> > &,
                                  const bool)
{
  // there should not be much to do
  // in 1d...
}


template<>
void
GridReordering<1>::invert_all_cells_of_negative_grid(const std::vector<Point<1> > &,
                                                     std::vector<CellData<1> > &)
{
  // nothing to be done in 1d
}

template<>
void
GridReordering<1,2>::reorder_cells (std::vector<CellData<1> > &,
                                    const bool)
{
  // there should not be much to do
  // in 1d...
}


template<>
void
GridReordering<1,2>::invert_all_cells_of_negative_grid(const std::vector<Point<2> > &,
                                                       std::vector<CellData<1> > &)
{
  // nothing to be done in 1d
}


template<>
void
GridReordering<1,3>::reorder_cells (std::vector<CellData<1> > &,
                                    const bool)
{
  // there should not be much to do
  // in 1d...
}


template<>
void
GridReordering<1,3>::invert_all_cells_of_negative_grid(const std::vector<Point<3> > &,
                                                       std::vector<CellData<1> > &)
{
  // nothing to be done in 1d
}


// anonymous namespace for internal helper functions
namespace
{
  /**
   * A set of functions that
   * reorder the data from the
   * "current" to the "classic"
   * format of vertex numbering of
   * cells and faces. These functions
   * do the reordering of their
   * arguments in-place.
   */
  void
  reorder_new_to_old_style (std::vector<CellData<2> > &cells)
  {
    for (unsigned int cell=0; cell<cells.size(); ++cell)
      std::swap(cells[cell].vertices[2], cells[cell].vertices[3]);
  }


  void
  reorder_new_to_old_style (std::vector<CellData<3> > &cells)
  {
    unsigned int tmp[GeometryInfo<3>::vertices_per_cell];
    for (unsigned int cell=0; cell<cells.size(); ++cell)
      {
        for (unsigned int i=0; i<GeometryInfo<3>::vertices_per_cell; ++i)
          tmp[i] = cells[cell].vertices[i];
        for (unsigned int i=0; i<GeometryInfo<3>::vertices_per_cell; ++i)
          cells[cell].vertices[i] = tmp[GeometryInfo<3>::ucd_to_deal[i]];
      }
  }


  /**
   * And now also in the opposite direction.
   */
  void
  reorder_old_to_new_style (std::vector<CellData<2> > &cells)
  {
    // just invert the permutation:
    reorder_new_to_old_style(cells);
  }


  void
  reorder_old_to_new_style (std::vector<CellData<3> > &cells)
  {
    // undo the ordering above
    unsigned int tmp[GeometryInfo<3>::vertices_per_cell];
    for (unsigned int cell=0; cell<cells.size(); ++cell)
      {
        for (unsigned int i=0; i<GeometryInfo<3>::vertices_per_cell; ++i)
          tmp[i] = cells[cell].vertices[i];
        for (unsigned int i=0; i<GeometryInfo<3>::vertices_per_cell; ++i)
          cells[cell].vertices[GeometryInfo<3>::ucd_to_deal[i]] = tmp[i];
      }
  }
}


template<>
void
GridReordering<2>::reorder_cells (std::vector<CellData<2> > &cells,
                                  const bool use_new_style_ordering)
{
  // if necessary, convert to old (compatibility) to new-style format
  if (!use_new_style_ordering)
    reorder_old_to_new_style(cells);

  // check if grids are already
  // consistent. if so, do
  // nothing. if not, then do the
  // reordering
  if (!internal::GridReordering2d::is_consistent (cells))
    internal::GridReordering2d::reorient(cells);


  // and convert back if necessary
  if (!use_new_style_ordering)
    reorder_new_to_old_style(cells);
}


template<>
void
GridReordering<2,3>::reorder_cells (std::vector<CellData<2> > &cells,
                                    const bool use_new_style_ordering)
{
  // if necessary, convert to old-style format
  if (use_new_style_ordering)
    reorder_new_to_old_style(cells);

  GridReordering<2>::reorder_cells(cells);


  // and convert back if necessary
  if (use_new_style_ordering)
    reorder_old_to_new_style(cells);
}



template<>
void
GridReordering<2>::invert_all_cells_of_negative_grid(const std::vector<Point<2> > &all_vertices,
                                                     std::vector<CellData<2> >    &cells)
{
  unsigned int vertices_lex[GeometryInfo<2>::vertices_per_cell];
  unsigned int n_negative_cells=0;
  for (unsigned int cell_no=0; cell_no<cells.size(); ++cell_no)
    {
      // GridTools::cell_measure
      // requires the vertices to be
      // in lexicographic ordering
      for (unsigned int i=0; i<GeometryInfo<2>::vertices_per_cell; ++i)
        vertices_lex[GeometryInfo<2>::ucd_to_deal[i]]=cells[cell_no].vertices[i];
      if (GridTools::cell_measure<2>(all_vertices, vertices_lex) < 0)
        {
          ++n_negative_cells;
          std::swap(cells[cell_no].vertices[1], cells[cell_no].vertices[3]);

          // check whether the
          // resulting cell is now ok.
          // if not, then the grid is
          // seriously broken and
          // should be sticked into the
          // bin
          for (unsigned int i=0; i<GeometryInfo<2>::vertices_per_cell; ++i)
            vertices_lex[GeometryInfo<2>::ucd_to_deal[i]]=cells[cell_no].vertices[i];
          AssertThrow(GridTools::cell_measure<2>(all_vertices, vertices_lex) > 0,
                      ExcInternalError());
        }
    }

  // We assume that all cells of a grid have
  // either positive or negative volumes but
  // not both mixed. Although above reordering
  // might work also on single cells, grids
  // with both kind of cells are very likely to
  // be broken. Check for this here.
  AssertThrow(n_negative_cells==0 || n_negative_cells==cells.size(),
              ExcMessage(std::string("This class assumes that either all cells have positive "
                                     "volume, or that all cells have been specified in an "
                                     "inverted vertex order so that their volume is negative. "
                                     "(In the latter case, this class automatically inverts "
                                     "every cell.) However, the mesh you have specified "
                                     "appears to have both cells with positive and cells with "
                                     "negative volume. You need to check your mesh which "
                                     "cells these are and how they got there.\n"
                                     "As a hint, of the total ")
                         + Utilities::to_string (cells.size())
                         + " cells in the mesh, "
                         + Utilities::to_string (n_negative_cells)
                         + " appear to have a negative volume."));
}



template<>
void
GridReordering<2,3>::invert_all_cells_of_negative_grid(const std::vector<Point<3> > &,
                                                       std::vector<CellData<2> > &)
{
  Assert(false, ExcNotImplemented());
}



namespace internal
{
  namespace GridReordering3d
  {
    DeclException1 (ExcGridOrientError,
                    char *,
                    <<  "Grid Orientation Error: " << arg1);

    const EdgeOrientation unoriented_edge = {'u'};
    const EdgeOrientation forward_edge    = {'f'};
    const EdgeOrientation backward_edge   = {'b'};


    inline
    bool
    EdgeOrientation::
    operator == (const EdgeOrientation &edge_orientation) const
    {
      Assert ((orientation == 'u') || (orientation == 'f') || (orientation == 'b'),
              ExcInternalError());
      return orientation == edge_orientation.orientation;
    }



    inline
    bool
    EdgeOrientation::
    operator != (const EdgeOrientation &edge_orientation) const
    {
      return ! (*this == edge_orientation);
    }



    namespace ElementInfo
    {
      /**
       * The numbers of the edges
       * coming into node i are
       * given by
       * edge_to_node[i][k] where
       * k=0,1,2.
       */
      static const unsigned int edge_to_node[8][3] =
      {
        {0,4,8},
        {0,5,9},
        {3,5,10},
        {3,4,11},
        {1,7,8},
        {1,6,9},
        {2,6,10},
        {2,7,11}
      };


      /**
       * The orientation of edge
       * coming into node i is
       * given by
       * edge_to_node_orient[i][k]
       * where k=0,1,2. 1 means the
       * given node is the start of
       * the edge -1 means the end
       * of the edge.
       */
      static const EdgeOrientation edge_to_node_orient[8][3] =
      {
        {forward_edge,  forward_edge,  forward_edge},
        {backward_edge, forward_edge,  forward_edge},
        {backward_edge, backward_edge, forward_edge},
        {forward_edge,  backward_edge, forward_edge},
        {forward_edge,  forward_edge,  backward_edge},
        {backward_edge, forward_edge,  backward_edge},
        {backward_edge, backward_edge, backward_edge},
        {forward_edge,  backward_edge, backward_edge}
      };

      /**
       * nodesonedge[i][0] is the
       * start node for edge i.
       * nodesonedge[i][1] is the
       * end node for edge i.
       */
      static const unsigned int nodes_on_edge[12][2] =
      {
        {0,1},
        {4,5},
        {7,6},
        {3,2},
        {0,3},
        {1,2},
        {5,6},
        {4,7},
        {0,4},
        {1,5},
        {2,6},
        {3,7}
      };
    }


    CheapEdge::CheapEdge (const unsigned int n0,
                          const unsigned int n1)
      :
      // sort the
      // entries so
      // that
      // node0<node1
      node0(std::min (n0, n1)),
      node1(std::max (n0, n1))
    {}



    bool CheapEdge::operator< (const CheapEdge &e2) const
    {
      if (node0 < e2.node0) return true;
      if (node0 > e2.node0) return false;
      if (node1 < e2.node1) return true;
      return false;
    }


    Edge::Edge (const unsigned int n0,
                const unsigned int n1)
      :
      orientation_flag (unoriented_edge),
      group (numbers::invalid_unsigned_int)
    {
      nodes[0] = n0;
      nodes[1] = n1;
    }



    Cell::Cell ()
    {
      for (unsigned int i=0; i<GeometryInfo<3>::lines_per_cell; ++i)
        {
          edges[i] = numbers::invalid_unsigned_int;
          local_orientation_flags[i] = forward_edge;
        }

      for (unsigned int i=0; i<GeometryInfo<3>::vertices_per_cell; ++i)
        nodes[i] = numbers::invalid_unsigned_int;

      waiting_to_be_processed = false;
    }



    Mesh::Mesh (const std::vector<CellData<3> > &incubes)
    {
      // copy the cells into our own
      // internal data format.
      const unsigned int numelems = incubes.size();
      for (unsigned int i=0; i<numelems; ++i)
        {
          Cell the_cell;
          std::copy (&incubes[i].vertices[0],
                     &incubes[i].vertices[GeometryInfo<3>::vertices_per_cell],
                     &the_cell.nodes[0]);

          cell_list.push_back(the_cell);
        }

      // then build edges and
      // connectivity
      build_connectivity ();
    }



    void
    Mesh::sanity_check () const
    {
      for (unsigned int i=0; i<cell_list.size(); ++i)
        for (unsigned int j=0; j<8; ++j)
          sanity_check_node (cell_list[i], j);
    }



    void
    Mesh::sanity_check_node (const Cell         &c,
                             const unsigned int local_node_num) const
    {
#ifdef DEBUG
      // check that every edge
      // coming into a node has the
      // same node value

      // Get the Local Node Numbers
      // of the incoming edges
      const unsigned int e0 = ElementInfo::edge_to_node[local_node_num][0];
      const unsigned int e1 = ElementInfo::edge_to_node[local_node_num][1];
      const unsigned int e2 = ElementInfo::edge_to_node[local_node_num][2];

      // Global Edge Numbers
      const unsigned int ge0 = c.edges[e0];
      const unsigned int ge1 = c.edges[e1];
      const unsigned int ge2 = c.edges[e2];

      const EdgeOrientation or0 = ElementInfo::edge_to_node_orient[local_node_num][0] ==
                                  c.local_orientation_flags[e0] ?
                                  forward_edge : backward_edge;
      const EdgeOrientation or1 = ElementInfo::edge_to_node_orient[local_node_num][1] ==
                                  c.local_orientation_flags[e1] ?
                                  forward_edge : backward_edge;
      const EdgeOrientation or2 = ElementInfo::edge_to_node_orient[local_node_num][2] ==
                                  c.local_orientation_flags[e2] ?
                                  forward_edge : backward_edge;

      // Make sure that edges agree
      // what the current node should
      // be.
      Assert ((edge_list[ge0].nodes[or0 == forward_edge ? 0 : 1] ==
               edge_list[ge1].nodes[or1 == forward_edge ? 0 : 1])
              &&
              (edge_list[ge1].nodes[or1 == forward_edge ? 0 : 1] ==
               edge_list[ge2].nodes[or2 == forward_edge ? 0 : 1]),
              ExcMessage ("This message does not satisfy the internal "
                          "consistency check"));
#else
      (void)c;
      (void)local_node_num;
#endif
    }



    // This is the guts of the matter...
    void Mesh::build_connectivity ()
    {
      const unsigned int n_cells = cell_list.size();

      unsigned int n_edges = 0;
      // Correctly build the edge
      // list
      {
        // edge_map stores the
        // edge_number associated
        // with a given CheapEdge
        std::map<CheapEdge,unsigned int> edge_map;
        unsigned int ctr = 0;
        for (unsigned int cur_cell_id = 0;
             cur_cell_id<n_cells;
             ++cur_cell_id)
          {
            // Get the local node
            // numbers on edge
            // edge_num
            const Cell &cur_cell = cell_list[cur_cell_id];

            for (unsigned short int edge_num = 0;
                 edge_num<12;
                 ++edge_num)
              {
                unsigned int gl_edge_num = 0;
                EdgeOrientation l_edge_orient = forward_edge;

                // Construct the
                // CheapEdge
                const unsigned int
                node0 = cur_cell.nodes[ElementInfo::nodes_on_edge[edge_num][0]],
                node1 = cur_cell.nodes[ElementInfo::nodes_on_edge[edge_num][1]];
                const CheapEdge cur_edge (node0, node1);

                if (edge_map.count(cur_edge) == 0)
                  // Edge not in map
                  {
                    // put edge in
                    // hash map with
                    // ctr value;
                    edge_map[cur_edge] = ctr;
                    gl_edge_num = ctr;

                    // put the edge
                    // into the
                    // global edge
                    // list
                    edge_list.push_back(Edge(node0,node1));
                    ctr++;
                  }
                else
                  {
                    // get edge_num
                    // from hash_map
                    gl_edge_num = edge_map[cur_edge];
                    if (edge_list[gl_edge_num].nodes[0] != node0)
                      l_edge_orient = backward_edge;
                  }
                // set edge number to
                // edgenum
                cell_list[cur_cell_id].edges[edge_num] = gl_edge_num;
                cell_list[cur_cell_id].local_orientation_flags[edge_num]
                  = l_edge_orient;
              }
          }
        n_edges = ctr;
      }

      // Count each of the edges.
      {
        std::vector<int> edge_count(n_edges,0);


        // Count every time an edge
        // occurs in a cube.
        for (unsigned int cur_cell_id=0; cur_cell_id<n_cells; ++cur_cell_id)
          for (unsigned short int edge_num = 0; edge_num<12; ++edge_num)
            ++edge_count[cell_list[cur_cell_id].edges[edge_num]];

        // So we now know how many
        // cubes contain a given
        // edge. Just need to store
        // the list of cubes in the
        // edge

        // Allocate the space for the
        // neighbor list
        for (unsigned int cur_edge_id=0; cur_edge_id<n_edges; ++cur_edge_id)
          edge_list[cur_edge_id].neighboring_cubes
          .resize (edge_count[cur_edge_id]);

        // Store the position of the
        // current neighbor in the
        // edge's neighbor list
        std::vector<int> cur_cell_edge_list_posn(n_edges,0);
        for (unsigned int cur_cell_id=0; cur_cell_id<n_cells; ++cur_cell_id)
          for (unsigned short int edge_num=0; edge_num<12; ++edge_num)
            {
              const unsigned int
              gl_edge_id = cell_list[cur_cell_id].edges[edge_num];
              Edge &cur_edge = edge_list[gl_edge_id];
              cur_edge.neighboring_cubes[cur_cell_edge_list_posn[gl_edge_id]]
                = cur_cell_id;
              cur_cell_edge_list_posn[gl_edge_id]++;
            }
      }
    }



    void
    Mesh::export_to_deal_format (std::vector<CellData<3> > &outcubes) const
    {
      Assert (outcubes.size() == cell_list.size(),
              ExcInternalError());

      // simply overwrite the output
      // array with the new
      // information
      for (unsigned int i=0; i<cell_list.size(); ++i)
        std::copy (&cell_list[i].nodes[0],
                   &cell_list[i].nodes[GeometryInfo<3>::vertices_per_cell],
                   &outcubes[i].vertices[0]);
    }



    Orienter::Orienter (const std::vector<CellData<3> > &incubes)
      :
      mesh (incubes),
      cur_posn (0),
      marker_cube (0),
      cur_edge_group  (0)
    {
      for (unsigned int i = 0; i<12; ++i)
        edge_orient_array[i] = false;
    }



    bool Orienter::orient_mesh (std::vector<CellData<3> > &incubes)
    {
      Orienter orienter (incubes);

      // First check that the mesh is
      // sensible
      orienter.mesh.sanity_check ();

      // Orient the mesh

      // if not successful, break here, else go
      // on
      if (!orienter.orient_edges ())
        return false;

      // Now we have a bunch of oriented
      // edges int the structure we only
      // have to turn the cubes so they
      // match the edge orientation.
      orienter.orient_cubes ();

      // Copy the elements from our
      // internal structure back into
      // their original location.
      orienter.mesh.export_to_deal_format (incubes);
      // reordering was successful
      return true;
    }

    /**
     * This assigns an orientation
     * to each edge so that every
     * cube is a rotated Deal.II
     * cube.
     */
    bool Orienter::orient_edges ()
    {
      // While there are still cubes
      // to orient
      while (get_next_unoriented_cube())
        // And there are edges in
        // the cube to orient
        while (orient_next_unoriented_edge())
          {
            // Make all the sides
            // in the current set
            // match
            orient_edges_in_current_cube();

            // Add the adjacent
            // cubes to the list
            // for processing
            get_adjacent_cubes();
            // Start working on
            // this list of cubes
            while (get_next_active_cube())
              {
                // Make sure the
                // Cube doesn't
                // have a
                // contradiction
                if (!cell_is_consistent(cur_posn))
                  return false;

                // If we needed to
                // orient any edges
                // in the current
                // cube then we may
                // have to process
                // the neighbor.
                if (orient_edges_in_current_cube())
                  get_adjacent_cubes();
              }

            // start the next sheet
            // (equivalence class
            // of edges)
            ++cur_edge_group;
          }
      return true;
    }



    bool Orienter::get_next_unoriented_cube ()
    {
      // The last cube in the list
      const unsigned int n_cubes = mesh.cell_list.size();
      // Keep shifting along the list
      // until we find a cube which
      // is not fully oriented or the
      // end.
      while ( (marker_cube<n_cubes) &&
              (is_oriented(marker_cube)) )
        ++marker_cube;
      cur_posn = marker_cube;
      // Return true if we now point
      // at a valid cube.
      return (cur_posn < n_cubes);
    }



    bool Orienter::is_oriented (const unsigned int cell_num) const
    {
      for (unsigned int i=0; i<12; ++i)
        if (mesh.edge_list[mesh.cell_list[cell_num].edges[i]].orientation_flag
            == unoriented_edge)
          return false;
      return true;
    }



    bool
    Orienter::cell_is_consistent(const unsigned int cell_num) const
    {

      const Cell &c = mesh.cell_list[cell_num];

      // Checks that all oriented
      // edges in the group are
      // oriented consistently.
      for (unsigned int group=0; group<3; ++group)
        {
          // When a nonzero
          // orientation is first
          // encountered in the group
          // it is stored in this
          EdgeOrientation value = unoriented_edge;
          // Loop over all parallel
          // edges
          for (unsigned int i=4*group; i<4*(group+1); ++i)
            {
              // If the edge has
              // orientation
              if ((c.local_orientation_flags[i] !=
                   unoriented_edge)
                  &&
                  (mesh.edge_list[c.edges[i]].orientation_flag !=
                   unoriented_edge))
                {
                  const EdgeOrientation this_edge_direction
                    = (c.local_orientation_flags[i]
                       == mesh.edge_list[c.edges[i]].orientation_flag  ?
                       forward_edge : backward_edge);

                  // If we haven't
                  // seen an oriented
                  // edge before,
                  // then store its
                  // value:
                  if (value == unoriented_edge)
                    value = this_edge_direction;
                  else
                    // If we have
                    // seen an
                    // oriented edge
                    // in this group
                    // we'd better
                    // have the same
                    // orientation.
                    if (value != this_edge_direction)
                      return false;
                }
            }
        }
      return true;
    }



    bool Orienter::orient_next_unoriented_edge ()
    {
      cur_posn = marker_cube;
      const Cell &c = mesh.cell_list[cur_posn];
      unsigned int edge = 0;

      // search for the unoriented
      // side
      while ((edge<12) &&
             (mesh.edge_list[c.edges[edge]].orientation_flag !=
              unoriented_edge))
        ++edge;

      // if we found none then return
      // false
      if (edge == 12)
        return false;

      // Which edge group we're in.
      const unsigned int edge_group = edge/4;

      // A sanity check that none of
      // the other edges in the group
      // have been oriented yet Each
      // of the edges in the group
      // should be un-oriented
      for (unsigned int j = edge_group*4; j<edge_group*4+4; ++j)
        Assert (mesh.edge_list[c.edges[j]].orientation_flag ==
                unoriented_edge,
                ExcGridOrientError("Tried to orient edge when other edges "
                                   "in group are already oriented!"));

      // Make the edge alignment
      // match that of the local
      // cube.
      mesh.edge_list[c.edges[edge]].orientation_flag
        = c.local_orientation_flags[edge];
      mesh.edge_list[c.edges[edge]].group = cur_edge_group;

      // Remember that we have oriented
      // this edge in the current cell.
      edge_orient_array[edge] = true;

      return true;
    }



    bool Orienter::orient_edges_in_current_cube ()
    {
      for (unsigned int edge_group=0; edge_group<3; ++edge_group)
        if (orient_edge_set_in_current_cube(edge_group) == true)
          return true;

      return false;
    }



    bool
    Orienter::orient_edge_set_in_current_cube (const unsigned int n)
    {
      const Cell &c = mesh.cell_list[cur_posn];

      // Check if any edge is
      // oriented
      unsigned int n_oriented = 0;
      EdgeOrientation glorient   = unoriented_edge;
      unsigned int edge_flags = 0;
      unsigned int cur_flag   = 1;
      for (unsigned int i = 4*n; i<4*(n+1); ++i, cur_flag<<=1)
        {
          if ((mesh.edge_list[c.edges[i]].orientation_flag !=
               unoriented_edge)
              &&
              (c.local_orientation_flags[i] !=
               unoriented_edge))
            {
              ++n_oriented;

              const EdgeOrientation orient
                = (mesh.edge_list[c.edges[i]].orientation_flag ==
                   c.local_orientation_flags[i] ?
                   forward_edge : backward_edge);

              if (glorient == unoriented_edge)
                glorient = orient;
              else
                AssertThrow(orient == glorient,
                            ExcGridOrientError("Attempted to Orient Misaligned cube"));
            }
          else
            edge_flags |= cur_flag;
        }

      // were any of the sides
      // oriented?  were they all
      // already oriented?
      if ((glorient == unoriented_edge) || (n_oriented == 4))
        return false;

      // If so orient all edges
      // consistently.
      cur_flag = 1;
      for (unsigned int i=4*n; i<4*(n+1); ++i, cur_flag<<=1)
        if ((edge_flags & cur_flag) != 0)
          {
            mesh.edge_list[c.edges[i]].orientation_flag
              = (c.local_orientation_flags[i] == glorient ?
                 forward_edge : backward_edge);

            mesh.edge_list[c.edges[i]].group = cur_edge_group;
            // Remember that we have oriented
            // this edge in the current cell.
            edge_orient_array[i] = true;
          }

      return true;
    }



    void Orienter::get_adjacent_cubes ()
    {
      const Cell &c = mesh.cell_list[cur_posn];
      for (unsigned int e=0; e<12; ++e)
        // Only need to add the adjacent
        // cubes for edges we recently
        // oriented
        if (edge_orient_array[e] == true)
          {
            const Edge &the_edge = mesh.edge_list[c.edges[e]];
            for (unsigned int local_cube_num = 0;
                 local_cube_num < the_edge.neighboring_cubes.size();
                 ++local_cube_num)
              {
                const unsigned int
                global_cell_num = the_edge.neighboring_cubes[local_cube_num];
                Cell &ncell = mesh.cell_list[global_cell_num];

                // If the cell is waiting to be
                // processed we dont want to add
                // it to the list a second time.
                if (!ncell.waiting_to_be_processed)
                  {
                    sheet_to_process.push_back(global_cell_num);
                    ncell.waiting_to_be_processed = true;
                  }
              }
          }
      // we're done with this cube so
      // clear its processing flags.
      for (unsigned int e=0; e<12; ++e)
        edge_orient_array[e] = false;

    }



    bool Orienter::get_next_active_cube ()
    {
      // Mark the curent Cube as
      // finished with.
      Cell &c = mesh.cell_list[cur_posn];
      c.waiting_to_be_processed = false;
      if (sheet_to_process.empty() == false)
        {
          cur_posn = sheet_to_process.back();
          sheet_to_process.pop_back();
          return true;
        }
      return false;
    }


    void Orienter::orient_cubes ()
    {
      // We assume that the mesh has
      // all edges oriented already.

      // This is a list of
      // permutations that take node
      // 0 to node i but only rotate
      // the cube.  (This set is far
      // from unique (there are 3 for
      // each node - for our
      // algorithm it doesn't matter
      // which of the three we use)
      static const unsigned int CubePermutations[8][8] =
      {
        {0,1,2,3,4,5,6,7},
        {1,2,3,0,5,6,7,4},
        {2,3,0,1,6,7,4,5},
        {3,0,1,2,7,4,5,6},
        {4,7,6,5,0,3,2,1},
        {5,4,7,6,1,0,3,2},
        {6,5,4,7,2,1,0,3},
        {7,6,5,4,3,2,1,0}
      };

      // So now we need to work out
      // which node needs to be
      // mapped to the zero node.
      // The trick is that the node
      // that should be the local
      // zero node has three edges
      // coming into it.
      for (unsigned int i=0; i<mesh.cell_list.size(); ++i)
        {
          Cell &the_cell = mesh.cell_list[i];

          // This stores whether the
          // global oriented edge
          // points in the same
          // direction as it's local
          // edge on the current
          // cube. (for each edge on
          // the curent cube)
          EdgeOrientation local_edge_orientation[12];
          for (unsigned int j = 0; j<12; ++j)
            {
              // get the global edge
              const Edge &the_edge = mesh.edge_list[the_cell.edges[j]];
              // All edges should be
              // oriented at this
              // stage..
              Assert (the_edge.orientation_flag != unoriented_edge,
                      ExcGridOrientError ("Unoriented edge encountered"));
              // calculate whether it
              // points the right way
              // or not
              local_edge_orientation[j] = (the_cell.local_orientation_flags[j] ==
                                           the_edge.orientation_flag ?
                                           forward_edge : backward_edge);
            }

          // Here the number of
          // incoming edges is
          // tallied for each node.
          unsigned int perm_num = numbers::invalid_unsigned_int;
          for (unsigned int node_num=0; node_num<8; ++node_num)
            {
              // The local edge
              // numbers coming into
              // the node
              const unsigned int e0 = ElementInfo::edge_to_node[node_num][0];
              const unsigned int e1 = ElementInfo::edge_to_node[node_num][1];
              const unsigned int e2 = ElementInfo::edge_to_node[node_num][2];

              // The local
              // orientation of the
              // edge coming into the
              // node.
              const EdgeOrientation sign0 = ElementInfo::edge_to_node_orient[node_num][0];
              const EdgeOrientation sign1 = ElementInfo::edge_to_node_orient[node_num][1];
              const EdgeOrientation sign2 = ElementInfo::edge_to_node_orient[node_num][2];

              // Add one to the total
              // for each edge
              // pointing in
              Assert (local_edge_orientation[e0] != unoriented_edge,
                      ExcInternalError());
              Assert (local_edge_orientation[e1] != unoriented_edge,
                      ExcInternalError());
              Assert (local_edge_orientation[e2] != unoriented_edge,
                      ExcInternalError());

              const unsigned int
              total  = (((local_edge_orientation[e0] == sign0) ? 1 : 0)
                        +((local_edge_orientation[e1] == sign1) ? 1 : 0)
                        +((local_edge_orientation[e2] == sign2) ? 1 : 0));

              if (total == 3)
                {
                  Assert (perm_num == numbers::invalid_unsigned_int,
                          ExcGridOrientError("More than one node with 3 incoming "
                                             "edges found in curent hex."));
                  perm_num = node_num;
                }
            }
          // We should now have a
          // valid permutation number
          Assert (perm_num != numbers::invalid_unsigned_int,
                  ExcGridOrientError("No node having 3 incoming edges found in curent hex."));

          // So use the appropriate
          // rotation to get the new
          // cube
          unsigned int temp[8];
          for (unsigned int v=0; v<8; ++v)
            temp[v] = the_cell.nodes[CubePermutations[perm_num][v]];
          for (unsigned int v=0; v<8; ++v)
            the_cell.nodes[v] = temp[v];
        }
    }
  } // namespace GridReordering3d
} // namespace internal



template<>
void
GridReordering<3>::reorder_cells (std::vector<CellData<3> > &cells,
                                  const bool use_new_style_ordering)
{
  Assert (cells.size() != 0,
          ExcMessage("List of elements to orient must have at least one cell"));

  // if necessary, convert to old-style format
  if (use_new_style_ordering)
    reorder_new_to_old_style(cells);

  // create a backup to use if GridReordering
  // was not successful
  std::vector<CellData<3> > backup=cells;

  // This does the real work
  const bool success=
    internal::GridReordering3d::Orienter::orient_mesh (cells);

  // if reordering was not successful use
  // original connectivity, otherwise do
  // nothing (i.e. use the reordered
  // connectivity)
  if (!success)
    cells=backup;

  // and convert back if necessary
  if (use_new_style_ordering)
    reorder_old_to_new_style(cells);
}



template<>
void
GridReordering<3>::invert_all_cells_of_negative_grid(
  const std::vector<Point<3> > &all_vertices,
  std::vector<CellData<3> > &cells)
{
  unsigned int vertices_lex[GeometryInfo<3>::vertices_per_cell];
  unsigned int n_negative_cells=0;
  for (unsigned int cell_no=0; cell_no<cells.size(); ++cell_no)
    {
      // GridTools::cell_measure
      // requires the vertices to be
      // in lexicographic ordering
      for (unsigned int i=0; i<GeometryInfo<3>::vertices_per_cell; ++i)
        vertices_lex[GeometryInfo<3>::ucd_to_deal[i]]=cells[cell_no].vertices[i];
      if (GridTools::cell_measure<3>(all_vertices, vertices_lex) < 0)
        {
          ++n_negative_cells;
          // reorder vertices: swap front and back face
          for (unsigned int i=0; i<4; ++i)
            std::swap(cells[cell_no].vertices[i], cells[cell_no].vertices[i+4]);

          // check whether the
          // resulting cell is now ok.
          // if not, then the grid is
          // seriously broken and
          // should be sticked into the
          // bin
          for (unsigned int i=0; i<GeometryInfo<3>::vertices_per_cell; ++i)
            vertices_lex[GeometryInfo<3>::ucd_to_deal[i]]=cells[cell_no].vertices[i];
          AssertThrow(GridTools::cell_measure<3>(all_vertices, vertices_lex) > 0,
                      ExcInternalError());
        }
    }

  // We assume that all cells of a
  // grid have either positive or
  // negative volumes but not both
  // mixed. Although above reordering
  // might work also on single cells,
  // grids with both kind of cells
  // are very likely to be
  // broken. Check for this here.
  AssertThrow(n_negative_cells==0 || n_negative_cells==cells.size(),
              ExcMessage("While sorting the cells that will be passed for "
                         "creating a Triangulation object, deal.II found that "
                         "some but not all cells have a negative volume. (If "
                         "all cells had a negative volume, they would simply "
                         "all have been inverted.) This usually happens in "
                         "hand-generated meshes if one accidentally uses an "
                         "incorrect convention for ordering the vertices in "
                         "one or more cells; in that case, you may want to "
                         "double check that you specified the vertex indices "
                         "in their correct order. If you are reading a mesh "
                         "that was created by a mesh generator, then this "
                         "exception indicates that some of the cells created "
                         "are so badly distorted that their volume becomes "
                         "negative; this commonly occurs at complex geometric "
                         "features, and you may see if the problem can be "
                         "fixed by playing with the parameters that control "
                         "mesh properties in your mesh generator, such as "
                         "the number of cells, the mesh density, etc."));
}


DEAL_II_NAMESPACE_CLOSE

