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
#include <deal.II/grid/grid_tools.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/timer.h>

#include <algorithm>
#include <set>
#include <iostream>
#include <fstream>
#include <functional>

DEAL_II_NAMESPACE_OPEN


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
    static const unsigned int n_other_parallel_edges = (1<<(dim-1)) - 1;
    static const unsigned int parallel_edges[GeometryInfo<dim>::lines_per_cell][n_other_parallel_edges];
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


  /**
   * A structure that store the index of a cell and, crucially, how a
   * given edge relates to this cell.
   */
  struct AdjacentCell
  {
    /**
     * Default constructor. Initialize the fields with invalid values.
     */
    AdjacentCell ()
      :
      cell_index (numbers::invalid_unsigned_int),
      edge_within_cell (numbers::invalid_unsigned_int)
    {}

    /**
     * Constructor. Initialize the fields with the given values.
     */
    AdjacentCell (const unsigned int cell_index,
                  const unsigned int edge_within_cell)
      :
      cell_index (cell_index),
      edge_within_cell (edge_within_cell)
    {}


    unsigned int cell_index;
    unsigned int edge_within_cell;
  };



  template <int dim> class AdjacentCells;

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
    typedef const AdjacentCell *const_iterator;

    /**
     * Add the given cell to the collection of cells adjacent to
     * the edge this object corresponds to. Since we are covering
     * the 2d case, the set of adjacent cells currently
     * represented by this object must have either zero or
     * one element already, since we can not add more than two
     * adjacent cells for each edge.
     */
    void push_back (const AdjacentCell &adjacent_cell)
    {
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
     * Return an iterator to the first valid cell stored as adjacent to the
     * edge represented by the current object.
     */
    const_iterator begin () const
    {
      return &adjacent_cells[0];
    }


    /**
     * Return an iterator to the element past the last valid cell stored
     * as adjacent to the edge represented by the current object.
     * @return
     */
    const_iterator end () const
    {
      // check whether the current object stores zero, one, or two
      // adjacent cells, and use this to point to the element past the
      // last valid one
      if (adjacent_cells[0].cell_index == numbers::invalid_unsigned_int)
        return &adjacent_cells[0];
      else if (adjacent_cells[1].cell_index == numbers::invalid_unsigned_int)
        return &adjacent_cells[0]+1;
      else
        return &adjacent_cells[0]+2;
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
     * Default constructor. Creates an invalid edge.
     */
    Edge ()
      :
      orientation_status (not_oriented)
    {
      for (unsigned int i=0; i<2; ++i)
        vertex_indices[i] = numbers::invalid_unsigned_int;
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
     * Construct a Cell object from a CellData object. Also take a
     * (sorted) list of edges and to point the edges of the current
     * object into this list of edges.
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



  template <int dim> class EdgeDeltaSet;

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
    typedef const unsigned int *const_iterator;

    /**
     * Default constructor. Initialize both slots as unused, corresponding
     * to an empty set.
     */
    EdgeDeltaSet ()
    {
      edge_indices[0] = edge_indices[1] = numbers::invalid_unsigned_int;
    }


    /**
     * Delete the elements of the set by marking both slots as unused.
     */
    void clear ()
    {
      edge_indices[0] = edge_indices[1] = numbers::invalid_unsigned_int;
    }

    /**
     * Insert one element into the set. This will fail if the set already
     * has two elements.
     */
    void insert (const unsigned int edge_index)
    {
      if (edge_indices[0] == numbers::invalid_unsigned_int)
        edge_indices[0] = edge_index;
      else
        {
          Assert (edge_indices[1] == numbers::invalid_unsigned_int,
                  ExcInternalError());
          edge_indices[1] = edge_index;
        }
    }


    /**
     * Return an iterator pointing to the first element of the set.
     */
    const_iterator begin () const
    {
      return &edge_indices[0];
    }


    /**
     * Return an iterator pointing to the element past the last used one.
     */
    const_iterator end () const
    {
      // check whether the current object stores zero, one, or two
      // indices, and use this to point to the element past the
      // last valid one
      if (edge_indices[0] == numbers::invalid_unsigned_int)
        return &edge_indices[0];
      else if (edge_indices[1] == numbers::invalid_unsigned_int)
        return &edge_indices[0]+1;
      else
        return &edge_indices[0]+2;
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
        edge_list.emplace_back (cells[i], l);

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
        cell_list.emplace_back (cells[i], edges);

        // then also inform the edges that they are adjacent
        // to the current cell, and where within this cell
        for (unsigned int l=0; l<GeometryInfo<dim>::lines_per_cell; ++l)
          edges[cell_list.back().edge_indices[l]].adjacent_cells.push_back (AdjacentCell (i, l));
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
  get_next_unoriented_cell(const std::vector<Cell<dim> > &cells,
                           const std::vector<Edge<dim> > &edges,
                           const unsigned int             current_cell)
  {
    for (unsigned int c=current_cell; c<cells.size(); ++c)
      for (unsigned int l=0; l<GeometryInfo<dim>::lines_per_cell; ++l)
        if (edges[cells[c].edge_indices[l]].orientation_status == Edge<dim>::not_oriented)
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
  orient_one_set_of_parallel_edges (const std::vector<Cell<dim> > &cells,
                                    std::vector<Edge<dim> >       &edges,
                                    const unsigned int             cell,
                                    const unsigned int             local_edge)
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
    if (edges[cells[cell].edge_indices[local_edge]].vertex_indices[0]
        ==
        cells[cell].vertex_indices[GeometryInfo<dim>::line_to_cell_vertices (local_edge, 0)])
      // orient initial edge *opposite* to the way it is in the cell
      // (see above for the reason)
      edges[cells[cell].edge_indices[local_edge]].orientation_status = (dim == 2 ?
          Edge<dim>::backward :
          Edge<dim>::forward);
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
        edges[cells[cell].edge_indices[local_edge]].orientation_status = (dim == 2 ?
            Edge<dim>::forward :
            Edge<dim>::backward);
      }

    // walk outward from the given edge as described in
    // the algorithm in the paper that documents all of
    // this
    //
    // note that in 2d, each of the Deltas can at most
    // contain two elements, whereas in 3d it can be arbitrarily many
    EdgeDeltaSet<dim> Delta_k;
    EdgeDeltaSet<dim> Delta_k_minus_1;
    Delta_k_minus_1.insert (cells[cell].edge_indices[local_edge]);

    while (Delta_k_minus_1.begin() != Delta_k_minus_1.end())   // while set is not empty
      {
        Delta_k.clear ();

        for (typename EdgeDeltaSet<dim>::const_iterator delta = Delta_k_minus_1.begin();
             delta != Delta_k_minus_1.end(); ++delta)
          {
            Assert (edges[*delta].orientation_status != Edge<dim>::not_oriented,
                    ExcInternalError());

            // now go through the cells adjacent to this edge
            for (typename AdjacentCells<dim>::const_iterator
                 adjacent_cell = edges[*delta].adjacent_cells.begin();
                 adjacent_cell != edges[*delta].adjacent_cells.end(); ++adjacent_cell)
              {
                const unsigned int K = adjacent_cell->cell_index;
                const unsigned int delta_is_edge_in_K = adjacent_cell->edge_within_cell;

                // figure out the direction of delta with respect to the cell K
                // (in the orientation in which the user has given it to us)
                const unsigned int first_edge_vertex
                  = (edges[*delta].orientation_status == Edge<dim>::forward
                     ?
                     edges[*delta].vertex_indices[0]
                     :
                     edges[*delta].vertex_indices[1]);
                const unsigned int first_edge_vertex_in_K
                  = cells[K].vertex_indices[GeometryInfo<dim>::line_to_cell_vertices(delta_is_edge_in_K, 0)];
                Assert (first_edge_vertex == first_edge_vertex_in_K
                        ||
                        first_edge_vertex == cells[K].vertex_indices[GeometryInfo<dim>::line_to_cell_vertices(delta_is_edge_in_K, 1)],
                        ExcInternalError());

                // now figure out which direction the each of the "opposite" edges
                // needs to be oriented into.
                for (unsigned int o_e=0; o_e<ParallelEdges<dim>::n_other_parallel_edges; ++o_e)
                  {
                    // get the index of the opposite edge and select which its first
                    // vertex needs to be based on how the current edge is oriented
                    // in the current cell
                    const unsigned int opposite_edge
                      = cells[K].edge_indices[ParallelEdges<dim>::parallel_edges[delta_is_edge_in_K][o_e]];
                    const unsigned int first_opposite_edge_vertex
                      =  cells[K].vertex_indices[GeometryInfo<dim>::line_to_cell_vertices(
                                                   ParallelEdges<dim>::parallel_edges[delta_is_edge_in_K][o_e],
                                                   (first_edge_vertex == first_edge_vertex_in_K
                                                    ?
                                                    0
                                                    :
                                                    1))];

                    // then determine the orientation of the edge based on
                    // whether the vertex we want to be the edge's first
                    // vertex is already the first vertex of the edge, or
                    // whether it points in the opposite direction
                    const typename Edge<dim>::OrientationStatus opposite_edge_orientation
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
                        // Delta_k
                        edges[opposite_edge].orientation_status = opposite_edge_orientation;
                        Delta_k.insert (opposite_edge);
                      }
                    else
                      {
                        // this opposite edge has already been oriented. it should be
                        // consistent with the current one in 2d, while in 3d it may in fact
                        // be mis-oriented, and in that case the mesh will not be
                        // orientable. indicate this by throwing an exception that we can
                        // catch further up; this has the advantage that we can propagate
                        // through a couple of functions without having to do error
                        // checking and without modifying the 'cells' array that the
                        // user gave us
                        if (dim == 2)
                          {
                            Assert (edges[opposite_edge].orientation_status == opposite_edge_orientation,
                                    ExcMeshNotOrientable());
                          }
                        else if (dim == 3)
                          {
                            if (edges[opposite_edge].orientation_status != opposite_edge_orientation)
                              throw ExcMeshNotOrientable ();
                          }
                        else
                          Assert (false, ExcNotImplemented());
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
  rotate_cell (const std::vector<Cell<dim> > &cell_list,
               const std::vector<Edge<dim> > &edge_list,
               const unsigned int             cell_index,
               std::vector<CellData<dim> >   &raw_cells)
  {
    // find the first vertex of the cell. this is the vertex where dim edges
    // originate, so for each of the edges record which the starting vertex is
    unsigned int starting_vertex_of_edge[GeometryInfo<dim>::lines_per_cell];
    for (unsigned int e=0; e<GeometryInfo<dim>::lines_per_cell; ++e)
      {
        Assert (edge_list[cell_list[cell_index].edge_indices[e]].orientation_status
                != Edge<dim>::not_oriented,
                ExcInternalError());
        if (edge_list[cell_list[cell_index].edge_indices[e]].orientation_status == Edge<dim>::forward)
          starting_vertex_of_edge[e] = edge_list[cell_list[cell_index].edge_indices[e]].vertex_indices[0];
        else
          starting_vertex_of_edge[e] = edge_list[cell_list[cell_index].edge_indices[e]].vertex_indices[1];
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
        if ((starting_vertex_of_edge[0] == starting_vertex_of_edge[2])
            ||
            (starting_vertex_of_edge[0] == starting_vertex_of_edge[3]))
          origin_vertex_of_cell = starting_vertex_of_edge[0];
        else if ((starting_vertex_of_edge[1] == starting_vertex_of_edge[2])
                 ||
                 (starting_vertex_of_edge[1] == starting_vertex_of_edge[3]))
          origin_vertex_of_cell = starting_vertex_of_edge[1];
        else
          Assert (false, ExcInternalError());

        break;
      }

      case 3:
      {
        // one could probably do something similar in 3d, but that seems
        // more complicated than one wants to write down. just go
        // through the list of possible starting vertices and check
        for (origin_vertex_of_cell=0;
             origin_vertex_of_cell<GeometryInfo<dim>::vertices_per_cell;
             ++origin_vertex_of_cell)
          if (std::count (&starting_vertex_of_edge[0],
                          &starting_vertex_of_edge[0]+GeometryInfo<dim>::lines_per_cell,
                          cell_list[cell_index].vertex_indices[origin_vertex_of_cell])
              == dim)
            break;
        Assert (origin_vertex_of_cell < GeometryInfo<dim>::vertices_per_cell,
                ExcInternalError());

        break;
      }

      default:
        Assert (false, ExcNotImplemented());
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
            raw_cells[cell_index].vertices[0] = raw_cells[cell_index].vertices[1];
            raw_cells[cell_index].vertices[1] = raw_cells[cell_index].vertices[3];
            raw_cells[cell_index].vertices[3] = raw_cells[cell_index].vertices[2];
            raw_cells[cell_index].vertices[2] = tmp;
          }
        break;
      }

      case 3:
      {
        // in 3d, the situation is a bit more complicated. from above, we
        // now know which vertex is at the origin (because 3 edges originate
        // from it), but that still leaves 3 possible rotations of the cube.
        // the important realization is that we can choose any of them:
        // in all 3 rotations, all edges originate from the one vertex,
        // and that fixes the directions of all 12 edges in the cube because
        // these 3 cover all 3 equivalence classes! consequently, we can
        // select an arbitrary one among the permutations -- for
        // example the following ones:
        static const unsigned int cube_permutations[8][8] =
        {
          {0,1,2,3,4,5,6,7},
          {1,5,3,7,0,4,2,6},
          {2,6,0,4,3,7,1,5},
          {3,2,1,0,7,6,5,4},
          {4,0,6,2,5,1,7,3},
          {5,4,7,6,1,0,3,2},
          {6,7,4,5,2,3,0,1},
          {7,3,5,1,6,2,4,0}
        };

        unsigned int temp_vertex_indices[GeometryInfo<dim>::vertices_per_cell];
        for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
          temp_vertex_indices[v]
            = raw_cells[cell_index].vertices[cube_permutations[origin_vertex_of_cell][v]];
        for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
          raw_cells[cell_index].vertices[v] = temp_vertex_indices[v];

        break;
      }

      default:
      {
        Assert (false, ExcNotImplemented());
      }
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
    unsigned int next_cell_with_unoriented_edge = 0;
    while ((next_cell_with_unoriented_edge = get_next_unoriented_cell(cell_list,
                                             edge_list,
                                             next_cell_with_unoriented_edge)) !=
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


  // overload of the function above for 1d -- there is nothing
  // to orient in that case
  void reorient (std::vector<CellData<1> > &)
  {}
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
  reorder_new_to_old_style (std::vector<CellData<1> > &)
  {}


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
  reorder_old_to_new_style (std::vector<CellData<1> > &)
  {}


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



template <int dim, int spacedim>
void
GridReordering<dim,spacedim>::reorder_cells (std::vector<CellData<dim> > &cells,
                                             const bool use_new_style_ordering)
{
  Assert (cells.size() != 0,
          ExcMessage("List of elements to orient must have at least one cell"));

  // there is nothing for us to do in 1d
  if (dim == 1)
    return;

  // if necessary, convert to new-style format
  if (use_new_style_ordering == false)
    reorder_old_to_new_style(cells);

  // check if grids are already consistent. if so, do
  // nothing. if not, then do the reordering
  if (!is_consistent (cells))
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

  // and convert back if necessary
  if (use_new_style_ordering == false)
    reorder_new_to_old_style(cells);
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
GridReordering<1,2>::invert_all_cells_of_negative_grid(const std::vector<Point<2> > &,
                                                       std::vector<CellData<1> > &)
{
  // nothing to be done in 1d
}



template<>
void
GridReordering<1,3>::invert_all_cells_of_negative_grid(const std::vector<Point<3> > &,
                                                       std::vector<CellData<1> > &)
{
  // nothing to be done in 1d
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



template<>
void
GridReordering<3>::invert_all_cells_of_negative_grid(const std::vector<Point<3> > &all_vertices,
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



/* ------------------------ explicit instantiations ------------------- */
template class GridReordering<1,1>;
template class GridReordering<1,2>;
template class GridReordering<1,3>;
template class GridReordering<2,2>;
template class GridReordering<2,3>;
template class GridReordering<3,3>;

DEAL_II_NAMESPACE_CLOSE
