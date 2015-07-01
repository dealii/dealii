// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2016 by the deal.II authors
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

#ifndef dealii__grid_reordering_internal_h
#define dealii__grid_reordering_internal_h


#include <deal.II/base/config.h>
#include <deal.II/grid/tria.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  /**
   * Implement the algorithm described in the documentation of the
   * GridReordering<3> class.
   *
   * @author Michael Anderson, 2003
   */
  namespace GridReordering3d
  {
    /**
     * A structure indicating the direction of an edge. In the implementation
     * file, we define three objects, <tt>unoriented_edge</tt>,
     * <tt>forward_edge</tt>, and <tt>backward_edge</tt>, that denote whether
     * an edge has already been oriented, whether it is in standard
     * orientation, or whether it has reverse direction. The state that each
     * of these objects encode is stored in the <tt>orientation</tt> member
     * variable -- we would really need only three such values, which we pick
     * in the implementation file, and make sure when we compare such objects
     * that only these three special values are actually used.
     *
     * The reason for this way of implementing things is as follows. Usually,
     * such a property would be implemented as an enum. However, in the
     * previous implementation, a signed integer was used with unoriented=0,
     * forward=+1, and backward=-1. A number of operations, such as equality
     * of ordered edges were mapped to checking whether the product of two
     * edge orientations equals +1. Such arithmetic isn't always portable and
     * sometimes flagged when using -ftrapv with gcc. Using this class instead
     * makes sure that there isn't going to be any arithmetic going on on edge
     * orientations, just comparisons for equality or inequality.
     *
     * @author Wolfgang Bangerth, 2005
     */
    struct EdgeOrientation
    {
      /**
       * A value indicating the orientation.
       */
      char orientation;

      /**
       * Comparison operator.
       */
      bool operator == (const EdgeOrientation &edge_orientation) const;

      /**
       * Comparison operator.
       */
      bool operator != (const EdgeOrientation &edge_orientation) const;
    };

    /**
     * During building the connectivity information we don't need all the
     * heavy duty information about edges that we will need later. So we can
     * save memory and time by using a light-weight class for edges. It stores
     * the two vertices, but no direction, so we make the optimization to
     * store the vertex number in sorted order to allow for easier comparison
     * of edge objects.
     */
    struct CheapEdge
    {
      /**
       * The first node
       */
      const unsigned int node0;

      /**
       * The second node
       */
      const unsigned int node1;

      /**
       * Constructor. Take the vertex numbers and store them sorted.
       */
      CheapEdge (const unsigned int n0,
                 const unsigned int n1);

      /**
       * Need a partial ordering for sorting algorithms
       */
      bool operator< (const CheapEdge &e2) const;
    };



    /**
     * A connectivity and orientation aware edge class.
     */
    struct Edge
    {
      /**
       * Simple constructor
       */
      Edge (const unsigned int n0,
            const unsigned int n1);

      /**
       * The IDs for the end nodes
       */
      unsigned int nodes[2];

      /**
       * Whether the edge has not already been oriented, points from node 0 to
       * node 1, or the reverse. The initial state of this flag is unoriented.
       */
      EdgeOrientation orientation_flag;

      /**
       * Used to determine which "sheet" or equivalence class of parallel
       * edges the edge falls in when oriented. numbers::invalid_unsigned_int
       * means not yet decided. This is also the default value after
       * construction. Each edge will later be assigned an index greater than
       * zero.
       */
      unsigned int group;

      /**
       * Indices of neighboring cubes.
       */
      std::vector<unsigned int> neighboring_cubes;
    };

    /**
     * A connectivity and orientation aware cell.
     *
     * The connectivity of the cell is not contained within. (This was for
     * flexibility in using deal.II's ordering of edges or the XDA format
     * etc.) For this information we need the ElemInfo class.
     *
     * One thing we do know is that the first four edges in the edge class are
     * parallel, as are the second four, and the third four.
     *
     * TODO: Need to move connectivity information out of cell and into edge.
     */
    struct Cell
    {
      /**
       * Default Constructor
       */
      Cell ();

      /**
       * The IDs for each of the edges.
       */
      unsigned int edges[GeometryInfo<3>::lines_per_cell];

      /**
       * The IDs for each of the nodes.
       */
      unsigned int nodes[GeometryInfo<3>::vertices_per_cell];

      /**
       * Which way do the edges point.  Whether node 0 of the edge is the base
       * of the edge in local element (1) or node 1 is the base (-1).
       */
      EdgeOrientation local_orientation_flags[GeometryInfo<3>::lines_per_cell];

      /**
       * An internal flag used to determine whether the cell is in the queue
       * of cells to be oriented in the current sheet.
       */
      bool waiting_to_be_processed;
    };


    /**
     * This holds all the pieces for orientation together.
     *
     * Contains lists of nodes, edges and cells.  As well as the information
     * about how they all connect together.
     */
    class Mesh
    {
    public:
      /**
       * Default Constructor
       */
      Mesh (const std::vector<CellData<3> > &incubes);

      /**
       * Export the data of this object to the deal.II format that the
       * Triangulation class wants as input.
       */
      void
      export_to_deal_format (std::vector<CellData<3> > &outcubes) const;

    private:
      /**
       * The list of edges
       */
      std::vector<Edge> edge_list;

      /**
       * The list of cells
       */
      std::vector<Cell> cell_list;

      /**
       * Check whether every cell in the mesh is sensible.
       */
      void sanity_check() const;

      /**
       * Given the cell list, build the edge list and all the connectivity
       * information and other stuff that we will need later.
       */
      void build_connectivity ();

      /**
       * Unimplemented private copy constructor to disable it.
       */
      Mesh (const Mesh &);

      /**
       * Unimplemented private assignment operator to disable it.
       */
      Mesh &operator=(const Mesh &);

      /**
       * Checks that each edge going into a node is correctly set up.
       */
      void sanity_check_node (const Cell        &cell,
                              const unsigned int local_node_num) const;

      /**
       * Let the orienter access out private fields.
       */
      friend class Orienter;
    };


    /**
     * The class that orients the edges of a triangulation in 3d. The member
     * variables basically only store the present state of the algorithm.
     */
    class Orienter
    {
    public:
      /**
       * Orient the given mesh. Creates an object of the present type and lets
       * that toil away at the task.
       *
       * This function is the single entry point to the functionality of this
       * class.
       *
       * Return, whether a consistent orientation of lines was possible for
       * the given mesh.
       */
      static
      bool
      orient_mesh (std::vector<CellData<3> > &incubes);

    private:
      /**
       * Internal representation of the given list of cells, including
       * connectivity information and the like.
       */
      Mesh mesh;

      /**
       * The cube we're looking at presently.
       */
      unsigned int cur_posn;

      /**
       * We have fully oriented all cubes before this one.
       */
      unsigned int marker_cube;

      /**
       * The index of the sheet or equivalence class we are presently
       * processing.
       */
      unsigned int cur_edge_group;

      /**
       * Indices of the cells to be processed within the present sheet. If a
       * cell is being processed presently, it is taken from this list.
       */
      std::vector<int> sheet_to_process;


      /**
       * Which edges of the current cell have been oriented during the current
       * iteration. Is reset when moving on to the next cube.
       */
      bool edge_orient_array[12];

      /**
       * Constructor. Take a list of cells and set up the internal data
       * structures of the mesh member variable.
       *
       * Since it is private, the only entry point of this class is the static
       * function orient_mesh().
       */
      Orienter (const std::vector<CellData<3> > &incubes);

      /**
       * Orient all the edges of a mesh.
       *
       * Return, whether this action was carried out successfully.
       */
      bool orient_edges ();

      /**
       * Given oriented edges, rotate the cubes so that the edges are in
       * standard direction.
       */
      void orient_cubes ();

      bool get_next_unoriented_cube ();

      /**
       * Return whether the cell with cell number @p cell_num is fully
       * oriented.
       */
      bool is_oriented (const unsigned int cell_num) const;

      bool orient_edges_in_current_cube ();
      bool orient_edge_set_in_current_cube (const unsigned int edge_set);
      bool orient_next_unoriented_edge ();

      /**
       * Return whether the cell is consistently oriented at present (i.e.
       * only considering those edges that are already oriented. This is a
       * sanity check that should be called from inside an assert macro.
       */
      bool cell_is_consistent (const unsigned int cell_num) const;


      void get_adjacent_cubes ();
      bool get_next_active_cube ();
    };
  }  // namespace GridReordering3d
}  // namespace internal


DEAL_II_NAMESPACE_CLOSE

#endif
