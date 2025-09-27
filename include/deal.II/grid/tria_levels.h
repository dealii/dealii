// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1998 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_tria_levels_h
#define dealii_tria_levels_h


#include <deal.II/base/config.h>

#include <deal.II/base/point.h>

#include <deal.II/grid/reference_cell.h>
#include <deal.II/grid/tria_objects.h>
#include <deal.II/grid/tria_objects_orientations.h>

#include <boost/serialization/utility.hpp>

#include <cstdint>
#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace TriangulationImplementation
  {
    /**
     * Store all information which belongs to one level of the multilevel
     * hierarchy.
     *
     * In TriaLevel, all cell data is stored which is not dependent on the
     * dimension, e.g. a field to store the refinement flag for the cells
     * (what a cell actually is is declared elsewhere), etc. See also
     * TriaObjects for non level-oriented data.
     *
     * There is another field, which may fit in here, namely the material data
     * (for cells) or the boundary indicators (for faces), but since we need
     * for a line or quad either boundary information or material data, we
     * store them with the lines and quads rather than with the common data.
     * Likewise, in 3d, we need boundary indicators for lines and quads (we
     * need to know how to refine a line if the two adjacent faces have
     * different boundary indicators), and material data for cells.
     */
    class TriaLevel
    {
    public:
      /**
       * Constructor.
       *
       * @param dim Dimension of the Triangulation.
       */
      TriaLevel(const unsigned int dim)
        : dim(dim)
        , cells(dim)
      {}

      /**
       * Default constructor (needed by Boost).
       */
      TriaLevel()
        : dim(numbers::invalid_unsigned_int)
        , cells(numbers::invalid_unsigned_int)
      {}

      /**
       * Dimension of the Triangulation.
       */
      unsigned int dim;

      /**
       * @p RefinementCase<dim>::Type flags for the cells to be refined with
       * or not (RefinementCase<dim>::no_refinement). The meaning what a cell
       * is, is dimension specific, therefore also the length of this vector
       * depends on the dimension: in one dimension, the length of this vector
       * equals the length of the @p lines vector, in two dimensions that of
       * the @p quads vector, etc.
       */
      std::vector<std::uint8_t> refine_flags;

      /**
       * @p IsotropicRefinementChoices direction for the tetrahedral cells to be
       * refined with
       * or choose the optimal way @p ReferenceCell::isotropic_refinement.
       */
      std::vector<std::uint8_t> refine_choice;

      /**
       * Same meaning as the one above, but specifies whether a cell must be
       * coarsened.
       */
      std::vector<bool> coarsen_flags;


      /**
       * An integer that, for every active cell, stores the how many-th active
       * cell this is. For non-active cells, this value is unused and set to
       * an invalid value.
       */
      std::vector<unsigned int> active_cell_indices;

      /**
       * Global cell index of each active cell.
       */
      std::vector<types::global_cell_index> global_active_cell_indices;

      /**
       * Global cell index of each cell on the given level.
       */
      std::vector<types::global_cell_index> global_level_cell_indices;

      /**
       * Levels and indices of the neighbors of the cells. Convention is, that
       * the neighbors of the cell with index @p i are stored in the fields
       * following $i*(2*real\_space\_dimension)$, e.g. in one spatial
       * dimension, the neighbors of cell 0 are stored in
       * <tt>neighbors[0]</tt> and <tt>neighbors[1]</tt>, the neighbors of
       * cell 1 are stored in <tt>neighbors[2]</tt> and <tt>neighbors[3]</tt>,
       * and so on.
       *
       * In neighbors, <tt>neighbors[i].first</tt> is the level, while
       * <tt>neighbors[i].second</tt> is the index of the neighbor.
       *
       * If a neighbor does not exist (cell is at the boundary),
       * <tt>level=index=-1</tt> is set.
       *
       * <em>Conventions:</em> The @p ith neighbor of a cell is the one which
       * shares the @p ith face (@p Line in 2d, @p Quad in 3d) of this cell.
       *
       * The neighbor of a cell has at most the same level as this cell, i.e.
       * it may or may not be refined.
       *
       * In one dimension, a neighbor may have any level less or equal the
       * level of this cell. If it has the same level, it may be refined an
       * arbitrary number of times, but the neighbor pointer still points to
       * the cell on the same level, while the neighbors of the children of
       * the neighbor may point to this cell or its children.
       *
       * In two and more dimensions, the neighbor is either on the same level
       * and refined (in which case its children have neighbor pointers to
       * this cell or its direct children), unrefined on the same level or one
       * level down (in which case its neighbor pointer points to the parent
       * cell of this cell).
       */
      std::vector<std::pair<int, int>> neighbors;

      /**
       * One integer per cell to store which subdomain it belongs to. This
       * field is most often used in parallel computations, where it denotes
       * which processor shall work on/owns the cells with a given subdomain
       * number.
       *
       * This number is only used on active cells.
       */
      std::vector<types::subdomain_id> subdomain_ids;

      /**
       * The subdomain id used on each level for parallel multigrid.
       *
       * In contrast to the subdomain_id, this number is also used on inactive
       * cells once the mesh has been partitioned also on the lower levels of
       * the multigrid hierarchy.
       */
      std::vector<types::subdomain_id> level_subdomain_ids;

      /**
       * One integer for every consecutive pair of cells to store which index
       * their parent has.
       *
       * (We store this information once for each pair of cells since every
       * refinement, isotropic or anisotropic, and in any space dimension,
       * always creates children in multiples of two, so there is no need to
       * store the parent index for every cell.)
       */
      std::vector<int> parents;

      /**
       * One bool per cell to indicate the direction of the normal true:  use
       * orientation from vertex false: revert the orientation. See
       * @ref GlossDirectionFlag.
       *
       * This is only used for codim==1 meshes.
       */
      std::vector<bool> direction_flags;

      /**
       * The object containing the data on lines and related functions
       */
      TriaObjects cells;

      /**
       * For edges, we enforce a standard convention that opposite
       * edges should be parallel. Now, that's enforceable in most
       * cases, and we have code that makes sure that if a mesh allows
       * this to happen, that we have this convention. We also know
       * that it is always possible to have opposite faces have
       * parallel normal vectors. (For both things, see the paper by
       * Agelek, Anderson, Bangerth, Barth in the ACM Transactions on
       * Mathematical Software mentioned in the documentation of the
       * GridTools::consistently_order_cells() function.)
       *
       * The problem is that we originally had another condition, namely that
       * faces 0, 2 and 4 have normals that point into the cell, while the
       * other faces have normals that point outward. It turns out that this
       * is not always possible. In effect, we have to store whether the
       * normal vector of each face of each cell follows this convention or
       * not. If this is so, then this variable stores a @p true value,
       * otherwise a @p false value.
       *
       * In effect, this field has <code>6*n_cells</code> elements, being the
       * number of cells times the six faces each has.
       *
       * @note This array is only used in dim == 2 or dim == 3: for dim == 1
       * meshes consist purely of lines which are always consistently oriented.
       */
      TriaObjectsOrientations face_orientations;

      /**
       * Reference cell type of each cell.
       *
       * @note Used only for dim=2 and dim=3.
       */
      std::vector<ReferenceCell> reference_cell;

      /**
       * A cache for the vertex indices of the cells (`structdim == dim`), in
       * order to more quickly retrieve these frequently accessed quantities.
       * For simplified addressing, the information is indexed by the maximum
       * number of vertices possible for a cell (quadrilateral/hexahedron).
       */
      std::vector<unsigned int> cell_vertex_indices_cache;

      /**
       * Determine an estimate for the memory consumption (in bytes) of this
       * object.
       */
      std::size_t
      memory_consumption() const;

      /**
       * Read or write the data of this object to or from a stream for the
       * purpose of serialization using the [BOOST serialization
       * library](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html).
       */
      template <class Archive>
      void
      serialize(Archive &ar, const unsigned int version);
    };


    template <class Archive>
    void
    TriaLevel::serialize(Archive &ar, const unsigned int)
    {
      ar &dim;

      ar &refine_flags &coarsen_flags;

      // do not serialize `active_cell_indices` and `vertex_indices_cache`
      // here. instead of storing them to the stream and re-reading them again
      // later, we just rebuild them in Triangulation::load()

      ar &neighbors;
      ar &subdomain_ids;
      ar &level_subdomain_ids;
      ar &parents;
      ar &direction_flags;
      ar &cells;
      ar &face_orientations;
      ar &reference_cell;
    }

  } // namespace TriangulationImplementation
} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
