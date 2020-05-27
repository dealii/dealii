// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_tria_levels_h
#define dealii_tria_levels_h


#include <deal.II/base/config.h>

#include <deal.II/base/point.h>

#include <deal.II/grid/tria_object.h>
#include <deal.II/grid/tria_objects.h>

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
     *
     * @author Wolfgang Bangerth, Guido Kanschat, 1998, 2007
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
       * Levels and indices of the neighbors of the cells. Convention is, that
       * the neighbors of the cell with index @p i are stored in the fields
       * following $i*(2*real\_space\_dimension)$, e.g. in one spatial
       * dimension, the neighbors of cell 0 are stored in
       * <tt>neighbors[0]</tt> and <tt>neighbors[1]</tt>, the neighbors of
       * cell 1 are stored in <tt>neighbors[2]</tt> and <tt>neighbors[3]</tt>,
       * and so on.
       *
       * In neighbors, <tt>neighbors[i].first</tt> is the level, while
       * <tt>neighbors[i].first</tt> is the index of the neighbor.
       *
       * If a neighbor does not exist (cell is at the boundary),
       * <tt>level=index=-1</tt> is set.
       *
       * <em>Conventions:</em> The @p ith neighbor of a cell is the one which
       * shares the @p ith face (@p Line in 2D, @p Quad in 3D) of this cell.
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
       * level down (in which case its neighbor pointer points to the mother
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
       * GridReordering class.)
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
       * @note Only needed for dim=3.
       */
      std::vector<unsigned char> face_orientations;

      /**
       * Reserve enough space to accommodate @p total_cells cells on this
       * level. Since there are no @p used flags on this level, you have to
       * give the total number of cells, not only the number of newly to
       * accommodate ones, like in the <tt>TriaLevel<N>::reserve_space</tt>
       * functions, with <tt>N>0</tt>.
       *
       * Since the number of neighbors per cell depends on the dimensions, you
       * have to pass that additionally.
       */

      void
      reserve_space(const unsigned int total_cells,
                    const unsigned int dimension,
                    const unsigned int space_dimension);

      /**
       * Check the memory consistency of the different containers. Should only
       * be called with the preprocessor flag @p DEBUG set. The function
       * should be called from the functions of the higher TriaLevel classes.
       */
      void
      monitor_memory(const unsigned int true_dimension) const;

      /**
       * Determine an estimate for the memory consumption (in bytes) of this
       * object.
       */
      std::size_t
      memory_consumption() const;

      /**
       * Read or write the data of this object to or from a stream for the
       * purpose of serialization
       */
      template <class Archive>
      void
      serialize(Archive &ar, const unsigned int version);

      /**
       * Exception
       */
      DeclException2(ExcMemoryInexact,
                     int,
                     int,
                     << "The containers have sizes " << arg1 << " and " << arg2
                     << ", which is not as expected.");
    };


    template <class Archive>
    void
    TriaLevel::serialize(Archive &ar, const unsigned int)
    {
      ar &dim;

      ar &refine_flags &coarsen_flags;

      // do not serialize 'active_cell_indices' here. instead of storing them
      // to the stream and re-reading them again later, we just rebuild them
      // in Triangulation::load()

      ar &neighbors;
      ar &subdomain_ids;
      ar &level_subdomain_ids;
      ar &parents;
      ar &direction_flags;
      ar &cells;

      if (dim == 3)
        ar &face_orientations;
    }

  } // namespace TriangulationImplementation
} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
