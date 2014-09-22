// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
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

#ifndef __deal2__hp_dof_level_h
#define __deal2__hp_dof_level_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>

#include <vector>


DEAL_II_NAMESPACE_OPEN

namespace hp
{
  template <int, int> class DoFHandler;
  template <int, int> class FECollection;
}


namespace internal
{
  namespace hp
  {
    namespace DoFHandler
    {
      struct Implementation;
    }
  }
  namespace DoFCellAccessor
  {
    struct Implementation;
  }
}


namespace internal
{
  namespace hp
  {
    /**
     * This is the class that stores the degrees of freedom on cells in a hp
     * hierarchy. Compared to faces and edges, the task here is simple since
     * each cell can only have a single active finite element index. Consequently,
     * all we need is one long array with DoF indices and one array of offsets
     * where each cell's indices start within the array of indices. This is in
     * contrast to the DoFObjects class where each face or edge may have more than
     * one associated finite element with corresponding degrees of freedom.
     *
     * The data stored here is represented by three arrays
     * - The @p active_fe_indices array stores for each cell which
     *   finite element is used on this cell. Since some cells are not active
     *   on the current level, some entries in this array may represent an
     *   invalid value.
     * - The @p dof_indices array stores for each active cell on the current
     *   level the dofs that associated with the <i>interior</i> of the cell,
     *   i.e., the @p dofs_per_line dofs associated with the line in 1d, and
     *   @p dofs_per_quad and @p dofs_per_hex in 2d and 3d. These numbers are
     *   in general smaller than @p dofs_per_cell.
     * - The @p dof_offsets array stores, for each cell, the starting point
     *   of the dof indices corresponding to this cell in the @p dof_indices
     *   array. This is analogous to how we store data in compressed row storage
     *   for sparse matrices. For cells that are not active on the current level,
     *   we store an invalid value for the starting index.
     *
     * <h3>Compression</h3>
     *
     * It is common for the indices stored in @p dof_indices for one cell to be
     * numbered consecutively. For example, using the standard numbering (without
     * renumbering DoFs), the quad dofs on the first cell of a mesh when using a
     * $Q_3$ element will be numbered <tt>12, 13, 14, 15</tt>. This allows for
     * compression if we only store the first entry and have some way to mark
     * the DoFs on this object as compressed. Here, compression means that we
     * know that subsequent DoF indices can be obtained from the previous ones
     * by just incrementing them by one -- in other words, we use a variant of doing
     * run-length encoding. The way to do this is that we
     * use positive FE indices for uncompressed sets of DoFs and if a set of indices
     * is compressed, then we instead store the FE index in binary complement (which
     * we can identify by looking at the sign bit when interpreting the number as a
     * signed one). There are two functions, compress_data() and uncompress_data()
     * that convert between the two possible representations.
     *
     * Note that compression is not always possible. For example, if one renumbered
     * the example above using DoFRenumbering::downstream with $(1,0)^T$ as
     * direction, then they would likely be numbered <tt>12, 14, 13, 15</tt>, which
     * can not be compressed using run-length encoding.
     */
    class DoFLevel
    {
    private:
      /**
       * The type in which we store the offsets into the dof_indices array.
       */
      typedef unsigned int offset_type;

      /**
       * The type in which we store the active FE index.
       */
      typedef unsigned short int active_fe_index_type;

      /**
       * A signed type that matches the type in which we store the active FE
       * index. We use this in computing binary complements.
       */
      typedef signed short int signed_active_fe_index_type;

      /**
       * Indices specifying the finite element of hp::FECollection to
       * use for the different cells on the current level. The vector
       * stores one element per cell since the active_fe_index is
       * unique for cells.
       *
       * If a cell is not active on the level corresponding to the
       * current object (i.e., it has children on higher levels) then
       * it does not have an associated fe index and we store
       * an invalid fe index marker instead.
       */
      std::vector<active_fe_index_type> active_fe_indices;

      /**
       * Store the start index for the degrees of freedom of each
       * object in the @p dof_indices array. If the cell corresponding to
       * a particular index in this array is not active on this level,
       * then we do not store any DoFs for it. In that case, the offset
       * we store here must be an invalid number and indeed we store
       * <code>(std::vector<types::global_dof_index>::size_type)(-1)</code>
       * for it.
       *
       * The type we store is then obviously the type the @p dof_indices array
       * uses for indexing.
       */
      std::vector<offset_type> dof_offsets;

      /**
       * Store the global indices of the degrees of freedom.
       * information. The dof_offsets field determines where each
       * (active) cell's data is stored.
       */
      std::vector<types::global_dof_index> dof_indices;

      /**
       * The offsets for each cell of the cache that holds all DoF indices.
       */
      std::vector<offset_type> cell_cache_offsets;

      /**
       * Cache for the DoF indices
       * on cells. The size of this
       * array equals the sum over all cells of selected_fe[active_fe_index[cell]].dofs_per_cell.
       */
      std::vector<types::global_dof_index> cell_dof_indices_cache;

    public:

      /**
       * Set the global index of
       * the @p local_index-th
       * degree of freedom located
       * on the object with number @p
       * obj_index to the value
       * given by @p global_index. The @p
       * dof_handler argument is
       * used to access the finite
       * element that is to be used
       * to compute the location
       * where this data is stored.
       *
       * The third argument, @p
       * fe_index, denotes which of
       * the finite elements
       * associated with this
       * object we shall
       * access. Refer to the
       * general documentation of
       * the internal::hp::DoFLevel
       * class template for more
       * information.
       */
      void
      set_dof_index (const unsigned int               obj_index,
                     const unsigned int               fe_index,
                     const unsigned int               local_index,
                     const types::global_dof_index    global_index);

      /**
       * Return the global index of
       * the @p local_index-th
       * degree of freedom located
       * on the object with number @p
       * obj_index. The @p
       * dof_handler argument is
       * used to access the finite
       * element that is to be used
       * to compute the location
       * where this data is stored.
       *
       * The third argument, @p
       * fe_index, denotes which of
       * the finite elements
       * associated with this
       * object we shall
       * access. Refer to the
       * general documentation of
       * the internal::hp::DoFLevel
       * class template for more
       * information.
       */
      types::global_dof_index
      get_dof_index (const unsigned int               obj_index,
                     const unsigned int               fe_index,
                     const unsigned int               local_index) const;

      /**
       * Return the fe_index of the
       * active finite element
       * on this object.
       */
      unsigned int
      active_fe_index (const unsigned int obj_index) const;

      /**
       * Check whether a given
       * finite element index is
       * used on the present
       * object or not.
       */
      bool
      fe_index_is_active (const unsigned int               obj_index,
                          const unsigned int               fe_index) const;

      /**
       * Set the fe_index of the
       * active finite element
       * on this object.
       */
      void
      set_active_fe_index (const unsigned int obj_index,
                           const unsigned int fe_index);

      /**
       * Return a pointer to the beginning of the DoF indices cache
       * for a given cell.
       *
       * @param obj_index The number of the cell we are looking at.
       * @param dofs_per_cell The number of DoFs per cell for this cell. This
       *   is not used for the hp case but necessary to keep the interface
       *   the same as for the non-hp case.
       * @return A pointer to the first DoF index for the current cell. The
       *   next dofs_per_cell indices are for the current cell.
       */
      const types::global_dof_index *
      get_cell_cache_start (const unsigned int obj_index,
                            const unsigned int dofs_per_cell) const;

      /**
       * Determine an estimate for the
       * memory consumption (in bytes)
       * of this object.
       */
      std::size_t memory_consumption () const;

    private:
      /**
       * Compress the arrays that store dof indices by using a variant
       * of run-length encoding. See the general documentation of this
       * class for more information.
       *
       * @param fe_collection The object that can tell us how many
       * degrees of freedom each of the finite elements has that we
       * store in this object.
       */
      template <int dim, int spacedim>
      void compress_data (const dealii::hp::FECollection<dim,spacedim> &fe_collection);

      /**
       * Uncompress the arrays that store dof indices by using a variant
       * of run-length encoding. See the general documentation of this
       * class for more information.
       *
       * @param fe_collection The object that can tell us how many
       * degrees of freedom each of the finite elements has that we
       * store in this object.
       */
      template <int dim, int spacedim>
      void uncompress_data (const dealii::hp::FECollection<dim,spacedim> &fe_collection);

      /**
       * Make hp::DoFHandler and its auxiliary class a friend since it
       * is the class that needs to create these data structures.
       */
      template <int, int> friend class dealii::hp::DoFHandler;
      friend struct dealii::internal::hp::DoFHandler::Implementation;
      friend struct dealii::internal::DoFCellAccessor::Implementation;
    };


    // -------------------- template functions --------------------------------

    inline
    types::global_dof_index
    DoFLevel::
    get_dof_index (const unsigned int                obj_index,
                   const unsigned int                fe_index,
                   const unsigned int                local_index) const
    {
      Assert (obj_index < dof_offsets.size(),
              ExcIndexRange (obj_index, 0, dof_offsets.size()));

      // make sure we are on an
      // object for which DoFs have
      // been allocated at all
      Assert (dof_offsets[obj_index] != (offset_type)(-1),
              ExcMessage ("You are trying to access degree of freedom "
                          "information for an object on which no such "
                          "information is available"));

      Assert (fe_index == active_fe_indices[obj_index],
              ExcMessage ("FE index does not match that of the present cell"));

      // see if the dof_indices array has been compressed for this
      // particular cell
      if ((signed_active_fe_index_type)active_fe_indices[obj_index]>=0)
        return dof_indices[dof_offsets[obj_index]+local_index];
      else
        return dof_indices[dof_offsets[obj_index]]+local_index;
    }



    inline
    void
    DoFLevel::
    set_dof_index (const unsigned int                obj_index,
                   const unsigned int                fe_index,
                   const unsigned int                local_index,
                   const types::global_dof_index     global_index)
    {
      Assert (obj_index < dof_offsets.size(),
              ExcIndexRange (obj_index, 0, dof_offsets.size()));

      // make sure we are on an
      // object for which DoFs have
      // been allocated at all
      Assert (dof_offsets[obj_index] != (offset_type)(-1),
              ExcMessage ("You are trying to access degree of freedom "
                          "information for an object on which no such "
                          "information is available"));
      Assert ((signed_active_fe_index_type)active_fe_indices[obj_index]>=0,
              ExcMessage ("This function can no longer be called after compressing the dof_indices array"));
      Assert (fe_index == active_fe_indices[obj_index],
              ExcMessage ("FE index does not match that of the present cell"));
      dof_indices[dof_offsets[obj_index]+local_index] = global_index;
    }



    inline
    unsigned int
    DoFLevel::
    active_fe_index (const unsigned int obj_index) const
    {
      Assert (obj_index < active_fe_indices.size(),
              ExcIndexRange (obj_index, 0, active_fe_indices.size()));

      if (((signed_active_fe_index_type)active_fe_indices[obj_index]) >= 0)
        return active_fe_indices[obj_index];
      else
        return (active_fe_index_type)~(signed_active_fe_index_type)active_fe_indices[obj_index];
    }



    inline
    bool
    DoFLevel::
    fe_index_is_active (const unsigned int                obj_index,
                        const unsigned int                fe_index) const
    {
      return (fe_index == active_fe_index(obj_index));
    }


    inline
    void
    DoFLevel::
    set_active_fe_index (const unsigned int obj_index,
                         const unsigned int fe_index)
    {
      Assert (obj_index < active_fe_indices.size(),
              ExcIndexRange (obj_index, 0, active_fe_indices.size()));

      active_fe_indices[obj_index] = fe_index;
    }



    inline
    const types::global_dof_index *
    DoFLevel::get_cell_cache_start (const unsigned int obj_index,
                                    const unsigned int dofs_per_cell) const
    {
      Assert (obj_index < cell_cache_offsets.size(),
              ExcInternalError());
      Assert (cell_cache_offsets[obj_index]+dofs_per_cell
              <=
              cell_dof_indices_cache.size(),
              ExcInternalError());

      return &cell_dof_indices_cache[cell_cache_offsets[obj_index]];
    }
  } // namespace hp

} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
