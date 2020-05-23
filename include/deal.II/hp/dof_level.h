// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2020 by the deal.II authors
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

#ifndef dealii_hp_dof_level_h
#define dealii_hp_dof_level_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <vector>


DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
namespace hp
{
  template <int, int>
  class DoFHandler;
  template <int, int>
  class FECollection;
} // namespace hp


namespace internal
{
  namespace hp
  {
    namespace DoFHandlerImplementation
    {
      struct Implementation;
    }
  } // namespace hp
  namespace DoFCellAccessorImplementation
  {
    struct Implementation;
  }
} // namespace internal
#endif


namespace internal
{
  namespace hp
  {
    /**
     * This is the class that stores the degrees of freedom on cells in a hp
     * hierarchy. Compared to faces and edges, the task here is simple since
     * each cell can only have a single active finite element index.
     * Consequently, all we need is one long array with DoF indices and one
     * array of offsets where each cell's indices start within the array of
     * indices. This is in contrast to the DoFObjects class where each face or
     * edge may have more than one associated finite element with
     * corresponding degrees of freedom.
     *
     * The data stored here is represented by three arrays - The @p
     * active_fe_indices array stores for each cell which finite element is
     * used on this cell. Since some cells are not active on the current
     * level, some entries in this array may represent an invalid value. - The
     * @p dof_indices array stores for each active cell on the current level
     * the dofs that associated with the <i>interior</i> of the cell, i.e.,
     * the @p dofs_per_line dofs associated with the line in 1d, and @p
     * dofs_per_quad and @p dofs_per_hex in 2d and 3d. These numbers are in
     * general smaller than @p dofs_per_cell. - The @p dof_offsets array
     * stores, for each cell, the starting point of the dof indices
     * corresponding to this cell in the @p dof_indices array. This is
     * analogous to how we store data in compressed row storage for sparse
     * matrices. For cells that are not active on the current level, we store
     * an invalid value for the starting index.
     *
     * <h3>Compression</h3>
     *
     * It is common for the indices stored in @p dof_indices for one cell to
     * be numbered consecutively. For example, using the standard numbering
     * (without renumbering DoFs), the quad dofs on the first cell of a mesh
     * when using a $Q_3$ element will be numbered <tt>12, 13, 14, 15</tt>.
     * This allows for compression if we only store the first entry and have
     * some way to mark the DoFs on this object as compressed. Here,
     * compression means that we know that subsequent DoF indices can be
     * obtained from the previous ones by just incrementing them by one -- in
     * other words, we use a variant of doing run-length encoding. The way to
     * do this is that we use positive FE indices for uncompressed sets of
     * DoFs and if a set of indices is compressed, then we instead store the
     * FE index in binary complement (which we can identify by looking at the
     * sign bit when interpreting the number as a signed one). There are two
     * functions, compress_data() and uncompress_data() that convert between
     * the two possible representations.
     *
     * Note that compression is not always possible. For example, if one
     * renumbered the example above using DoFRenumbering::downstream with
     * $(1,0)^T$ as direction, then they would likely be numbered <tt>12, 14,
     * 13, 15</tt>, which can not be compressed using run-length encoding.
     */
    class DoFLevel
    {
    private:
      /**
       * The type in which we store the offsets into the dof_indices array.
       */
      using offset_type = unsigned int;

      /**
       * The type in which we store the active FE index.
       */
      using active_fe_index_type = unsigned short int;

      /**
       * A signed type that matches the type in which we store the active FE
       * index. We use this in computing binary complements.
       */
      using signed_active_fe_index_type = signed short int;

      /**
       * Invalid active_fe_index which will be used as a default value to
       * determine whether a future_fe_index has been set or not.
       */
      static const active_fe_index_type invalid_active_fe_index;

      /**
       * Given an active_fe_index, return whether the corresponding
       * set of DoF indices are compressed. See the general documentation
       * of this class for a description of when this is the case.
       */
      static bool
      is_compressed_entry(const active_fe_index_type active_fe_index);

      /**
       * Given an active_fe_index (either corresponding to an uncompressed
       * or compressed state), return the active_fe_index that corresponds
       * to the respectively other state. See the general documentation
       * of this class for a description of how compression is indicated.
       */
      static active_fe_index_type
      get_toggled_compression_state(const active_fe_index_type active_fe_index);

      /**
       * Indices specifying the finite element of hp::FECollection to use for
       * the different cells on the current level. The vector stores one
       * element per cell since the active_fe_index is unique for cells.
       *
       * Each active_fe_index will be preset to zero, i.e. the default finite
       * element. However, only active cells are eligible to have a finite
       * element assigned, which will be verified by corresponding accessor
       * functions of the DoFCellAccessor class.
       */
      std::vector<active_fe_index_type> active_fe_indices;

      /**
       * Indices specifying the finite element of hp::FECollection that the cell
       * will be assigned to after the triangulation changed (due to refinement
       * and/or coarsening, or repartitioning). The vector stores one element
       * per cell since the future_fe_index is unique for cells.
       *
       * Each future_fe_index will be preset to invalid_active_fe_index,
       * indicating that the active_fe_index will remain unchanged. However,
       * only active cells are eligible to have a finite element assigned, which
       * will be verified by corresponding accessor functions of the
       * DoFCellAccessor class.
       */
      std::vector<active_fe_index_type> future_fe_indices;

      /**
       * Store the start index for the degrees of freedom of each object in
       * the @p dof_indices array. If the cell corresponding to a particular
       * index in this array is not active on this level, then we do not store
       * any DoFs for it. In that case, the offset we store here must be an
       * invalid number and indeed we store
       * <code>(std::vector<types::global_dof_index>::size_type)(-1)</code>
       * for it.
       *
       * The type we store is then obviously the type the @p dof_indices array
       * uses for indexing.
       */
      std::vector<offset_type> dof_offsets;

      /**
       * Store the global indices of the degrees of freedom. information. The
       * dof_offsets field determines where each (active) cell's data is
       * stored.
       */
      std::vector<types::global_dof_index> dof_indices;

      /**
       * The offsets for each cell of the cache that holds all DoF indices.
       */
      std::vector<offset_type> cell_cache_offsets;

      /**
       * Cache for the DoF indices on cells. The size of this array equals the
       * sum over all cells of
       * selected_fe[active_fe_index[cell]].dofs_per_cell.
       */
      std::vector<types::global_dof_index> cell_dof_indices_cache;

    public:
      /**
       * Set the global index of the @p local_index-th degree of freedom
       * located on the object with number @p obj_index to the value given by
       * @p global_index. The @p dof_handler argument is used to access the
       * finite element that is to be used to compute the location where this
       * data is stored.
       *
       * The third argument, @p fe_index, denotes which of the finite elements
       * associated with this object we shall access. Refer to the general
       * documentation of the internal::hp::DoFLevel class template for more
       * information.
       */
      void
      set_dof_index(const unsigned int            obj_index,
                    const unsigned int            fe_index,
                    const unsigned int            local_index,
                    const types::global_dof_index global_index);

      /**
       * Return the global index of the @p local_index-th degree of freedom
       * located on the object with number @p obj_index. The @p dof_handler
       * argument is used to access the finite element that is to be used to
       * compute the location where this data is stored.
       *
       * The third argument, @p fe_index, denotes which of the finite elements
       * associated with this object we shall access. Refer to the general
       * documentation of the internal::hp::DoFLevel class template for more
       * information.
       */
      types::global_dof_index
      get_dof_index(const unsigned int obj_index,
                    const unsigned int fe_index,
                    const unsigned int local_index) const;

      /**
       * Return the fe_index of the active finite element on this object.
       */
      unsigned int
      active_fe_index(const unsigned int obj_index) const;

      /**
       * Check whether a given finite element index is used on the present
       * object or not.
       */
      bool
      fe_index_is_active(const unsigned int obj_index,
                         const unsigned int fe_index) const;

      /**
       * Set the fe_index of the active finite element on this object.
       */
      void
      set_active_fe_index(const unsigned int obj_index,
                          const unsigned int fe_index);

      /**
       * Return the fe_index of the future finite element on this object. If no
       * future_fe_index has been specified, return the active_fe_index instead.
       */
      unsigned int
      future_fe_index(const unsigned int obj_index) const;

      /**
       * Set the fe_index of the future finite element on this object.
       */
      void
      set_future_fe_index(const unsigned int obj_index,
                          const unsigned int fe_index);

      /**
       * Return whether a future fe index has been set on this object.
       */
      bool
      future_fe_index_set(const unsigned int obj_index) const;

      /**
       * Revoke the future finite element assigned to this object.
       */
      void
      clear_future_fe_index(const unsigned int obj_index);

      /**
       * Return a pointer to the beginning of the DoF indices cache for a
       * given cell.
       *
       * @param obj_index The number of the cell we are looking at.
       * @param dofs_per_cell The number of DoFs per cell for this cell. This
       * is not used for the hp case but necessary to keep the interface the
       * same as for the non-hp case.
       * @return A pointer to the first DoF index for the current cell. The
       * next dofs_per_cell indices are for the current cell.
       */
      const types::global_dof_index *
      get_cell_cache_start(const unsigned int obj_index,
                           const unsigned int dofs_per_cell) const;

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

    private:
      /**
       * Compress the arrays that store dof indices by using a variant of run-
       * length encoding. See the general documentation of this class for more
       * information.
       *
       * @param fe_collection The object that can tell us how many degrees of
       * freedom each of the finite elements has that we store in this object.
       */
      template <int dim, int spacedim>
      void
      compress_data(
        const dealii::hp::FECollection<dim, spacedim> &fe_collection);

      /**
       * Uncompress the arrays that store dof indices by using a variant of
       * run-length encoding. See the general documentation of this class for
       * more information.
       *
       * @param fe_collection The object that can tell us how many degrees of
       * freedom each of the finite elements has that we store in this object.
       */
      template <int dim, int spacedim>
      void
      uncompress_data(
        const dealii::hp::FECollection<dim, spacedim> &fe_collection);


      /**
       * Restore the active fe indices stored by the current object to
       * an uncompressed state. At the same time, do not touch any of
       * the other data structures. This leaves the active_fe_indices
       * array out-of-synch from the other member variables, and the
       * function can consequently only be used when the remaining
       * data structures are going to be rebuilt next. This is
       * specifically the case in hp::DoFHandler::distribute_dofs()
       * where we throw away all data in the DoFLevels objects *except
       * for the active_fe_indices* array. This function therefore
       * simply makes sure that that one array is in uncompressed
       * format so that we can use the information about active
       * fe_indices for all cells that were used in the previous mesh
       * refinement cycle (or the previous time distribute_dofs() was
       * called) without having to care about any of the other data
       * fields.
       */
      void
      normalize_active_fe_indices();


      // Make hp::DoFHandler and its auxiliary class a friend since it is the
      // class that needs to create these data structures.
      template <int, int>
      friend class dealii::hp::DoFHandler;
      friend struct dealii::internal::hp::DoFHandlerImplementation::
        Implementation;
      friend struct dealii::internal::DoFCellAccessorImplementation::
        Implementation;
    };


    // -------------------- template functions --------------------------------


    inline bool
    DoFLevel::is_compressed_entry(const active_fe_index_type active_fe_index)
    {
      return (static_cast<signed_active_fe_index_type>(active_fe_index) < 0);
    }



    inline DoFLevel::active_fe_index_type
    DoFLevel::get_toggled_compression_state(
      const active_fe_index_type active_fe_index)
    {
      // convert the active_fe_index into a signed type, flip all
      // bits, and get the unsigned representation back
      return static_cast<active_fe_index_type>(
        ~(static_cast<signed_active_fe_index_type>(active_fe_index)));
    }



    inline types::global_dof_index
    DoFLevel::get_dof_index(const unsigned int obj_index,
                            const unsigned int fe_index,
                            const unsigned int local_index) const
    {
      (void)fe_index;
      AssertIndexRange(obj_index, dof_offsets.size());

      // make sure we are on an object for which DoFs have been
      // allocated at all
      Assert(dof_offsets[obj_index] != static_cast<offset_type>(-1),
             ExcMessage("You are trying to access degree of freedom "
                        "information for an object on which no such "
                        "information is available"));

      Assert(fe_index ==
               (is_compressed_entry(active_fe_indices[obj_index]) == false ?
                  active_fe_indices[obj_index] :
                  get_toggled_compression_state(active_fe_indices[obj_index])),
             ExcMessage("FE index does not match that of the present cell"));

      // see if the dof_indices array has been compressed for this
      // particular cell
      if (is_compressed_entry(active_fe_indices[obj_index]) == false)
        return dof_indices[dof_offsets[obj_index] + local_index];
      else
        return dof_indices[dof_offsets[obj_index]] + local_index;
    }



    inline void
    DoFLevel::set_dof_index(const unsigned int            obj_index,
                            const unsigned int            fe_index,
                            const unsigned int            local_index,
                            const types::global_dof_index global_index)
    {
      (void)fe_index;
      AssertIndexRange(obj_index, dof_offsets.size());

      // make sure we are on an
      // object for which DoFs have
      // been allocated at all
      Assert(dof_offsets[obj_index] != static_cast<offset_type>(-1),
             ExcMessage("You are trying to access degree of freedom "
                        "information for an object on which no such "
                        "information is available"));
      Assert(
        is_compressed_entry(active_fe_indices[obj_index]) == false,
        ExcMessage(
          "This function can no longer be called after compressing the dof_indices array"));
      Assert(fe_index == active_fe_indices[obj_index],
             ExcMessage("FE index does not match that of the present cell"));
      dof_indices[dof_offsets[obj_index] + local_index] = global_index;
    }



    inline unsigned int
    DoFLevel::active_fe_index(const unsigned int obj_index) const
    {
      AssertIndexRange(obj_index, active_fe_indices.size());

      if (is_compressed_entry(active_fe_indices[obj_index]) == false)
        return active_fe_indices[obj_index];
      else
        return get_toggled_compression_state(active_fe_indices[obj_index]);
    }



    inline bool
    DoFLevel::fe_index_is_active(const unsigned int obj_index,
                                 const unsigned int fe_index) const
    {
      return (fe_index == active_fe_index(obj_index));
    }



    inline void
    DoFLevel::set_active_fe_index(const unsigned int obj_index,
                                  const unsigned int fe_index)
    {
      AssertIndexRange(obj_index, active_fe_indices.size());

      // check whether the given fe_index is within the range of
      // values that we interpret as "not compressed". if not, then
      // the index is so large that we cannot accept it. (but this
      // will not likely happen because it requires someone using an
      // FECollection that has more than 32k entries.)
      Assert(is_compressed_entry(fe_index) == false,
             ExcMessage(
               "You are using an active_fe_index that is larger than an "
               "internal limitation for these objects. Try to work with "
               "hp::FECollection objects that have a more modest size."));
      Assert(fe_index != invalid_active_fe_index,
             ExcMessage(
               "You are using an active_fe_index that is reserved for "
               "internal purposes for these objects. Try to work with "
               "hp::FECollection objects that have a more modest size."));

      active_fe_indices[obj_index] = fe_index;
    }



    inline unsigned int
    DoFLevel::future_fe_index(const unsigned int obj_index) const
    {
      AssertIndexRange(obj_index, future_fe_indices.size());

      if (future_fe_index_set(obj_index))
        return future_fe_indices[obj_index];

      return active_fe_index(obj_index);
    }



    inline void
    DoFLevel::set_future_fe_index(const unsigned int obj_index,
                                  const unsigned int fe_index)
    {
      AssertIndexRange(obj_index, future_fe_indices.size());

      // check whether the given fe_index is within the range of
      // values that we interpret as "not compressed". if not, then
      // the index is so large that we cannot accept it. (but this
      // will not likely happen because it requires someone using an
      // FECollection that has more than 32k entries.)
      Assert(is_compressed_entry(fe_index) == false,
             ExcMessage(
               "You are using a future_fe_index that is larger than an "
               "internal limitation for these objects. Try to work with "
               "hp::FECollection objects that have a more modest size."));
      Assert(fe_index != invalid_active_fe_index,
             ExcMessage(
               "You are using a future_fe_index that is reserved for "
               "internal purposes for these objects. Try to work with "
               "hp::FECollection objects that have a more modest size."));

      future_fe_indices[obj_index] = fe_index;
    }



    inline bool
    DoFLevel::future_fe_index_set(const unsigned int obj_index) const
    {
      AssertIndexRange(obj_index, future_fe_indices.size());

      return (future_fe_indices[obj_index] != invalid_active_fe_index);
    }



    inline void
    DoFLevel::clear_future_fe_index(const unsigned int obj_index)
    {
      AssertIndexRange(obj_index, future_fe_indices.size());

      future_fe_indices[obj_index] = invalid_active_fe_index;
    }



    inline const types::global_dof_index *
    DoFLevel::get_cell_cache_start(const unsigned int obj_index,
                                   const unsigned int dofs_per_cell) const
    {
      (void)dofs_per_cell;
      Assert(
        (obj_index < cell_cache_offsets.size()) &&
          (cell_cache_offsets[obj_index] + dofs_per_cell <=
           cell_dof_indices_cache.size()),
        ExcMessage(
          "You are trying to access an element of the cache that stores "
          "the indices of all degrees of freedom that live on one cell. "
          "However, this element does not exist. Did you forget to call "
          "DoFHandler::distribute_dofs(), or did you forget to call it "
          "again after changing the active_fe_index of one of the cells?"));

      return &cell_dof_indices_cache[cell_cache_offsets[obj_index]];
    }



    template <class Archive>
    inline void
    DoFLevel::serialize(Archive &ar, const unsigned int)
    {
      ar & this->active_fe_indices;
      ar & this->cell_cache_offsets;
      ar & this->cell_dof_indices_cache;
      ar & this->dof_indices;
      ar & this->dof_offsets;
      ar & this->future_fe_indices;
    }
  } // namespace hp

} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
