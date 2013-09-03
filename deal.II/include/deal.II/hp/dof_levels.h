// ---------------------------------------------------------------------
// $Id$
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

#ifndef __deal2__hp_dof_levels_h
#define __deal2__hp_dof_levels_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace hp
{
  template <int, int> class DoFHandler;
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
     */
    template <int dim>
    class DoFLevel
    {
    private:
      /**
       *  Indices specifying the finite
       *  element of hp::FECollection to use
       *  for the different cells on the current level. The
       *  meaning what a cell is, is
       *  dimension specific, therefore also
       *  the length of this vector depends
       *  on the dimension: in one dimension,
       *  the length of this vector equals
       *  the length of the @p lines vector,
       *  in two dimensions that of the @p
       *  quads vector, etc. The vector stores one element per cell
       *  since the actiev_fe_index is unique for cells.
       */
      std::vector<unsigned int> active_fe_indices;

      /**
       * Store the start index for
       * the degrees of freedom of each
       * object in the @p dofs array.
       *
       * The type we store is then obviously the type the @p dofs array
       * uses for indexing.
       */
      std::vector<std::vector<types::global_dof_index>::size_type> dof_offsets;

      /**
       * Store the global indices of
       * the degrees of freedom. See
       * DoFLevel() for detailed
       * information.
       */
      std::vector<types::global_dof_index> dofs;

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
       * Determine an estimate for the
       * memory consumption (in bytes)
       * of this object.
       */
      std::size_t memory_consumption () const;

      /**
       * Make hp::DoFHandler and its auxiliary class a friend since it
       * is the class that needs to create these data structures.
       */
      template <int, int> friend class dealii::hp::DoFHandler;
      friend class dealii::internal::hp::DoFHandler::Implementation;
    };


    // -------------------- template functions --------------------------------

    template <int dim>
    inline
    types::global_dof_index
    DoFLevel<dim>::
    get_dof_index (const unsigned int                obj_index,
                   const unsigned int                fe_index,
                   const unsigned int                local_index) const
    {
      Assert (obj_index < dof_offsets.size(),
              ExcIndexRange (obj_index, 0, dof_offsets.size()));

      // make sure we are on an
      // object for which DoFs have
      // been allocated at all
      Assert (dof_offsets[obj_index] != numbers::invalid_dof_index,
              ExcMessage ("You are trying to access degree of freedom "
                          "information for an object on which no such "
                          "information is available"));

      Assert (fe_index == active_fe_indices[obj_index],
	      ExcMessage ("FE index does not match that of the present cell"));
      return dofs[dof_offsets[obj_index]+local_index];
    }



    template <int dim>
    inline
    void
    DoFLevel<dim>::
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
      Assert (dof_offsets[obj_index] != numbers::invalid_dof_index,
              ExcMessage ("You are trying to access degree of freedom "
                          "information for an object on which no such "
                          "information is available"));

      Assert (fe_index == active_fe_indices[obj_index],
	      ExcMessage ("FE index does not match that of the present cell"));
      dofs[dof_offsets[obj_index]+local_index] = global_index;
    }



    template <int dim>
    inline
    unsigned int
    DoFLevel<dim>::
    active_fe_index (const unsigned int obj_index) const
    {
      Assert (obj_index < active_fe_indices.size(),
              ExcIndexRange (obj_index, 0, active_fe_indices.size()));

      return active_fe_indices[obj_index];
    }



    template <int dim>
    inline
    bool
    DoFLevel<dim>::
    fe_index_is_active (const unsigned int                obj_index,
                        const unsigned int                fe_index) const
    {
      const unsigned int invalid_fe_index = numbers::invalid_unsigned_int;
      Assert ((fe_index != invalid_fe_index),
              ExcMessage ("You need to specify a FE index when working "
                          "with hp DoFHandlers"));

      // make sure we are on an
      // object for which DoFs have
      // been allocated at all
      Assert (dof_offsets[obj_index] != numbers::invalid_dof_index,
              ExcMessage ("You are trying to access degree of freedom "
                          "information for an object on which no such "
                          "information is available"));

      Assert (obj_index < active_fe_indices.size(),
	      ExcInternalError());
      return (fe_index == active_fe_indices[obj_index]);
    }


    template <int dim>
    inline
    void
    DoFLevel<dim>::
    set_active_fe_index (const unsigned int obj_index,
			 const unsigned int fe_index)
    {
      Assert (obj_index < active_fe_indices.size(),
              ExcIndexRange (obj_index, 0, active_fe_indices.size()));

      active_fe_indices[obj_index] = fe_index;
    }

  } // namespace hp

} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
