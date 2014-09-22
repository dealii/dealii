// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2013 by the deal.II authors
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

#ifndef __deal2__hp_dof_faces_h
#define __deal2__hp_dof_faces_h


#include <deal.II/base/config.h>
#include <deal.II/hp/fe_collection.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace hp
{
  template <int dim, int spacedim>
  class FECollection;
}


namespace internal
{
  namespace hp
  {
    /**
     * Store the indices of the degrees of freedom which are located on
     * objects of dimension @p structdim < dim, i.e., for faces or edges
     * of cells. This is opposed to the internal::hp::DoFLevels class
     * that stores the DoF indices on cells.
     *
     * The things we store here is very similar to what is stored in the
     * internal::DoFHandler::DoFObjects classes (see there for more
     * information, in particular on the layout of the class hierarchy,
     * and the use of file names).
     *
     * <h4>Offset computations</h4>
     *
     * For hp methods, not all cells may use the same finite element, and
     * it is consequently more complicated to determine where the DoF
     * indices for a given line, quad, or hex are stored. As described in
     * the documentation of the internal::DoFHandler::DoFLevel class, we
     * can compute the location of the first line DoF, for example, by
     * calculating the offset as <code>line_index *
     * dof_handler.get_fe().dofs_per_line</code>. This of course doesn't
     * work any more if different lines may have different numbers of
     * degrees of freedom associated with them. Consequently, rather than
     * using this simple multiplication, the dofs array has an associated
     * array dof_offsets. The data corresponding to a
     * line then starts at index <code>line_dof_offsets[line_index]</code>
     * within the <code>line_dofs</code> array.
     *
     *
     * <h4>Multiple data sets per object</h4>
     *
     * If two adjacent cells use different finite elements, then
     * the face that they share needs to store DoF indices for both
     * involved finite elements. While faces therefore have to have at
     * most two sets of DoF indices, it is easy to see that edges and
     * vertices can have as many sets of DoF indices associated with them
     * as there are adjacent cells.
     *
     * Consequently, for objects that have a lower dimensionality than
     * cells, we have to store a map from the finite element index to the
     * set of DoF indices associated. Since real sets are typically very
     * inefficient to store, and since most of the time we expect the
     * number of individual keys to be small (frequently, adjacent cells
     * will have the same finite element, and only a single entry will
     * exist in the map), what we do is instead to store a linked list. In
     * this format, the first entry starting at position
     * <code>lines.dofs[lines.dof_offsets[line_index]]</code> will denote
     * the finite element index of the set of DoF indices following; after
     * this set, we will store the finite element index of the second set
     * followed by the corresponding DoF indices; and so on. Finally, when
     * all finite element indices adjacent to this object have been
     * covered, we write a -1 to indicate the end of the list.
     *
     * Access to this kind of data, as well as the distinction between
     * cells and objects of lower dimensionality are encoded in the
     * accessor functions, DoFObjects::set_dof_index() and
     * DoFLevel::get_dof_index(). They are able to traverse this
     * list and pick out or set a DoF index given the finite element index
     * and its location within the set of DoFs corresponding to this
     * finite element.
     *
     *
     * @ingroup hp
     * @author Tobias Leicht, 2006
     */
    template <int structdim>
    class DoFIndicesOnFacesOrEdges
    {
    public:
      /**
       * Store the start index for
       * the degrees of freedom of each
       * object in the @p dofs array.
       *
       * The type we store is then obviously the type the @p dofs array
       * uses for indexing.
       */
      std::vector<unsigned int> dof_offsets;

      /**
       * Store the global indices of
       * the degrees of freedom. See
       * DoFLevel() for detailed
       * information.
       */
      std::vector<types::global_dof_index> dofs;

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
      template <int dim, int spacedim>
      void
      set_dof_index (const dealii::hp::DoFHandler<dim,spacedim> &dof_handler,
                     const unsigned int               obj_index,
                     const unsigned int               fe_index,
                     const unsigned int               local_index,
                     const types::global_dof_index    global_index,
                     const unsigned int               obj_level);

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
      template <int dim, int spacedim>
      types::global_dof_index
      get_dof_index (const dealii::hp::DoFHandler<dim,spacedim> &dof_handler,
                     const unsigned int               obj_index,
                     const unsigned int               fe_index,
                     const unsigned int               local_index,
                     const unsigned int               obj_level) const;

      /**
       * Return the number of
       * finite elements that are
       * active on a given
       * object. If this is a cell,
       * the answer is of course
       * one. If it is a face, the
       * answer may be one or two,
       * depending on whether the
       * two adjacent cells use the
       * same finite element or
       * not. If it is an edge in
       * 3d, the possible return
       * value may be one or any
       * other value larger than
       * that.
       *
       * If the object is not part
       * of an active cell, then no
       * degrees of freedom have
       * been distributed and zero
       * is returned.
       */
      template <int dim, int spacedim>
      unsigned int
      n_active_fe_indices (const dealii::hp::DoFHandler<dim,spacedim> &dof_handler,
                           const unsigned int               obj_index) const;

      /**
       * Return the fe_index of the
       * n-th active finite element
       * on this object.
       */
      template <int dim, int spacedim>
      types::global_dof_index
      nth_active_fe_index (const dealii::hp::DoFHandler<dim,spacedim> &dof_handler,
                           const unsigned int               obj_level,
                           const unsigned int               obj_index,
                           const unsigned int               n) const;

      /**
       * Check whether a given
       * finite element index is
       * used on the present
       * object or not.
       */
      template <int dim, int spacedim>
      bool
      fe_index_is_active (const dealii::hp::DoFHandler<dim,spacedim> &dof_handler,
                          const unsigned int               obj_index,
                          const unsigned int               fe_index,
                          const unsigned int               obj_level) const;

      /**
       * Determine an estimate for the
       * memory consumption (in bytes)
       * of this object.
       */
      std::size_t memory_consumption () const;
    };



    /**
     * These classes are similar to the internal::hp::DoFLevel classes. We here store
     * information that is associated with faces, rather than cells, as this information is
     * independent of the hierarchical structure of cells, which are organized in levels. In 2D
     * we store information on degrees of freedom located on lines whereas in 3D we store
     * information on drefrees of freedom located on quads and lines. In 1D we do nothing, as
     * the faces of lines are vertices which are treated separately.
     *
     * Apart from the internal::hp::DoFObjects object containing the data to store
     * (degree of freedom indices) and all the access-functionality to this data, we do not
     * store any data or provide any functionality. However, we do implement a function to
     * determine an estimate of the memory consumption of the contained
     * internal::hp::DoFObjects object(s).
     *
     * The data contained isn't usually directly accessed. Rather, except for some access from
     * the DoFHandler class, access is usually through the DoFAccessor::set_dof_index() and
     * DoFAccessor::dof_index() functions or similar functions of derived classes that in turn
     * access the member variables using the DoFHandler::get_dof_index() and corresponding
     * setter functions. Knowledge of the actual data format is therefore encapsulated to the
     * present hierarchy of classes as well as the ::DoFHandler class.
     *
     * @ingroup dofs
     * @author Tobias Leicht, 2006
     */
    template<int dim>
    class DoFIndicesOnFaces;


    /**
     * Store the indices of degrees of freedom on faces in 1D. As these would be vertices, which
     * are treated separately, don't do anything.
     *
     * @ingroup hp
     * @author Tobias Leicht, 2006
     */
    template<>
    class DoFIndicesOnFaces<1>
    {
    public:
      /**
       * Determine an estimate for the
       * memory consumption (in bytes)
       * of this object.
       */
      std::size_t memory_consumption () const;
    };

    /**
     * Store the indices of degrees of freedom on faces in 2D, which are lines.
     *
     * @ingroup hp
     * @author Tobias Leicht, 2006
     */
    template<>
    class DoFIndicesOnFaces<2>
    {
    public:
      /**
       * Indices of DoFs on the lines that bound cells.
       */
      internal::hp::DoFIndicesOnFacesOrEdges<1> lines;

      /**
       * Determine an estimate for the
       * memory consumption (in bytes)
       * of this object.
       */
      std::size_t memory_consumption () const;
    };

    /**
     * Store the indices of degrees of freedom on faces in 3D, which are
     * quads, additionally also on lines.
     *
     * @ingroup hp
     * @author Tobias Leicht, 2006
     */
    template<>
    class DoFIndicesOnFaces<3>
    {
    public:
      /**
       * Indices of DoFs on the lines that form the edges of cells.
       */
      internal::hp::DoFIndicesOnFacesOrEdges<1> lines;

      /**
       * Indices of DoFs on the quads that bound cells.
       */
      internal::hp::DoFIndicesOnFacesOrEdges<2> quads;

      /**
       * Determine an estimate for the
       * memory consumption (in bytes)
       * of this object.
       */
      std::size_t memory_consumption () const;
    };


    // --------------------- inline and template functions ------------------

    template <int structdim>
    template <int dim, int spacedim>
    inline
    types::global_dof_index
    DoFIndicesOnFacesOrEdges<structdim>::
    get_dof_index (const dealii::hp::DoFHandler<dim,spacedim> &dof_handler,
                   const unsigned int                obj_index,
                   const unsigned int                fe_index,
                   const unsigned int                local_index,
                   const unsigned int                obj_level) const
    {
      Assert ((fe_index != dealii::hp::DoFHandler<dim,spacedim>::default_fe_index),
              ExcMessage ("You need to specify a FE index when working "
                          "with hp DoFHandlers"));
      Assert (&dof_handler != 0,
              ExcMessage ("No DoFHandler is specified for this iterator"));
      Assert (&dof_handler.get_fe() != 0,
              ExcMessage ("No finite element collection is associated with "
                          "this DoFHandler"));
      Assert (fe_index < dof_handler.get_fe().size(),
              ExcIndexRange (fe_index, 0, dof_handler.get_fe().size()));
      Assert (local_index <
              dof_handler.get_fe()[fe_index].template n_dofs_per_object<structdim>(),
              ExcIndexRange(local_index, 0,
                            dof_handler.get_fe()[fe_index]
                            .template n_dofs_per_object<structdim>()));
      Assert (obj_index < dof_offsets.size(),
              ExcIndexRange (obj_index, 0, dof_offsets.size()));

      // make sure we are on an
      // object for which DoFs have
      // been allocated at all
      Assert (dof_offsets[obj_index] != numbers::invalid_unsigned_int,
              ExcMessage ("You are trying to access degree of freedom "
                          "information for an object on which no such "
                          "information is available"));

      Assert (structdim<dim, ExcMessage ("This object can not be used for cells."));

      // there may be multiple finite elements associated with
      // this object. hop along the list of index sets until we
      // find the one with the correct fe_index, and then poke
      // into that part. trigger an exception if we can't find a
      // set for this particular fe_index
      const types::global_dof_index starting_offset = dof_offsets[obj_index];
      const types::global_dof_index *pointer        = &dofs[starting_offset];
      while (true)
        {
          Assert (*pointer != numbers::invalid_dof_index,
                  ExcInternalError());
          if (*pointer == fe_index)
            return *(pointer + 1 + local_index);
          else
            pointer += static_cast<types::global_dof_index>(
                         dof_handler.get_fe()[*pointer]
                         .template n_dofs_per_object<structdim>() + 1);
        }
    }



    template <int structdim>
    template <int dim, int spacedim>
    inline
    void
    DoFIndicesOnFacesOrEdges<structdim>::
    set_dof_index (const dealii::hp::DoFHandler<dim,spacedim> &dof_handler,
                   const unsigned int                obj_index,
                   const unsigned int                fe_index,
                   const unsigned int                local_index,
                   const types::global_dof_index     global_index,
                   const unsigned int                obj_level)
    {
      Assert ((fe_index != dealii::hp::DoFHandler<dim,spacedim>::default_fe_index),
              ExcMessage ("You need to specify a FE index when working "
                          "with hp DoFHandlers"));
      Assert (&dof_handler != 0,
              ExcMessage ("No DoFHandler is specified for this iterator"));
      Assert (&dof_handler.get_fe() != 0,
              ExcMessage ("No finite element collection is associated with "
                          "this DoFHandler"));
      Assert (fe_index < dof_handler.get_fe().size(),
              ExcIndexRange (fe_index, 0, dof_handler.get_fe().size()));
      Assert (local_index <
              dof_handler.get_fe()[fe_index].template n_dofs_per_object<structdim>(),
              ExcIndexRange(local_index, 0,
                            dof_handler.get_fe()[fe_index]
                            .template n_dofs_per_object<structdim>()));
      Assert (obj_index < dof_offsets.size(),
              ExcIndexRange (obj_index, 0, dof_offsets.size()));

      // make sure we are on an
      // object for which DoFs have
      // been allocated at all
      Assert (dof_offsets[obj_index] != numbers::invalid_unsigned_int,
              ExcMessage ("You are trying to access degree of freedom "
                          "information for an object on which no such "
                          "information is available"));

      Assert (structdim<dim, ExcMessage ("This object can not be used for cells."));

      // there may be multiple finite elements associated with
      // this object.  hop along the list of index sets until we
      // find the one with the correct fe_index, and then poke
      // into that part. trigger an exception if we can't find a
      // set for this particular fe_index
      const types::global_dof_index starting_offset = dof_offsets[obj_index];
      types::global_dof_index      *pointer         = &dofs[starting_offset];
      while (true)
        {
          Assert (*pointer != numbers::invalid_dof_index,
                  ExcInternalError());
          if (*pointer == fe_index)
            {
              *(pointer + 1 + local_index) = global_index;
              return;
            }
          else
            pointer += dof_handler.get_fe()[*pointer]
                       .template n_dofs_per_object<structdim>() + 1;
        }
    }



    template <int structdim>
    template <int dim, int spacedim>
    inline
    unsigned int
    DoFIndicesOnFacesOrEdges<structdim>::
    n_active_fe_indices (const dealii::hp::DoFHandler<dim,spacedim> &dof_handler,
                         const unsigned int                obj_index) const
    {
      Assert (&dof_handler != 0,
              ExcMessage ("No DoFHandler is specified for this iterator"));
      Assert (&dof_handler.get_fe() != 0,
              ExcMessage ("No finite element collection is associated with "
                          "this DoFHandler"));
      Assert (obj_index < dof_offsets.size(),
              ExcIndexRange (obj_index, 0, dof_offsets.size()));

      // make sure we are on an
      // object for which DoFs have
      // been allocated at all
      if (dof_offsets[obj_index] == numbers::invalid_unsigned_int)
        return 0;

      Assert (structdim<dim, ExcMessage ("This object can not be used for cells."));

      // there may be multiple finite elements associated with this
      // object. hop along the list of index sets until we find the
      // one with the correct fe_index, and then poke into that
      // part. trigger an exception if we can't find a set for this
      // particular fe_index
      const unsigned int starting_offset = dof_offsets[obj_index];
      const types::global_dof_index *pointer        = &dofs[starting_offset];
      unsigned int counter = 0;
      while (true)
        {
          if (*pointer == numbers::invalid_dof_index)
            // end of list reached
            return counter;
          else
            {
              ++counter;
              pointer += dof_handler.get_fe()[*pointer]
                         .template n_dofs_per_object<structdim>() + 1;
            }
        }
    }



    template <int structdim>
    template <int dim, int spacedim>
    inline
    types::global_dof_index
    DoFIndicesOnFacesOrEdges<structdim>::
    nth_active_fe_index (const dealii::hp::DoFHandler<dim,spacedim> &dof_handler,
                         const unsigned int                obj_level,
                         const unsigned int                obj_index,
                         const unsigned int                n) const
    {
      Assert (&dof_handler != 0,
              ExcMessage ("No DoFHandler is specified for this iterator"));
      Assert (&dof_handler.get_fe() != 0,
              ExcMessage ("No finite element collection is associated with "
                          "this DoFHandler"));
      Assert (obj_index < dof_offsets.size(),
              ExcIndexRange (obj_index, 0, dof_offsets.size()));

      // make sure we are on an
      // object for which DoFs have
      // been allocated at all
      Assert (dof_offsets[obj_index] != numbers::invalid_unsigned_int,
              ExcMessage ("You are trying to access degree of freedom "
                          "information for an object on which no such "
                          "information is available"));

      Assert (structdim<dim, ExcMessage ("This object can not be used for cells."));

      Assert (n < n_active_fe_indices(dof_handler, obj_index),
              ExcIndexRange (n, 0,
                             n_active_fe_indices(dof_handler, obj_index)));

      // there may be multiple finite elements associated with
      // this object. hop along the list of index sets until we
      // find the one with the correct fe_index, and then poke
      // into that part. trigger an exception if we can't find a
      // set for this particular fe_index
      const unsigned int starting_offset = dof_offsets[obj_index];
      const types::global_dof_index *pointer = &dofs[starting_offset];
      unsigned int counter = 0;
      while (true)
        {
          Assert (*pointer != numbers::invalid_dof_index,
                  ExcInternalError());

          const unsigned int fe_index = *pointer;

          Assert (fe_index < dof_handler.get_fe().size(),
                  ExcInternalError());

          if (counter == n)
            return fe_index;

          ++counter;
          pointer += dof_handler.get_fe()[fe_index]
                     .template n_dofs_per_object<structdim>() + 1;
        }
    }



    template <int structdim>
    template <int dim, int spacedim>
    inline
    bool
    DoFIndicesOnFacesOrEdges<structdim>::
    fe_index_is_active (const dealii::hp::DoFHandler<dim,spacedim> &dof_handler,
                        const unsigned int                obj_index,
                        const unsigned int                fe_index,
                        const unsigned int                obj_level) const
    {
      Assert (&dof_handler != 0,
              ExcMessage ("No DoFHandler is specified for this iterator"));
      Assert (&dof_handler.get_fe() != 0,
              ExcMessage ("No finite element collection is associated with "
                          "this DoFHandler"));
      Assert (obj_index < dof_offsets.size(),
              ExcIndexRange (obj_index, 0, static_cast<unsigned int>(dof_offsets.size())));
      Assert ((fe_index != dealii::hp::DoFHandler<dim,spacedim>::default_fe_index),
              ExcMessage ("You need to specify a FE index when working "
                          "with hp DoFHandlers"));
      Assert (fe_index < dof_handler.get_fe().size(),
              ExcIndexRange (fe_index, 0, dof_handler.get_fe().size()));

      // make sure we are on an
      // object for which DoFs have
      // been allocated at all
      Assert (dof_offsets[obj_index] != numbers::invalid_unsigned_int,
              ExcMessage ("You are trying to access degree of freedom "
                          "information for an object on which no such "
                          "information is available"));

      Assert (structdim<dim, ExcMessage ("This object can not be used for cells."));

      // there may be multiple finite elements associated with
      // this object. hop along the list of index sets until we
      // find the one with the correct fe_index, and then poke
      // into that part. trigger an exception if we can't find a
      // set for this particular fe_index
      const types::global_dof_index starting_offset = dof_offsets[obj_index];
      const types::global_dof_index *pointer = &dofs[starting_offset];
      while (true)
        {
          if (*pointer == numbers::invalid_dof_index)
            // end of list reached
            return false;
          else if (*pointer == fe_index)
            return true;
          else
            pointer += static_cast<types::global_dof_index>(
                         dof_handler.get_fe()[*pointer]
                         .template n_dofs_per_object<structdim>()+1);
        }
    }


  } // namespace hp

} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
