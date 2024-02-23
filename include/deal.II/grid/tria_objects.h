// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_tria_objects_h
#define dealii_tria_objects_h

#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/geometry_info.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
class Triangulation;
template <class Accessor>
class TriaRawIterator;
template <int, int, int>
class TriaAccessor;
#endif

namespace internal
{
  namespace TriangulationImplementation
  {
    /**
     * General template for information belonging to the geometrical objects
     * of a triangulation, i.e. lines, quads, hexahedra...  Apart from the
     * vector of objects additional information is included, namely vectors
     * indicating the children, the used-status, user-flags, material-ids..
     *
     * Objects of these classes are included in the TriaLevel and TriaFaces
     * classes.
     */
    class TriaObjects
    {
    public:
      /**
       * Constructor resetting some data.
       */
      TriaObjects();

      /**
       * Constructor for a specific dimension.
       */
      TriaObjects(const unsigned int structdim);

      unsigned int structdim;

      /**
       * Vector of the objects belonging to this level. The index of the
       * object equals the index in this container.
       */
      std::vector<int> cells;

      /**
       * Return number of geometric objects stored by this class.
       */
      unsigned int
      n_objects() const;

      /**
       * Return a view on the indices of the objects that bound the @p
       * index-th object stored by the current object. For example, if
       * the current object stores cells, then this function returns
       * the equivalent of an array containing the indices of the
       * faces that bound the @p index-th cell.
       */
      ArrayView<int>
      get_bounding_object_indices(const unsigned int index);

      /**
       * Index of the even children of an object. Since when objects are
       * refined, all children are created at the same time, they are appended
       * to the list at least in pairs after each other. We therefore only
       * store the index of the even children, the uneven follow immediately
       * afterwards.
       *
       * If an object has no children, -1 is stored in this list. An object is
       * called active if it has no children. The function
       * TriaAccessorBase::has_children() tests for this.
       */
      std::vector<int> children;

      /**
       * Store the refinement case each of the cells is refined with. This
       * vector might be replaced by vector<vector<bool> > (dim, vector<bool>
       * (n_cells)) which is more memory efficient.
       */
      std::vector<std::uint8_t> refinement_cases;

      /**
       * Vector storing whether an object is used in the @p cells vector.
       *
       * Since it is difficult to delete elements in a @p vector, when an
       * element is not needed any more (e.g. after derefinement), it is not
       * deleted from the list, but rather the according @p used flag is set
       * to @p false.
       */
      std::vector<bool> used;

      /**
       * Make available a field for user data, one bit per object. This field
       * is usually used when an operation runs over all cells and needs
       * information whether another cell (e.g. a neighbor) has already been
       * processed.
       *
       * You can clear all used flags using
       * Triangulation::clear_user_flags().
       */
      std::vector<bool> user_flags;


      /**
       * We use this union to store boundary and material data. Because only
       * one out of these two is actually needed here, we use an union.
       */
      struct BoundaryOrMaterialId
      {
        union
        {
          types::boundary_id boundary_id;
          types::material_id material_id;
        };


        /**
         * Default constructor.
         */
        BoundaryOrMaterialId();

        /**
         * Return the size of objects of this kind.
         */
        static std::size_t
        memory_consumption();

        /**
         * Read or write the data of this object to or from a stream for the
         * purpose of serialization using the [BOOST serialization
         * library](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html).
         */
        template <class Archive>
        void
        serialize(Archive &ar, const unsigned int version);
      };

      /**
       * Store boundary and material data. For example, in one dimension, this
       * field stores the material id of a line, which is a number between 0
       * and numbers::invalid_material_id-1. In more than one dimension, lines
       * have no material id, but they may be at the boundary; then, we store
       * the boundary indicator in this field, which denotes to which part of
       * the boundary this line belongs and which boundary conditions hold on
       * this part. The boundary indicator also is a number between zero and
       * numbers::internal_face_boundary_id-1; the id
       * numbers::internal_face_boundary_id is reserved for lines in the
       * interior and may be used to check whether a line is at the boundary
       * or not, which otherwise is not possible if you don't know which cell
       * it belongs to.
       */
      std::vector<BoundaryOrMaterialId> boundary_or_material_id;

      /**
       * Store manifold ids. This field stores the manifold id of each object,
       * which is a number between 0 and numbers::flat_manifold_id-1.
       */
      std::vector<types::manifold_id> manifold_id;

      /**
       * Return an iterator to the next free slot for a single object. This
       * function is only used by Triangulation::execute_refinement()
       * in 3d.
       *
       * @warning Interestingly, this function is not used for 1d or 2d
       * triangulations, where it seems the authors of the refinement function
       * insist on reimplementing its contents.
       *
       * @todo This function is not instantiated for the codim-one case
       */
      template <int structdim, int dim, int spacedim>
      dealii::TriaRawIterator<dealii::TriaAccessor<structdim, dim, spacedim>>
      next_free_single_object(const Triangulation<dim, spacedim> &tria);

      /**
       * Return an iterator to the next free slot for a pair of objects. This
       * function is only used by Triangulation::execute_refinement()
       * in 3d.
       *
       * @warning Interestingly, this function is not used for 1d or 2d
       * triangulations, where it seems the authors of the refinement function
       * insist on reimplementing its contents.
       *
       * @todo This function is not instantiated for the codim-one case
       */
      template <int structdim, int dim, int spacedim>
      dealii::TriaRawIterator<dealii::TriaAccessor<structdim, dim, spacedim>>
      next_free_pair_object(const Triangulation<dim, spacedim> &tria);

      /**
       * Return an iterator to the next free slot for a pair of hexes. Only
       * implemented for <code>G=Hexahedron</code>.
       */
      template <int dim, int spacedim>
      typename Triangulation<dim, spacedim>::raw_hex_iterator
      next_free_hex(const Triangulation<dim, spacedim> &tria,
                    const unsigned int                  level);

      /**
       * Access to user pointers.
       */
      void *&
      user_pointer(const unsigned int i);

      /**
       * Read-only access to user pointers.
       */
      const void *
      user_pointer(const unsigned int i) const;

      /**
       * Access to user indices.
       */
      unsigned int &
      user_index(const unsigned int i);

      /**
       * Read-only access to user pointers.
       */
      unsigned int
      user_index(const unsigned int i) const;

      /**
       * Reset user data to zero.
       */
      void
      clear_user_data(const unsigned int i);

      /**
       * Clear all user pointers or indices and reset their type, such that
       * the next access may be either or.
       */
      void
      clear_user_data();

      /**
       * Clear all user flags.
       */
      void
      clear_user_flags();

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

      /**
       * Triangulation objects can either access a user pointer or a
       * user index. What you tried to do is trying to access one of those
       * after using the other.
       *
       * @ingroup Exceptions
       */
      DeclException0(ExcPointerIndexClash);

      /**
       * Counter for next_free_single_* functions
       */
      unsigned int next_free_single;

      /**
       * Counter for next_free_pair_* functions
       */
      unsigned int next_free_pair;

      /**
       * Bool flag for next_free_single_* functions
       */
      bool reverse_order_next_free_single;

      /**
       * The data type storing user pointers or user indices.
       */
      struct UserData
      {
        union
        {
          /// The entry used as user
          /// pointer.
          void *p;
          /// The entry used as user
          /// index.
          unsigned int i;
        };

        /**
         * Default constructor.
         */
        UserData()
        {
          p = nullptr;
        }

        /**
         * Write the data of this object to a stream for the purpose of
         * serialization using the [BOOST serialization
         * library](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html).
         */
        template <class Archive>
        void
        serialize(Archive &ar, const unsigned int version);
      };

      /**
       * Enum describing the possible types of userdata.
       */
      enum UserDataType
      {
        /// No userdata used yet.
        data_unknown,
        /// UserData contains pointers.
        data_pointer,
        /// UserData contains indices.
        data_index
      };


      /**
       * Pointer which is not used by the library but may be accessed and set
       * by the user to handle data local to a line/quad/etc.
       */
      std::vector<UserData> user_data;

      /**
       * In order to avoid confusion between user pointers and indices, this
       * enum is set by the first function accessing either and subsequent
       * access will not be allowed to change the type of data accessed.
       */
      mutable UserDataType user_data_type;
    };


    //----------------------------------------------------------------------//

    inline unsigned int
    TriaObjects::n_objects() const
    {
      // ensure that sizes are consistent, and then return one that
      // corresponds to the number of objects
      AssertDimension(cells.size(), manifold_id.size() * 2 * this->structdim);
      return manifold_id.size();
    }



    inline ArrayView<int>
    TriaObjects::get_bounding_object_indices(const unsigned int index)
    {
      // assume that each cell has the same number of faces
      const unsigned int faces_per_cell = 2 * this->structdim;
      return ArrayView<int>(cells.data() + index * faces_per_cell,
                            faces_per_cell);
    }



    inline TriaObjects::BoundaryOrMaterialId::BoundaryOrMaterialId()
    {
      material_id = numbers::invalid_material_id;
    }



    inline std::size_t
    TriaObjects::BoundaryOrMaterialId::memory_consumption()
    {
      return sizeof(BoundaryOrMaterialId);
    }



    template <class Archive>
    void
    TriaObjects::BoundaryOrMaterialId::serialize(Archive &ar,
                                                 const unsigned int /*version*/)
    {
      // serialize this
      // structure by
      // writing and
      // reading the larger
      // of the two values,
      // in order to make
      // sure we get all
      // bits
      if (sizeof(material_id) > sizeof(boundary_id))
        ar &material_id;
      else
        ar &boundary_id;
    }


    inline void *&
    TriaObjects::user_pointer(const unsigned int i)
    {
      Assert(user_data_type == data_unknown || user_data_type == data_pointer,
             ExcPointerIndexClash());
      user_data_type = data_pointer;

      AssertIndexRange(i, user_data.size());
      return user_data[i].p;
    }


    inline const void *
    TriaObjects::user_pointer(const unsigned int i) const
    {
      Assert(user_data_type == data_unknown || user_data_type == data_pointer,
             ExcPointerIndexClash());
      user_data_type = data_pointer;

      AssertIndexRange(i, user_data.size());
      return user_data[i].p;
    }


    inline unsigned int &
    TriaObjects::user_index(const unsigned int i)
    {
      Assert(user_data_type == data_unknown || user_data_type == data_index,
             ExcPointerIndexClash());
      user_data_type = data_index;

      AssertIndexRange(i, user_data.size());
      return user_data[i].i;
    }


    inline void
    TriaObjects::clear_user_data(const unsigned int i)
    {
      AssertIndexRange(i, user_data.size());
      user_data[i].i = 0;
    }


    inline TriaObjects::TriaObjects()
      : structdim(numbers::invalid_unsigned_int)
      , next_free_single(numbers::invalid_unsigned_int)
      , next_free_pair(numbers::invalid_unsigned_int)
      , reverse_order_next_free_single(false)
      , user_data_type(data_unknown)
    {}


    inline TriaObjects::TriaObjects(const unsigned int structdim)
      : structdim(structdim)
      , next_free_single(numbers::invalid_unsigned_int)
      , next_free_pair(numbers::invalid_unsigned_int)
      , reverse_order_next_free_single(false)
      , user_data_type(data_unknown)
    {}


    inline unsigned int
    TriaObjects::user_index(const unsigned int i) const
    {
      Assert(user_data_type == data_unknown || user_data_type == data_index,
             ExcPointerIndexClash());
      user_data_type = data_index;

      AssertIndexRange(i, user_data.size());
      return user_data[i].i;
    }


    inline void
    TriaObjects::clear_user_data()
    {
      user_data_type = data_unknown;
      for (auto &data : user_data)
        data.p = nullptr;
    }


    inline void
    TriaObjects::clear_user_flags()
    {
      user_flags.assign(user_flags.size(), false);
    }


    template <class Archive>
    void
    TriaObjects::UserData::serialize(Archive &ar, const unsigned int)
    {
      // serialize this as an integer
      ar &i;
    }



    template <class Archive>
    void
    TriaObjects::serialize(Archive &ar, const unsigned int)
    {
      ar                                   &structdim;
      ar &cells                            &children;
      ar                                   &refinement_cases;
      ar                                   &used;
      ar                                   &user_flags;
      ar                                   &boundary_or_material_id;
      ar                                   &manifold_id;
      ar &next_free_single &next_free_pair &reverse_order_next_free_single;
      ar &user_data                        &user_data_type;
    }


    //----------------------------------------------------------------------//

    template <int structdim_, int dim, int spacedim>
    dealii::TriaRawIterator<dealii::TriaAccessor<structdim_, dim, spacedim>>
    TriaObjects::next_free_single_object(
      const Triangulation<dim, spacedim> &tria)
    {
      // TODO: Think of a way to ensure that we are using the correct
      // triangulation, i.e. the one containing *this.

      AssertDimension(structdim_, this->structdim);

      int pos = next_free_single, last = used.size() - 1;
      if (!reverse_order_next_free_single)
        {
          // first sweep forward, only use really single slots, do not use
          // pair slots
          for (; pos < last; ++pos)
            if (!used[pos])
              if (used[++pos])
                {
                  // this was a single slot
                  pos -= 1;
                  break;
                }
          if (pos >= last)
            {
              reverse_order_next_free_single = true;
              next_free_single               = used.size() - 1;
              pos                            = used.size() - 1;
            }
          else
            next_free_single = pos + 1;
        }

      if (reverse_order_next_free_single)
        {
          // second sweep, use all slots, even
          // in pairs
          for (; pos >= 0; --pos)
            if (!used[pos])
              break;
          if (pos > 0)
            next_free_single = pos - 1;
          else
            // no valid single object anymore
            return dealii::TriaRawIterator<
              dealii::TriaAccessor<structdim_, dim, spacedim>>(&tria, -1, -1);
        }

      return dealii::TriaRawIterator<
        dealii::TriaAccessor<structdim_, dim, spacedim>>(&tria, 0, pos);
    }



    template <int structdim_, int dim, int spacedim>
    dealii::TriaRawIterator<dealii::TriaAccessor<structdim_, dim, spacedim>>
    TriaObjects::next_free_pair_object(const Triangulation<dim, spacedim> &tria)
    {
      // TODO: Think of a way to ensure that we are using the correct
      // triangulation, i.e. the one containing *this.

      AssertDimension(structdim_, this->structdim);

      int pos = next_free_pair, last = used.size() - 1;
      for (; pos < last; ++pos)
        if (!used[pos])
          if (!used[++pos])
            {
              // this was a pair slot
              pos -= 1;
              break;
            }
      if (pos >= last)
        // no free slot
        return dealii::TriaRawIterator<
          dealii::TriaAccessor<structdim_, dim, spacedim>>(&tria, -1, -1);
      else
        next_free_pair = pos + 2;

      return dealii::TriaRawIterator<
        dealii::TriaAccessor<structdim_, dim, spacedim>>(&tria, 0, pos);
    }
  } // namespace TriangulationImplementation
} // namespace internal



DEAL_II_NAMESPACE_CLOSE

#endif
