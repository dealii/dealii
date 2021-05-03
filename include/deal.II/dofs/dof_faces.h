// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2018 by the deal.II authors
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

#ifndef dealii_dof_faces_h
#define dealii_dof_faces_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <deal.II/dofs/dof_objects.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  /**
   * A namespace for internal data structures of the DoFHandler group of
   * classes.
   *
   * @ingroup dofs
   */
  namespace DoFHandlerImplementation
  {
    /**
     *
     * <h4>DoFFaces</h4>
     *
     * These classes are similar to the DoFLevel classes. We here store
     * information that is associated with faces, rather than cells, as this
     * information is independent of the hierarchical structure of cells,
     * which are organized in levels. In 2D we store information on degrees of
     * freedom located on lines whereas in 3D we store information on degrees
     * of freedom located on quads and lines. In 1D we do nothing, as the
     * faces of lines are vertices which are treated separately.
     *
     * Apart from the DoFObjects object containing the data to store (degree
     * of freedom indices) we do not store any data or provide any
     * functionality. However, we do implement a function to determine an
     * estimate of the memory consumption of the contained DoFObjects
     * object(s).
     *
     * The data contained isn't usually directly accessed. Rather, except for
     * some access from the DoFHandler class, access is usually through the
     * DoFAccessor::set_dof_index() and DoFAccessor::dof_index() functions or
     * similar functions of derived classes that in turn access the member
     * variables using the DoFHandler::get_dof_index() and corresponding
     * setter functions. Knowledge of the actual data format is therefore
     * encapsulated to the present hierarchy of classes as well as the
     * dealii::DoFHandler class.
     */
    template <int dim>
    class DoFFaces
    {
    public:
      /**
       * Constructor. This constructor is deleted to prevent the use of this
       * template, as only the specializations should be used
       */
      DoFFaces() = delete;
    };

    /**
     * Store the indices of degrees of freedom on faces in 1D. As these would
     * be vertices, which are treated separately, don't do anything.
     */
    template <>
    class DoFFaces<1>
    {
    public:
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
       *
       */
      template <class Archive>
      void
      serialize(Archive &ar, const unsigned int version);
    };

    /**
     * Store the indices of degrees of freedom on faces in 2D, which are
     * lines.
     */
    template <>
    class DoFFaces<2>
    {
    public:
      /**
       * The object containing the data of DoFs on lines.
       */
      internal::DoFHandlerImplementation::DoFObjects<1> lines;

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

    /**
     * Store the indices of degrees of freedom on faces in 3D, which are
     * quads, additionally also on lines.
     */
    template <>
    class DoFFaces<3>
    {
    public:
      /**
       * The object containing the data of DoFs on lines.
       */
      internal::DoFHandlerImplementation::DoFObjects<1> lines;

      /**
       * The object containing the data of DoFs on quads.
       */
      internal::DoFHandlerImplementation::DoFObjects<2> quads;

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
    DoFFaces<1>::serialize(Archive &, const unsigned int)
    {}


    template <class Archive>
    void
    DoFFaces<2>::serialize(Archive &ar, const unsigned int)
    {
      ar &lines;
    }


    template <class Archive>
    void
    DoFFaces<3>::serialize(Archive &ar, const unsigned int)
    {
      ar &lines &quads;
    }

  } // namespace DoFHandlerImplementation
} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
