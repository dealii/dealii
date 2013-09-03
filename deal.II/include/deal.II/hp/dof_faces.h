// ---------------------------------------------------------------------
// $Id$
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
#include <deal.II/dofs/dof_levels.h>
#include <deal.II/hp/dof_objects.h>
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
    class DoFFaces;


    /**
     * Store the indices of degrees of freedom on faces in 1D. As these would be vertices, which
     * are treated separately, don't do anything.
     *
     * @ingroup hp
     * @author Tobias Leicht, 2006
     */
    template<>
    class DoFFaces<1>
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
    class DoFFaces<2>
    {
    public:
      /**
       * Indices of DoFs on the lines that bound cells.
       */
      internal::hp::DoFObjects<1> lines;

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
    class DoFFaces<3>
    {
    public:
      /**
       * Indices of DoFs on the lines that form the edges of cells.
       */
      internal::hp::DoFObjects<1> lines;

      /**
       * Indices of DoFs on the quads that bound cells.
       */
      internal::hp::DoFObjects<2> quads;

      /**
       * Determine an estimate for the
       * memory consumption (in bytes)
       * of this object.
       */
      std::size_t memory_consumption () const;
    };

  } // namespace hp

} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
