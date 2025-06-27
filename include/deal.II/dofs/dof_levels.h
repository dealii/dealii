// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1998 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_dof_levels_h
#define dealii_dof_levels_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/types.h>

#include <deal.II/dofs/dof_objects.h>

#include <vector>


DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace DoFHandlerImplementation
  {
    /**
     * Structure for storing degree of freedom information for cells,
     * organized by levels.
     *
     * Note that vertices use a storage scheme that is entirely separate from
     * the one used for the cells. The indices of degrees of freedom located
     * on vertices are therefore not stored here, but rather in member
     * variables of the DoFHandler class.
     *
     * The indices of degrees of freedom located on lower dimensional objects,
     * i.e. on lines for 2d and on quads and lines for 3d are treated
     * similarly than that on cells. However, these geometrical objects, which
     * are called faces as a generalization, are not organised in a
     * hierarchical structure of levels. Therefore, the degrees of freedom
     * located on these objects are stored in separate classes, namely the
     * <tt>DoFFaces</tt> classes.
     *
     * Access to this object is usually through the
     * DoFAccessor::set_dof_index() and DoFAccessor::dof_index() functions or
     * similar functions of derived classes that in turn access the member
     * variables using the DoFHandler::get_dof_index() and corresponding
     * setter functions. Knowledge of the actual data format is therefore
     * encapsulated to the present hierarchy of classes as well as the
     * DoFHandler class.
     */
    template <int dim>
    class DoFLevel
    {
    public:
      /**
       * Cache for the DoF indices on cells. The size of this array equals the
       * number of cells on a given level times selected_fe.n_dofs_per_cell().
       */
      std::vector<types::global_dof_index> cell_dof_indices_cache;

      /**
       * The object containing dof-indices and related access-functions
       */
      DoFObjects<dim> dof_object;

      /**
       * Return a pointer to the beginning of the DoF indices cache for a
       * given cell.
       *
       * @param obj_index The number of the cell we are looking at.
       * @param dofs_per_cell The number of DoFs per cell for this cell.
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
       * purpose of serialization using the [BOOST serialization
       * library](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html).
       */
      template <class Archive>
      void
      serialize(Archive &ar, const unsigned int version);
    };



    template <int dim>
    inline const types::global_dof_index *
    DoFLevel<dim>::get_cell_cache_start(const unsigned int obj_index,
                                        const unsigned int dofs_per_cell) const
    {
      Assert(obj_index * dofs_per_cell + dofs_per_cell <=
               cell_dof_indices_cache.size(),
             ExcInternalError());

      return cell_dof_indices_cache.data() + (obj_index * dofs_per_cell);
    }



    template <int dim>
    inline std::size_t
    DoFLevel<dim>::memory_consumption() const
    {
      return (MemoryConsumption::memory_consumption(cell_dof_indices_cache) +
              MemoryConsumption::memory_consumption(dof_object));
    }


    template <int dim>
    template <class Archive>
    inline void
    DoFLevel<dim>::serialize(Archive &ar, const unsigned int)
    {
      ar &cell_dof_indices_cache;
      ar &dof_object;
    }
  } // namespace DoFHandlerImplementation
} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
