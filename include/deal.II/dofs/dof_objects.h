// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_dof_objects_h
#define dealii_dof_objects_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/types.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
class DoFHandler;
#endif

namespace internal
{
  namespace DoFHandlerImplementation
  {
#ifndef DOXYGEN
    template <int>
    class DoFLevel;
    template <int>
    class DoFFaces;
#endif

    /**
     * Store the indices of the degrees of freedom which are located on
     * objects of dimension @p dim.
     *
     * <h3>Information for all DoFObjects classes</h3>
     *
     * The DoFObjects classes store the global indices of the degrees of
     * freedom for each cell on a certain level. The global index or number of
     * a degree of freedom is the zero-based index of the according value in
     * the solution vector and the row and column index in the global matrix
     * or the multigrid matrix for this level. These indices refer to the
     * unconstrained vectors and matrices, where we have not taken account of
     * the constraints introduced by hanging nodes.
     *
     * Since vertices are not associated with a particular level, the indices
     * associated with vertices are not stored in the DoFObjects classes but
     * rather in the DoFHandler::vertex_dofs array.
     *
     * The DoFObjects classes are not used directly, but objects of these
     * classes are included in the DoFLevel and DoFFaces classes.
     *
     * @ingroup dofs
     */
    template <int dim>
    class DoFObjects
    {
    public:
      /**
       * Store the global indices of the degrees of freedom.
       */
      std::vector<types::global_dof_index> dofs;

      /**
       * Return the global index of the @p local_index-th degree of freedom
       * located on the object with number @p obj_index. The @p dof_handler
       * argument is used to access the finite element that is to be used to
       * compute the location where this data is stored.
       *
       * The third argument, @p fe_index, must equal zero. It is otherwise
       * unused, but we retain the argument so that we can use the same
       * interface for non-hp- and hp-finite element methods, in effect making
       * it possible to share the DoFAccessor class hierarchy between hp- and
       * non-hp-classes.
       */
      template <int dh_dim, int spacedim>
      types::global_dof_index &
      access_dof_index(const DoFHandler<dh_dim, spacedim> &dof_handler,
                       const unsigned int                  obj_index,
                       const types::fe_index               fe_index,
                       const unsigned int                  local_index);

      /**
       * Return the value 1. The meaning of this function becomes clear by
       * looking at what the corresponding functions in the classes
       * internal::hp::DoFObjects
       */
      template <int dh_dim, int spacedim>
      unsigned int
      n_active_fe_indices(const DoFHandler<dh_dim, spacedim> &dof_handler,
                          const types::global_dof_index       index) const;

      /**
       * Similar to the function above. Assert that the given index is zero,
       * and then return true.
       */
      template <int dh_dim, int spacedim>
      bool
      fe_index_is_active(const DoFHandler<dh_dim, spacedim> &dof_handler,
                         const types::global_dof_index       index,
                         const types::fe_index               fe_index) const;

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

      // Declare the classes that store levels and faces of DoFs friends so
      // that they can resize arrays.
      template <int>
      friend class DoFLevel;
      template <int>
      friend class DoFFaces;
    };


    // --------------------- template and inline functions ------------------

    template <int dim>
    template <int dh_dim, int spacedim>
    inline unsigned int
    DoFObjects<dim>::n_active_fe_indices(const DoFHandler<dh_dim, spacedim> &,
                                         const types::global_dof_index) const
    {
      return 1;
    }



    template <int dim>
    template <int dh_dim, int spacedim>
    inline bool
    DoFObjects<dim>::fe_index_is_active(const DoFHandler<dh_dim, spacedim> &,
                                        const types::global_dof_index,
                                        const types::fe_index fe_index) const
    {
      (void)fe_index;
      Assert((fe_index == DoFHandler<dh_dim, spacedim>::default_fe_index),
             ExcMessage("Only zero fe_index values are allowed for "
                        "non-hp-DoFHandlers."));
      return true;
    }



    template <int dim>
    template <int dh_dim, int spacedim>
    inline types::global_dof_index &
    DoFObjects<dim>::access_dof_index(
      const DoFHandler<dh_dim, spacedim> &dof_handler,
      const unsigned int                  obj_index,
      const types::fe_index               fe_index,
      const unsigned int                  local_index)
    {
      (void)fe_index;
      Assert((fe_index == DoFHandler<dh_dim, spacedim>::default_fe_index),
             ExcMessage("Only the default FE index is allowed for DoFHandler "
                        "objects without hp capability"));
      AssertIndexRange(local_index,
                       dof_handler.get_fe().template n_dofs_per_object<dim>());
      AssertIndexRange(
        obj_index * dof_handler.get_fe().template n_dofs_per_object<dim>() +
          local_index,
        dofs.size());

      return dofs[obj_index *
                    dof_handler.get_fe().template n_dofs_per_object<dim>() +
                  local_index];
    }


    template <int dim>
    template <class Archive>
    void
    DoFObjects<dim>::serialize(Archive &ar, const unsigned int)
    {
      ar &dofs;
    }

  } // namespace DoFHandlerImplementation
} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
