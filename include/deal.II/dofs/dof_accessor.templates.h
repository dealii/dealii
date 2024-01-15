// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2023 by the deal.II authors
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

#ifndef dealii_dof_accessor_templates_h
#define dealii_dof_accessor_templates_h


#include <deal.II/base/config.h>

#include <deal.II/base/types.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_faces.h>
#include <deal.II/dofs/dof_levels.h>

#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_iterator.templates.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/read_write_vector.h>

#include <boost/container/small_vector.hpp>

#include <limits>
#include <type_traits>
#include <vector>

DEAL_II_NAMESPACE_OPEN


/*------------------------- Functions: DoFAccessor ---------------------------*/


template <int structdim, int dim, int spacedim, bool level_dof_access>
inline DoFAccessor<structdim, dim, spacedim, level_dof_access>::DoFAccessor()
{
  Assert(false, ExcInvalidObject());
}



template <int structdim, int dim, int spacedim, bool level_dof_access>
inline DoFAccessor<structdim, dim, spacedim, level_dof_access>::DoFAccessor(
  const Triangulation<dim, spacedim> *tria,
  const int                           level,
  const int                           index,
  const DoFHandler<dim, spacedim>    *dof_handler)
  : dealii::internal::DoFAccessorImplementation::
      Inheritance<structdim, dim, spacedim>::BaseClass(tria, level, index)
  , dof_handler(const_cast<DoFHandler<dim, spacedim> *>(dof_handler))
{
  Assert(
    tria == nullptr || &dof_handler->get_triangulation() == tria,
    ExcMessage(
      "You can't create a DoF accessor in which the DoFHandler object "
      "uses a different triangulation than the one you pass as argument."));
}



template <int structdim, int dim, int spacedim, bool level_dof_access>
template <int structdim2, int dim2, int spacedim2>
DoFAccessor<structdim, dim, spacedim, level_dof_access>::DoFAccessor(
  const InvalidAccessor<structdim2, dim2, spacedim2> &)
{
  Assert(false, ExcInvalidObject());
}



template <int structdim, int dim, int spacedim, bool level_dof_access>
template <int structdim2, int dim2, int spacedim2, bool level_dof_access2>
inline DoFAccessor<structdim, dim, spacedim, level_dof_access>::DoFAccessor(
  const DoFAccessor<structdim2, dim2, spacedim2, level_dof_access2> &other)
  : BaseClass(other)
  , dof_handler(nullptr)
{
  Assert(false,
         ExcMessage(
           "You are trying to assign iterators that are incompatible. "
           "The reason for incompatibility is that they refer to objects of "
           "different dimensionality (e.g., assigning a line iterator "
           "to a quad iterator)."));
}



template <int structdim, int dim, int spacedim, bool level_dof_access>
template <bool level_dof_access2>
inline DoFAccessor<structdim, dim, spacedim, level_dof_access>::DoFAccessor(
  const DoFAccessor<structdim, dim, spacedim, level_dof_access2> &other)
  : BaseClass(other)
  , dof_handler(const_cast<DoFHandler<dim, spacedim> *>(other.dof_handler))
{}



template <int structdim, int dim, int spacedim, bool level_dof_access>
inline void
DoFAccessor<structdim, dim, spacedim, level_dof_access>::set_dof_handler(
  DoFHandler<dim, spacedim> *dh)
{
  Assert(dh != nullptr, ExcInvalidObject());
  this->dof_handler = dh;
}



template <int structdim, int dim, int spacedim, bool level_dof_access>
inline const DoFHandler<dim, spacedim> &
DoFAccessor<structdim, dim, spacedim, level_dof_access>::get_dof_handler() const
{
  Assert(this->dof_handler != nullptr, ExcInvalidObject());
  return *this->dof_handler;
}



template <int structdim, int dim, int spacedim, bool level_dof_access>
inline void
DoFAccessor<structdim, dim, spacedim, level_dof_access>::copy_from(
  const TriaAccessorBase<structdim, dim, spacedim> &da)
{
  Assert(this->dof_handler != nullptr, ExcInvalidObject());
  BaseClass::copy_from(da);
}



template <int structdim, int dim, int spacedim, bool level_dof_access>
template <bool level_dof_access2>
inline void
DoFAccessor<structdim, dim, spacedim, level_dof_access>::copy_from(
  const DoFAccessor<structdim, dim, spacedim, level_dof_access2> &a)
{
  BaseClass::copy_from(a);
  this->dof_handler = a.dof_handler;
}



template <int structdim, int dim, int spacedim, bool level_dof_access>
template <int structdim2, int dim2, int spacedim2, bool level_dof_access2>
inline bool
DoFAccessor<structdim, dim, spacedim, level_dof_access>::operator==(
  const DoFAccessor<structdim2, dim2, spacedim2, level_dof_access2> &a) const
{
  Assert(structdim == structdim2, ExcCantCompareIterators());
  Assert(this->dof_handler == a.dof_handler, ExcCantCompareIterators());
  return (BaseClass::operator==(a));
}



template <int structdim, int dim, int spacedim, bool level_dof_access>
template <int structdim2, int dim2, int spacedim2, bool level_dof_access2>
inline bool
DoFAccessor<structdim, dim, spacedim, level_dof_access>::operator!=(
  const DoFAccessor<structdim2, dim2, spacedim2, level_dof_access2> &a) const
{
  Assert(structdim == structdim2, ExcCantCompareIterators());
  Assert(this->dof_handler == a.dof_handler, ExcCantCompareIterators());
  return (BaseClass::operator!=(a));
}



template <int structdim, int dim, int spacedim, bool level_dof_access>
inline TriaIterator<DoFAccessor<structdim, dim, spacedim, level_dof_access>>
DoFAccessor<structdim, dim, spacedim, level_dof_access>::child(
  const unsigned int i) const
{
  Assert(static_cast<unsigned int>(this->present_level) <
           this->dof_handler->object_dof_indices.size(),
         ExcMessage("DoFHandler not initialized"));

  TriaIterator<TriaAccessor<structdim, dim, spacedim>> t =
    TriaAccessor<structdim, dim, spacedim>::child(i);

  TriaIterator<DoFAccessor<structdim, dim, spacedim, level_dof_access>> q(
    *t, this->dof_handler);
  return q;
}


namespace internal
{
  namespace DoFAccessorImplementation
  {
    /**
     * Convert an FE index that might contain the right value but also
     * numbers::invalid_fe_index to a right value if needed/possible.
     */
    template <int structdim, int dim, int spacedim, bool level_dof_access>
    types::fe_index
    get_fe_index_or_default(
      const DoFAccessor<structdim, dim, spacedim, level_dof_access> &cell,
      const types::fe_index                                          fe_index)
    {
      if (cell.get_dof_handler().has_hp_capabilities() == false)
        {
          // No hp enabled, and the argument is at its default value -> we
          // can translate to the default active fe index
          Assert(
            (fe_index == numbers::invalid_fe_index) ||
              (fe_index == DoFHandler<dim, spacedim>::default_fe_index),
            ExcMessage(
              "It is not possible to specify a FE index if no hp support is used!"));

          return DoFHandler<dim, spacedim>::default_fe_index;
        }
      else
        {
          // Otherwise: If anything other than the default is provided by
          // the caller, then we should take just that. As an exception, if
          // we are on a cell (rather than a face/edge/vertex), then we know
          // that there is only one active fe index on this cell and we can
          // use that:
          if ((dim == structdim) && (fe_index == numbers::invalid_fe_index))
            {
              AssertDimension(cell.n_active_fe_indices(), 1);

              return cell.nth_active_fe_index(0);
            }

          Assert((fe_index != numbers::invalid_fe_index),
                 ExcMessage(
                   "You need to specify a FE index if hp support is used!"));

          return fe_index;
        }
    }

    /**
     * A class like the one with same name in tria.cc. See there for more
     * information.
     */
    struct Implementation
    {
      /**
       * In several applications of DoFAccessor::get_dof_values(), we want to
       * extract some indices without having to allocate memory. We do this by
       * setting a boost small_vector with 27 elements on the stack, and only
       * allocate when we exceed 27. The number 27 is heuristic and allows up
       * to quadratic shape functions on scalar problems in 3d, or linear
       * shape functions on systems (elasticity).
       */
      using dof_index_vector_type =
        boost::container::small_vector<dealii::types::global_dof_index, 27>;

      /**
       * Process the @p local_index-th degree of freedom corresponding to the
       * finite element specified by @p fe_index on the vertex with global
       * number @p vertex_index to @p global_index.
       *
       * The template argument `structdim` indicates the
       * dimensionality of the object on which we seek to know the DoF
       * index. For example, if `structdim==0`, then we are looking to
       * get a DoF index on a vertex of the indicated cell.
       */
      template <int dim,
                int spacedim,
                int structdim,
                typename GlobalIndexType,
                typename DoFPProcessor>
      static void
      process_dof_index(const DoFHandler<dim, spacedim> &dof_handler,
                        const unsigned int               obj_level,
                        const unsigned int               obj_index,
                        const types::fe_index            fe_index,
                        const unsigned int               local_index,
                        const std::integral_constant<int, structdim> &,
                        GlobalIndexType     &global_index,
                        const DoFPProcessor &process)
      {
        Assert(structdim == dim || obj_level == 0, ExcNotImplemented());

        // 1) no hp used -> fe_index == 0
        if (dof_handler.hp_capability_enabled == false)
          {
            AssertDimension(fe_index,
                            (DoFHandler<dim, spacedim>::default_fe_index));

            process(
              dof_handler.object_dof_indices
                [obj_level][structdim]
                [dof_handler.object_dof_ptr[obj_level][structdim][obj_index] +
                 local_index],
              global_index);

            return;
          }

        // 2) cell and hp is used -> there is only one fe_index
        if (structdim == dim)
          {
            process(
              dof_handler.object_dof_indices
                [obj_level][structdim]
                [dof_handler.object_dof_ptr[obj_level][structdim][obj_index] +
                 local_index],
              global_index);
            return;
          }

        // 3) general entity and hp is used
        AssertIndexRange(obj_level, dof_handler.object_dof_indices.size());
        AssertIndexRange(structdim,
                         dof_handler.object_dof_indices[obj_level].size());

        Assert(dof_handler.hp_capability_enabled, ExcInternalError());

        AssertIndexRange(structdim, dof_handler.hp_object_fe_ptr.size());
        AssertIndexRange(obj_index,
                         dof_handler.hp_object_fe_ptr[structdim].size());

        const auto ptr =
          std::find(dof_handler.hp_object_fe_indices[structdim].begin() +
                      dof_handler.hp_object_fe_ptr[structdim][obj_index],
                    dof_handler.hp_object_fe_indices[structdim].begin() +
                      dof_handler.hp_object_fe_ptr[structdim][obj_index + 1],
                    fe_index);

        Assert(ptr != dof_handler.hp_object_fe_indices[structdim].begin() +
                        dof_handler.hp_object_fe_ptr[structdim][obj_index + 1],
               ExcMessage(
                 "You are requesting an active FE index that is not assigned "
                 "to any of the cells connected to this entity."));

        const types::fe_index fe_index_ =
          std::distance(dof_handler.hp_object_fe_indices[structdim].begin() +
                          dof_handler.hp_object_fe_ptr[structdim][obj_index],
                        ptr);

        AssertIndexRange(
          dof_handler.hp_capability_enabled ?
            (dof_handler.hp_object_fe_ptr[structdim][obj_index] + fe_index_) :
            obj_index,
          dof_handler.object_dof_ptr[obj_level][structdim].size());

        AssertIndexRange(
          dof_handler.object_dof_ptr
              [obj_level][structdim]
              [dof_handler.hp_capability_enabled ?
                 (dof_handler.hp_object_fe_ptr[structdim][obj_index] +
                  fe_index_) :
                 obj_index] +
            local_index,
          dof_handler.object_dof_indices[obj_level][structdim].size());

        process(dof_handler.object_dof_indices
                  [obj_level][structdim]
                  [dof_handler.object_dof_ptr
                     [obj_level][structdim]
                     [dof_handler.hp_capability_enabled ?
                        (dof_handler.hp_object_fe_ptr[structdim][obj_index] +
                         fe_index_) :
                        obj_index] +
                   local_index],
                global_index);
      }

      /**
       * Determine start index and number of dofs of object in global data
       * structure.
       *
       * The template argument `structdim` indicates the
       * dimensionality of the object on which we seek to know the DoF
       * index. For example, if `structdim==0`, then we are looking to
       * get a DoF index on a vertex of the indicated cell.
       */
      template <int dim, int spacedim, int structdim>
      static std::pair<unsigned int, unsigned int>
      process_object_range(const DoFHandler<dim, spacedim> &dof_handler,
                           const unsigned int               obj_level,
                           const unsigned int               obj_index,
                           const types::fe_index            fe_index,
                           const std::integral_constant<int, structdim> &)
      {
        Assert(structdim == dim || obj_level == 0, ExcNotImplemented());

        // determine range of dofs in global data structure
        // 1) cell
        if (structdim == dim)
          {
            const unsigned int ptr_0 =
              dof_handler.object_dof_ptr[obj_level][structdim][obj_index];
            const unsigned int length =
              dof_handler.get_fe(fe_index).template n_dofs_per_object<dim>(0);

            return {ptr_0, length};
          }

        // 2) hp is not used -> fe_index == 0
        if (dof_handler.hp_capability_enabled == false)
          {
            AssertDimension(fe_index,
                            (DoFHandler<dim, spacedim>::default_fe_index));

            const unsigned int ptr_0 =
              dof_handler.object_dof_ptr[obj_level][structdim][obj_index];
            const unsigned int length =
              dof_handler.object_dof_ptr[obj_level][structdim][obj_index + 1] -
              ptr_0;

            return {ptr_0, length};
          }

        // 3) hp is used
        AssertIndexRange(obj_level, dof_handler.object_dof_indices.size());
        AssertIndexRange(structdim,
                         dof_handler.object_dof_indices[obj_level].size());

        AssertIndexRange(structdim, dof_handler.hp_object_fe_ptr.size());
        AssertIndexRange(obj_index,
                         dof_handler.hp_object_fe_ptr[structdim].size());

        const auto fe_index_local_ptr =
          std::find(dof_handler.hp_object_fe_indices[structdim].begin() +
                      dof_handler.hp_object_fe_ptr[structdim][obj_index],
                    dof_handler.hp_object_fe_indices[structdim].begin() +
                      dof_handler.hp_object_fe_ptr[structdim][obj_index + 1],
                    fe_index);

        Assert(fe_index_local_ptr !=
                 dof_handler.hp_object_fe_indices[structdim].begin() +
                   dof_handler.hp_object_fe_ptr[structdim][obj_index + 1],
               ExcMessage("Call distribute_dofs() first."));

        const types::fe_index fe_index_local =
          std::distance(dof_handler.hp_object_fe_indices[structdim].begin() +
                          dof_handler.hp_object_fe_ptr[structdim][obj_index],
                        fe_index_local_ptr);

        AssertIndexRange(
          dof_handler.hp_object_fe_ptr[structdim][obj_index] + fe_index_local,
          dof_handler.object_dof_ptr[obj_level][structdim].size());

        const unsigned int ptr_0 =
          dof_handler
            .object_dof_ptr[obj_level][structdim]
                           [dof_handler.hp_object_fe_ptr[structdim][obj_index] +
                            fe_index_local];
        const unsigned int ptr_1 =
          dof_handler
            .object_dof_ptr[obj_level][structdim]
                           [dof_handler.hp_object_fe_ptr[structdim][obj_index] +
                            fe_index_local + 1];

        return {ptr_0, ptr_1 - ptr_0};
      }

      template <int dim, int spacedim, int structdim, bool level_dof_access>
      static std::pair<unsigned int, unsigned int>
      process_object_range(
        const dealii::DoFAccessor<structdim, dim, spacedim, level_dof_access>
                              accessor,
        const types::fe_index fe_index)
      {
        return process_object_range(accessor.get_dof_handler(),
                                    accessor.level(),
                                    accessor.index(),
                                    fe_index,
                                    std::integral_constant<int, structdim>());
      }

      template <int dim, int spacedim, int structdim>
      static std::pair<unsigned int, unsigned int>
      process_object_range(dealii::DoFInvalidAccessor<structdim, dim, spacedim>,
                           const unsigned int)
      {
        Assert(false, ExcInternalError());

        return {0, 0};
      }



      /**
       * Process all dofs of an object.
       *
       * The template argument `structdim` indicates the
       * dimensionality of the object on which we seek to know the DoF
       * index. For example, if `structdim==0`, then we are looking to
       * get a DoF index on a vertex of the indicated cell.
       */
      template <int dim,
                int spacedim,
                int structdim,
                typename DoFProcessor,
                typename DoFMapping>
      static DEAL_II_ALWAYS_INLINE void
      process_object(const DoFHandler<dim, spacedim>              &dof_handler,
                     const unsigned int                            obj_level,
                     const unsigned int                            obj_index,
                     const types::fe_index                         fe_index,
                     const DoFMapping                             &mapping,
                     const std::integral_constant<int, structdim> &dd,
                     types::global_dof_index *&dof_indices_ptr,
                     const DoFProcessor       &process)
      {
        Assert(structdim == dim || obj_level == 0, ExcNotImplemented());

        // determine range of dofs in global data structure
        const auto range =
          process_object_range(dof_handler, obj_level, obj_index, fe_index, dd);
        if (range.second == 0)
          return;

        std::vector<types::global_dof_index> &object_dof_indices =
          dof_handler
            .object_dof_indices[structdim < dim ? 0 : obj_level][structdim];
        AssertIndexRange(range.first, object_dof_indices.size());
        types::global_dof_index *DEAL_II_RESTRICT stored_indices =
          object_dof_indices.data() + range.first;

        // process dofs
        for (unsigned int i = 0; i < range.second; ++i, ++dof_indices_ptr)
          process(
            stored_indices[(structdim == 0 || structdim == dim) ? i :
                                                                  mapping(i)],
            dof_indices_ptr);
      }



      /**
       * Set the @p local_index-th degree of freedom corresponding to the
       * finite element specified by @p fe_index on the vertex with global
       * number @p vertex_index to @p global_index.
       *
       * The template argument `structdim` indicates the
       * dimensionality of the object on which we seek to know the DoF
       * index. For example, if `structdim==0`, then we are looking to
       * get a DoF index on a vertex of the indicated cell.
       */
      template <int dim, int spacedim, int structdim>
      static void
      set_dof_index(const DoFHandler<dim, spacedim>              &dof_handler,
                    const unsigned int                            obj_level,
                    const unsigned int                            obj_index,
                    const types::fe_index                         fe_index,
                    const unsigned int                            local_index,
                    const std::integral_constant<int, structdim> &dd,
                    const types::global_dof_index                 global_index)
      {
        process_dof_index(dof_handler,
                          obj_level,
                          obj_index,
                          fe_index,
                          local_index,
                          dd,
                          global_index,
                          [](auto &ptr, const auto &value) { ptr = value; });
      }


      /**
       * Get the @p local_index-th degree of freedom corresponding to the
       * finite element specified by @p fe_index on the vertex with global
       * number @p vertex_index to @p global_index.
       *
       * The template argument `structdim` indicates the
       * dimensionality of the object on which we seek to know the DoF
       * index. For example, if `structdim==0`, then we are looking to
       * get a DoF index on a vertex of the indicated cell.
       */
      template <int dim, int spacedim, int structdim>
      static types::global_dof_index
      get_dof_index(const DoFHandler<dim, spacedim>              &dof_handler,
                    const unsigned int                            obj_level,
                    const unsigned int                            obj_index,
                    const types::fe_index                         fe_index,
                    const unsigned int                            local_index,
                    const std::integral_constant<int, structdim> &dd)
      {
        types::global_dof_index global_index;
        process_dof_index(dof_handler,
                          obj_level,
                          obj_index,
                          fe_index,
                          local_index,
                          dd,
                          global_index,
                          [](const auto &ptr, auto &value) { value = ptr; });
        return global_index;
      }


      template <int dim, int spacedim>
      static types::global_dof_index &
      mg_vertex_dof_index(DoFHandler<dim, spacedim> &dof_handler,
                          const int                  level,
                          const unsigned int         vertex_index,
                          const unsigned int         i)
      {
        Assert(dof_handler.hp_capability_enabled == false,
               ExcMessage(
                 "DoFHandler in hp-mode does not implement multilevel DoFs."));

        return dof_handler.mg_vertex_dofs[vertex_index].access_index(
          level, i, dof_handler.get_fe().n_dofs_per_vertex());
      }



      /**
       * Return the number of different finite elements that are active on a
       * given object such as a vertex, line, or cell.
       *
       * The template argument `structdim` indicates the
       * dimensionality of the object on which we seek to know the DoF
       * index. For example, if `structdim==0`, then we are looking to
       * get a DoF index on a vertex of the indicated cell.
       */
      template <int dim, int spacedim, int structdim>
      static unsigned int
      n_active_fe_indices(const DoFHandler<dim, spacedim> &dof_handler,
                          const unsigned int               obj_level,
                          const unsigned int               obj_index,
                          const std::integral_constant<int, structdim> &)
      {
        (void)obj_level;

        Assert(structdim == dim || obj_level == 0, ExcNotImplemented());

        // 1) no hp used -> fe_index == 0
        if (dof_handler.hp_capability_enabled == false)
          return 1;

        // 2) cell and hp is used -> there is only one fe_index
        if (structdim == dim)
          return 1;

        // 3) general entity and hp is used
        AssertIndexRange(structdim, dof_handler.hp_object_fe_ptr.size());
        AssertIndexRange(obj_index + 1,
                         dof_handler.hp_object_fe_ptr[structdim].size());

        return dof_handler.hp_object_fe_ptr[structdim][obj_index + 1] -
               dof_handler.hp_object_fe_ptr[structdim][obj_index];
      }



      /**
       * Return the FE index of the n-th finite element active on a given
       * object such as a vertex, line, or cell.
       *
       * The template argument `structdim` indicates the
       * dimensionality of the object on which we seek to know the DoF
       * index. For example, if `structdim==0`, then we are looking to
       * get a DoF index on a vertex of the indicated cell.
       */
      template <int dim, int spacedim, int structdim>
      static types::fe_index
      nth_active_fe_index(const DoFHandler<dim, spacedim> &dof_handler,
                          const unsigned int               obj_level,
                          const unsigned int               obj_index,
                          const unsigned int               local_index,
                          const std::integral_constant<int, structdim> &)
      {
        Assert(structdim == dim || obj_level == 0, ExcNotImplemented());

        // for cells only one active FE index available
        Assert(((structdim == dim) &&
                (local_index != DoFHandler<dim, spacedim>::default_fe_index)) ==
                 false,
               ExcNotImplemented());

        // 1) no hp used -> fe_index == 0
        if (dof_handler.hp_capability_enabled == false)
          return DoFHandler<dim, spacedim>::default_fe_index;

        // 2) cell and hp is used -> there is only one fe_index
        if (structdim == dim)
          return dof_handler.hp_cell_active_fe_indices[obj_level][obj_index];

        // 3) general entity and hp is used
        AssertIndexRange(structdim, dof_handler.hp_object_fe_indices.size());
        AssertIndexRange(structdim, dof_handler.hp_object_fe_ptr.size());
        AssertIndexRange(obj_index,
                         dof_handler.hp_object_fe_ptr[structdim].size());
        AssertIndexRange(dof_handler.hp_object_fe_ptr[structdim][obj_index] +
                           local_index,
                         dof_handler.hp_object_fe_indices[structdim].size());

        return dof_handler.hp_object_fe_indices
          [structdim]
          [dof_handler.hp_object_fe_ptr[structdim][obj_index] + local_index];
      }



      /**
       * Returns all active FE indices on a given object such as a
       * vertex, line, or cell.
       *
       * The size of the returned set equals the number of finite elements that
       * are active on this vertex.
       *
       * The template argument `structdim` indicates the
       * dimensionality of the object on which we seek to know the DoF
       * index. For example, if `structdim==0`, then we are looking to
       * get a DoF index on a vertex of the indicated cell.
       */
      template <int dim, int spacedim, int structdim>
      static std::set<types::fe_index>
      get_active_fe_indices(const DoFHandler<dim, spacedim> &dof_handler,
                            const unsigned int               obj_level,
                            const unsigned int               obj_index,
                            const std::integral_constant<int, structdim> &t)
      {
        Assert(structdim == dim || obj_level == 0, ExcNotImplemented());

        // 1) no hp used -> fe_index == 0
        if (dof_handler.hp_capability_enabled == false)
          return {DoFHandler<dim, spacedim>::default_fe_index};

        // 2) cell and hp is used -> there is only one fe_index
        if (structdim == dim)
          return {dof_handler.hp_cell_active_fe_indices[obj_level][obj_index]};

        // 3) general entity and hp is used
        std::set<types::fe_index> active_fe_indices;
        for (unsigned int i = 0;
             i < n_active_fe_indices(dof_handler, obj_level, obj_index, t);
             ++i)
          active_fe_indices.insert(
            nth_active_fe_index(dof_handler, obj_level, obj_index, i, t));
        return active_fe_indices;
      }



      template <int dim, int spacedim, int structdim>
      static bool
      fe_index_is_active(const DoFHandler<dim, spacedim> &dof_handler,
                         const unsigned int               obj_level,
                         const unsigned int               obj_index,
                         const types::fe_index            fe_index,
                         const std::integral_constant<int, structdim> &)
      {
        Assert(structdim == dim || obj_level == 0, ExcNotImplemented());

        // 1) no hp used -> fe_index == 0
        if (dof_handler.hp_capability_enabled == false)
          return (fe_index == DoFHandler<dim, spacedim>::default_fe_index);

        // 2) cell and hp is used -> there is only one fe_index
        if (structdim == dim)
          return dof_handler.hp_cell_active_fe_indices[obj_level][obj_index] ==
                 fe_index;

        // 3) general entity and hp is used
        return std::find(
                 dof_handler.hp_object_fe_indices[structdim].begin() +
                   dof_handler.hp_object_fe_ptr[structdim][obj_index],
                 dof_handler.hp_object_fe_indices[structdim].begin() +
                   dof_handler.hp_object_fe_ptr[structdim][obj_index + 1],
                 fe_index) !=
               (dof_handler.hp_object_fe_indices[structdim].begin() +
                dof_handler.hp_object_fe_ptr[structdim][obj_index + 1]);
      }



      template <typename InputVector, typename ForwardIterator>
      static void
      extract_subvector_to(const InputVector             &values,
                           const types::global_dof_index *cache,
                           const types::global_dof_index *cache_end,
                           ForwardIterator                local_values_begin)
      {
        values.extract_subvector_to(cache, cache_end, local_values_begin);
      }



#ifdef DEAL_II_WITH_TRILINOS
      static std::vector<unsigned int>
      sort_indices(const types::global_dof_index *v_begin,
                   const types::global_dof_index *v_end)
      {
        // initialize original index locations
        std::vector<unsigned int> idx(v_end - v_begin);
        std::iota(idx.begin(), idx.end(), 0u);

        // sort indices based on comparing values in v
        std::sort(idx.begin(),
                  idx.end(),
                  [&v_begin](unsigned int i1, unsigned int i2) {
                    return *(v_begin + i1) < *(v_begin + i2);
                  });

        return idx;
      }



#  ifdef DEAL_II_TRILINOS_WITH_TPETRA
      template <typename ForwardIterator, typename Number>
      static void
      extract_subvector_to(
        const LinearAlgebra::TpetraWrappers::Vector<Number> &values,
        const types::global_dof_index                       *cache_begin,
        const types::global_dof_index                       *cache_end,
        ForwardIterator                                      local_values_begin)
      {
        std::vector<unsigned int> sorted_indices_pos =
          sort_indices(cache_begin, cache_end);
        const unsigned int cache_size = cache_end - cache_begin;
        std::vector<types::global_dof_index> cache_indices(cache_size);
        for (unsigned int i = 0; i < cache_size; ++i)
          cache_indices[i] = *(cache_begin + sorted_indices_pos[i]);

        IndexSet index_set(cache_indices.back() + 1);
        index_set.add_indices(cache_indices.begin(), cache_indices.end());
        index_set.compress();
        LinearAlgebra::ReadWriteVector<Number> read_write_vector(index_set);
        read_write_vector.import_elements(values, VectorOperation::insert);

        // Copy the elements from read_write_vector and reorder them.
        for (unsigned int i = 0; i < cache_size; ++i, ++local_values_begin)
          *local_values_begin = read_write_vector[sorted_indices_pos[i]];
      }
#  endif



      template <typename ForwardIterator>
      static void
      extract_subvector_to(const LinearAlgebra::EpetraWrappers::Vector &values,
                           const types::global_dof_index *cache_begin,
                           const types::global_dof_index *cache_end,
                           ForwardIterator                local_values_begin)
      {
        std::vector<unsigned int> sorted_indices_pos =
          sort_indices(cache_begin, cache_end);
        const unsigned int cache_size = cache_end - cache_begin;
        std::vector<types::global_dof_index> cache_indices(cache_size);
        for (unsigned int i = 0; i < cache_size; ++i)
          cache_indices[i] = *(cache_begin + sorted_indices_pos[i]);

        IndexSet index_set(cache_indices.back() + 1);
        index_set.add_indices(cache_indices.begin(), cache_indices.end());
        index_set.compress();
        LinearAlgebra::ReadWriteVector<double> read_write_vector(index_set);
        read_write_vector.import_elements(values, VectorOperation::insert);

        // Copy the elements from read_write_vector and reorder them.
        for (unsigned int i = 0; i < cache_size; ++i, ++local_values_begin)
          *local_values_begin = read_write_vector[sorted_indices_pos[i]];
      }
#endif

      /**
       * Loop over all degrees of freedom of the object described by the
       * provided @p accessor and @p fe_index and count them.
       */
      template <int dim, int spacedim, bool level_dof_access, int structdim>
      static unsigned int
      n_dof_indices(
        const dealii::DoFAccessor<structdim, dim, spacedim, level_dof_access>
                             &accessor,
        const types::fe_index fe_index_,
        const bool            count_level_dofs)
      {
        // note: we cannot rely on the template parameter level_dof_access here,
        // since the function get_mg_dof_indices()/set_mg_dof_indices() can be
        // called even if level_dof_access==false.
        if (count_level_dofs)
          {
            const auto &fe = accessor.get_fe(fe_index_);

            const unsigned int                                   //
              dofs_per_vertex = fe.n_dofs_per_vertex(),          //
              dofs_per_line   = fe.n_dofs_per_line(),            //
              dofs_per_quad   = fe.n_dofs_per_quad(0 /*dummy*/), //
              dofs_per_hex    = fe.n_dofs_per_hex();             //

            unsigned int index = 0;

            // 1) VERTEX dofs
            index += dofs_per_vertex * accessor.n_vertices();

            // 2) LINE dofs
            if (structdim == 2 || structdim == 3)
              index += dofs_per_line * accessor.n_lines();

            // 3) FACE dofs
            if (structdim == 3)
              index += dofs_per_quad * accessor.n_faces();

            // 4) INNER dofs
            const unsigned int interior_dofs =
              structdim == 1 ? dofs_per_line :
                               (structdim == 2 ? dofs_per_quad : dofs_per_hex);

            index += interior_dofs;

            return index;
          }
        else
          {
            const auto fe_index =
              internal::DoFAccessorImplementation::get_fe_index_or_default(
                accessor, fe_index_);

            unsigned int index = 0;

            // 1) VERTEX dofs
            for (const auto vertex : accessor.vertex_indices())
              index += process_object_range(accessor.get_dof_handler(),
                                            0,
                                            accessor.vertex_index(vertex),
                                            fe_index,
                                            std::integral_constant<int, 0>())
                         .second;

            // 2) LINE dofs
            if (structdim == 2 || structdim == 3)
              for (const auto line : accessor.line_indices())
                index +=
                  process_object_range(*accessor.line(line), fe_index).second;

            // 3) FACE dofs
            if (structdim == 3)
              for (const auto face : accessor.face_indices())
                index +=
                  process_object_range(*accessor.quad(face), fe_index).second;

            // 4) INNER dofs
            index += process_object_range(accessor, fe_index).second;

            return index;
          }
      }



      // The next few internal helper functions are needed to support various
      // DoFIndicesType kinds, e.g. actual vectors of DoFIndices or empty
      // types that we use when we only want to work on the internally stored
      // DoFs and never extract any number.
      template <typename ArrayType>
      static unsigned int
      get_array_length(const ArrayType &array)
      {
        return array.size();
      }

      static unsigned int
      get_array_length(const std::tuple<> &)
      {
        return 0;
      }

      template <typename ArrayType>
      static types::global_dof_index *
      get_array_ptr(const ArrayType &array)
      {
        return const_cast<types::global_dof_index *>(array.data());
      }

      static types::global_dof_index *
      get_array_ptr(const std::tuple<> &)
      {
        return nullptr;
      }



      /**
       * Loop over all degrees of freedom of the object described by the
       * provided @p accessor and @p fe_index and perform the static functions
       * provided by DoFOperation (set/get) on these.
       */
      template <int  dim,
                int  spacedim,
                bool level_dof_access,
                int  structdim,
                typename DoFIndicesType,
                typename DoFOperation,
                typename DoFProcessor>
      static void
      process_dof_indices(
        const dealii::DoFAccessor<structdim, dim, spacedim, level_dof_access>
                             &accessor,
        const DoFIndicesType &const_dof_indices,
        const types::fe_index fe_index_,
        const DoFOperation   &dof_operation,
        const DoFProcessor   &dof_processor,
        const bool            count_level_dofs)
      {
        const types::fe_index fe_index =
          internal::DoFAccessorImplementation::get_fe_index_or_default(
            accessor, fe_index_);

        // we cannot rely on the template parameter level_dof_access here, since
        // the function get_mg_dof_indices()/set_mg_dof_indices() can be called
        // even if level_dof_access==false.
        (void)count_level_dofs;

        const auto &fe = accessor.get_fe(fe_index);

        // we want to pass in rvalue 'std::tuple<>' types as `DoFIndicesType`,
        // but we need non-const references for std::vector<> types, so get in
        // a const reference here and immediately cast the constness away -
        // note that any use of the dereferenced invalid type will result in a
        // segfault
        types::global_dof_index *dof_indices_ptr =
          get_array_ptr(const_dof_indices);
        types::global_dof_index *end_dof_indices =
          dof_indices_ptr + get_array_length(const_dof_indices);

        // 1) VERTEX dofs, only step into the functions if we actually have
        // DoFs on them
        if (fe.n_dofs_per_vertex() > 0)
          for (const auto vertex : accessor.vertex_indices())
            dof_operation.process_vertex_dofs(*accessor.dof_handler,
                                              accessor.vertex_index(vertex),
                                              fe_index,
                                              dof_indices_ptr,
                                              dof_processor);

        // 2) copy dof numbers from the LINE. for lines with the wrong
        // orientation (which might occur in 3d), we have already made sure that
        // we're ok by picking the correct vertices (this happens automatically
        // in the vertex() function). however, if the line is in wrong
        // orientation, we look at it in flipped orientation and we will have to
        // adjust the shape function indices that we see to correspond to the
        // correct (face/cell-local) ordering.
        if ((structdim == 2 || structdim == 3) && fe.n_dofs_per_line() > 0)
          {
            const auto line_indices = internal::TriaAccessorImplementation::
              Implementation::get_line_indices_of_cell(accessor);
            const auto line_orientations =
              internal::TriaAccessorImplementation::Implementation::
                get_line_orientations_of_cell(accessor);

            for (const auto line : accessor.line_indices())
              {
                const bool line_orientation = line_orientations[line];
                if (line_orientation)
                  dof_operation.process_dofs(
                    accessor.get_dof_handler(),
                    0,
                    line_indices[line],
                    fe_index,
                    [](const auto d) { return d; },
                    std::integral_constant<int, 1>(),
                    dof_indices_ptr,
                    dof_processor);
                else
                  dof_operation.process_dofs(
                    accessor.get_dof_handler(),
                    0,
                    line_indices[line],
                    fe_index,
                    [&fe, line_orientation](const auto d) {
                      return fe.adjust_line_dof_index_for_line_orientation(
                        d, line_orientation);
                    },
                    std::integral_constant<int, 1>(),
                    dof_indices_ptr,
                    dof_processor);
              }
          }

        // 3) copy dof numbers from the FACE. for faces with the wrong
        // orientation, we have already made sure that we're ok by picking the
        // correct lines and vertices (this happens automatically in the line()
        // and vertex() functions). however, if the face is in wrong
        // orientation, we look at it in flipped orientation and we will have to
        // adjust the shape function indices that we see to correspond to the
        // correct (cell-local) ordering. The same applies, if the face_rotation
        // or face_orientation is non-standard
        if (structdim == 3 && fe.max_dofs_per_quad() > 0)
          for (const auto face_no : accessor.face_indices())
            {
              const auto combined_orientation = TriaAccessorImplementation::
                Implementation::combined_face_orientation(accessor, face_no);
              const unsigned int quad_index = accessor.quad_index(face_no);
              if (combined_orientation ==
                  ReferenceCell::default_combined_face_orientation())
                dof_operation.process_dofs(
                  accessor.get_dof_handler(),
                  0,
                  quad_index,
                  fe_index,
                  [](const auto d) { return d; },
                  std::integral_constant<int, 2>(),
                  dof_indices_ptr,
                  dof_processor);
              else
                dof_operation.process_dofs(
                  accessor.get_dof_handler(),
                  0,
                  quad_index,
                  fe_index,
                  [&](const auto d) {
                    return fe.adjust_quad_dof_index_for_face_orientation(
                      d, face_no, combined_orientation);
                  },
                  std::integral_constant<int, 2>(),
                  dof_indices_ptr,
                  dof_processor);
            }

        // 4) INNER dofs - here we need to make sure that the shortcut to not
        // run the function does not miss the faces of wedge and pyramid
        // elements where n_dofs_per_object might not return the largest
        // possible value
        if (((dim == 3 && structdim == 2) ?
               fe.max_dofs_per_quad() :
               fe.template n_dofs_per_object<structdim>()) > 0)
          dof_operation.process_dofs(
            accessor.get_dof_handler(),
            accessor.level(),
            accessor.index(),
            fe_index,
            [&](const auto d) { return d; },
            std::integral_constant<int, structdim>(),
            dof_indices_ptr,
            dof_processor);

        AssertDimension(n_dof_indices(accessor, fe_index, count_level_dofs),
                        dof_indices_ptr - get_array_ptr(const_dof_indices));

        // PM: This is a part that should not be reached since it indicates that
        // an object (and/or its subobjects) is not active. However,
        // unfortunately this function is called by
        // DoFTools::set_periodicity_constraints() indirectly by
        // get_dof_indices() also for artificial faces to determine if a face
        // is artificial.
        types::global_dof_index invalid_index = numbers::invalid_dof_index;
        for (; dof_indices_ptr < end_dof_indices; ++dof_indices_ptr)
          dof_processor(invalid_index, dof_indices_ptr);
      }



      /**
       * An internal struct encapsulating the task of getting (vertex)
       * DoF indices.
       */
      template <int dim, int spacedim>
      struct DoFIndexProcessor
      {
        /**
         * Return vertex DoF indices.
         */
        template <typename DoFProcessor>
        DEAL_II_ALWAYS_INLINE void
        process_vertex_dofs(DoFHandler<dim, spacedim> &dof_handler,
                            const unsigned int         vertex_index,
                            const types::fe_index      fe_index,
                            types::global_dof_index  *&dof_indices_ptr,
                            const DoFProcessor        &dof_processor) const
        {
          process_object(
            dof_handler,
            0,
            vertex_index,
            fe_index,
            [](const auto d) {
              Assert(false, ExcInternalError());
              return d;
            },
            std::integral_constant<int, 0>(),
            dof_indices_ptr,
            dof_processor);
        }

        /**
         * Return DoF indices for lines, quads, and inner degrees of freedom.
         */
        template <int structdim, typename DoFMapping, typename DoFProcessor>
        DEAL_II_ALWAYS_INLINE void
        process_dofs(const DoFHandler<dim, spacedim> &dof_handler,
                     const unsigned int               obj_level,
                     const unsigned int               obj_index,
                     const types::fe_index            fe_index,
                     const DoFMapping                &mapping,
                     const std::integral_constant<int, structdim>,
                     types::global_dof_index *&dof_indices_ptr,
                     const DoFProcessor       &dof_processor) const
        {
          process_object(
            dof_handler,
            obj_level,
            obj_index,
            fe_index,
            mapping,
            std::integral_constant<int, std::min(structdim, dim)>(),
            dof_indices_ptr,
            dof_processor);
        }
      };



      /**
       * An internal struct encapsulating the task of getting level (vertex)
       * DoF indices.
       */
      template <int dim, int spacedim>
      struct MGDoFIndexProcessor
      {
        /**
         * Constructor.
         */
        MGDoFIndexProcessor(const unsigned int level)
          : level(level)
        {}

        /**
         * Return vertex DoF indices.
         */
        template <typename DoFProcessor>
        DEAL_II_ALWAYS_INLINE void
        process_vertex_dofs(DoFHandler<dim, spacedim> &dof_handler,
                            const unsigned int         vertex_index,
                            const types::fe_index,
                            types::global_dof_index *&dof_indices_ptr,
                            const DoFProcessor       &dof_processor) const
        {
          const unsigned int n_indices =
            dof_handler.get_fe(0).template n_dofs_per_object<0>();
          types::global_dof_index *stored_indices =
            &dof_handler.mg_vertex_dofs[vertex_index].access_index(level,
                                                                   0,
                                                                   n_indices);
          for (unsigned int d = 0; d < n_indices; ++d, ++dof_indices_ptr)
            dof_processor(stored_indices[d], dof_indices_ptr);
        }

        /**
         * Return DoF indices for lines, quads, and inner degrees of freedom.
         */
        template <int structdim, typename DoFMapping, typename DoFProcessor>
        DEAL_II_ALWAYS_INLINE void
        process_dofs(const DoFHandler<dim, spacedim> &dof_handler,
                     const unsigned int,
                     const unsigned int    obj_index,
                     const types::fe_index fe_index,
                     const DoFMapping     &mapping,
                     const std::integral_constant<int, structdim>,
                     types::global_dof_index *&dof_indices_ptr,
                     const DoFProcessor       &dof_processor) const
        {
          const unsigned int n_indices =
            dof_handler.get_fe(0).template n_dofs_per_object<structdim>();
          types::global_dof_index *stored_indices = &get_mg_dof_index(
            dof_handler,
            dof_handler.mg_levels[level],
            dof_handler.mg_faces,
            obj_index,
            fe_index,
            0,
            std::integral_constant<int, std::min(structdim, dim)>());
          for (unsigned int d = 0; d < n_indices; ++d, ++dof_indices_ptr)
            dof_processor(stored_indices[structdim < dim ? mapping(d) : d],
                          dof_indices_ptr);
        }

      private:
        const unsigned int level;
      };



      template <int dim, int spacedim, bool level_dof_access, int structdim>
      static void
      get_dof_indices(
        const dealii::DoFAccessor<structdim, dim, spacedim, level_dof_access>
                                             &accessor,
        std::vector<types::global_dof_index> &dof_indices,
        const types::fe_index                 fe_index)
      {
        process_dof_indices(
          accessor,
          dof_indices,
          fe_index,
          DoFIndexProcessor<dim, spacedim>(),
          [](auto stored_index, auto dof_ptr) { *dof_ptr = stored_index; },
          false);
      }



      template <int dim, int spacedim, bool level_dof_access, int structdim>
      static void
      set_dof_indices(
        const dealii::DoFAccessor<structdim, dim, spacedim, level_dof_access>
                                                   &accessor,
        const std::vector<types::global_dof_index> &dof_indices,
        const types::fe_index                       fe_index)
      {
        // Note: this function is as general as `get_dof_indices()`. This
        // assert is placed here since it is currently only used by the
        // function DoFCellAccessor::set_dof_indices(), which is called by
        // internal::DoFHandlerImplementation::Policy::Implementation::distribute_dofs().
        // In the case of new use cases, this assert can be removed.
        Assert(
          dim == structdim,
          ExcMessage(
            "This function is intended to be used for DoFCellAccessor, i.e., "
            "dimension == structdim."));

        process_dof_indices(
          accessor,
          dof_indices,
          fe_index,
          DoFIndexProcessor<dim, spacedim>(),
          [](auto &stored_index, auto dof_ptr) { stored_index = *dof_ptr; },
          false);
      }



      template <int dim, int spacedim, bool level_dof_access, int structdim>
      static void
      get_mg_dof_indices(
        const dealii::DoFAccessor<structdim, dim, spacedim, level_dof_access>
                                             &accessor,
        const int                             level,
        std::vector<types::global_dof_index> &dof_indices,
        const types::fe_index                 fe_index)
      {
        Assert((fe_index == DoFHandler<dim, spacedim>::default_fe_index),
               ExcMessage("MG DoF indices cannot be queried in hp case"));
        process_dof_indices(
          accessor,
          dof_indices,
          fe_index,
          MGDoFIndexProcessor<dim, spacedim>(level),
          [](auto stored_index, auto dof_ptr) { *dof_ptr = stored_index; },
          true);
      }



      template <int dim, int spacedim, bool level_dof_access, int structdim>
      static void
      set_mg_dof_indices(
        const dealii::DoFAccessor<structdim, dim, spacedim, level_dof_access>
                                                   &accessor,
        const int                                   level,
        const std::vector<types::global_dof_index> &dof_indices,
        const types::fe_index                       fe_index)
      {
        Assert((fe_index == DoFHandler<dim, spacedim>::default_fe_index),
               ExcMessage("MG DoF indices cannot be queried in hp case"));

        // Note: this function is as general as `get_mg_dof_indices()`. This
        // assert is placed here since it is currently only used by the
        // function DoFCellAccessor::set_mg_dof_indices(), which is called by
        // internal::DoFHandlerImplementation::Policy::Implementation::distribute_mg_dofs().
        // In the case of new use cases, this assert can be removed.
        Assert(dim == structdim,
               ExcMessage("This function is intended to be used for "
                          "DoFCellAccessor, i.e., dimension == structdim."));

        process_dof_indices(
          accessor,
          dof_indices,
          fe_index,
          MGDoFIndexProcessor<dim, spacedim>(level),
          [](auto &stored_index, auto dof_ptr) { stored_index = *dof_ptr; },
          true);
      }



      template <int dim, int spacedim>
      static types::global_dof_index &
      get_mg_dof_index(
        const DoFHandler<dim, spacedim> &dof_handler,
        const std::unique_ptr<internal::DoFHandlerImplementation::DoFLevel<dim>>
          &mg_level,
        const std::unique_ptr<internal::DoFHandlerImplementation::DoFFaces<dim>>
          &,
        const unsigned int    obj_index,
        const types::fe_index fe_index,
        const unsigned int    local_index,
        const std::integral_constant<int, dim>)
      {
        Assert(dof_handler.hp_capability_enabled == false,
               (typename DoFHandler<dim, spacedim>::ExcNotImplementedWithHP()));

        return mg_level->dof_object.access_dof_index(
          static_cast<const DoFHandler<dim, spacedim> &>(dof_handler),
          obj_index,
          fe_index,
          local_index);
      }



      template <int dim, int spacedim, std::enable_if_t<(dim > 1), int> = 0>
      static types::global_dof_index &
      get_mg_dof_index(
        const DoFHandler<dim, spacedim> &dof_handler,
        const std::unique_ptr<internal::DoFHandlerImplementation::DoFLevel<dim>>
          &,
        const std::unique_ptr<internal::DoFHandlerImplementation::DoFFaces<dim>>
                             &mg_faces,
        const unsigned int    obj_index,
        const types::fe_index fe_index,
        const unsigned int    local_index,
        const std::integral_constant<int, 1>)
      {
        return mg_faces->lines.access_dof_index(
          static_cast<const DoFHandler<dim, spacedim> &>(dof_handler),
          obj_index,
          fe_index,
          local_index);
      }



      template <int spacedim>
      static types::global_dof_index &
      get_mg_dof_index(
        const DoFHandler<3, spacedim> &dof_handler,
        const std::unique_ptr<internal::DoFHandlerImplementation::DoFLevel<3>>
          &,
        const std::unique_ptr<internal::DoFHandlerImplementation::DoFFaces<3>>
                             &mg_faces,
        const unsigned int    obj_index,
        const types::fe_index fe_index,
        const unsigned int    local_index,
        const std::integral_constant<int, 2>)
      {
        Assert(dof_handler.hp_capability_enabled == false,
               (typename DoFHandler<3, spacedim>::ExcNotImplementedWithHP()));
        return mg_faces->quads.access_dof_index(
          static_cast<const DoFHandler<3, spacedim> &>(dof_handler),
          obj_index,
          fe_index,
          local_index);
      }
    };



    template <int dim, int spacedim, bool level_dof_access>
    void
    get_cell_dof_indices(
      const dealii::DoFCellAccessor<dim, spacedim, level_dof_access> &accessor,
      Implementation::dof_index_vector_type &dof_indices,
      const unsigned int                     fe_index);
  } // namespace DoFAccessorImplementation
} // namespace internal



template <int structdim, int dim, int spacedim, bool level_dof_access>
inline types::global_dof_index
DoFAccessor<structdim, dim, spacedim, level_dof_access>::dof_index(
  const unsigned int    i,
  const types::fe_index fe_index_) const
{
  const auto fe_index =
    internal::DoFAccessorImplementation::get_fe_index_or_default(*this,
                                                                 fe_index_);

  // access the respective DoF
  return dealii::internal::DoFAccessorImplementation::Implementation::
    get_dof_index(*this->dof_handler,
                  this->level(),
                  this->index(),
                  fe_index,
                  i,
                  std::integral_constant<int, structdim>());
}


template <int structdim, int dim, int spacedim, bool level_dof_access>
inline types::global_dof_index
DoFAccessor<structdim, dim, spacedim, level_dof_access>::mg_dof_index(
  const int          level,
  const unsigned int i) const
{
  return internal::DoFAccessorImplementation::Implementation::get_mg_dof_index(
    *this->dof_handler,
    this->dof_handler->mg_levels[level],
    this->dof_handler->mg_faces,
    this->index(),
    0,
    i,
    std::integral_constant<int, structdim>());
}


template <int structdim, int dim, int spacedim, bool level_dof_access>
inline void
DoFAccessor<structdim, dim, spacedim, level_dof_access>::set_dof_index(
  const unsigned int            i,
  const types::global_dof_index index,
  const types::fe_index         fe_index_) const
{
  const auto fe_index =
    internal::DoFAccessorImplementation::get_fe_index_or_default(*this,
                                                                 fe_index_);

  // access the respective DoF
  dealii::internal::DoFAccessorImplementation::Implementation::set_dof_index(
    *this->dof_handler,
    this->level(),
    this->index(),
    fe_index,
    i,
    std::integral_constant<int, structdim>(),
    index);
}



template <int structdim, int dim, int spacedim, bool level_dof_access>
inline void
DoFAccessor<structdim, dim, spacedim, level_dof_access>::set_mg_dof_index(
  const int                     level,
  const unsigned int            i,
  const types::global_dof_index index) const
{
  internal::DoFAccessorImplementation::Implementation::get_mg_dof_index(
    *this->dof_handler,
    this->dof_handler->mg_levels[level],
    this->dof_handler->mg_faces,
    this->index(),
    0,
    i,
    std::integral_constant<int, structdim>()) = index;
}



template <int structdim, int dim, int spacedim, bool level_dof_access>
inline unsigned int
DoFAccessor<structdim, dim, spacedim, level_dof_access>::n_active_fe_indices()
  const
{
  // access the respective DoF
  return dealii::internal::DoFAccessorImplementation::Implementation::
    n_active_fe_indices(*this->dof_handler,
                        this->level(),
                        this->index(),
                        std::integral_constant<int, structdim>());
}



template <int structdim, int dim, int spacedim, bool level_dof_access>
inline types::fe_index
DoFAccessor<structdim, dim, spacedim, level_dof_access>::nth_active_fe_index(
  const unsigned int n) const
{
  // access the respective DoF
  return dealii::internal::DoFAccessorImplementation::Implementation::
    nth_active_fe_index(*this->dof_handler,
                        this->level(),
                        this->index(),
                        n,
                        std::integral_constant<int, structdim>());
}



template <int structdim, int dim, int spacedim, bool level_dof_access>
inline std::set<types::fe_index>
DoFAccessor<structdim, dim, spacedim, level_dof_access>::get_active_fe_indices()
  const
{
  std::set<types::fe_index> active_fe_indices;
  for (unsigned int i = 0; i < n_active_fe_indices(); ++i)
    active_fe_indices.insert(nth_active_fe_index(i));
  return active_fe_indices;
}



template <int structdim, int dim, int spacedim, bool level_dof_access>
inline bool
DoFAccessor<structdim, dim, spacedim, level_dof_access>::fe_index_is_active(
  const types::fe_index fe_index) const
{
  // access the respective DoF
  return dealii::internal::DoFAccessorImplementation::Implementation::
    fe_index_is_active(*this->dof_handler,
                       this->level(),
                       this->index(),
                       fe_index,
                       std::integral_constant<int, structdim>());
}



template <int structdim, int dim, int spacedim, bool level_dof_access>
inline types::global_dof_index
DoFAccessor<structdim, dim, spacedim, level_dof_access>::vertex_dof_index(
  const unsigned int    vertex,
  const unsigned int    i,
  const types::fe_index fe_index_) const
{
  const types::fe_index fe_index =
    (((this->dof_handler->hp_capability_enabled == false) &&
      (fe_index_ == numbers::invalid_fe_index)) ?
       // No hp enabled, and the argument is at its default value -> we
       // can translate to the default active fe index
       DoFHandler<dim, spacedim>::default_fe_index :
       // Otherwise: If anything other than the default is provided by
       // the caller, then we should take just that. As an exception, if
       // we are on a cell (rather than a face/edge/vertex), then we know
       // that there is only one active fe index on this cell and we can
       // use that:
       ((dim == structdim) && (fe_index_ == numbers::invalid_fe_index) ?
          this->nth_active_fe_index(0) :
          fe_index_));

  return dealii::internal::DoFAccessorImplementation::Implementation::
    get_dof_index(*this->dof_handler,
                  0,
                  this->vertex_index(vertex),
                  fe_index,
                  i,
                  std::integral_constant<int, 0>());
}


template <int structdim, int dim, int spacedim, bool level_dof_access>
inline types::global_dof_index
DoFAccessor<structdim, dim, spacedim, level_dof_access>::mg_vertex_dof_index(
  const int             level,
  const unsigned int    vertex,
  const unsigned int    i,
  const types::fe_index fe_index_) const
{
  const auto fe_index =
    internal::DoFAccessorImplementation::get_fe_index_or_default(*this,
                                                                 fe_index_);
  (void)fe_index;
  Assert(this->dof_handler != nullptr, ExcInvalidObject());
  Assert(this->dof_handler->mg_vertex_dofs.size() > 0,
         ExcMessage("Multigrid DoF indices can only be accessed after "
                    "DoFHandler::distribute_mg_dofs() has been called!"));
  AssertIndexRange(vertex, this->n_vertices());
  AssertIndexRange(i, this->dof_handler->get_fe(fe_index).n_dofs_per_vertex());

  Assert(dof_handler->hp_capability_enabled == false,
         ExcMessage(
           "DoFHandler in hp-mode does not implement multilevel DoFs."));

  return this->dof_handler->mg_vertex_dofs[this->vertex_index(vertex)]
    .access_index(level, i, this->dof_handler->get_fe().n_dofs_per_vertex());
}



template <int structdim, int dim, int spacedim, bool level_dof_access>
inline void
DoFAccessor<structdim, dim, spacedim, level_dof_access>::
  set_mg_vertex_dof_index(const int                     level,
                          const unsigned int            vertex,
                          const unsigned int            i,
                          const types::global_dof_index index,
                          const types::fe_index         fe_index_) const
{
  const auto fe_index =
    internal::DoFAccessorImplementation::get_fe_index_or_default(*this,
                                                                 fe_index_);
  (void)fe_index;
  Assert(this->dof_handler != nullptr, ExcInvalidObject());
  AssertIndexRange(vertex, this->n_vertices());
  AssertIndexRange(i, this->dof_handler->get_fe(fe_index).n_dofs_per_vertex());

  Assert(dof_handler->hp_capability_enabled == false,
         ExcMessage(
           "DoFHandler in hp-mode does not implement multilevel DoFs."));

  this->dof_handler->mg_vertex_dofs[this->vertex_index(vertex)].access_index(
    level, i, this->dof_handler->get_fe().n_dofs_per_vertex()) = index;
}



template <int structdim, int dim, int spacedim, bool level_dof_access>
inline const FiniteElement<dim, spacedim> &
DoFAccessor<structdim, dim, spacedim, level_dof_access>::get_fe(
  const types::fe_index fe_index) const
{
  Assert(fe_index_is_active(fe_index) == true,
         ExcMessage("This function can only be called for active FE indices"));

  return this->dof_handler->get_fe(fe_index);
}



template <int structdim, int dim, int spacedim, bool level_dof_access>
inline typename dealii::internal::DoFHandlerImplementation::
  Iterators<dim, spacedim, level_dof_access>::line_iterator
  DoFAccessor<structdim, dim, spacedim, level_dof_access>::line(
    const unsigned int i) const
{
  // if we are asking for a particular line and this object refers to
  // a line, then the only valid index is i==0 and we should return
  // *this
  if (structdim == 1)
    {
      Assert(i == 0,
             ExcMessage("You can only ask for line zero if the "
                        "current object is a line itself."));
      return typename dealii::internal::DoFHandlerImplementation::
        Iterators<dim, spacedim, level_dof_access>::cell_iterator(
          &this->get_triangulation(),
          this->level(),
          this->index(),
          &this->get_dof_handler());
    }

  // otherwise we need to be in structdim>=2
  Assert(structdim > 1, ExcImpossibleInDim(structdim));
  Assert(dim > 1, ExcImpossibleInDim(dim));

  // checking of 'i' happens in line_index(i)
  return typename dealii::internal::DoFHandlerImplementation::
    Iterators<dim, spacedim, level_dof_access>::line_iterator(
      this->tria,
      0, // only sub-objects are allowed, which have no level
      this->line_index(i),
      this->dof_handler);
}


template <int structdim, int dim, int spacedim, bool level_dof_access>
inline typename dealii::internal::DoFHandlerImplementation::
  Iterators<dim, spacedim, level_dof_access>::quad_iterator
  DoFAccessor<structdim, dim, spacedim, level_dof_access>::quad(
    const unsigned int i) const
{
  // if we are asking for a
  // particular quad and this object
  // refers to a quad, then the only
  // valid index is i==0 and we
  // should return *this
  if (structdim == 2)
    {
      Assert(i == 0,
             ExcMessage("You can only ask for quad zero if the "
                        "current object is a quad itself."));
      return typename dealii::internal::DoFHandlerImplementation::
        Iterators<dim, spacedim>::cell_iterator(&this->get_triangulation(),
                                                this->level(),
                                                this->index(),
                                                &this->get_dof_handler());
    }

  // otherwise we need to be in structdim>=3
  Assert(structdim > 2, ExcImpossibleInDim(structdim));
  Assert(dim > 2, ExcImpossibleInDim(dim));

  // checking of 'i' happens in quad_index(i)
  return typename dealii::internal::DoFHandlerImplementation::
    Iterators<dim, spacedim, level_dof_access>::quad_iterator(
      this->tria,
      0, // only sub-objects are allowed, which have no level
      this->quad_index(i),
      this->dof_handler);
}


/*----------------- Functions: DoFAccessor<0,1,spacedim> --------------------*/


template <int spacedim, bool level_dof_access>
inline DoFAccessor<0, 1, spacedim, level_dof_access>::DoFAccessor()
{
  Assert(false, ExcInvalidObject());
}



template <int spacedim, bool level_dof_access>
inline DoFAccessor<0, 1, spacedim, level_dof_access>::DoFAccessor(
  const Triangulation<1, spacedim>                       *tria,
  const typename TriaAccessor<0, 1, spacedim>::VertexKind vertex_kind,
  const unsigned int                                      vertex_index,
  const DoFHandler<1, spacedim>                          *dof_handler)
  : BaseClass(tria, vertex_kind, vertex_index)
  , dof_handler(const_cast<DoFHandler<1, spacedim> *>(dof_handler))
{}



template <int spacedim, bool level_dof_access>
inline DoFAccessor<0, 1, spacedim, level_dof_access>::DoFAccessor(
  const Triangulation<1, spacedim> *tria,
  const int                         level,
  const int                         index,
  const DoFHandler<1, spacedim>    *dof_handler)
  // This is the constructor signature for "ordinary" (non-vertex)
  // accessors and we shouldn't be calling it altogether. But it is also
  // the constructor that the default-constructor of TriaRawIterator
  // calls when default-constructing an iterator object. If so, this
  // happens with level==-2 and index==-2, and this is the only case we
  // would like to support. We do this by just forwarding to the
  // other constructor of this class, and then asserting the condition
  // on level and index.
  : DoFAccessor<0, 1, spacedim, level_dof_access>(
      tria,
      TriaAccessor<0, 1, spacedim>::interior_vertex,
      0U,
      dof_handler)
{
  (void)level;
  (void)index;
  Assert((tria == nullptr) && (level == -2) && (index == -2) &&
           (dof_handler == nullptr),
         ExcMessage(
           "This constructor can not be called for face iterators in 1d, "
           "except to default-construct iterator objects."));
}



template <int spacedim, bool level_dof_access>
template <int structdim2, int dim2, int spacedim2>
DoFAccessor<0, 1, spacedim, level_dof_access>::DoFAccessor(
  const InvalidAccessor<structdim2, dim2, spacedim2> &)
{
  Assert(false, ExcInvalidObject());
}



template <int spacedim, bool level_dof_access>
template <int structdim2, int dim2, int spacedim2, bool level_dof_access2>
inline DoFAccessor<0, 1, spacedim, level_dof_access>::DoFAccessor(
  const DoFAccessor<structdim2, dim2, spacedim2, level_dof_access2> &)
{
  Assert(false, ExcInvalidObject());
}



template <int spacedim, bool level_dof_access>
inline void
DoFAccessor<0, 1, spacedim, level_dof_access>::set_dof_handler(
  DoFHandler<1, spacedim> *dh)
{
  Assert(dh != nullptr, ExcInvalidObject());
  this->dof_handler = dh;
}



template <int spacedim, bool level_dof_access>
inline void
DoFAccessor<0, 1, spacedim, level_dof_access>::set_dof_index(
  const unsigned int /*i*/,
  const types::global_dof_index /*index*/,
  const types::fe_index /*fe_index*/) const
{
  Assert(false, ExcNotImplemented());
}



template <int spacedim, bool level_dof_access>
inline const DoFHandler<1, spacedim> &
DoFAccessor<0, 1, spacedim, level_dof_access>::get_dof_handler() const
{
  return *this->dof_handler;
}



template <int spacedim, bool level_dof_access>
inline void
DoFAccessor<0, 1, spacedim, level_dof_access>::get_dof_indices(
  std::vector<types::global_dof_index> &dof_indices,
  const types::fe_index                 fe_index_) const
{
  const auto fe_index =
    internal::DoFAccessorImplementation::get_fe_index_or_default(*this,
                                                                 fe_index_);

  for (unsigned int i = 0; i < dof_indices.size(); ++i)
    dof_indices[i] = dealii::internal::DoFAccessorImplementation::
      Implementation::get_dof_index(*dof_handler,
                                    0,
                                    this->global_vertex_index,
                                    fe_index,
                                    i,
                                    std::integral_constant<int, 0>());
}



template <int spacedim, bool level_dof_access>
inline void
DoFAccessor<0, 1, spacedim, level_dof_access>::get_mg_dof_indices(
  const int                             level,
  std::vector<types::global_dof_index> &dof_indices,
  const types::fe_index                 fe_index_) const
{
  const auto fe_index =
    internal::DoFAccessorImplementation::get_fe_index_or_default(*this,
                                                                 fe_index_);
  (void)fe_index;
  AssertDimension(fe_index, (DoFHandler<1, spacedim>::default_fe_index));

  for (unsigned int i = 0; i < dof_indices.size(); ++i)
    dof_indices[i] =
      dealii::internal::DoFAccessorImplementation::Implementation::
        mg_vertex_dof_index(*dof_handler, level, this->global_vertex_index, i);
}



template <int spacedim, bool level_dof_access>
inline types::global_dof_index
DoFAccessor<0, 1, spacedim, level_dof_access>::vertex_dof_index(
  const unsigned int    vertex,
  const unsigned int    i,
  const types::fe_index fe_index_) const
{
  const auto fe_index =
    internal::DoFAccessorImplementation::get_fe_index_or_default(*this,
                                                                 fe_index_);

  (void)vertex;
  AssertIndexRange(vertex, 1);
  return dealii::internal::DoFAccessorImplementation::Implementation::
    get_dof_index(*dof_handler,
                  0,
                  this->global_vertex_index,
                  fe_index,
                  i,
                  std::integral_constant<int, 0>());
}



template <int spacedim, bool level_dof_access>
inline types::global_dof_index
DoFAccessor<0, 1, spacedim, level_dof_access>::dof_index(
  const unsigned int    i,
  const types::fe_index fe_index_) const
{
  const auto fe_index =
    internal::DoFAccessorImplementation::get_fe_index_or_default(*this,
                                                                 fe_index_);

  return dealii::internal::DoFAccessorImplementation::Implementation::
    get_dof_index(*this->dof_handler,
                  0,
                  this->vertex_index(0),
                  fe_index,
                  i,
                  std::integral_constant<int, 0>());
}



template <int spacedim, bool level_dof_access>
inline unsigned int
DoFAccessor<0, 1, spacedim, level_dof_access>::n_active_fe_indices() const
{
  return 1;
}



template <int spacedim, bool level_dof_access>
inline types::fe_index
DoFAccessor<0, 1, spacedim, level_dof_access>::nth_active_fe_index(
  const unsigned int /*n*/) const
{
  return 0;
}



template <int spacedim, bool level_dof_access>
inline bool
DoFAccessor<0, 1, spacedim, level_dof_access>::fe_index_is_active(
  const types::fe_index /*fe_index*/) const
{
  Assert(false, ExcNotImplemented());
  return false;
}



template <int spacedim, bool level_dof_access>
inline const FiniteElement<1, spacedim> &
DoFAccessor<0, 1, spacedim, level_dof_access>::get_fe(
  const types::fe_index fe_index) const
{
  Assert(this->dof_handler != nullptr, ExcInvalidObject());
  return dof_handler->get_fe(fe_index);
}



template <int spacedim, bool level_dof_access>
inline void
DoFAccessor<0, 1, spacedim, level_dof_access>::copy_from(
  const TriaAccessorBase<0, 1, spacedim> &da)
{
  Assert(this->dof_handler != nullptr, ExcInvalidObject());
  BaseClass::copy_from(da);
}



template <int spacedim, bool level_dof_access>
template <bool level_dof_access2>
inline void
DoFAccessor<0, 1, spacedim, level_dof_access>::copy_from(
  const DoFAccessor<0, 1, spacedim, level_dof_access2> &a)
{
  BaseClass::copy_from(a);
  set_dof_handler(a.dof_handler);
}



template <int spacedim, bool level_dof_access>
inline TriaIterator<DoFAccessor<0, 1, spacedim, level_dof_access>>
DoFAccessor<0, 1, spacedim, level_dof_access>::child(
  const unsigned int /*i*/) const
{
  return TriaIterator<DoFAccessor<0, 1, spacedim, level_dof_access>>();
}



template <int spacedim, bool level_dof_access>
inline typename dealii::internal::DoFHandlerImplementation::
  Iterators<1, spacedim, level_dof_access>::line_iterator
  DoFAccessor<0, 1, spacedim, level_dof_access>::line(
    const unsigned int /*c*/) const
{
  Assert(false, ExcNotImplemented());
  return typename dealii::internal::DoFHandlerImplementation::
    Iterators<1, spacedim, level_dof_access>::line_iterator();
}



template <int spacedim, bool level_dof_access>
inline typename dealii::internal::DoFHandlerImplementation::
  Iterators<1, spacedim, level_dof_access>::quad_iterator
  DoFAccessor<0, 1, spacedim, level_dof_access>::quad(
    const unsigned int /*c*/) const
{
  Assert(false, ExcNotImplemented());
  return typename dealii::internal::DoFHandlerImplementation::
    Iterators<1, spacedim, level_dof_access>::quad_iterator();
}



template <int spacedim, bool level_dof_access>
template <int structdim2, int dim2, int spacedim2, bool level_dof_access2>
inline bool
DoFAccessor<0, 1, spacedim, level_dof_access>::operator==(
  const DoFAccessor<structdim2, dim2, spacedim2, level_dof_access2> &a) const
{
  Assert(structdim2 == 0, ExcCantCompareIterators());
  Assert(this->dof_handler == a.dof_handler, ExcCantCompareIterators());
  return (BaseClass::operator==(a));
}



template <int spacedim, bool level_dof_access>
template <int structdim2, int dim2, int spacedim2, bool level_dof_access2>
inline bool
DoFAccessor<0, 1, spacedim, level_dof_access>::operator!=(
  const DoFAccessor<structdim2, dim2, spacedim2, level_dof_access2> &a) const
{
  Assert(structdim2 == 0, ExcCantCompareIterators());
  Assert(this->dof_handler == a.dof_handler, ExcCantCompareIterators());
  return (BaseClass::operator!=(a));
}



/*------------------------- Functions: DoFCellAccessor -----------------------*/


namespace internal
{
  namespace DoFCellAccessorImplementation
  {
    /**
     * A class with the same purpose as the similarly named class of the
     * Triangulation class. See there for more information.
     */
    struct Implementation
    {
      /**
       * Do what the active_fe_index function in the parent class is supposed to
       * do.
       */
      template <int dim, int spacedim, bool level_dof_access>
      static types::fe_index
      active_fe_index(
        const DoFCellAccessor<dim, spacedim, level_dof_access> &accessor)
      {
        if (accessor.dof_handler->hp_capability_enabled == false)
          return DoFHandler<dim, spacedim>::default_fe_index;

        Assert(accessor.dof_handler != nullptr,
               (typename std::decay_t<decltype(accessor)>::ExcInvalidObject()));
        Assert(static_cast<unsigned int>(accessor.level()) <
                 accessor.dof_handler->hp_cell_future_fe_indices.size(),
               ExcMessage("DoFHandler not initialized"));

        return accessor.dof_handler
          ->hp_cell_active_fe_indices[accessor.level()][accessor.present_index];
      }



      /**
       * Do what the set_active_fe_index function in the parent class is
       * supposed to do.
       */
      template <int dim, int spacedim, bool level_dof_access>
      static void
      set_active_fe_index(
        const DoFCellAccessor<dim, spacedim, level_dof_access> &accessor,
        const types::fe_index                                   i)
      {
        if (accessor.dof_handler->hp_capability_enabled == false)
          {
            AssertDimension(i, (DoFHandler<dim, spacedim>::default_fe_index));
            return;
          }

        Assert(accessor.dof_handler != nullptr,
               (typename std::decay_t<decltype(accessor)>::ExcInvalidObject()));
        Assert(static_cast<unsigned int>(accessor.level()) <
                 accessor.dof_handler->hp_cell_future_fe_indices.size(),
               ExcMessage("DoFHandler not initialized"));
        Assert(i != numbers::invalid_fe_index,
               ExcMessage("Invalid finite element index."));

        accessor.dof_handler
          ->hp_cell_active_fe_indices[accessor.level()]
                                     [accessor.present_index] = i;
      }



      /**
       * Do what the future_fe_index function in the parent class is supposed to
       * do.
       */
      template <int dim, int spacedim, bool level_dof_access>
      static types::fe_index
      future_fe_index(
        const DoFCellAccessor<dim, spacedim, level_dof_access> &accessor)
      {
        if (accessor.dof_handler->hp_capability_enabled == false)
          return DoFHandler<dim, spacedim>::default_fe_index;

        Assert(accessor.dof_handler != nullptr,
               (typename std::decay_t<decltype(accessor)>::ExcInvalidObject()));
        Assert(static_cast<unsigned int>(accessor.level()) <
                 accessor.dof_handler->hp_cell_future_fe_indices.size(),
               ExcMessage("DoFHandler not initialized"));

        if (future_fe_index_set(accessor))
          return accessor.dof_handler
            ->hp_cell_future_fe_indices[accessor.level()]
                                       [accessor.present_index];
        else
          return accessor.dof_handler
            ->hp_cell_active_fe_indices[accessor.level()]
                                       [accessor.present_index];
      }


      /**
       * Do what the set_future_fe_index function in the parent class is
       * supposed to do.
       */
      template <int dim, int spacedim, bool level_dof_access>
      static void
      set_future_fe_index(
        const DoFCellAccessor<dim, spacedim, level_dof_access> &accessor,
        const types::fe_index                                   i)
      {
        if (accessor.dof_handler->hp_capability_enabled == false)
          {
            AssertDimension(i, (DoFHandler<dim, spacedim>::default_fe_index));
            return;
          }

        Assert(accessor.dof_handler != nullptr,
               (typename std::decay_t<decltype(accessor)>::ExcInvalidObject()));
        Assert(static_cast<unsigned int>(accessor.level()) <
                 accessor.dof_handler->hp_cell_future_fe_indices.size(),
               ExcMessage("DoFHandler not initialized"));
        Assert(i != numbers::invalid_fe_index,
               ExcMessage("Invalid finite element index."));

        accessor.dof_handler
          ->hp_cell_future_fe_indices[accessor.level()]
                                     [accessor.present_index] = i;
      }



      /**
       * Do what the future_fe_index_set function in the parent class is
       * supposed to do.
       */
      template <int dim, int spacedim, bool level_dof_access>
      static bool
      future_fe_index_set(
        const DoFCellAccessor<dim, spacedim, level_dof_access> &accessor)
      {
        if (accessor.dof_handler->hp_capability_enabled == false)
          return false;

        Assert(accessor.dof_handler != nullptr,
               (typename std::decay_t<decltype(accessor)>::ExcInvalidObject()));
        Assert(static_cast<unsigned int>(accessor.level()) <
                 accessor.dof_handler->hp_cell_future_fe_indices.size(),
               ExcMessage("DoFHandler not initialized"));

        return accessor.dof_handler
                 ->hp_cell_future_fe_indices[accessor.level()]
                                            [accessor.present_index] !=
               numbers::invalid_fe_index;
      }



      /**
       * Do what the clear_fe_index function in the parent class is supposed to
       * do.
       */
      template <int dim, int spacedim, bool level_dof_access>
      static void
      clear_future_fe_index(
        const DoFCellAccessor<dim, spacedim, level_dof_access> &accessor)
      {
        if (accessor.dof_handler->hp_capability_enabled == false)
          return;

        Assert(accessor.dof_handler != nullptr,
               (typename std::decay_t<decltype(accessor)>::ExcInvalidObject()));
        Assert(static_cast<unsigned int>(accessor.level()) <
                 accessor.dof_handler->hp_cell_future_fe_indices.size(),
               ExcMessage("DoFHandler not initialized"));

        accessor.dof_handler
          ->hp_cell_future_fe_indices[accessor.level()]
                                     [accessor.present_index] =
          numbers::invalid_fe_index;
      }
    };
  } // namespace DoFCellAccessorImplementation
} // namespace internal



template <int dimension_, int space_dimension_, bool level_dof_access>
inline DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::
  DoFCellAccessor(const Triangulation<dimension_, space_dimension_> *tria,
                  const int                                          level,
                  const int                                          index,
                  const AccessorData                                *local_data)
  : DoFAccessor<dimension_, dimension_, space_dimension_, level_dof_access>(
      tria,
      level,
      index,
      local_data)
{}



template <int dimension_, int space_dimension_, bool level_dof_access>
template <int structdim2, int dim2, int spacedim2>
inline DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::
  DoFCellAccessor(const InvalidAccessor<structdim2, dim2, spacedim2> &)
{
  Assert(false, typename BaseClass::ExcInvalidObject());
}



template <int dimension_, int space_dimension_, bool level_dof_access>
template <int structdim2, int dim2, int spacedim2, bool level_dof_access2>
inline DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::
  DoFCellAccessor(
    const DoFAccessor<structdim2, dim2, spacedim2, level_dof_access2> &other)
  : BaseClass(other)
{}



template <int dimension_, int space_dimension_, bool level_dof_access>
inline TriaIterator<
  DoFCellAccessor<dimension_, space_dimension_, level_dof_access>>
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::neighbor(
  const unsigned int i) const
{
  TriaIterator<DoFCellAccessor<dimension_, space_dimension_, level_dof_access>>
    q(this->tria,
      this->neighbor_level(i),
      this->neighbor_index(i),
      this->dof_handler);

#ifdef DEBUG
  if (q.state() != IteratorState::past_the_end)
    Assert(q->used(), ExcInternalError());
#endif
  return q;
}



template <int dimension_, int space_dimension_, bool level_dof_access>
inline TriaIterator<
  DoFCellAccessor<dimension_, space_dimension_, level_dof_access>>
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::child(
  const unsigned int i) const
{
  TriaIterator<DoFCellAccessor<dimension_, space_dimension_, level_dof_access>>
    q(this->tria, this->level() + 1, this->child_index(i), this->dof_handler);

#ifdef DEBUG
  if (q.state() != IteratorState::past_the_end)
    Assert(q->used(), ExcInternalError());
#endif
  return q;
}



template <int dimension_, int space_dimension_, bool level_dof_access>
inline boost::container::small_vector<
  TriaIterator<DoFCellAccessor<dimension_, space_dimension_, level_dof_access>>,
  GeometryInfo<dimension_>::max_children_per_cell>
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::
  child_iterators() const
{
  boost::container::small_vector<
    TriaIterator<
      DoFCellAccessor<dimension_, space_dimension_, level_dof_access>>,
    GeometryInfo<dimension_>::max_children_per_cell>
    child_iterators(this->n_children());

  for (unsigned int i = 0; i < this->n_children(); ++i)
    child_iterators[i] = this->child(i);

  return child_iterators;
}



template <int dimension_, int space_dimension_, bool level_dof_access>
inline TriaIterator<
  DoFCellAccessor<dimension_, space_dimension_, level_dof_access>>
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::parent() const
{
  TriaIterator<DoFCellAccessor<dimension_, space_dimension_, level_dof_access>>
    q(this->tria, this->level() - 1, this->parent_index(), this->dof_handler);

  return q;
}



namespace internal
{
  namespace DoFCellAccessorImplementation
  {
    template <int dim, int spacedim, bool level_dof_access>
    inline TriaIterator<
      dealii::DoFAccessor<dim - 1, dim, spacedim, level_dof_access>>
    get_face(
      const dealii::DoFCellAccessor<dim, spacedim, level_dof_access> &cell,
      const unsigned int                                              i,
      const std::integral_constant<int, 1>)
    {
      dealii::DoFAccessor<0, dim, spacedim, level_dof_access> a(
        &cell.get_triangulation(),
        ((i == 0) && cell.at_boundary(0) ?
           dealii::TriaAccessor<0, 1, spacedim>::left_vertex :
           ((i == 1) && cell.at_boundary(1) ?
              dealii::TriaAccessor<0, 1, spacedim>::right_vertex :
              dealii::TriaAccessor<0, 1, spacedim>::interior_vertex)),
        cell.vertex_index(i),
        &cell.get_dof_handler());
      return dealii::TriaIterator<
        dealii::DoFAccessor<0, dim, spacedim, level_dof_access>>(a);
    }


    template <int dim, int spacedim, bool level_dof_access>
    inline TriaIterator<
      dealii::DoFAccessor<dim - 1, dim, spacedim, level_dof_access>>
    get_face(
      const dealii::DoFCellAccessor<dim, spacedim, level_dof_access> &cell,
      const unsigned int                                              i,
      const std::integral_constant<int, 2>)
    {
      return cell.line(i);
    }


    template <int dim, int spacedim, bool level_dof_access>
    inline TriaIterator<
      dealii::DoFAccessor<dim - 1, dim, spacedim, level_dof_access>>
    get_face(
      const dealii::DoFCellAccessor<dim, spacedim, level_dof_access> &cell,
      const unsigned int                                              i,
      const std::integral_constant<int, 3>)
    {
      return cell.quad(i);
    }
  } // namespace DoFCellAccessorImplementation
} // namespace internal



template <int dimension_, int space_dimension_, bool level_dof_access>
inline typename DoFCellAccessor<dimension_,
                                space_dimension_,
                                level_dof_access>::face_iterator
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::face(
  const unsigned int i) const
{
  AssertIndexRange(i, this->n_faces());

  return dealii::internal::DoFCellAccessorImplementation::get_face(
    *this, i, std::integral_constant<int, dimension_>());
}



template <int dimension_, int space_dimension_, bool level_dof_access>
inline boost::container::small_vector<
  typename DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::
    face_iterator,
  GeometryInfo<dimension_>::faces_per_cell>
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::
  face_iterators() const
{
  boost::container::small_vector<
    typename DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::
      face_iterator,
    GeometryInfo<dimension_>::faces_per_cell>
    face_iterators(this->n_faces());

  for (const unsigned int i : this->face_indices())
    face_iterators[i] =
      dealii::internal::DoFCellAccessorImplementation::get_face(
        *this, i, std::integral_constant<int, dimension_>());

  return face_iterators;
}



template <int dimension_, int space_dimension_, bool level_dof_access>
inline void
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::
  get_active_or_mg_dof_indices(
    std::vector<types::global_dof_index> &dof_indices) const
{
  if (level_dof_access)
    get_mg_dof_indices(dof_indices);
  else
    get_dof_indices(dof_indices);
}



template <int dimension_, int space_dimension_, bool level_dof_access>
template <class InputVector, typename number>
inline void
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::get_dof_values(
  const InputVector &values,
  Vector<number>    &local_values) const
{
  get_dof_values(values, local_values.begin(), local_values.end());
}



template <int dimension_, int space_dimension_, bool level_dof_access>
template <typename Number, typename ForwardIterator>
inline void
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::get_dof_values(
  const ReadVector<Number> &values,
  ForwardIterator           local_values_begin,
  ForwardIterator           local_values_end) const
{
  (void)local_values_end;
  Assert(this->is_artificial() == false,
         ExcMessage("Can't ask for DoF indices on artificial cells."));
  Assert(this->is_active(), ExcMessage("Cell must be active."));
  Assert(this->dof_handler != nullptr, typename BaseClass::ExcInvalidObject());

  Assert(static_cast<unsigned int>(local_values_end - local_values_begin) ==
           this->get_fe().n_dofs_per_cell(),
         typename DoFCellAccessor::ExcVectorDoesNotMatch());
  Assert(values.size() == this->get_dof_handler().n_dofs(),
         typename DoFCellAccessor::ExcVectorDoesNotMatch());

  internal::DoFAccessorImplementation::Implementation::dof_index_vector_type
    dof_indices(this->get_fe().n_dofs_per_cell());
  internal::DoFAccessorImplementation::get_cell_dof_indices(
    *this, dof_indices, this->active_fe_index());

  boost::container::small_vector<Number, 27> values_temp(local_values_end -
                                                         local_values_begin);
  auto view = make_array_view(values_temp.begin(), values_temp.end());
  values.extract_subvector_to(make_array_view(dof_indices.begin(),
                                              dof_indices.end()),
                              view);
  using view_type = std::remove_reference_t<decltype(*local_values_begin)>;
  ArrayView<view_type> values_view2(&*local_values_begin,
                                    local_values_end - local_values_begin);
  std::copy(values_temp.begin(), values_temp.end(), values_view2.begin());
}



template <int dimension_, int space_dimension_, bool level_dof_access>
template <class InputVector, typename ForwardIterator>
inline void
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::get_dof_values(
  const AffineConstraints<typename InputVector::value_type> &constraints,
  const InputVector                                         &values,
  ForwardIterator                                            local_values_begin,
  ForwardIterator local_values_end) const
{
  Assert(this->is_artificial() == false,
         ExcMessage("Can't ask for DoF indices on artificial cells."));
  Assert(this->is_active(), ExcMessage("Cell must be active."));

  Assert(static_cast<unsigned int>(local_values_end - local_values_begin) ==
           this->get_fe().n_dofs_per_cell(),
         typename DoFCellAccessor::ExcVectorDoesNotMatch());
  Assert(values.size() == this->get_dof_handler().n_dofs(),
         typename DoFCellAccessor::ExcVectorDoesNotMatch());


  internal::DoFAccessorImplementation::Implementation::dof_index_vector_type
    dof_indices(this->get_fe().n_dofs_per_cell());
  internal::DoFAccessorImplementation::get_cell_dof_indices(
    *this, dof_indices, this->active_fe_index());

  constraints.get_dof_values(values,
                             dof_indices.data(),
                             local_values_begin,
                             local_values_end);
}



template <int dimension_, int space_dimension_, bool level_dof_access>
template <class OutputVector, typename number>
inline void
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::set_dof_values(
  const Vector<number> &local_values,
  OutputVector         &values) const
{
  Assert(this->is_artificial() == false,
         ExcMessage("Can't ask for DoF indices on artificial cells."));
  Assert(this->is_active(), ExcMessage("Cell must be active."));

  Assert(static_cast<unsigned int>(local_values.size()) ==
           this->get_fe().n_dofs_per_cell(),
         typename DoFCellAccessor::ExcVectorDoesNotMatch());
  Assert(values.size() == this->get_dof_handler().n_dofs(),
         typename DoFCellAccessor::ExcVectorDoesNotMatch());


  Assert(this->dof_handler != nullptr, typename BaseClass::ExcInvalidObject());
  internal::DoFAccessorImplementation::Implementation::dof_index_vector_type
    dof_indices(this->get_fe().n_dofs_per_cell());
  internal::DoFAccessorImplementation::get_cell_dof_indices(
    *this, dof_indices, this->active_fe_index());

  for (unsigned int i = 0; i < this->get_fe().n_dofs_per_cell(); ++i)
    internal::ElementAccess<OutputVector>::set(local_values(i),
                                               dof_indices[i],
                                               values);
}



template <int dimension_, int space_dimension_, bool level_dof_access>
inline const FiniteElement<dimension_, space_dimension_> &
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::get_fe() const
{
  Assert(this->dof_handler != nullptr, typename BaseClass::ExcInvalidObject());
  Assert((this->dof_handler->hp_capability_enabled == false) ||
           this->is_active(),
         ExcMessage(
           "For DoFHandler objects in hp-mode, finite elements are only "
           "associated with active cells. Consequently, you can not ask "
           "for the active finite element on cells with children."));

  const auto &fe = this->dof_handler->get_fe(active_fe_index());

  Assert(this->reference_cell() == fe.reference_cell(),
         internal::ExcNonMatchingReferenceCellTypes(this->reference_cell(),
                                                    fe.reference_cell()));

  return fe;
}



template <int dimension_, int space_dimension_, bool level_dof_access>
inline types::fe_index
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::
  active_fe_index() const
{
  Assert((this->dof_handler->hp_capability_enabled == false) ||
           this->is_active(),
         ExcMessage(
           "You can not ask for the active FE index on a cell that has "
           "children because no degrees of freedom are assigned "
           "to this cell and, consequently, no finite element "
           "is associated with it."));
  Assert((this->dof_handler->hp_capability_enabled == false) ||
           (this->is_locally_owned() || this->is_ghost()),
         ExcMessage("You can only query active FE index information on cells "
                    "that are either locally owned or (after distributing "
                    "degrees of freedom) are ghost cells."));

  return dealii::internal::DoFCellAccessorImplementation::Implementation::
    active_fe_index(*this);
}



template <int dimension_, int space_dimension_, bool level_dof_access>
inline void
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::
  set_active_fe_index(const types::fe_index i) const
{
  Assert((this->dof_handler->hp_capability_enabled == false) ||
           this->is_active(),
         ExcMessage("You can not set the active FE index on a cell that has "
                    "children because no degrees of freedom will be assigned "
                    "to this cell."));

  Assert((this->dof_handler->hp_capability_enabled == false) ||
           this->is_locally_owned(),
         ExcMessage("You can only set active FE index information on cells "
                    "that are locally owned. On ghost cells, this information "
                    "will automatically be propagated from the owning process "
                    "of that cell, and there is no information at all on "
                    "artificial cells."));

  dealii::internal::DoFCellAccessorImplementation::Implementation::
    set_active_fe_index(*this, i);
}



template <int dimension_, int space_dimension_, bool level_dof_access>
inline const FiniteElement<dimension_, space_dimension_> &
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::get_future_fe()
  const
{
  Assert(this->dof_handler != nullptr, typename BaseClass::ExcInvalidObject());
  Assert((this->dof_handler->hp_capability_enabled == false) ||
           this->is_active(),
         ExcMessage(
           "For DoFHandler objects in hp-mode, finite elements are only "
           "associated with active cells. Consequently, you can not ask "
           "for the future finite element on cells with children."));

  return this->dof_handler->get_fe(future_fe_index());
}



template <int dimension_, int space_dimension_, bool level_dof_access>
inline types::fe_index
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::
  future_fe_index() const
{
  Assert((this->dof_handler->hp_capability_enabled == false) ||
           (this->has_children() == false),
         ExcMessage(
           "You can not ask for the future FE index on a cell that has "
           "children because no degrees of freedom are assigned "
           "to this cell and, consequently, no finite element "
           "is associated with it."));
  Assert((this->dof_handler->hp_capability_enabled == false) ||
           (this->is_locally_owned()),
         ExcMessage("You can only query future FE index information on cells "
                    "that are locally owned."));

  return dealii::internal::DoFCellAccessorImplementation::Implementation::
    future_fe_index(*this);
}



template <int dimension_, int space_dimension_, bool level_dof_access>
inline void
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::
  set_future_fe_index(const types::fe_index i) const
{
  Assert((this->dof_handler->hp_capability_enabled == false) ||
           (this->has_children() == false),
         ExcMessage("You can not set the future FE index on a cell that has "
                    "children because no degrees of freedom will be assigned "
                    "to this cell."));

  Assert((this->dof_handler->hp_capability_enabled == false) ||
           this->is_locally_owned(),
         ExcMessage("You can only set future FE index information on cells "
                    "that are locally owned."));

  dealii::internal::DoFCellAccessorImplementation::Implementation::
    set_future_fe_index(*this, i);
}



template <int dimension_, int space_dimension_, bool level_dof_access>
inline bool
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::
  future_fe_index_set() const
{
  Assert((this->dof_handler->hp_capability_enabled == false) ||
           (this->has_children() == false),
         ExcMessage(
           "You can not ask for the future FE index on a cell that has "
           "children because no degrees of freedom are assigned "
           "to this cell and, consequently, no finite element "
           "is associated with it."));
  Assert((this->dof_handler->hp_capability_enabled == false) ||
           (this->is_locally_owned()),
         ExcMessage("You can only query future FE index information on cells "
                    "that are locally owned."));

  return dealii::internal::DoFCellAccessorImplementation::Implementation::
    future_fe_index_set(*this);
}



template <int dimension_, int space_dimension_, bool level_dof_access>
inline void
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::
  clear_future_fe_index() const
{
  Assert((this->dof_handler->hp_capability_enabled == false) ||
           (this->has_children() == false),
         ExcMessage(
           "You can not ask for the future FE index on a cell that has "
           "children because no degrees of freedom are assigned "
           "to this cell and, consequently, no finite element "
           "is associated with it."));
  Assert((this->dof_handler->hp_capability_enabled == false) ||
           (this->is_locally_owned()),
         ExcMessage("You can only query future FE index information on cells "
                    "that are locally owned."));

  dealii::internal::DoFCellAccessorImplementation::Implementation::
    clear_future_fe_index(*this);
}



template <int dimension_, int space_dimension_, bool level_dof_access>
template <typename number, typename OutputVector>
inline void
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::
  distribute_local_to_global(const Vector<number> &local_source,
                             OutputVector         &global_destination) const
{
  this->distribute_local_to_global(local_source.begin(),
                                   local_source.end(),
                                   global_destination);
}



template <int dimension_, int space_dimension_, bool level_dof_access>
template <typename ForwardIterator, typename OutputVector>
inline void
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::
  distribute_local_to_global(ForwardIterator local_source_begin,
                             ForwardIterator local_source_end,
                             OutputVector   &global_destination) const
{
  Assert(this->dof_handler != nullptr,
         (typename std::decay_t<decltype(*this)>::ExcInvalidObject()));
  Assert(static_cast<unsigned int>(local_source_end - local_source_begin) ==
           this->get_fe().n_dofs_per_cell(),
         (typename std::decay_t<decltype(*this)>::ExcVectorDoesNotMatch()));
  Assert(this->dof_handler->n_dofs() == global_destination.size(),
         (typename std::decay_t<decltype(*this)>::ExcVectorDoesNotMatch()));

  Assert(!this->has_children(), ExcMessage("Cell must be active"));

  Assert(
    internal::ArrayViewHelper::is_contiguous(local_source_begin,
                                             local_source_end),
    ExcMessage(
      "This function can not be called with iterator types that do not point to contiguous memory."));

  const unsigned int n_dofs = local_source_end - local_source_begin;

  internal::DoFAccessorImplementation::Implementation::dof_index_vector_type
    dof_indices(n_dofs);
  internal::DoFAccessorImplementation::get_cell_dof_indices(
    *this, dof_indices, this->active_fe_index());

  // distribute cell vector
  global_destination.add(n_dofs, dof_indices.data(), &(*local_source_begin));
}



template <int dimension_, int space_dimension_, bool level_dof_access>
template <typename ForwardIterator, typename OutputVector>
inline void
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::
  distribute_local_to_global(
    const AffineConstraints<typename OutputVector::value_type> &constraints,
    ForwardIterator local_source_begin,
    ForwardIterator local_source_end,
    OutputVector   &global_destination) const
{
  Assert(this->dof_handler != nullptr,
         (typename std::decay_t<decltype(*this)>::ExcInvalidObject()));
  Assert(local_source_end - local_source_begin ==
           this->get_fe().n_dofs_per_cell(),
         (typename std::decay_t<decltype(*this)>::ExcVectorDoesNotMatch()));
  Assert(this->dof_handler->n_dofs() == global_destination.size(),
         (typename std::decay_t<decltype(*this)>::ExcVectorDoesNotMatch()));

  Assert(!this->has_children(), ExcMessage("Cell must be active."));

  internal::DoFAccessorImplementation::Implementation::dof_index_vector_type
    dof_indices(this->get_fe().n_dofs_per_cell());
  internal::DoFAccessorImplementation::get_cell_dof_indices(
    *this, dof_indices, this->active_fe_index());

  // distribute cell vector
  constraints.distribute_local_to_global(local_source_begin,
                                         local_source_end,
                                         dof_indices.data(),
                                         global_destination);
}



template <int dimension_, int space_dimension_, bool level_dof_access>
template <typename number, typename OutputMatrix>
inline void
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::
  distribute_local_to_global(const FullMatrix<number> &local_source,
                             OutputMatrix             &global_destination) const
{
  Assert(this->dof_handler != nullptr,
         (typename std::decay_t<decltype(*this)>::ExcInvalidObject()));
  Assert(local_source.m() == this->get_fe().n_dofs_per_cell(),
         (typename std::decay_t<decltype(*this)>::ExcMatrixDoesNotMatch()));
  Assert(local_source.n() == this->get_fe().n_dofs_per_cell(),
         (typename std::decay_t<decltype(*this)>::ExcMatrixDoesNotMatch()));
  Assert(this->dof_handler->n_dofs() == global_destination.m(),
         (typename std::decay_t<decltype(*this)>::ExcMatrixDoesNotMatch()));
  Assert(this->dof_handler->n_dofs() == global_destination.n(),
         (typename std::decay_t<decltype(*this)>::ExcMatrixDoesNotMatch()));

  Assert(!this->has_children(), ExcMessage("Cell must be active."));

  const unsigned int n_dofs = local_source.m();

  internal::DoFAccessorImplementation::Implementation::dof_index_vector_type
    dof_indices(n_dofs);
  internal::DoFAccessorImplementation::get_cell_dof_indices(
    *this, dof_indices, this->active_fe_index());

  // distribute cell matrix
  for (unsigned int i = 0; i < n_dofs; ++i)
    global_destination.add(dof_indices[i],
                           n_dofs,
                           dof_indices.data(),
                           &local_source(i, 0));
}



template <int dimension_, int space_dimension_, bool level_dof_access>
template <typename number, typename OutputMatrix, typename OutputVector>
inline void
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::
  distribute_local_to_global(const FullMatrix<number> &local_matrix,
                             const Vector<number>     &local_vector,
                             OutputMatrix             &global_matrix,
                             OutputVector             &global_vector) const
{
  Assert(this->dof_handler != nullptr,
         (typename std::decay_t<decltype(*this)>::ExcInvalidObject()));
  Assert(local_matrix.m() == this->get_fe().n_dofs_per_cell(),
         (typename std::decay_t<decltype(*this)>::ExcMatrixDoesNotMatch()));
  Assert(local_matrix.n() == this->get_fe().n_dofs_per_cell(),
         (typename std::decay_t<decltype(*this)>::ExcVectorDoesNotMatch()));
  Assert(this->dof_handler->n_dofs() == global_matrix.m(),
         (typename std::decay_t<decltype(*this)>::ExcMatrixDoesNotMatch()));
  Assert(this->dof_handler->n_dofs() == global_matrix.n(),
         (typename std::decay_t<decltype(*this)>::ExcMatrixDoesNotMatch()));
  Assert(local_vector.size() == this->get_fe().n_dofs_per_cell(),
         (typename std::decay_t<decltype(*this)>::ExcVectorDoesNotMatch()));
  Assert(this->dof_handler->n_dofs() == global_vector.size(),
         (typename std::decay_t<decltype(*this)>::ExcVectorDoesNotMatch()));

  Assert(!this->has_children(), ExcMessage("Cell must be active."));

  const unsigned int n_dofs = this->get_fe().n_dofs_per_cell();
  internal::DoFAccessorImplementation::Implementation::dof_index_vector_type
    dof_indices(n_dofs);
  internal::DoFAccessorImplementation::get_cell_dof_indices(
    *this, dof_indices, this->active_fe_index());

  // distribute cell matrices
  for (unsigned int i = 0; i < n_dofs; ++i)
    {
      global_matrix.add(dof_indices[i],
                        n_dofs,
                        dof_indices.data(),
                        &local_matrix(i, 0));
      global_vector(dof_indices[i]) += local_vector(i);
    }
}



DEAL_II_NAMESPACE_CLOSE

#endif
