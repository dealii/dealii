// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2020 by the deal.II authors
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

#include <deal.II/hp/dof_faces.h>
#include <deal.II/hp/dof_level.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/read_write_vector.h>

#include <limits>
#include <type_traits>
#include <vector>

DEAL_II_NAMESPACE_OPEN


/*------------------------- Functions: DoFAccessor ---------------------------*/


template <int structdim, typename DoFHandlerType, bool level_dof_access>
inline DoFAccessor<structdim, DoFHandlerType, level_dof_access>::DoFAccessor()
{
  Assert(false, ExcInvalidObject());
}



template <int structdim, typename DoFHandlerType, bool level_dof_access>
inline DoFAccessor<structdim, DoFHandlerType, level_dof_access>::DoFAccessor(
  const Triangulation<DoFHandlerType::dimension,
                      DoFHandlerType::space_dimension> *tria,
  const int                                             level,
  const int                                             index,
  const DoFHandlerType *                                dof_handler)
  : dealii::internal::DoFAccessorImplementation::Inheritance<
      structdim,
      DoFHandlerType::dimension,
      DoFHandlerType::space_dimension>::BaseClass(tria, level, index)
  , dof_handler(const_cast<DoFHandlerType *>(dof_handler))
{
  Assert(
    tria == nullptr || &dof_handler->get_triangulation() == tria,
    ExcMessage(
      "You can't create a DoF accessor in which the DoFHandler object "
      "uses a different triangulation than the one you pass as argument."));
}



template <int structdim, typename DoFHandlerType, bool level_dof_access>
template <int structdim2, int dim2, int spacedim2>
DoFAccessor<structdim, DoFHandlerType, level_dof_access>::DoFAccessor(
  const InvalidAccessor<structdim2, dim2, spacedim2> &)
{
  Assert(false, ExcInvalidObject());
}



template <int structdim, typename DoFHandlerType, bool level_dof_access>
template <int dim2, class DoFHandlerType2, bool level_dof_access2>
inline DoFAccessor<structdim, DoFHandlerType, level_dof_access>::DoFAccessor(
  const DoFAccessor<dim2, DoFHandlerType2, level_dof_access2> &other)
  : BaseClass(other)
  , dof_handler(nullptr)
{
  Assert(false,
         ExcMessage(
           "You are trying to assign iterators that are incompatible. "
           "Reasons for incompatibility are that they point to different "
           "types of DoFHandlers (e.g., dealii::DoFHandler and "
           "dealii::hp::DoFHandler) or that they refer to objects of "
           "different dimensionality (e.g., assigning a line iterator "
           "to a quad iterator)."));
}



template <int structdim, typename DoFHandlerType, bool level_dof_access>
template <bool level_dof_access2>
inline DoFAccessor<structdim, DoFHandlerType, level_dof_access>::DoFAccessor(
  const DoFAccessor<structdim, DoFHandlerType, level_dof_access2> &other)
  : BaseClass(other)
  , dof_handler(const_cast<DoFHandlerType *>(other.dof_handler))
{}



template <int structdim, typename DoFHandlerType, bool level_dof_access>
inline void
DoFAccessor<structdim, DoFHandlerType, level_dof_access>::set_dof_handler(
  DoFHandlerType *dh)
{
  Assert(dh != nullptr, ExcInvalidObject());
  this->dof_handler = dh;
}



template <int structdim, typename DoFHandlerType, bool level_dof_access>
inline const DoFHandlerType &
DoFAccessor<structdim, DoFHandlerType, level_dof_access>::get_dof_handler()
  const
{
  Assert(this->dof_handler != nullptr, ExcInvalidObject());
  return *this->dof_handler;
}



template <int structdim, typename DoFHandlerType, bool level_dof_access>
inline void
DoFAccessor<structdim, DoFHandlerType, level_dof_access>::copy_from(
  const TriaAccessorBase<structdim,
                         DoFHandlerType::dimension,
                         DoFHandlerType::space_dimension> &da)
{
  Assert(this->dof_handler != nullptr, ExcInvalidObject());
  BaseClass::copy_from(da);
}



template <int structdim, typename DoFHandlerType, bool level_dof_access>
template <bool level_dof_access2>
inline void
DoFAccessor<structdim, DoFHandlerType, level_dof_access>::copy_from(
  const DoFAccessor<structdim, DoFHandlerType, level_dof_access2> &a)
{
  BaseClass::copy_from(a);
  this->dof_handler = a.dof_handler;
}



template <int structdim, typename DoFHandlerType, bool level_dof_access>
template <int dim2, class DoFHandlerType2, bool level_dof_access2>
inline bool
DoFAccessor<structdim, DoFHandlerType, level_dof_access>::
operator==(const DoFAccessor<dim2, DoFHandlerType2, level_dof_access2> &a) const
{
  Assert(structdim == dim2, ExcCantCompareIterators());
  Assert(this->dof_handler == a.dof_handler, ExcCantCompareIterators());
  return (BaseClass::operator==(a));
}



template <int structdim, typename DoFHandlerType, bool level_dof_access>
template <int dim2, class DoFHandlerType2, bool level_dof_access2>
inline bool
DoFAccessor<structdim, DoFHandlerType, level_dof_access>::
operator!=(const DoFAccessor<dim2, DoFHandlerType2, level_dof_access2> &a) const
{
  Assert(structdim == dim2, ExcCantCompareIterators());
  Assert(this->dof_handler == a.dof_handler, ExcCantCompareIterators());
  return (BaseClass::operator!=(a));
}



template <int structdim, typename DoFHandlerType, bool level_dof_access>
inline TriaIterator<DoFAccessor<structdim, DoFHandlerType, level_dof_access>>
DoFAccessor<structdim, DoFHandlerType, level_dof_access>::child(
  const unsigned int i) const
{
  Assert(static_cast<unsigned int>(this->level()) <
           this->dof_handler->levels.size(),
         ExcMessage("DoFHandler not initialized"));

  TriaIterator<TriaAccessor<structdim,
                            DoFHandlerType::dimension,
                            DoFHandlerType::space_dimension>>
    t = TriaAccessor<structdim,
                     DoFHandlerType::dimension,
                     DoFHandlerType::space_dimension>::child(i);

  TriaIterator<DoFAccessor<structdim, DoFHandlerType, level_dof_access>> q(
    *t, this->dof_handler);
  return q;
}


namespace internal
{
  namespace DoFAccessorImplementation
  {
    /**
     * A class like the one with same name in tria.cc. See there for more
     * information.
     */
    struct Implementation
    {
      /**
       * Implementations of the get_dof_index/set_dof_index functions.
       */
      template <int spacedim>
      static types::global_dof_index
      get_dof_index(const dealii::DoFHandler<1, spacedim> &dof_handler,
                    const unsigned int                     obj_level,
                    const unsigned int                     obj_index,
                    const unsigned int                     fe_index,
                    const unsigned int                     local_index,
                    std::integral_constant<int, 1>)
      {
        return dof_handler.levels[obj_level]->dof_object.get_dof_index(
          dof_handler, obj_index, fe_index, local_index);
      }


      template <int spacedim>
      static void
      set_dof_index(const dealii::DoFHandler<1, spacedim> &dof_handler,
                    const unsigned int                     obj_level,
                    const unsigned int                     obj_index,
                    const unsigned int                     fe_index,
                    const unsigned int                     local_index,
                    std::integral_constant<int, 1>,
                    const types::global_dof_index global_index)
      {
        dof_handler.levels[obj_level]->dof_object.set_dof_index(
          dof_handler, obj_index, fe_index, local_index, global_index);
      }


      template <int spacedim>
      static types::global_dof_index
      get_dof_index(const dealii::DoFHandler<2, spacedim> &dof_handler,
                    const unsigned int                     obj_level,
                    const unsigned int                     obj_index,
                    const unsigned int                     fe_index,
                    const unsigned int                     local_index,
                    std::integral_constant<int, 1>)
      {
        (void)obj_level;
        // faces have no levels
        Assert(obj_level == 0, ExcInternalError());
        return dof_handler.faces->lines.get_dof_index(dof_handler,
                                                      obj_index,
                                                      fe_index,
                                                      local_index);
      }


      template <int spacedim>
      static void
      set_dof_index(const dealii::DoFHandler<2, spacedim> &dof_handler,
                    const unsigned int                     obj_level,
                    const unsigned int                     obj_index,
                    const unsigned int                     fe_index,
                    const unsigned int                     local_index,
                    std::integral_constant<int, 1>,
                    const types::global_dof_index global_index)
      {
        (void)obj_level;
        // faces have no levels
        Assert(obj_level == 0, ExcInternalError());
        dof_handler.faces->lines.set_dof_index(
          dof_handler, obj_index, fe_index, local_index, global_index);
      }


      template <int spacedim>
      static types::global_dof_index
      get_dof_index(const dealii::DoFHandler<2, spacedim> &dof_handler,
                    const unsigned int                     obj_level,
                    const unsigned int                     obj_index,
                    const unsigned int                     fe_index,
                    const unsigned int                     local_index,
                    std::integral_constant<int, 2>)
      {
        return dof_handler.levels[obj_level]->dof_object.get_dof_index(
          dof_handler, obj_index, fe_index, local_index);
      }


      template <int spacedim>
      static void
      set_dof_index(const dealii::DoFHandler<2, spacedim> &dof_handler,
                    const unsigned int                     obj_level,
                    const unsigned int                     obj_index,
                    const unsigned int                     fe_index,
                    const unsigned int                     local_index,
                    std::integral_constant<int, 2>,
                    const types::global_dof_index global_index)
      {
        dof_handler.levels[obj_level]->dof_object.set_dof_index(
          dof_handler, obj_index, fe_index, local_index, global_index);
      }


      template <int spacedim>
      static types::global_dof_index
      get_dof_index(const dealii::DoFHandler<3, spacedim> &dof_handler,
                    const unsigned int                     obj_level,
                    const unsigned int                     obj_index,
                    const unsigned int                     fe_index,
                    const unsigned int                     local_index,
                    std::integral_constant<int, 1>)
      {
        (void)obj_level;
        // faces have no levels
        Assert(obj_level == 0, ExcInternalError());
        return dof_handler.faces->lines.get_dof_index(dof_handler,
                                                      obj_index,
                                                      fe_index,
                                                      local_index);
      }


      template <int spacedim>
      static void
      set_dof_index(const dealii::DoFHandler<3, spacedim> &dof_handler,
                    const unsigned int                     obj_level,
                    const unsigned int                     obj_index,
                    const unsigned int                     fe_index,
                    const unsigned int                     local_index,
                    std::integral_constant<int, 1>,
                    const types::global_dof_index global_index)
      {
        (void)obj_level;
        // faces have no levels
        Assert(obj_level == 0, ExcInternalError());
        dof_handler.faces->lines.set_dof_index(
          dof_handler, obj_index, fe_index, local_index, global_index);
      }



      template <int spacedim>
      static types::global_dof_index
      get_dof_index(const dealii::DoFHandler<3, spacedim> &dof_handler,
                    const unsigned int                     obj_level,
                    const unsigned int                     obj_index,
                    const unsigned int                     fe_index,
                    const unsigned int                     local_index,
                    std::integral_constant<int, 2>)
      {
        (void)obj_level;
        // faces have no levels
        Assert(obj_level == 0, ExcInternalError());
        return dof_handler.faces->quads.get_dof_index(dof_handler,
                                                      obj_index,
                                                      fe_index,
                                                      local_index);
      }


      template <int spacedim>
      static void
      set_dof_index(const dealii::DoFHandler<3, spacedim> &dof_handler,
                    const unsigned int                     obj_level,
                    const unsigned int                     obj_index,
                    const unsigned int                     fe_index,
                    const unsigned int                     local_index,
                    std::integral_constant<int, 2>,
                    const types::global_dof_index global_index)
      {
        (void)obj_level;
        // faces have no levels
        Assert(obj_level == 0, ExcInternalError());
        dof_handler.faces->quads.set_dof_index(
          dof_handler, obj_index, fe_index, local_index, global_index);
      }



      template <int spacedim>
      static types::global_dof_index
      get_dof_index(const dealii::DoFHandler<3, spacedim> &dof_handler,
                    const unsigned int                     obj_level,
                    const unsigned int                     obj_index,
                    const unsigned int                     fe_index,
                    const unsigned int                     local_index,
                    std::integral_constant<int, 3>)
      {
        return dof_handler.levels[obj_level]->dof_object.get_dof_index(
          dof_handler, obj_index, fe_index, local_index);
      }


      template <int spacedim>
      static void
      set_dof_index(const dealii::DoFHandler<3, spacedim> &dof_handler,
                    const unsigned int                     obj_level,
                    const unsigned int                     obj_index,
                    const unsigned int                     fe_index,
                    const unsigned int                     local_index,
                    std::integral_constant<int, 3>,
                    const types::global_dof_index global_index)
      {
        dof_handler.levels[obj_level]->dof_object.set_dof_index(
          dof_handler, obj_index, fe_index, local_index, global_index);
      }


      template <int spacedim>
      static types::global_dof_index
      get_dof_index(const dealii::hp::DoFHandler<1, spacedim> &dof_handler,
                    const unsigned int                         obj_level,
                    const unsigned int                         obj_index,
                    const unsigned int                         fe_index,
                    const unsigned int                         local_index,
                    const std::integral_constant<int, 1> &)
      {
        return dof_handler.levels[obj_level]->get_dof_index(obj_index,
                                                            fe_index,
                                                            local_index);
      }


      template <int spacedim>
      static void
      set_dof_index(const dealii::hp::DoFHandler<1, spacedim> &dof_handler,
                    const unsigned int                         obj_level,
                    const unsigned int                         obj_index,
                    const unsigned int                         fe_index,
                    const unsigned int                         local_index,
                    const std::integral_constant<int, 1> &,
                    const types::global_dof_index global_index)
      {
        dof_handler.levels[obj_level]->set_dof_index(obj_index,
                                                     fe_index,
                                                     local_index,
                                                     global_index);
      }


      template <int spacedim>
      static types::global_dof_index
      get_dof_index(const dealii::hp::DoFHandler<2, spacedim> &dof_handler,
                    const unsigned int                         obj_level,
                    const unsigned int                         obj_index,
                    const unsigned int                         fe_index,
                    const unsigned int                         local_index,
                    const std::integral_constant<int, 1> &)
      {
        return dof_handler.faces->lines.get_dof_index(
          dof_handler, obj_index, fe_index, local_index, obj_level);
      }


      template <int spacedim>
      static void
      set_dof_index(const dealii::hp::DoFHandler<2, spacedim> &dof_handler,
                    const unsigned int                         obj_level,
                    const unsigned int                         obj_index,
                    const unsigned int                         fe_index,
                    const unsigned int                         local_index,
                    const std::integral_constant<int, 1> &,
                    const types::global_dof_index global_index)
      {
        dof_handler.faces->lines.set_dof_index(dof_handler,
                                               obj_index,
                                               fe_index,
                                               local_index,
                                               global_index,
                                               obj_level);
      }


      template <int spacedim>
      static types::global_dof_index
      get_dof_index(const dealii::hp::DoFHandler<2, spacedim> &dof_handler,
                    const unsigned int                         obj_level,
                    const unsigned int                         obj_index,
                    const unsigned int                         fe_index,
                    const unsigned int                         local_index,
                    const std::integral_constant<int, 2> &)
      {
        return dof_handler.levels[obj_level]->get_dof_index(obj_index,
                                                            fe_index,
                                                            local_index);
      }


      template <int spacedim>
      static void
      set_dof_index(const dealii::hp::DoFHandler<2, spacedim> &dof_handler,
                    const unsigned int                         obj_level,
                    const unsigned int                         obj_index,
                    const unsigned int                         fe_index,
                    const unsigned int                         local_index,
                    const std::integral_constant<int, 2> &,
                    const types::global_dof_index global_index)
      {
        dof_handler.levels[obj_level]->set_dof_index(obj_index,
                                                     fe_index,
                                                     local_index,
                                                     global_index);
      }


      template <int spacedim>
      static types::global_dof_index
      get_dof_index(const dealii::hp::DoFHandler<3, spacedim> &dof_handler,
                    const unsigned int                         obj_level,
                    const unsigned int                         obj_index,
                    const unsigned int                         fe_index,
                    const unsigned int                         local_index,
                    const std::integral_constant<int, 1> &)
      {
        return dof_handler.faces->lines.get_dof_index(
          dof_handler, obj_index, fe_index, local_index, obj_level);
      }


      template <int spacedim>
      static void
      set_dof_index(const dealii::hp::DoFHandler<3, spacedim> &dof_handler,
                    const unsigned int                         obj_level,
                    const unsigned int                         obj_index,
                    const unsigned int                         fe_index,
                    const unsigned int                         local_index,
                    const std::integral_constant<int, 1> &,
                    const types::global_dof_index global_index)
      {
        dof_handler.faces->lines.set_dof_index(dof_handler,
                                               obj_index,
                                               fe_index,
                                               local_index,
                                               global_index,
                                               obj_level);
      }


      template <int spacedim>
      static types::global_dof_index
      get_dof_index(const dealii::hp::DoFHandler<3, spacedim> &dof_handler,
                    const unsigned int                         obj_level,
                    const unsigned int                         obj_index,
                    const unsigned int                         fe_index,
                    const unsigned int                         local_index,
                    const std::integral_constant<int, 2> &)
      {
        return dof_handler.faces->quads.get_dof_index(
          dof_handler, obj_index, fe_index, local_index, obj_level);
      }


      template <int spacedim>
      static void
      set_dof_index(const dealii::hp::DoFHandler<3, spacedim> &dof_handler,
                    const unsigned int                         obj_level,
                    const unsigned int                         obj_index,
                    const unsigned int                         fe_index,
                    const unsigned int                         local_index,
                    const std::integral_constant<int, 2> &,
                    const types::global_dof_index global_index)
      {
        dof_handler.faces->quads.set_dof_index(dof_handler,
                                               obj_index,
                                               fe_index,
                                               local_index,
                                               global_index,
                                               obj_level);
      }


      template <int spacedim>
      static types::global_dof_index
      get_dof_index(const dealii::hp::DoFHandler<3, spacedim> &dof_handler,
                    const unsigned int                         obj_level,
                    const unsigned int                         obj_index,
                    const unsigned int                         fe_index,
                    const unsigned int                         local_index,
                    const std::integral_constant<int, 3> &)
      {
        return dof_handler.levels[obj_level]->get_dof_index(obj_index,
                                                            fe_index,
                                                            local_index);
      }


      template <int spacedim>
      static void
      set_dof_index(const dealii::hp::DoFHandler<3, spacedim> &dof_handler,
                    const unsigned int                         obj_level,
                    const unsigned int                         obj_index,
                    const unsigned int                         fe_index,
                    const unsigned int                         local_index,
                    const std::integral_constant<int, 3> &,
                    const types::global_dof_index global_index)
      {
        dof_handler.levels[obj_level]->set_dof_index(obj_index,
                                                     fe_index,
                                                     local_index,
                                                     global_index);
      }


      template <int dim, int spacedim>
      static types::global_dof_index
      mg_vertex_dof_index(const dealii::DoFHandler<dim, spacedim> &dof_handler,
                          const int                                level,
                          const unsigned int                       vertex_index,
                          const unsigned int                       i)
      {
        return dof_handler.mg_vertex_dofs[vertex_index].get_index(
          level, i, dof_handler.get_fe().dofs_per_vertex);
      }


      template <int dim, int spacedim>
      static types::global_dof_index
      mg_vertex_dof_index(const dealii::hp::DoFHandler<dim, spacedim> &,
                          const int,
                          const unsigned int,
                          const unsigned int)
      {
        Assert(false,
               ExcMessage(
                 "hp::DoFHandler does not implement multilevel DoFs."));
        return numbers::invalid_dof_index;
      }


      template <int dim, int spacedim>
      static void
      set_mg_vertex_dof_index(dealii::DoFHandler<dim, spacedim> &dof_handler,
                              const int                          level,
                              const unsigned int                 vertex_index,
                              const unsigned int                 i,
                              types::global_dof_index            index)
      {
        return dof_handler.mg_vertex_dofs[vertex_index].set_index(
          level, i, dof_handler.get_fe().dofs_per_vertex, index);
      }


      template <int dim, int spacedim>
      static void
      set_mg_vertex_dof_index(dealii::hp::DoFHandler<dim, spacedim> &,
                              const int,
                              const unsigned int,
                              const unsigned int,
                              types::global_dof_index)
      {
        Assert(false,
               ExcMessage(
                 "hp::DoFHandler does not implement multilevel DoFs."));
      }



      template <int structdim, int dim, int spacedim>
      static bool
      fe_index_is_active(const dealii::DoFHandler<dim, spacedim> &,
                         const unsigned int,
                         const unsigned int,
                         const unsigned int fe_index,
                         const std::integral_constant<int, structdim> &)
      {
        return (fe_index == 0);
      }



      template <int structdim, int dim, int spacedim>
      static unsigned int
      n_active_fe_indices(const dealii::DoFHandler<dim, spacedim> &dof_handler,
                          const unsigned int                       obj_level,
                          const unsigned int                       obj_index,
                          const std::integral_constant<int, structdim> &)
      {
        (void)dof_handler;
        (void)obj_level;
        (void)obj_index;
        // check that the object we look
        // at is in fact active. the
        // problem is that we have
        // templatized on the
        // dimensionality of the object,
        // so it may be a cell, a face,
        // or a line. we have a bit of
        // trouble doing this all in the
        // generic case, so only check if
        // it is either a cell or a
        // line. the only case this
        // leaves out is faces in 3d --
        // let's hope that this never is
        // a problem
        Assert((dim == structdim ?
                  TriaRawIterator<dealii::CellAccessor<dim, spacedim>>(
                    &dof_handler.get_triangulation(), obj_level, obj_index)
                    ->used() :
                  (structdim == 1 ?
                     typename internal::TriangulationImplementation::
                       Iterators<dim, spacedim>::raw_line_iterator(
                         &dof_handler.get_triangulation(), obj_level, obj_index)
                         ->used() :
                     true)) == true,
               ExcMessage("This cell is not active and therefore can't be "
                          "queried for its active FE indices"));
        return 1;
      }



      template <int structdim, int dim, int spacedim>
      static unsigned int
      nth_active_fe_index(const dealii::DoFHandler<dim, spacedim> &dof_handler,
                          const unsigned int                       obj_level,
                          const unsigned int                       obj_index,
                          const unsigned int                       n,
                          const std::integral_constant<int, structdim> &)
      {
        (void)dof_handler;
        (void)obj_level;
        (void)obj_index;
        (void)n;
        // check that the object we look
        // at is in fact active. the
        // problem is that we have
        // templatized on the
        // dimensionality of the object,
        // so it may be a cell, a face,
        // or a line. we have a bit of
        // trouble doing this all in the
        // generic case, so only check if
        // it is either a cell or a
        // line. the only case this
        // leaves out is faces in 3d --
        // let's hope that this never is
        // a problem
        Assert((dim == structdim ?
                  TriaRawIterator<dealii::CellAccessor<dim, spacedim>>(
                    &dof_handler.get_triangulation(), obj_level, obj_index)
                    ->used() :
                  (structdim == 1 ?
                     typename internal::TriangulationImplementation::
                       Iterators<dim, spacedim>::raw_line_iterator(
                         &dof_handler.get_triangulation(), obj_level, obj_index)
                         ->used() :
                     true)) == true,
               ExcMessage("This cell is not active and therefore can't be "
                          "queried for its active FE indices"));
        AssertIndexRange(n, 1);

        return dealii::DoFHandler<dim, spacedim>::default_fe_index;
      }


      template <int spacedim>
      static bool
      fe_index_is_active(const dealii::hp::DoFHandler<1, spacedim> &dof_handler,
                         const unsigned int                         obj_level,
                         const unsigned int                         obj_index,
                         const unsigned int                         fe_index,
                         const std::integral_constant<int, 1> &)
      {
        return dof_handler.levels[obj_level]->fe_index_is_active(obj_index,
                                                                 fe_index);
      }


      template <int spacedim>
      static unsigned int
      n_active_fe_indices(const dealii::hp::DoFHandler<1, spacedim> &,
                          const unsigned int /*obj_level*/,
                          const unsigned int /*obj_index*/,
                          const std::integral_constant<int, 1> &)
      {
        // on a cell, the number of active elements is one
        return 1;
      }



      template <int spacedim>
      static unsigned int
      nth_active_fe_index(
        const dealii::hp::DoFHandler<1, spacedim> &dof_handler,
        const unsigned int                         obj_level,
        const unsigned int                         obj_index,
        const unsigned int                         n,
        const std::integral_constant<int, 1> &)
      {
        (void)n;
        Assert(n == 0,
               ExcMessage("On cells, there can only be one active FE index"));
        return dof_handler.levels[obj_level]->active_fe_index(obj_index);
      }


      template <int spacedim>
      static bool
      fe_index_is_active(const dealii::hp::DoFHandler<2, spacedim> &dof_handler,
                         const unsigned int                         obj_level,
                         const unsigned int                         obj_index,
                         const unsigned int                         fe_index,
                         const std::integral_constant<int, 1> &)
      {
        return dof_handler.faces->lines.fe_index_is_active(dof_handler,
                                                           obj_index,
                                                           fe_index,
                                                           obj_level);
      }


      template <int spacedim>
      static unsigned int
      n_active_fe_indices(
        const dealii::hp::DoFHandler<2, spacedim> &dof_handler,
        const unsigned int,
        const unsigned int obj_index,
        const std::integral_constant<int, 1> &)
      {
        return dof_handler.faces->lines.n_active_fe_indices(dof_handler,
                                                            obj_index);
      }


      template <int spacedim>
      static unsigned int
      nth_active_fe_index(
        const dealii::hp::DoFHandler<2, spacedim> &dof_handler,
        const unsigned int                         obj_level,
        const unsigned int                         obj_index,
        const unsigned int                         n,
        const std::integral_constant<int, 1> &)
      {
        return dof_handler.faces->lines.nth_active_fe_index(dof_handler,
                                                            obj_level,
                                                            obj_index,
                                                            n);
      }



      template <int spacedim>
      static bool
      fe_index_is_active(const dealii::hp::DoFHandler<2, spacedim> &dof_handler,
                         const unsigned int                         obj_level,
                         const unsigned int                         obj_index,
                         const unsigned int                         fe_index,
                         const std::integral_constant<int, 2> &)
      {
        return dof_handler.levels[obj_level]->fe_index_is_active(obj_index,
                                                                 fe_index);
      }


      template <int spacedim>
      static unsigned int
      n_active_fe_indices(const dealii::hp::DoFHandler<2, spacedim> &,
                          const unsigned int /*obj_level*/,
                          const unsigned int /*obj_index*/,
                          const std::integral_constant<int, 2> &)
      {
        // on a cell, the number of active elements is one
        return 1;
      }



      template <int spacedim>
      static unsigned int
      nth_active_fe_index(
        const dealii::hp::DoFHandler<2, spacedim> &dof_handler,
        const unsigned int                         obj_level,
        const unsigned int                         obj_index,
        const unsigned int                         n,
        const std::integral_constant<int, 2> &)
      {
        (void)n;
        Assert(n == 0,
               ExcMessage("On cells, there can only be one active FE index"));
        return dof_handler.levels[obj_level]->active_fe_index(obj_index);
      }



      template <int spacedim>
      static bool
      fe_index_is_active(const dealii::hp::DoFHandler<3, spacedim> &dof_handler,
                         const unsigned int                         obj_level,
                         const unsigned int                         obj_index,
                         const unsigned int                         fe_index,
                         const std::integral_constant<int, 1> &)
      {
        return dof_handler.faces->lines.fe_index_is_active(dof_handler,
                                                           obj_index,
                                                           fe_index,
                                                           obj_level);
      }


      template <int spacedim>
      static unsigned int
      n_active_fe_indices(
        const dealii::hp::DoFHandler<3, spacedim> &dof_handler,
        const unsigned int,
        const unsigned int obj_index,
        const std::integral_constant<int, 1> &)
      {
        return dof_handler.faces->lines.n_active_fe_indices(dof_handler,
                                                            obj_index);
      }



      template <int spacedim>
      static unsigned int
      nth_active_fe_index(
        const dealii::hp::DoFHandler<3, spacedim> &dof_handler,
        const unsigned int                         obj_level,
        const unsigned int                         obj_index,
        const unsigned int                         n,
        const std::integral_constant<int, 1> &)
      {
        return dof_handler.faces->lines.nth_active_fe_index(dof_handler,
                                                            obj_level,
                                                            obj_index,
                                                            n);
      }



      template <int spacedim>
      static bool
      fe_index_is_active(const dealii::hp::DoFHandler<3, spacedim> &dof_handler,
                         const unsigned int                         obj_level,
                         const unsigned int                         obj_index,
                         const unsigned int                         fe_index,
                         const std::integral_constant<int, 2> &)
      {
        return dof_handler.faces->quads.fe_index_is_active(dof_handler,
                                                           obj_index,
                                                           fe_index,
                                                           obj_level);
      }

      template <int spacedim>
      static bool
      fe_index_is_active(const dealii::hp::DoFHandler<3, spacedim> &dof_handler,
                         const unsigned int                         obj_level,
                         const unsigned int                         obj_index,
                         const unsigned int                         fe_index,
                         const std::integral_constant<int, 3> &)
      {
        return dof_handler.levels[obj_level]->fe_index_is_active(obj_index,
                                                                 fe_index);
      }


      template <int spacedim>
      static unsigned int
      n_active_fe_indices(
        const dealii::hp::DoFHandler<3, spacedim> &dof_handler,
        const unsigned int,
        const unsigned int obj_index,
        const std::integral_constant<int, 2> &)
      {
        return dof_handler.faces->quads.n_active_fe_indices(dof_handler,
                                                            obj_index);
      }



      template <int spacedim>
      static unsigned int
      nth_active_fe_index(
        const dealii::hp::DoFHandler<3, spacedim> &dof_handler,
        const unsigned int                         obj_level,
        const unsigned int                         obj_index,
        const unsigned int                         n,
        const std::integral_constant<int, 2> &)
      {
        return dof_handler.faces->quads.nth_active_fe_index(dof_handler,
                                                            obj_level,
                                                            obj_index,
                                                            n);
      }



      template <int spacedim>
      static unsigned int
      n_active_fe_indices(const dealii::hp::DoFHandler<3, spacedim> &,
                          const unsigned int /*obj_level*/,
                          const unsigned int /*obj_index*/,
                          const std::integral_constant<int, 3> &)
      {
        // on a cell, the number of active elements is one
        return 1;
      }



      template <int spacedim>
      static unsigned int
      nth_active_fe_index(
        const dealii::hp::DoFHandler<3, spacedim> &dof_handler,
        const unsigned int                         obj_level,
        const unsigned int                         obj_index,
        const unsigned int                         n,
        const std::integral_constant<int, 3> &)
      {
        (void)n;
        Assert(n == 0,
               ExcMessage("On cells, there can only be one active FE index"));
        return dof_handler.levels[obj_level]->active_fe_index(obj_index);
      }

      /**
       * Set the @p local_index-th degree of freedom corresponding to the
       * finite element specified by @p fe_index on the vertex with global
       * number @p vertex_index to @p global_index.
       */
      template <int dim, int spacedim>
      static void
      set_vertex_dof_index(dealii::DoFHandler<dim, spacedim> &dof_handler,
                           const unsigned int                 vertex_index,
                           const unsigned int                 fe_index,
                           const unsigned int                 local_index,
                           const types::global_dof_index      global_index)
      {
        (void)fe_index;
        Assert(
          (fe_index == dealii::DoFHandler<dim, spacedim>::default_fe_index),
          ExcMessage(
            "Only the default FE index is allowed for non-hp DoFHandler objects"));
        AssertIndexRange(local_index, dof_handler.get_fe().dofs_per_vertex);

        dof_handler
          .vertex_dofs[vertex_index * dof_handler.get_fe().dofs_per_vertex +
                       local_index] = global_index;
      }


      template <int dim, int spacedim>
      static void
      set_vertex_dof_index(dealii::hp::DoFHandler<dim, spacedim> &dof_handler,
                           const unsigned int                     vertex_index,
                           const unsigned int                     fe_index,
                           const unsigned int                     local_index,
                           const types::global_dof_index          global_index)
      {
        Assert((fe_index !=
                dealii::hp::DoFHandler<dim, spacedim>::default_fe_index),
               ExcMessage("You need to specify a FE index when working "
                          "with hp DoFHandlers"));
        Assert(dof_handler.fe_collection.size() > 0,
               ExcMessage("No finite element collection is associated with "
                          "this DoFHandler"));
        AssertIndexRange(local_index,
                         dof_handler.get_fe(fe_index).dofs_per_vertex);
        Assert(fe_index < dof_handler.fe_collection.size(), ExcInternalError());
        Assert(dof_handler.vertex_dof_offsets[vertex_index] !=
                 numbers::invalid_unsigned_int,
               ExcMessage(
                 "This vertex is unused and has no DoFs associated with it"));

        // hop along the list of index
        // sets until we find the one
        // with the correct fe_index, and
        // then poke into that
        // part. trigger an exception if
        // we can't find a set for this
        // particular fe_index
        const unsigned int starting_offset =
          dof_handler.vertex_dof_offsets[vertex_index];
        types::global_dof_index *pointer =
          &dof_handler.vertex_dofs[starting_offset];
        while (true)
          {
            Assert(pointer <= &dof_handler.vertex_dofs.back(),
                   ExcInternalError());

            // a fe index is always small
            Assert((*pointer) < std::numeric_limits<unsigned int>::max(),
                   ExcInternalError());
            const types::global_dof_index this_fe_index = *pointer;

            Assert(this_fe_index != numbers::invalid_dof_index,
                   ExcInternalError());
            Assert(this_fe_index < dof_handler.fe_collection.size(),
                   ExcInternalError());

            if (this_fe_index == fe_index)
              {
                *(pointer + 1 + local_index) = global_index;
                return;
              }
            else
              pointer += static_cast<types::global_dof_index>(
                dof_handler.get_fe(this_fe_index).dofs_per_vertex + 1);
          }
      }


      /**
       * Get the @p local_index-th degree of freedom corresponding to the
       * finite element specified by @p fe_index on the vertex with global
       * number @p vertex_index to @p global_index.
       */

      template <int dim, int spacedim>
      static types::global_dof_index
      get_vertex_dof_index(const dealii::DoFHandler<dim, spacedim> &dof_handler,
                           const unsigned int vertex_index,
                           const unsigned int fe_index,
                           const unsigned int local_index)
      {
        (void)fe_index;
        Assert(
          (fe_index == dealii::DoFHandler<dim, spacedim>::default_fe_index),
          ExcMessage(
            "Only the default FE index is allowed for non-hp DoFHandler objects"));
        AssertIndexRange(local_index, dof_handler.get_fe().dofs_per_vertex);

        return dof_handler
          .vertex_dofs[vertex_index * dof_handler.get_fe().dofs_per_vertex +
                       local_index];
      }


      template <int dim, int spacedim>
      static types::global_dof_index
      get_vertex_dof_index(
        const dealii::hp::DoFHandler<dim, spacedim> &dof_handler,
        const unsigned int                           vertex_index,
        const unsigned int                           fe_index,
        const unsigned int                           local_index)
      {
        Assert((fe_index !=
                dealii::hp::DoFHandler<dim, spacedim>::default_fe_index),
               ExcMessage("You need to specify a FE index when working "
                          "with hp DoFHandlers"));
        Assert(dof_handler.fe_collection.size() > 0,
               ExcMessage("No finite element collection is associated with "
                          "this DoFHandler"));
        AssertIndexRange(local_index,
                         dof_handler.get_fe(fe_index).dofs_per_vertex);
        AssertIndexRange(vertex_index, dof_handler.vertex_dof_offsets.size());
        Assert(dof_handler.vertex_dof_offsets[vertex_index] !=
                 numbers::invalid_unsigned_int,
               ExcMessage(
                 "This vertex is unused and has no DoFs associated with it"));

        // hop along the list of index
        // sets until we find the one
        // with the correct fe_index, and
        // then poke into that
        // part. trigger an exception if
        // we can't find a set for this
        // particular fe_index
        const unsigned int starting_offset =
          dof_handler.vertex_dof_offsets[vertex_index];
        const types::global_dof_index *pointer =
          &dof_handler.vertex_dofs[starting_offset];
        while (true)
          {
            Assert(pointer <= &dof_handler.vertex_dofs.back(),
                   ExcInternalError());

            Assert((*pointer) <
                     std::numeric_limits<types::global_dof_index>::max(),
                   ExcInternalError());
            const types::global_dof_index this_fe_index = *pointer;

            Assert(this_fe_index != numbers::invalid_dof_index,
                   ExcInternalError());
            Assert(this_fe_index < dof_handler.fe_collection.size(),
                   ExcInternalError());

            if (this_fe_index == fe_index)
              return *(pointer + 1 + local_index);
            else
              pointer += static_cast<types::global_dof_index>(
                dof_handler.get_fe(this_fe_index).dofs_per_vertex + 1);
          }
      }


      /**
       * Return the number of different finite elements that are active on a
       * given vertex.
       */
      template <int dim, int spacedim>
      static unsigned int
      n_active_vertex_fe_indices(
        const dealii::hp::DoFHandler<dim, spacedim> &dof_handler,
        const unsigned int                           vertex_index)
      {
        Assert(dof_handler.fe_collection.size() > 0,
               ExcMessage("No finite element collection is associated with "
                          "this DoFHandler"));

        // if this vertex is unused, return 0
        if (dof_handler.vertex_dof_offsets[vertex_index] ==
            numbers::invalid_unsigned_int)
          return 0;

        // hop along the list of index
        // sets and count the number of
        // hops
        const unsigned int starting_offset =
          dof_handler.vertex_dof_offsets[vertex_index];
        const types::global_dof_index *pointer =
          &dof_handler.vertex_dofs[starting_offset];

        Assert(*pointer != numbers::invalid_dof_index, ExcInternalError());

        unsigned int counter = 0;
        while (true)
          {
            Assert(pointer <= &dof_handler.vertex_dofs.back(),
                   ExcInternalError());

            const types::global_dof_index this_fe_index = *pointer;

            if (this_fe_index == numbers::invalid_dof_index)
              return counter;
            else
              {
                pointer += static_cast<types::global_dof_index>(
                  dof_handler.get_fe(this_fe_index).dofs_per_vertex + 1);
                ++counter;
              }
          }
      }



      /**
       * Return the fe index of the n-th finite element active on a given
       * vertex.
       */
      template <int dim, int spacedim>
      static unsigned int
      nth_active_vertex_fe_index(
        const dealii::hp::DoFHandler<dim, spacedim> &dof_handler,
        const unsigned int                           vertex_index,
        const unsigned int                           n)
      {
        Assert(dof_handler.fe_collection.size() > 0,
               ExcMessage("No finite element collection is associated with "
                          "this DoFHandler"));
        Assert(n < n_active_vertex_fe_indices(dof_handler, vertex_index),
               ExcIndexRange(
                 n, 0, n_active_vertex_fe_indices(dof_handler, vertex_index)));
        // make sure we don't ask on
        // unused vertices
        Assert(dof_handler.vertex_dof_offsets[vertex_index] !=
                 numbers::invalid_unsigned_int,
               ExcInternalError());

        // hop along the list of index
        // sets and count the number of
        // hops
        const unsigned int starting_offset =
          dof_handler.vertex_dof_offsets[vertex_index];
        const types::global_dof_index *pointer =
          &dof_handler.vertex_dofs[starting_offset];

        Assert(*pointer != numbers::invalid_dof_index, ExcInternalError());

        unsigned int counter = 0;
        while (true)
          {
            Assert(pointer <= &dof_handler.vertex_dofs.back(),
                   ExcInternalError());

            Assert((*pointer) < std::numeric_limits<unsigned int>::max(),
                   ExcInternalError());
            const types::global_dof_index this_fe_index = *pointer;

            Assert(this_fe_index < dof_handler.fe_collection.size(),
                   ExcInternalError());

            if (counter == n)
              return this_fe_index;

            Assert(this_fe_index != numbers::invalid_dof_index,
                   ExcInternalError());

            pointer += static_cast<types::global_dof_index>(
              dof_handler.get_fe(this_fe_index).dofs_per_vertex + 1);
            ++counter;
          }
      }



      /**
       * Returns all active fe indices on a given vertex.
       *
       * The size of the returned set equals the number of finite elements that
       * are active on this vertex.
       */
      template <int dim, int spacedim>
      static std::set<unsigned int>
      get_active_vertex_fe_indices(
        const dealii::hp::DoFHandler<dim, spacedim> &dof_handler,
        const unsigned int                           vertex_index)
      {
        std::set<unsigned int> active_fe_indices;
        for (unsigned int i = 0;
             i < n_active_vertex_fe_indices(dof_handler, vertex_index);
             ++i)
          active_fe_indices.insert(
            nth_active_vertex_fe_index(dof_handler, vertex_index, i));
        return active_fe_indices;
      }



      /**
       * Return whether a particular finite element index is active on the
       * specified vertex.
       */
      template <int dim, int spacedim>
      static bool
      fe_is_active_on_vertex(
        const dealii::hp::DoFHandler<dim, spacedim> &dof_handler,
        const unsigned int                           vertex_index,
        const unsigned int                           fe_index)
      {
        Assert((fe_index !=
                dealii::hp::DoFHandler<dim, spacedim>::default_fe_index),
               ExcMessage("You need to specify a FE index when working "
                          "with hp DoFHandlers"));
        Assert(dof_handler.fe_collection.size() > 0,
               ExcMessage("No finite element collection is associated with "
                          "this DoFHandler"));
        Assert(fe_index < dof_handler.fe_collection.size(), ExcInternalError());

        // make sure we don't ask on
        // unused vertices
        Assert(dof_handler.vertex_dof_offsets[vertex_index] !=
                 numbers::invalid_unsigned_int,
               ExcInternalError());

        // hop along the list of index
        // sets and see whether we find
        // the given index
        const unsigned int starting_offset =
          dof_handler.vertex_dof_offsets[vertex_index];
        const types::global_dof_index *pointer =
          &dof_handler.vertex_dofs[starting_offset];

        Assert(*pointer != numbers::invalid_dof_index, ExcInternalError());

        while (true)
          {
            Assert(pointer <= &dof_handler.vertex_dofs.back(),
                   ExcInternalError());

            Assert((*pointer) <
                     std::numeric_limits<types::global_dof_index>::max(),
                   ExcInternalError());
            const types::global_dof_index this_fe_index = *pointer;

            Assert(this_fe_index < dof_handler.fe_collection.size(),
                   ExcInternalError());

            if (this_fe_index == numbers::invalid_dof_index)
              return false;
            else if (this_fe_index == fe_index)
              return true;
            else
              pointer += dof_handler.get_fe(this_fe_index).dofs_per_vertex + 1;
          }
      }

      template <typename DoFHandlerType, bool level_dof_access>
      static void
      set_mg_dof_indices(
        const dealii::DoFAccessor<1, DoFHandlerType, level_dof_access>
          &                                         accessor,
        const int                                   level,
        const std::vector<types::global_dof_index> &dof_indices,
        const unsigned int                          fe_index)
      {
        const FiniteElement<DoFHandlerType::dimension,
                            DoFHandlerType::space_dimension> &fe =
          accessor.get_dof_handler().get_fe(fe_index);
        std::vector<types::global_dof_index>::const_iterator next =
          dof_indices.begin();

        for (const unsigned int vertex : GeometryInfo<1>::vertex_indices())
          for (unsigned int dof = 0; dof < fe.dofs_per_vertex; ++dof)
            accessor.set_mg_vertex_dof_index(
              level, vertex, dof, *next++, fe_index);

        for (unsigned int dof = 0; dof < fe.dofs_per_line; ++dof)
          accessor.set_mg_dof_index(level, dof, *next++);

        Assert(next == dof_indices.end(), ExcInternalError());
      }



      template <typename DoFHandlerType, bool level_dof_access>
      static void set_mg_dof_indices(
        dealii::DoFAccessor<2, DoFHandlerType, level_dof_access> &accessor,
        const int                                                 level,
        const std::vector<types::global_dof_index> &              dof_indices,
        const unsigned int                                        fe_index)
      {
        const FiniteElement<DoFHandlerType::dimension,
                            DoFHandlerType::space_dimension> &fe =
          accessor.get_dof_handler().get_fe(fe_index);
        std::vector<types::global_dof_index>::const_iterator next =
          dof_indices.begin();

        for (const unsigned int vertex : GeometryInfo<2>::vertex_indices())
          for (unsigned int dof = 0; dof < fe.dofs_per_vertex; ++dof)
            accessor.set_mg_vertex_dof_index(
              level, vertex, dof, *next++, fe_index);

        for (unsigned int line = 0; line < GeometryInfo<2>::lines_per_cell;
             ++line)
          for (unsigned int dof = 0; dof < fe.dofs_per_line; ++dof)
            accessor.line(line)->set_mg_dof_index(level, dof, *next++);

        for (unsigned int dof = 0; dof < fe.dofs_per_quad; ++dof)
          accessor.set_mg_dof_index(level, dof, *next++);

        Assert(next == dof_indices.end(), ExcInternalError());
      }



      template <typename DoFHandlerType, bool level_dof_access>
      static void
      set_mg_dof_indices(
        const dealii::DoFAccessor<3, DoFHandlerType, level_dof_access>
          &                                         accessor,
        const int                                   level,
        const std::vector<types::global_dof_index> &dof_indices,
        const unsigned int                          fe_index)
      {
        const FiniteElement<DoFHandlerType::dimension,
                            DoFHandlerType::space_dimension> &fe =
          accessor.get_dof_handler().get_fe(fe_index);
        std::vector<types::global_dof_index>::const_iterator next =
          dof_indices.begin();

        for (const unsigned int vertex : GeometryInfo<3>::vertex_indices())
          for (unsigned int dof = 0; dof < fe.dofs_per_vertex; ++dof)
            accessor.set_mg_vertex_dof_index(
              level, vertex, dof, *next++, fe_index);

        for (unsigned int line = 0; line < GeometryInfo<3>::lines_per_cell;
             ++line)
          for (unsigned int dof = 0; dof < fe.dofs_per_line; ++dof)
            accessor.line(line)->set_mg_dof_index(
              level,
              fe.adjust_line_dof_index_for_line_orientation(
                dof, accessor.line_orientation(line)),
              *next++);

        for (unsigned int quad = 0; quad < GeometryInfo<3>::quads_per_cell;
             ++quad)
          for (unsigned int dof = 0; dof < fe.dofs_per_quad; ++dof)
            accessor.quad(quad)->set_mg_dof_index(
              level,
              fe.adjust_quad_dof_index_for_face_orientation(
                dof,
                accessor.face_orientation(quad),
                accessor.face_flip(quad),
                accessor.face_rotation(quad)),
              *next++);

        for (unsigned int dof = 0; dof < fe.dofs_per_hex; ++dof)
          accessor.set_mg_dof_index(level, dof, *next++);

        Assert(next == dof_indices.end(), ExcInternalError());
      }



      template <typename InputVector, typename ForwardIterator>
      static void
      extract_subvector_to(const InputVector &            values,
                           const types::global_dof_index *cache,
                           const types::global_dof_index *cache_end,
                           ForwardIterator                local_values_begin)
      {
        values.extract_subvector_to(cache, cache_end, local_values_begin);
      }



#if defined(DEAL_II_WITH_TRILINOS) && defined(DEAL_II_WITH_MPI)
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
        const types::global_dof_index *                      cache_begin,
        const types::global_dof_index *                      cache_end,
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
        read_write_vector.import(values, VectorOperation::insert);

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
        read_write_vector.import(values, VectorOperation::insert);

        // Copy the elements from read_write_vector and reorder them.
        for (unsigned int i = 0; i < cache_size; ++i, ++local_values_begin)
          *local_values_begin = read_write_vector[sorted_indices_pos[i]];
      }
#endif



      /**
       * Loop over all degrees of freedom of the object described by the
       * provided @p accessor and @p fe_index and perform the static functions
       * provided by DoFOperation (set/get) on these.
       */
      template <typename DoFHandlerType,
                bool level_dof_access,
                int  structdim,
                typename DoFIndicesType,
                typename DoFOperation>
      static void
      process_dof_indices(
        const dealii::DoFAccessor<structdim, DoFHandlerType, level_dof_access>
          &                accessor,
        DoFIndicesType &   dof_indices,
        const unsigned int fe_index,
        const DoFOperation &)
      {
        const unsigned int                                             //
          dofs_per_vertex = accessor.get_fe(fe_index).dofs_per_vertex, //
          dofs_per_line   = accessor.get_fe(fe_index).dofs_per_line,   //
          dofs_per_quad   = accessor.get_fe(fe_index).dofs_per_quad,   //
          dofs_per_hex    = accessor.get_fe(fe_index).dofs_per_hex;    //

        const unsigned int inner_dofs =
          structdim == 1 ? dofs_per_line :
                           (structdim == 2 ? dofs_per_quad : dofs_per_hex);

        unsigned int index = 0;

        // 1) VERTEX dofs
        for (const unsigned int vertex :
             GeometryInfo<structdim>::vertex_indices())
          for (unsigned int d = 0; d < dofs_per_vertex; ++d, ++index)
            DoFOperation::process_vertex_dof(
              accessor, vertex, d, dof_indices[index], fe_index);

        // 2) copy dof numbers from the LINE. for lines with the wrong
        // orientation (which might occur in 3d), we have already made sure that
        // we're ok by picking the correct vertices (this happens automatically
        // in the vertex() function). however, if the line is in wrong
        // orientation, we look at it in flipped orientation and we will have to
        // adjust the shape function indices that we see to correspond to the
        // correct (face/cell-local) ordering.
        if (structdim == 2 || structdim == 3)
          for (unsigned int line = 0;
               line < GeometryInfo<structdim>::lines_per_cell;
               ++line)
            for (unsigned int d = 0; d < dofs_per_line; ++d, ++index)
              DoFOperation::process_dof(
                *accessor.line(line),
                accessor.get_fe(fe_index)
                  .adjust_line_dof_index_for_line_orientation(
                    d, accessor.line_orientation(line)),
                dof_indices[index],
                fe_index);

        // 3) copy dof numbers from the FACE. for faces with the wrong
        // orientation, we have already made sure that we're ok by picking the
        // correct lines and vertices (this happens automatically in the line()
        // and vertex() functions). however, if the face is in wrong
        // orientation, we look at it in flipped orientation and we will have to
        // adjust the shape function indices that we see to correspond to the
        // correct (cell-local) ordering. The same applies, if the face_rotation
        // or face_orientation is non-standard
        if (structdim == 3)
          for (unsigned int quad = 0;
               quad < GeometryInfo<structdim>::quads_per_cell;
               ++quad)
            for (unsigned int d = 0; d < dofs_per_quad; ++d, ++index)
              DoFOperation::process_dof(
                *accessor.quad(quad),
                accessor.get_fe(fe_index)
                  .adjust_quad_dof_index_for_face_orientation(
                    d,
                    accessor.face_orientation(quad),
                    accessor.face_flip(quad),
                    accessor.face_rotation(quad)),
                dof_indices[index],
                fe_index);

        // 4) INNER dofs
        for (unsigned int d = 0; d < inner_dofs; ++d, ++index)
          DoFOperation::process_dof(accessor, d, dof_indices[index], fe_index);

        AssertDimension(dof_indices.size(), index);
      }



      /**
       * An internal struct encapsulating the task of getting (vertex)
       * DoF indices.
       */
      template <typename DoFHandlerType, bool level_dof_access, int structdim>
      struct DoFIndexGetter
      {
        /**
         * Return vertex DoF index.
         */
        static DEAL_II_ALWAYS_INLINE void
        process_vertex_dof(
          const dealii::DoFAccessor<structdim, DoFHandlerType, level_dof_access>
            &                      accessor,
          const unsigned int       vertex,
          const unsigned int       d,
          types::global_dof_index &index_value,
          const unsigned int       fe_index)
        {
          index_value = accessor.vertex_dof_index(vertex, d, fe_index);
        }

        /**
         * Return DoF index for lines, quads, and inner degrees of freedom.
         */
        template <int structdim_>
        static DEAL_II_ALWAYS_INLINE void
        process_dof(const dealii::DoFAccessor<structdim_,
                                              DoFHandlerType,
                                              level_dof_access> &accessor,
                    const unsigned int                           d,
                    types::global_dof_index &                    index_value,
                    const unsigned int                           fe_index)
        {
          index_value = accessor.dof_index(d, fe_index);
        }

        /**
         * Fallback for DoFInvalidAccessor.
         */
        template <int structdim_>
        static DEAL_II_ALWAYS_INLINE void
        process_dof(
          const dealii::DoFInvalidAccessor<structdim_,
                                           DoFHandlerType::dimension,
                                           DoFHandlerType::space_dimension> &,
          const unsigned int,
          types::global_dof_index &,
          const unsigned int)
        {
          Assert(false, ExcInternalError());
        }
      };



      /**
       * An internal struct encapsulating the task of setting (vertex)
       * DoF indices.
       */
      template <typename DoFHandlerType, bool level_dof_access, int structdim>
      struct DoFIndexSetter
      {
        /**
         * Set vertex DoF index.
         */
        static DEAL_II_ALWAYS_INLINE void
        process_vertex_dof(
          const dealii::DoFAccessor<structdim, DoFHandlerType, level_dof_access>
            &                            accessor,
          const unsigned int             vertex,
          const unsigned int             d,
          const types::global_dof_index &index_value,
          const unsigned int             fe_index)
        {
          accessor.set_vertex_dof_index(vertex, d, index_value, fe_index);
        }

        /**
         * Set DoF index for lines, quads, and inner degrees of freedom.
         */
        template <int structdim_>
        static DEAL_II_ALWAYS_INLINE void
        process_dof(const dealii::DoFAccessor<structdim_,
                                              DoFHandlerType,
                                              level_dof_access> &accessor,
                    const unsigned int                           d,
                    const types::global_dof_index &              index_value,
                    const unsigned int                           fe_index)
        {
          accessor.set_dof_index(d, index_value, fe_index);
        }

        /**
         * Fallback for DoFInvalidAccessor.
         */
        template <int structdim_>
        static DEAL_II_ALWAYS_INLINE void
        process_dof(
          const dealii::DoFInvalidAccessor<structdim_,
                                           DoFHandlerType::dimension,
                                           DoFHandlerType::space_dimension> &,
          const unsigned int,
          const types::global_dof_index &,
          const unsigned int)
        {
          Assert(false, ExcInternalError());
        }
      };



      template <typename DoFHandlerType, bool level_dof_access, int structdim>
      static void
      get_dof_indices(
        const dealii::DoFAccessor<structdim, DoFHandlerType, level_dof_access>
          &                                   accessor,
        std::vector<types::global_dof_index> &dof_indices,
        const unsigned int                    fe_index)
      {
        process_dof_indices(
          accessor,
          dof_indices,
          fe_index,
          DoFIndexGetter<DoFHandlerType, level_dof_access, structdim>());
      }



      template <typename DoFHandlerType, bool level_dof_access, int structdim>
      static void
      set_dof_indices(
        const dealii::DoFAccessor<structdim, DoFHandlerType, level_dof_access>
          &                                         accessor,
        const std::vector<types::global_dof_index> &dof_indices,
        const unsigned int                          fe_index)
      {
        // Note: this function is as general as `get_dof_indices()`. This
        // assert is placed here since it is currently only used by the
        // function DoFCellAccessor::set_dof_indices(), which is called by
        // internal::DoFHandlerImplementation::Policy::Implementation::distribute_dofs().
        // In the case of new use cases, this assert can be removed.
        Assert(
          DoFHandlerType::dimension == structdim,
          ExcMessage(
            "This function is intended to be used for DoFCellAccessor, i.e., dimension == structdim."));

        process_dof_indices(
          accessor,
          dof_indices,
          fe_index,
          DoFIndexSetter<DoFHandlerType, level_dof_access, structdim>());
      }
    };
  } // namespace DoFAccessorImplementation
} // namespace internal



template <int dim, typename DoFHandlerType, bool level_dof_access>
inline types::global_dof_index
DoFAccessor<dim, DoFHandlerType, level_dof_access>::dof_index(
  const unsigned int i,
  const unsigned int fe_index) const
{
  // access the respective DoF
  return dealii::internal::DoFAccessorImplementation::Implementation::
    get_dof_index(*this->dof_handler,
                  this->level(),
                  this->present_index,
                  fe_index,
                  i,
                  std::integral_constant<int, dim>());
}


template <int structdim, typename DoFHandlerType, bool level_dof_access>
inline types::global_dof_index
DoFAccessor<structdim, DoFHandlerType, level_dof_access>::mg_dof_index(
  const int          level,
  const unsigned int i) const
{
  return this->dof_handler->template get_dof_index<structdim>(
    level, this->present_index, 0, i);
}


template <int dim, typename DoFHandlerType, bool level_dof_access>
inline void
DoFAccessor<dim, DoFHandlerType, level_dof_access>::set_dof_index(
  const unsigned int            i,
  const types::global_dof_index index,
  const unsigned int            fe_index) const
{
  // access the respective DoF
  dealii::internal::DoFAccessorImplementation::Implementation::set_dof_index(
    *this->dof_handler,
    this->level(),
    this->present_index,
    fe_index,
    i,
    std::integral_constant<int, dim>(),
    index);
}



template <int dim, typename DoFHandlerType, bool level_dof_access>
inline unsigned int
DoFAccessor<dim, DoFHandlerType, level_dof_access>::n_active_fe_indices() const
{
  // access the respective DoF
  return dealii::internal::DoFAccessorImplementation::Implementation::
    n_active_fe_indices(*this->dof_handler,
                        this->level(),
                        this->present_index,
                        std::integral_constant<int, dim>());
}



template <int dim, typename DoFHandlerType, bool level_dof_access>
inline unsigned int
DoFAccessor<dim, DoFHandlerType, level_dof_access>::nth_active_fe_index(
  const unsigned int n) const
{
  // access the respective DoF
  return dealii::internal::DoFAccessorImplementation::Implementation::
    nth_active_fe_index(*this->dof_handler,
                        this->level(),
                        this->present_index,
                        n,
                        std::integral_constant<int, dim>());
}



template <int dim, typename DoFHandlerType, bool level_dof_access>
inline std::set<unsigned int>
DoFAccessor<dim, DoFHandlerType, level_dof_access>::get_active_fe_indices()
  const
{
  std::set<unsigned int> active_fe_indices;
  for (unsigned int i = 0; i < n_active_fe_indices(); ++i)
    active_fe_indices.insert(nth_active_fe_index(i));
  return active_fe_indices;
}



template <int dim, typename DoFHandlerType, bool level_dof_access>
inline bool
DoFAccessor<dim, DoFHandlerType, level_dof_access>::fe_index_is_active(
  const unsigned int fe_index) const
{
  // access the respective DoF
  return dealii::internal::DoFAccessorImplementation::Implementation::
    fe_index_is_active(*this->dof_handler,
                       this->level(),
                       this->present_index,
                       fe_index,
                       std::integral_constant<int, dim>());
}



template <int structdim, typename DoFHandlerType, bool level_dof_access>
inline types::global_dof_index
DoFAccessor<structdim, DoFHandlerType, level_dof_access>::vertex_dof_index(
  const unsigned int vertex,
  const unsigned int i,
  const unsigned int fe_index) const
{
  return dealii::internal::DoFAccessorImplementation::Implementation::
    get_vertex_dof_index(*this->dof_handler,
                         this->vertex_index(vertex),
                         fe_index,
                         i);
}


template <int structdim, typename DoFHandlerType, bool level_dof_access>
inline types::global_dof_index
DoFAccessor<structdim, DoFHandlerType, level_dof_access>::mg_vertex_dof_index(
  const int          level,
  const unsigned int vertex,
  const unsigned int i,
  const unsigned int fe_index) const
{
  (void)fe_index;
  Assert(this->dof_handler != nullptr, ExcInvalidObject());
  AssertIndexRange(vertex, GeometryInfo<structdim>::vertices_per_cell);
  AssertIndexRange(i, this->dof_handler->get_fe(fe_index).dofs_per_vertex);

  return dealii::internal::DoFAccessorImplementation::Implementation::
    mg_vertex_dof_index(*this->dof_handler,
                        level,
                        this->vertex_index(vertex),
                        i);
}


template <int structdim, typename DoFHandlerType, bool level_dof_access>
inline void
DoFAccessor<structdim, DoFHandlerType, level_dof_access>::set_vertex_dof_index(
  const unsigned int            vertex,
  const unsigned int            i,
  const types::global_dof_index index,
  const unsigned int            fe_index) const
{
  dealii::internal::DoFAccessorImplementation::Implementation::
    set_vertex_dof_index(
      *this->dof_handler, this->vertex_index(vertex), fe_index, i, index);
}


template <int structdim, typename DoFHandlerType, bool level_dof_access>
inline void
DoFAccessor<structdim, DoFHandlerType, level_dof_access>::
  set_mg_vertex_dof_index(const int                     level,
                          const unsigned int            vertex,
                          const unsigned int            i,
                          const types::global_dof_index index,
                          const unsigned int            fe_index) const
{
  (void)fe_index;
  Assert(this->dof_handler != nullptr, ExcInvalidObject());
  AssertIndexRange(vertex, GeometryInfo<structdim>::vertices_per_cell);
  AssertIndexRange(i, this->dof_handler->get_fe(fe_index).dofs_per_vertex);

  return dealii::internal::DoFAccessorImplementation::Implementation::
    set_mg_vertex_dof_index(
      *this->dof_handler, level, this->vertex_index(vertex), i, index);
}


template <int structdim, typename DoFHandlerType, bool level_dof_access>
inline void
DoFAccessor<structdim, DoFHandlerType, level_dof_access>::set_mg_dof_index(
  const int                     level,
  const unsigned int            i,
  const types::global_dof_index index) const
{
  this->dof_handler->template set_dof_index<structdim>(
    level, this->present_index, 0, i, index);
}



template <int dim, typename DoFHandlerType, bool level_dof_access>
inline const FiniteElement<DoFHandlerType::dimension,
                           DoFHandlerType::space_dimension> &
DoFAccessor<dim, DoFHandlerType, level_dof_access>::get_fe(
  const unsigned int fe_index) const
{
  Assert(fe_index_is_active(fe_index) == true,
         ExcMessage("This function can only be called for active fe indices"));

  return this->dof_handler->get_fe(fe_index);
}



namespace internal
{
  namespace DoFAccessorImplementation
  {
    template <typename DoFHandlerType, bool level_dof_access>
    void
    get_mg_dof_indices(
      const dealii::DoFAccessor<1, DoFHandlerType, level_dof_access> &accessor,
      const int                                                       level,
      std::vector<types::global_dof_index> &dof_indices,
      const unsigned int                    fe_index)
    {
      const DoFHandlerType &handler = accessor.get_dof_handler();

      const FiniteElement<DoFHandlerType::dimension,
                          DoFHandlerType::space_dimension> &fe =
        handler.get_fe(fe_index);
      std::vector<types::global_dof_index>::iterator next = dof_indices.begin();

      for (const unsigned int vertex : GeometryInfo<1>::vertex_indices())
        for (unsigned int dof = 0; dof < fe.dofs_per_vertex; ++dof)
          *next++ = accessor.mg_vertex_dof_index(level, vertex, dof);

      for (unsigned int dof = 0; dof < fe.dofs_per_line; ++dof)
        *next++ = accessor.mg_dof_index(level, dof);

      Assert(next == dof_indices.end(), ExcInternalError());
    }



    template <typename DoFHandlerType, bool level_dof_access>
    void
    get_mg_dof_indices(
      const dealii::DoFAccessor<2, DoFHandlerType, level_dof_access> &accessor,
      const int                                                       level,
      std::vector<types::global_dof_index> &dof_indices,
      const unsigned int                    fe_index)
    {
      const DoFHandlerType &handler = accessor.get_dof_handler();

      const FiniteElement<DoFHandlerType::dimension,
                          DoFHandlerType::space_dimension> &fe =
        handler.get_fe(fe_index);
      std::vector<types::global_dof_index>::iterator next = dof_indices.begin();

      for (const unsigned int vertex : GeometryInfo<2>::vertex_indices())
        for (unsigned int dof = 0; dof < fe.dofs_per_vertex; ++dof)
          *next++ = accessor.mg_vertex_dof_index(level, vertex, dof);

      for (unsigned int line = 0; line < GeometryInfo<2>::lines_per_cell;
           ++line)
        for (unsigned int dof = 0; dof < fe.dofs_per_line; ++dof)
          *next++ = accessor.line(line)->mg_dof_index(level, dof);

      for (unsigned int dof = 0; dof < fe.dofs_per_quad; ++dof)
        *next++ = accessor.mg_dof_index(level, dof);

      Assert(next == dof_indices.end(), ExcInternalError());
    }



    template <typename DoFHandlerType, bool level_dof_access>
    void
    get_mg_dof_indices(
      const dealii::DoFAccessor<3, DoFHandlerType, level_dof_access> &accessor,
      const int                                                       level,
      std::vector<types::global_dof_index> &dof_indices,
      const unsigned int                    fe_index)
    {
      const DoFHandlerType &handler = accessor.get_dof_handler();

      const FiniteElement<DoFHandlerType::dimension,
                          DoFHandlerType::space_dimension> &fe =
        handler.get_fe(fe_index);
      std::vector<types::global_dof_index>::iterator next = dof_indices.begin();

      for (const unsigned int vertex : GeometryInfo<3>::vertex_indices())
        for (unsigned int dof = 0; dof < fe.dofs_per_vertex; ++dof)
          *next++ = accessor.mg_vertex_dof_index(level, vertex, dof);

      for (unsigned int line = 0; line < GeometryInfo<3>::lines_per_cell;
           ++line)
        for (unsigned int dof = 0; dof < fe.dofs_per_line; ++dof)
          *next++ = accessor.line(line)->mg_dof_index(
            level,
            accessor.get_fe(fe_index)
              .adjust_line_dof_index_for_line_orientation(
                dof, accessor.line_orientation(line)));

      for (unsigned int quad = 0; quad < GeometryInfo<3>::quads_per_cell;
           ++quad)
        for (unsigned int dof = 0; dof < fe.dofs_per_quad; ++dof)
          *next++ = accessor.quad(quad)->mg_dof_index(
            level,
            accessor.get_fe(fe_index)
              .adjust_quad_dof_index_for_face_orientation(
                dof,
                accessor.face_orientation(quad),
                accessor.face_flip(quad),
                accessor.face_rotation(quad)));

      for (unsigned int dof = 0; dof < fe.dofs_per_hex; ++dof)
        *next++ = accessor.mg_dof_index(level, dof);

      Assert(next == dof_indices.end(), ExcInternalError());
    }


  } // namespace DoFAccessorImplementation
} // namespace internal


template <int structdim, typename DoFHandlerType, bool level_dof_access>
inline void
DoFAccessor<structdim, DoFHandlerType, level_dof_access>::get_dof_indices(
  std::vector<types::global_dof_index> &dof_indices,
  const unsigned int                    fe_index) const
{
  Assert(this->dof_handler != nullptr, ExcInvalidObject());
  Assert(static_cast<unsigned int>(this->level()) <
           this->dof_handler->levels.size(),
         ExcMessage(
           "The DoFHandler to which this accessor points has not "
           "been initialized, i.e., it doesn't appear that DoF indices "
           "have been distributed on it."));

  switch (structdim)
    {
      case 1:
        Assert(dof_indices.size() ==
                 (GeometryInfo<1>::vertices_per_cell *
                    this->dof_handler->get_fe(fe_index).dofs_per_vertex +
                  this->dof_handler->get_fe(fe_index).dofs_per_line),
               ExcVectorDoesNotMatch());
        break;
      case 2:
        Assert(dof_indices.size() ==
                 (GeometryInfo<2>::vertices_per_cell *
                    this->dof_handler->get_fe(fe_index).dofs_per_vertex +
                  GeometryInfo<2>::lines_per_cell *
                    this->dof_handler->get_fe(fe_index).dofs_per_line +
                  this->dof_handler->get_fe(fe_index).dofs_per_quad),
               ExcVectorDoesNotMatch());
        break;
      case 3:
        Assert(dof_indices.size() ==
                 (GeometryInfo<3>::vertices_per_cell *
                    this->dof_handler->get_fe(fe_index).dofs_per_vertex +
                  GeometryInfo<3>::lines_per_cell *
                    this->dof_handler->get_fe(fe_index).dofs_per_line +
                  GeometryInfo<3>::faces_per_cell *
                    this->dof_handler->get_fe(fe_index).dofs_per_quad +
                  this->dof_handler->get_fe(fe_index).dofs_per_hex),
               ExcVectorDoesNotMatch());
        break;
      default:
        Assert(false, ExcNotImplemented());
    }


  // this function really only makes
  // sense if either a) there are
  // degrees of freedom defined on
  // the present object, or b) the
  // object is non-active objects but
  // all degrees of freedom are
  // located on vertices, since
  // otherwise there are degrees of
  // freedom on sub-objects which are
  // not allocated for this
  // non-active thing
  Assert(this->fe_index_is_active(fe_index) ||
           (this->dof_handler->get_fe(fe_index).dofs_per_cell ==
            GeometryInfo<structdim>::vertices_per_cell *
              this->dof_handler->get_fe(fe_index).dofs_per_vertex),
         ExcInternalError());

  // now do the actual work
  dealii::internal::DoFAccessorImplementation::Implementation::get_dof_indices(
    *this, dof_indices, fe_index);
}



template <int structdim, typename DoFHandlerType, bool level_dof_access>
inline void
DoFAccessor<structdim, DoFHandlerType, level_dof_access>::get_mg_dof_indices(
  const int                             level,
  std::vector<types::global_dof_index> &dof_indices,
  const unsigned int                    fe_index) const
{
  Assert(this->dof_handler != nullptr, ExcInvalidObject());

  switch (structdim)
    {
      case 1:
        {
          Assert(dof_indices.size() ==
                   GeometryInfo<1>::vertices_per_cell *
                       this->dof_handler->get_fe(fe_index).dofs_per_vertex +
                     this->dof_handler->get_fe(fe_index).dofs_per_line,
                 ExcVectorDoesNotMatch());
          break;
        }

      case 2:
        {
          Assert(dof_indices.size() ==
                   GeometryInfo<2>::vertices_per_cell *
                       this->dof_handler->get_fe(fe_index).dofs_per_vertex +
                     GeometryInfo<2>::lines_per_cell *
                       this->dof_handler->get_fe(fe_index).dofs_per_line +
                     this->dof_handler->get_fe(fe_index).dofs_per_quad,
                 ExcVectorDoesNotMatch());
          break;
        }

      case 3:
        {
          Assert(dof_indices.size() ==
                   GeometryInfo<3>::vertices_per_cell *
                       this->dof_handler->get_fe(fe_index).dofs_per_vertex +
                     GeometryInfo<3>::lines_per_cell *
                       this->dof_handler->get_fe(fe_index).dofs_per_line +
                     GeometryInfo<3>::faces_per_cell *
                       this->dof_handler->get_fe(fe_index).dofs_per_quad +
                     this->dof_handler->get_fe(fe_index).dofs_per_hex,
                 ExcVectorDoesNotMatch());
          break;
        }

      default:
        Assert(false, ExcNotImplemented());
    }

  internal::DoFAccessorImplementation::get_mg_dof_indices(*this,
                                                          level,
                                                          dof_indices,
                                                          fe_index);
}


template <int structdim, typename DoFHandlerType, bool level_dof_access>
inline void
DoFAccessor<structdim, DoFHandlerType, level_dof_access>::set_mg_dof_indices(
  const int                                   level,
  const std::vector<types::global_dof_index> &dof_indices,
  const unsigned int                          fe_index)
{
  Assert(this->dof_handler != nullptr, ExcInvalidObject());

  switch (structdim)
    {
      case 1:
        {
          Assert(dof_indices.size() ==
                   GeometryInfo<1>::vertices_per_cell *
                       this->dof_handler->get_fe(fe_index).dofs_per_vertex +
                     this->dof_handler->get_fe(fe_index).dofs_per_line,
                 ExcVectorDoesNotMatch());
          break;
        }

      case 2:
        {
          Assert(dof_indices.size() ==
                   GeometryInfo<2>::vertices_per_cell *
                       this->dof_handler->get_fe(fe_index).dofs_per_vertex +
                     GeometryInfo<2>::lines_per_cell *
                       this->dof_handler->get_fe(fe_index).dofs_per_line +
                     this->dof_handler->get_fe(fe_index).dofs_per_quad,
                 ExcVectorDoesNotMatch());
          break;
        }

      case 3:
        {
          Assert(dof_indices.size() ==
                   GeometryInfo<3>::vertices_per_cell *
                       this->dof_handler->get_fe(fe_index).dofs_per_vertex +
                     GeometryInfo<3>::lines_per_cell *
                       this->dof_handler->get_fe(fe_index).dofs_per_line +
                     GeometryInfo<3>::faces_per_cell *
                       this->dof_handler->get_fe(fe_index).dofs_per_quad +
                     this->dof_handler->get_fe(fe_index).dofs_per_hex,
                 ExcVectorDoesNotMatch());
          break;
        }

      default:
        Assert(false, ExcNotImplemented());
    }

  internal::DoFAccessorImplementation::Implementation::set_mg_dof_indices(
    *this, level, dof_indices, fe_index);
}


template <int structdim, typename DoFHandlerType, bool level_dof_access>
inline typename dealii::internal::DoFHandlerImplementation::
  Iterators<DoFHandlerType, level_dof_access>::line_iterator
  DoFAccessor<structdim, DoFHandlerType, level_dof_access>::line(
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
      return typename dealii::internal::DoFHandlerImplementation::Iterators<
        DoFHandlerType,
        level_dof_access>::cell_iterator(&this->get_triangulation(),
                                         this->level(),
                                         this->index(),
                                         &this->get_dof_handler());
    }

  // otherwise we need to be in structdim>=2
  Assert(structdim > 1, ExcImpossibleInDim(structdim));
  Assert(DoFHandlerType::dimension > 1,
         ExcImpossibleInDim(DoFHandlerType::dimension));

  // checking of 'i' happens in line_index(i)
  return typename dealii::internal::DoFHandlerImplementation::
    Iterators<DoFHandlerType, level_dof_access>::line_iterator(
      this->tria,
      0, // only sub-objects are allowed, which have no level
      this->line_index(i),
      this->dof_handler);
}


template <int structdim, typename DoFHandlerType, bool level_dof_access>
inline typename dealii::internal::DoFHandlerImplementation::
  Iterators<DoFHandlerType, level_dof_access>::quad_iterator
  DoFAccessor<structdim, DoFHandlerType, level_dof_access>::quad(
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
      return typename dealii::internal::DoFHandlerImplementation::Iterators<
        DoFHandlerType>::cell_iterator(&this->get_triangulation(),
                                       this->level(),
                                       this->index(),
                                       &this->get_dof_handler());
    }

  // otherwise we need to be in structdim>=3
  Assert(structdim > 2, ExcImpossibleInDim(structdim));
  Assert(DoFHandlerType::dimension > 2,
         ExcImpossibleInDim(DoFHandlerType::dimension));

  // checking of 'i' happens in quad_index(i)
  return typename dealii::internal::DoFHandlerImplementation::
    Iterators<DoFHandlerType, level_dof_access>::quad_iterator(
      this->tria,
      0, // only sub-objects are allowed, which have no level
      this->quad_index(i),
      this->dof_handler);
}


/*----------------- Functions: DoFAccessor<0,1,spacedim> --------------------*/


template <template <int, int> class DoFHandlerType,
          int  spacedim,
          bool level_dof_access>
inline DoFAccessor<0, DoFHandlerType<1, spacedim>, level_dof_access>::
  DoFAccessor()
{
  Assert(false, ExcInvalidObject());
}



template <template <int, int> class DoFHandlerType,
          int  spacedim,
          bool level_dof_access>
inline DoFAccessor<0, DoFHandlerType<1, spacedim>, level_dof_access>::
  DoFAccessor(
    const Triangulation<1, spacedim> *                      tria,
    const typename TriaAccessor<0, 1, spacedim>::VertexKind vertex_kind,
    const unsigned int                                      vertex_index,
    const DoFHandlerType<1, spacedim> *                     dof_handler)
  : BaseClass(tria, vertex_kind, vertex_index)
  , dof_handler(const_cast<DoFHandlerType<1, spacedim> *>(dof_handler))
{}



template <template <int, int> class DoFHandlerType,
          int  spacedim,
          bool level_dof_access>
inline DoFAccessor<0, DoFHandlerType<1, spacedim>, level_dof_access>::
  DoFAccessor(const Triangulation<1, spacedim> *,
              const int,
              const int,
              const DoFHandlerType<1, spacedim> *)
  : dof_handler(nullptr)
{
  Assert(false,
         ExcMessage(
           "This constructor can not be called for face iterators in 1d."));
}



template <template <int, int> class DoFHandlerType,
          int  spacedim,
          bool level_dof_access>
template <int structdim2, int dim2, int spacedim2>
DoFAccessor<0, DoFHandlerType<1, spacedim>, level_dof_access>::DoFAccessor(
  const InvalidAccessor<structdim2, dim2, spacedim2> &)
{
  Assert(false, ExcInvalidObject());
}



template <template <int, int> class DoFHandlerType,
          int  spacedim,
          bool level_dof_access>
template <int dim2, class DoFHandlerType2, bool level_dof_access2>
inline DoFAccessor<0, DoFHandlerType<1, spacedim>, level_dof_access>::
  DoFAccessor(const DoFAccessor<dim2, DoFHandlerType2, level_dof_access2> &)
{
  Assert(false, ExcInvalidObject());
}



template <template <int, int> class DoFHandlerType,
          int  spacedim,
          bool level_dof_access>
inline void DoFAccessor<0, DoFHandlerType<1, spacedim>, level_dof_access>::
  set_dof_handler(DoFHandlerType<1, spacedim> *dh)
{
  Assert(dh != nullptr, ExcInvalidObject());
  this->dof_handler = dh;
}



template <template <int, int> class DoFHandlerType,
          int  spacedim,
          bool level_dof_access>
inline void
DoFAccessor<0, DoFHandlerType<1, spacedim>, level_dof_access>::set_dof_index(
  const unsigned int /*i*/,
  const types::global_dof_index /*index*/,
  const unsigned int /*fe_index*/) const
{
  Assert(false, ExcNotImplemented());
}



template <template <int, int> class DoFHandlerType,
          int  spacedim,
          bool level_dof_access>
inline void
DoFAccessor<0, DoFHandlerType<1, spacedim>, level_dof_access>::
  set_vertex_dof_index(const unsigned int /*vertex*/,
                       const unsigned int /*i*/,
                       const types::global_dof_index /*index*/,
                       const unsigned int /*fe_index*/) const
{
  Assert(false, ExcNotImplemented());
}



template <template <int, int> class DoFHandlerType,
          int  spacedim,
          bool level_dof_access>
inline const DoFHandlerType<1, spacedim> &
DoFAccessor<0, DoFHandlerType<1, spacedim>, level_dof_access>::get_dof_handler()
  const
{
  return *this->dof_handler;
}



template <template <int, int> class DoFHandlerType,
          int  spacedim,
          bool level_dof_access>
inline void
DoFAccessor<0, DoFHandlerType<1, spacedim>, level_dof_access>::get_dof_indices(
  std::vector<types::global_dof_index> &dof_indices,
  const unsigned int                    fe_index) const
{
  for (unsigned int i = 0; i < dof_indices.size(); ++i)
    dof_indices[i] = dealii::internal::DoFAccessorImplementation::
      Implementation::get_vertex_dof_index(*dof_handler,
                                           this->global_vertex_index,
                                           fe_index,
                                           i);
}



template <template <int, int> class DoFHandlerType,
          int  spacedim,
          bool level_dof_access>
inline void
DoFAccessor<0, DoFHandlerType<1, spacedim>, level_dof_access>::
  get_mg_dof_indices(const int                             level,
                     std::vector<types::global_dof_index> &dof_indices,
                     const unsigned int                    fe_index) const
{
  AssertThrow(fe_index == 0, ExcMessage("Unknown triangulation!"));

  for (unsigned int i = 0; i < dof_indices.size(); ++i)
    dof_indices[i] =
      dealii::internal::DoFAccessorImplementation::Implementation::
        mg_vertex_dof_index(*dof_handler, level, this->global_vertex_index, i);
}



template <template <int, int> class DoFHandlerType,
          int  spacedim,
          bool level_dof_access>
inline types::global_dof_index
DoFAccessor<0, DoFHandlerType<1, spacedim>, level_dof_access>::vertex_dof_index(
  const unsigned int vertex,
  const unsigned int i,
  const unsigned int fe_index) const
{
  (void)vertex;
  AssertIndexRange(vertex, 1);
  return dealii::internal::DoFAccessorImplementation::Implementation::
    get_vertex_dof_index(*dof_handler, this->global_vertex_index, fe_index, i);
}



template <template <int, int> class DoFHandlerType,
          int  spacedim,
          bool level_dof_access>
inline types::global_dof_index
DoFAccessor<0, DoFHandlerType<1, spacedim>, level_dof_access>::dof_index(
  const unsigned int i,
  const unsigned int fe_index) const
{
  return dealii::internal::DoFAccessorImplementation::Implementation::
    get_vertex_dof_index(*this->dof_handler,
                         this->vertex_index(0),
                         fe_index,
                         i);
}



template <template <int, int> class DoFHandlerType,
          int  spacedim,
          bool level_dof_access>
inline unsigned int
DoFAccessor<0, DoFHandlerType<1, spacedim>, level_dof_access>::
  n_active_fe_indices() const
{
  Assert((std::is_same<DoFHandlerType<1, spacedim>,
                       dealii::DoFHandler<1, spacedim>>::value == true),
         ExcNotImplemented());

  return 1;
}



template <template <int, int> class DoFHandlerType,
          int  spacedim,
          bool level_dof_access>
inline unsigned int
DoFAccessor<0, DoFHandlerType<1, spacedim>, level_dof_access>::
  nth_active_fe_index(const unsigned int /*n*/) const
{
  Assert((std::is_same<DoFHandlerType<1, spacedim>,
                       dealii::DoFHandler<1, spacedim>>::value == true),
         ExcNotImplemented());

  return 0;
}



template <template <int, int> class DoFHandlerType,
          int  spacedim,
          bool level_dof_access>
inline bool
DoFAccessor<0, DoFHandlerType<1, spacedim>, level_dof_access>::
  fe_index_is_active(const unsigned int /*fe_index*/) const
{
  Assert(false, ExcNotImplemented());
  return false;
}



template <template <int, int> class DoFHandlerType,
          int  spacedim,
          bool level_dof_access>
inline const FiniteElement<DoFHandlerType<1, spacedim>::dimension,
                           DoFHandlerType<1, spacedim>::space_dimension> &
DoFAccessor<0, DoFHandlerType<1, spacedim>, level_dof_access>::get_fe(
  const unsigned int fe_index) const
{
  Assert(this->dof_handler != nullptr, ExcInvalidObject());
  return dof_handler->get_fe(fe_index);
}



template <template <int, int> class DoFHandlerType,
          int  spacedim,
          bool level_dof_access>
inline void
DoFAccessor<0, DoFHandlerType<1, spacedim>, level_dof_access>::copy_from(
  const TriaAccessorBase<0, 1, spacedim> &da)
{
  Assert(this->dof_handler != nullptr, ExcInvalidObject());
  BaseClass::copy_from(da);
}



template <template <int, int> class DoFHandlerType,
          int  spacedim,
          bool level_dof_access>
template <bool level_dof_access2>
inline void
DoFAccessor<0, DoFHandlerType<1, spacedim>, level_dof_access>::copy_from(
  const DoFAccessor<0, DoFHandlerType<1, spacedim>, level_dof_access2> &a)
{
  BaseClass::copy_from(a);
  set_dof_handler(a.dof_handler);
}



template <template <int, int> class DoFHandlerType,
          int  spacedim,
          bool level_dof_access>
inline TriaIterator<
  DoFAccessor<0, DoFHandlerType<1, spacedim>, level_dof_access>>
DoFAccessor<0, DoFHandlerType<1, spacedim>, level_dof_access>::child(
  const unsigned int /*i*/) const
{
  return TriaIterator<
    DoFAccessor<0, DoFHandlerType<1, spacedim>, level_dof_access>>();
}



template <template <int, int> class DoFHandlerType,
          int  spacedim,
          bool level_dof_access>
inline typename dealii::internal::DoFHandlerImplementation::
  Iterators<DoFHandlerType<1, spacedim>, level_dof_access>::line_iterator
  DoFAccessor<0, DoFHandlerType<1, spacedim>, level_dof_access>::line(
    const unsigned int /*c*/) const
{
  Assert(false, ExcNotImplemented());
  return typename dealii::internal::DoFHandlerImplementation::
    Iterators<DoFHandlerType<1, spacedim>, level_dof_access>::line_iterator();
}



template <template <int, int> class DoFHandlerType,
          int  spacedim,
          bool level_dof_access>
inline typename dealii::internal::DoFHandlerImplementation::
  Iterators<DoFHandlerType<1, spacedim>, level_dof_access>::quad_iterator
  DoFAccessor<0, DoFHandlerType<1, spacedim>, level_dof_access>::quad(
    const unsigned int /*c*/) const
{
  Assert(false, ExcNotImplemented());
  return typename dealii::internal::DoFHandlerImplementation::
    Iterators<DoFHandlerType<1, spacedim>, level_dof_access>::quad_iterator();
}



template <template <int, int> class DoFHandlerType,
          int  spacedim,
          bool level_dof_access>
template <int dim2, class DoFHandlerType2, bool level_dof_access2>
inline bool
DoFAccessor<0, DoFHandlerType<1, spacedim>, level_dof_access>::
operator==(const DoFAccessor<dim2, DoFHandlerType2, level_dof_access2> &a) const
{
  Assert(dim2 == 0, ExcCantCompareIterators());
  Assert(this->dof_handler == a.dof_handler, ExcCantCompareIterators());
  return (BaseClass::operator==(a));
}



template <template <int, int> class DoFHandlerType,
          int  spacedim,
          bool level_dof_access>
template <int dim2, class DoFHandlerType2, bool level_dof_access2>
inline bool
DoFAccessor<0, DoFHandlerType<1, spacedim>, level_dof_access>::
operator!=(const DoFAccessor<dim2, DoFHandlerType2, level_dof_access2> &a) const
{
  Assert(dim2 == 0, ExcCantCompareIterators());
  Assert(this->dof_handler == a.dof_handler, ExcCantCompareIterators());
  return (BaseClass::operator!=(a));
}



/*------------------------- Functions: DoFCellAccessor -----------------------*/


namespace internal
{
  namespace DoFCellAccessorImplementation
  {
    // make sure we refer to class
    // dealii::DoFCellAccessor, not
    // namespace
    // dealii::internal::DoFCellAccessor
    using dealii::DoFCellAccessor;
    using dealii::DoFHandler;

    /**
     * A class with the same purpose as the similarly named class of the
     * Triangulation class. See there for more information.
     */
    struct Implementation
    {
      /**
       * Implement the updating of the cache.
       */
      template <typename DoFHandlerType, bool level_dof_access>
      static void
      update_cell_dof_indices_cache(
        const DoFCellAccessor<DoFHandlerType, level_dof_access> &accessor)
      {
        // caches are only for cells with DoFs, i.e., for active ones and not
        // FE_Nothing
        if (accessor.has_children())
          return;
        const unsigned int dofs_per_cell = accessor.get_fe().dofs_per_cell;
        if (dofs_per_cell == 0)
          return;

        // call the get_dof_indices() function of DoFAccessor, which goes
        // through all the parts of the cell to get the indices by hand. the
        // corresponding function of DoFCellAccessor can then later use the
        // cache
        std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
        static_cast<const dealii::DoFAccessor<DoFHandlerType::dimension,
                                              DoFHandlerType,
                                              level_dof_access> &>(accessor)
          .get_dof_indices(dof_indices, accessor.active_fe_index());

        types::global_dof_index *next_dof_index =
          const_cast<types::global_dof_index *>(
            accessor.dof_handler->levels[accessor.present_level]
              ->get_cell_cache_start(accessor.present_index, dofs_per_cell));

        for (unsigned int i = 0; i < dofs_per_cell; ++i, ++next_dof_index)
          *next_dof_index = dof_indices[i];
      }



      /**
       * Do what the active_fe_index function in the parent class is supposed to
       * do.
       */
      template <int dim, int spacedim, bool level_dof_access>
      static unsigned int
      active_fe_index(
        const DoFCellAccessor<DoFHandler<dim, spacedim>, level_dof_access> &)
      {
        // ::DoFHandler only supports a single active fe with index zero
        return 0;
      }



      template <int dim, int spacedim, bool level_dof_access>
      static unsigned int
      active_fe_index(
        const DoFCellAccessor<dealii::hp::DoFHandler<dim, spacedim>,
                              level_dof_access> &accessor)
      {
        Assert(
          accessor.dof_handler != nullptr,
          (typename std::decay<decltype(accessor)>::type::ExcInvalidObject()));
        Assert(static_cast<unsigned int>(accessor.level()) <
                 accessor.dof_handler->levels.size(),
               ExcMessage("DoFHandler not initialized"));

        return accessor.dof_handler->levels[accessor.level()]->active_fe_index(
          accessor.present_index);
      }



      /**
       * Do what the set_active_fe_index function in the parent class is
       * supposed to do.
       */
      template <int dim, int spacedim, bool level_dof_access>
      static void
      set_active_fe_index(const DoFCellAccessor<DoFHandler<dim, spacedim>,
                                                level_dof_access> &accessor,
                          const unsigned int                       i)
      {
        (void)accessor;
        (void)i;
        // ::DoFHandler only supports a single active fe with index zero
        Assert(
          i == 0,
          (typename std::decay<decltype(accessor)>::type::ExcInvalidObject()));
      }



      template <int dim, int spacedim, bool level_dof_access>
      static void
      set_active_fe_index(
        const DoFCellAccessor<dealii::hp::DoFHandler<dim, spacedim>,
                              level_dof_access> &accessor,
        const unsigned int                       i)
      {
        Assert(
          accessor.dof_handler != nullptr,
          (typename std::decay<decltype(accessor)>::type::ExcInvalidObject()));
        Assert(static_cast<unsigned int>(accessor.level()) <
                 accessor.dof_handler->levels.size(),
               ExcMessage("DoFHandler not initialized"));

        accessor.dof_handler->levels[accessor.level()]->set_active_fe_index(
          accessor.present_index, i);
      }



      /**
       * Do what the future_fe_index function in the parent class is supposed to
       * do.
       */
      template <int dim, int spacedim, bool level_dof_access>
      static unsigned int
      future_fe_index(
        const DoFCellAccessor<DoFHandler<dim, spacedim>, level_dof_access> &)
      {
        // ::DoFHandler only supports a single active fe with index zero
        return 0;
      }



      template <int dim, int spacedim, bool level_dof_access>
      static unsigned int
      future_fe_index(
        const DoFCellAccessor<dealii::hp::DoFHandler<dim, spacedim>,
                              level_dof_access> &accessor)
      {
        Assert(
          accessor.dof_handler != nullptr,
          (typename std::decay<decltype(accessor)>::type::ExcInvalidObject()));
        Assert(static_cast<unsigned int>(accessor.level()) <
                 accessor.dof_handler->levels.size(),
               ExcMessage("DoFHandler not initialized"));

        return accessor.dof_handler->levels[accessor.level()]->future_fe_index(
          accessor.present_index);
      }


      /**
       * Do what the set_future_fe_index function in the parent class is
       * supposed to do.
       */
      template <int dim, int spacedim, bool level_dof_access>
      static void
      set_future_fe_index(const DoFCellAccessor<DoFHandler<dim, spacedim>,
                                                level_dof_access> &accessor,
                          const unsigned int                       i)
      {
        (void)accessor;
        (void)i;
        // ::DoFHandler only supports a single active fe with index zero
        Assert(
          i == 0,
          (typename std::decay<decltype(accessor)>::type::ExcInvalidObject()));
      }



      template <int dim, int spacedim, bool level_dof_access>
      static void
      set_future_fe_index(
        const DoFCellAccessor<dealii::hp::DoFHandler<dim, spacedim>,
                              level_dof_access> &accessor,
        const unsigned int                       i)
      {
        Assert(
          accessor.dof_handler != nullptr,
          (typename std::decay<decltype(accessor)>::type::ExcInvalidObject()));
        Assert(static_cast<unsigned int>(accessor.level()) <
                 accessor.dof_handler->levels.size(),
               ExcMessage("DoFHandler not initialized"));

        accessor.dof_handler->levels[accessor.level()]->set_future_fe_index(
          accessor.present_index, i);
      }



      /**
       * Do what the future_fe_index_set function in the parent class is
       * supposed to do.
       */
      template <int dim, int spacedim, bool level_dof_access>
      static bool
      future_fe_index_set(
        const DoFCellAccessor<dealii::DoFHandler<dim, spacedim>,
                              level_dof_access> &)
      {
        // ::DoFHandler only supports a single active fe with index zero
        return false;
      }



      template <int dim, int spacedim, bool level_dof_access>
      static bool
      future_fe_index_set(
        const DoFCellAccessor<dealii::hp::DoFHandler<dim, spacedim>,
                              level_dof_access> &accessor)
      {
        Assert(
          accessor.dof_handler != nullptr,
          (typename std::decay<decltype(accessor)>::type::ExcInvalidObject()));
        Assert(static_cast<unsigned int>(accessor.level()) <
                 accessor.dof_handler->levels.size(),
               ExcMessage("DoFHandler not initialized"));

        return accessor.dof_handler->levels[accessor.level()]
          ->future_fe_index_set(accessor.present_index);
      }



      /**
       * Do what the clear_fe_index function in the parent class is supposed to
       * do.
       */
      template <int dim, int spacedim, bool level_dof_access>
      static void
      clear_future_fe_index(
        const DoFCellAccessor<dealii::DoFHandler<dim, spacedim>,
                              level_dof_access> &)
      {
        // ::DoFHandler only supports a single active fe with index zero
      }



      template <int dim, int spacedim, bool level_dof_access>
      static void
      clear_future_fe_index(
        const DoFCellAccessor<dealii::hp::DoFHandler<dim, spacedim>,
                              level_dof_access> &accessor)
      {
        Assert(
          accessor.dof_handler != nullptr,
          (typename std::decay<decltype(accessor)>::type::ExcInvalidObject()));
        Assert(static_cast<unsigned int>(accessor.level()) <
                 accessor.dof_handler->levels.size(),
               ExcMessage("DoFHandler not initialized"));

        accessor.dof_handler->levels[accessor.level()]->clear_future_fe_index(
          accessor.present_index);
      }
    };
  } // namespace DoFCellAccessorImplementation
} // namespace internal


template <typename DoFHandlerType, bool level_dof_access>
inline DoFCellAccessor<DoFHandlerType, level_dof_access>::DoFCellAccessor(
  const Triangulation<DoFHandlerType::dimension,
                      DoFHandlerType::space_dimension> *tria,
  const int                                             level,
  const int                                             index,
  const AccessorData *                                  local_data)
  : DoFAccessor<DoFHandlerType::dimension, DoFHandlerType, level_dof_access>(
      tria,
      level,
      index,
      local_data)
{}


template <typename DoFHandlerType, bool level_dof_access>
template <int structdim2, int dim2, int spacedim2>
inline DoFCellAccessor<DoFHandlerType, level_dof_access>::DoFCellAccessor(
  const InvalidAccessor<structdim2, dim2, spacedim2> &)
{
  Assert(false, typename BaseClass::ExcInvalidObject());
}



template <typename DoFHandlerType, bool level_dof_access>
template <int dim2, class DoFHandlerType2, bool level_dof_access2>
inline DoFCellAccessor<DoFHandlerType, level_dof_access>::DoFCellAccessor(
  const DoFAccessor<dim2, DoFHandlerType2, level_dof_access2> &other)
  : BaseClass(other)
{}


template <typename DoFHandlerType, bool level_dof_access>
inline TriaIterator<DoFCellAccessor<DoFHandlerType, level_dof_access>>
DoFCellAccessor<DoFHandlerType, level_dof_access>::neighbor(
  const unsigned int i) const
{
  TriaIterator<DoFCellAccessor<DoFHandlerType, level_dof_access>> q(
    this->tria,
    this->neighbor_level(i),
    this->neighbor_index(i),
    this->dof_handler);

#ifdef DEBUG
  if (q.state() != IteratorState::past_the_end)
    Assert(q->used(), ExcInternalError());
#endif
  return q;
}


template <typename DoFHandlerType, bool level_dof_access>
inline TriaIterator<DoFCellAccessor<DoFHandlerType, level_dof_access>>
DoFCellAccessor<DoFHandlerType, level_dof_access>::child(
  const unsigned int i) const
{
  TriaIterator<DoFCellAccessor<DoFHandlerType, level_dof_access>> q(
    this->tria, this->level() + 1, this->child_index(i), this->dof_handler);

#ifdef DEBUG
  if (q.state() != IteratorState::past_the_end)
    Assert(q->used(), ExcInternalError());
#endif
  return q;
}


template <typename DoFHandlerType, bool level_dof_access>
inline TriaIterator<DoFCellAccessor<DoFHandlerType, level_dof_access>>
DoFCellAccessor<DoFHandlerType, level_dof_access>::parent() const
{
  TriaIterator<DoFCellAccessor<DoFHandlerType, level_dof_access>> q(
    this->tria, this->level() - 1, this->parent_index(), this->dof_handler);

  return q;
}


namespace internal
{
  namespace DoFCellAccessorImplementation
  {
    template <typename DoFHandlerType, bool level_dof_access>
    inline TriaIterator<dealii::DoFAccessor<DoFHandlerType::dimension - 1,
                                            DoFHandlerType,
                                            level_dof_access>>
    get_face(
      const dealii::DoFCellAccessor<DoFHandlerType, level_dof_access> &cell,
      const unsigned int                                               i,
      const std::integral_constant<int, 1>)
    {
      dealii::DoFAccessor<0, DoFHandlerType, level_dof_access> a(
        &cell.get_triangulation(),
        ((i == 0) && cell.at_boundary(0) ?
           dealii::TriaAccessor<0, 1, DoFHandlerType::space_dimension>::
             left_vertex :
           ((i == 1) && cell.at_boundary(1) ?
              dealii::TriaAccessor<0, 1, DoFHandlerType::space_dimension>::
                right_vertex :
              dealii::TriaAccessor<0, 1, DoFHandlerType::space_dimension>::
                interior_vertex)),
        cell.vertex_index(i),
        &cell.get_dof_handler());
      return dealii::TriaIterator<
        dealii::DoFAccessor<0, DoFHandlerType, level_dof_access>>(a);
    }


    template <typename DoFHandlerType, bool level_dof_access>
    inline TriaIterator<dealii::DoFAccessor<DoFHandlerType::dimension - 1,
                                            DoFHandlerType,
                                            level_dof_access>>
    get_face(
      const dealii::DoFCellAccessor<DoFHandlerType, level_dof_access> &cell,
      const unsigned int                                               i,
      const std::integral_constant<int, 2>)
    {
      return cell.line(i);
    }


    template <typename DoFHandlerType, bool level_dof_access>
    inline TriaIterator<dealii::DoFAccessor<DoFHandlerType::dimension - 1,
                                            DoFHandlerType,
                                            level_dof_access>>
    get_face(
      const dealii::DoFCellAccessor<DoFHandlerType, level_dof_access> &cell,
      const unsigned int                                               i,
      const std::integral_constant<int, 3>)
    {
      return cell.quad(i);
    }
  } // namespace DoFCellAccessorImplementation
} // namespace internal


template <typename DoFHandlerType, bool level_dof_access>
inline typename DoFCellAccessor<DoFHandlerType, level_dof_access>::face_iterator
DoFCellAccessor<DoFHandlerType, level_dof_access>::face(
  const unsigned int i) const
{
  AssertIndexRange(i, GeometryInfo<dim>::faces_per_cell);

  const unsigned int dim = DoFHandlerType::dimension;
  return dealii::internal::DoFCellAccessorImplementation::get_face(
    *this, i, std::integral_constant<int, dim>());
}



template <typename DoFHandlerType, bool level_dof_access>
inline std::array<
  typename DoFCellAccessor<DoFHandlerType, level_dof_access>::face_iterator,
  GeometryInfo<DoFHandlerType::dimension>::faces_per_cell>
DoFCellAccessor<DoFHandlerType, level_dof_access>::face_iterators() const
{
  std::array<
    typename DoFCellAccessor<DoFHandlerType, level_dof_access>::face_iterator,
    GeometryInfo<dim>::faces_per_cell>
    face_iterators;

  const unsigned int dim = DoFHandlerType::dimension;
  for (unsigned int i : GeometryInfo<dim>::face_indices())
    face_iterators[i] =
      dealii::internal::DoFCellAccessorImplementation::get_face(
        *this, i, std::integral_constant<int, dim>());

  return face_iterators;
}



template <typename DoFHandlerType, bool level_dof_access>
inline void
DoFCellAccessor<DoFHandlerType, level_dof_access>::get_dof_indices(
  std::vector<types::global_dof_index> &dof_indices) const
{
  Assert(this->is_active(),
         ExcMessage("get_dof_indices() only works on active cells."));
  Assert(this->is_artificial() == false,
         ExcMessage("Can't ask for DoF indices on artificial cells."));
  AssertDimension(dof_indices.size(), this->get_fe().dofs_per_cell);

  const auto dofs_per_cell = this->get_fe().dofs_per_cell;
  if (dofs_per_cell > 0)
    {
      const types::global_dof_index *cache =
        this->dof_handler->levels[this->present_level]->get_cell_cache_start(
          this->present_index, dofs_per_cell);
      for (unsigned int i = 0; i < dofs_per_cell; ++i, ++cache)
        dof_indices[i] = *cache;
    }
}



template <typename DoFHandlerType, bool level_dof_access>
inline void
DoFCellAccessor<DoFHandlerType, level_dof_access>::get_mg_dof_indices(
  std::vector<types::global_dof_index> &dof_indices) const
{
  DoFAccessor<dim, DoFHandlerType, level_dof_access>::get_mg_dof_indices(
    this->level(), dof_indices);
}



template <typename DoFHandlerType, bool level_dof_access>
inline void
DoFCellAccessor<DoFHandlerType, level_dof_access>::set_mg_dof_indices(
  const std::vector<types::global_dof_index> &dof_indices)
{
  DoFAccessor<dim, DoFHandlerType, level_dof_access>::set_mg_dof_indices(
    this->level(), dof_indices);
}



template <typename DoFHandlerType, bool level_dof_access>
inline void
DoFCellAccessor<DoFHandlerType, level_dof_access>::get_active_or_mg_dof_indices(
  std::vector<types::global_dof_index> &dof_indices) const
{
  if (level_dof_access)
    get_mg_dof_indices(dof_indices);
  else
    get_dof_indices(dof_indices);
}



template <typename DoFHandlerType, bool level_dof_access>
template <class InputVector, typename number>
inline void
DoFCellAccessor<DoFHandlerType, level_dof_access>::get_dof_values(
  const InputVector &values,
  Vector<number> &   local_values) const
{
  get_dof_values(values, local_values.begin(), local_values.end());
}



template <typename DoFHandlerType, bool level_dof_access>
template <class InputVector, typename ForwardIterator>
inline void
DoFCellAccessor<DoFHandlerType, level_dof_access>::get_dof_values(
  const InputVector &values,
  ForwardIterator    local_values_begin,
  ForwardIterator    local_values_end) const
{
  (void)local_values_end;
  Assert(this->is_artificial() == false,
         ExcMessage("Can't ask for DoF indices on artificial cells."));
  Assert(this->is_active(), ExcMessage("Cell must be active."));
  Assert(this->dof_handler != nullptr, typename BaseClass::ExcInvalidObject());

  Assert(static_cast<unsigned int>(local_values_end - local_values_begin) ==
           this->get_fe().dofs_per_cell,
         typename DoFCellAccessor::ExcVectorDoesNotMatch());
  Assert(values.size() == this->get_dof_handler().n_dofs(),
         typename DoFCellAccessor::ExcVectorDoesNotMatch());

  const types::global_dof_index *cache =
    this->dof_handler->levels[this->present_level]->get_cell_cache_start(
      this->present_index, this->get_fe().dofs_per_cell);
  dealii::internal::DoFAccessorImplementation::Implementation::
    extract_subvector_to(values,
                         cache,
                         cache + this->get_fe().dofs_per_cell,
                         local_values_begin);
}



template <typename DoFHandlerType, bool level_dof_access>
template <class InputVector, typename ForwardIterator>
inline void
DoFCellAccessor<DoFHandlerType, level_dof_access>::get_dof_values(
  const AffineConstraints<typename InputVector::value_type> &constraints,
  const InputVector &                                        values,
  ForwardIterator                                            local_values_begin,
  ForwardIterator local_values_end) const
{
  Assert(this->is_artificial() == false,
         ExcMessage("Can't ask for DoF indices on artificial cells."));
  Assert(this->is_active(), ExcMessage("Cell must be active."));

  Assert(static_cast<unsigned int>(local_values_end - local_values_begin) ==
           this->get_fe().dofs_per_cell,
         typename DoFCellAccessor::ExcVectorDoesNotMatch());
  Assert(values.size() == this->get_dof_handler().n_dofs(),
         typename DoFCellAccessor::ExcVectorDoesNotMatch());


  const types::global_dof_index *cache =
    this->dof_handler->levels[this->present_level]->get_cell_cache_start(
      this->present_index, this->get_fe().dofs_per_cell);

  constraints.get_dof_values(values,
                             *cache,
                             local_values_begin,
                             local_values_end);
}



template <typename DoFHandlerType, bool level_dof_access>
template <class OutputVector, typename number>
inline void
DoFCellAccessor<DoFHandlerType, level_dof_access>::set_dof_values(
  const Vector<number> &local_values,
  OutputVector &        values) const
{
  Assert(this->is_artificial() == false,
         ExcMessage("Can't ask for DoF indices on artificial cells."));
  Assert(this->is_active(), ExcMessage("Cell must be active."));

  Assert(static_cast<unsigned int>(local_values.size()) ==
           this->get_fe().dofs_per_cell,
         typename DoFCellAccessor::ExcVectorDoesNotMatch());
  Assert(values.size() == this->get_dof_handler().n_dofs(),
         typename DoFCellAccessor::ExcVectorDoesNotMatch());


  Assert(this->dof_handler != nullptr, typename BaseClass::ExcInvalidObject());
  const types::global_dof_index *cache =
    this->dof_handler->levels[this->present_level]->get_cell_cache_start(
      this->present_index, this->get_fe().dofs_per_cell);

  for (unsigned int i = 0; i < this->get_fe().dofs_per_cell; ++i, ++cache)
    internal::ElementAccess<OutputVector>::set(local_values(i), *cache, values);
}



template <typename DoFHandlerType, bool level_dof_access>
inline const FiniteElement<DoFHandlerType::dimension,
                           DoFHandlerType::space_dimension> &
DoFCellAccessor<DoFHandlerType, level_dof_access>::get_fe() const
{
  Assert(this->dof_handler != nullptr, typename BaseClass::ExcInvalidObject());
  Assert(
    (dynamic_cast<const dealii::DoFHandler<DoFHandlerType::dimension,
                                           DoFHandlerType::space_dimension> *>(
       this->dof_handler) != nullptr) ||
      this->is_active(),
    ExcMessage("In hp::DoFHandler objects, finite elements are only associated "
               "with active cells. Consequently, you can not ask for the "
               "active finite element on cells with children."));

  return this->dof_handler->get_fe(active_fe_index());
}



template <typename DoFHandlerType, bool level_dof_access>
inline unsigned int
DoFCellAccessor<DoFHandlerType, level_dof_access>::active_fe_index() const
{
  Assert(
    (dynamic_cast<const dealii::DoFHandler<DoFHandlerType::dimension,
                                           DoFHandlerType::space_dimension> *>(
       this->dof_handler) != nullptr) ||
      this->is_active(),
    ExcMessage("You can not ask for the active_fe_index on a cell that has "
               "children because no degrees of freedom are assigned "
               "to this cell and, consequently, no finite element "
               "is associated with it."));
  Assert(
    (dynamic_cast<const dealii::DoFHandler<DoFHandlerType::dimension,
                                           DoFHandlerType::space_dimension> *>(
       this->dof_handler) != nullptr) ||
      (this->is_locally_owned() || this->is_ghost()),
    ExcMessage("You can only query active_fe_index information on cells "
               "that are either locally owned or (after distributing "
               "degrees of freedom) are ghost cells."));

  return dealii::internal::DoFCellAccessorImplementation::Implementation::
    active_fe_index(*this);
}



template <typename DoFHandlerType, bool level_dof_access>
inline void
DoFCellAccessor<DoFHandlerType, level_dof_access>::set_active_fe_index(
  const unsigned int i) const
{
  Assert(
    (dynamic_cast<const dealii::DoFHandler<DoFHandlerType::dimension,
                                           DoFHandlerType::space_dimension> *>(
       this->dof_handler) != nullptr) ||
      this->is_active(),
    ExcMessage("You can not set the active_fe_index on a cell that has "
               "children because no degrees of freedom will be assigned "
               "to this cell."));

  Assert(
    (dynamic_cast<const dealii::DoFHandler<DoFHandlerType::dimension,
                                           DoFHandlerType::space_dimension> *>(
       this->dof_handler) != nullptr) ||
      this->is_locally_owned(),
    ExcMessage("You can only set active_fe_index information on cells "
               "that are locally owned. On ghost cells, this information "
               "will automatically be propagated from the owning process "
               "of that cell, and there is no information at all on "
               "artificial cells."));

  dealii::internal::DoFCellAccessorImplementation::Implementation::
    set_active_fe_index(*this, i);
}



template <typename DoFHandlerType, bool level_dof_access>
inline const FiniteElement<DoFHandlerType::dimension,
                           DoFHandlerType::space_dimension> &
DoFCellAccessor<DoFHandlerType, level_dof_access>::get_future_fe() const
{
  Assert(this->dof_handler != nullptr, typename BaseClass::ExcInvalidObject());
  Assert(
    (dynamic_cast<const dealii::DoFHandler<DoFHandlerType::dimension,
                                           DoFHandlerType::space_dimension> *>(
       this->dof_handler) != nullptr) ||
      this->is_active(),
    ExcMessage("In hp::DoFHandler objects, finite elements are only associated "
               "with active cells. Consequently, you can not ask for the "
               "future finite element on cells with children."));

  return this->dof_handler->get_fe(future_fe_index());
}



template <typename DoFHandlerType, bool level_dof_access>
inline unsigned int
DoFCellAccessor<DoFHandlerType, level_dof_access>::future_fe_index() const
{
  Assert(
    (dynamic_cast<const dealii::DoFHandler<DoFHandlerType::dimension,
                                           DoFHandlerType::space_dimension> *>(
       this->dof_handler) != nullptr) ||
      (this->has_children() == false),
    ExcMessage("You can not ask for the future_fe_index on a cell that has "
               "children because no degrees of freedom are assigned "
               "to this cell and, consequently, no finite element "
               "is associated with it."));
  Assert(
    (dynamic_cast<const dealii::DoFHandler<DoFHandlerType::dimension,
                                           DoFHandlerType::space_dimension> *>(
       this->dof_handler) != nullptr) ||
      (this->is_locally_owned()),
    ExcMessage("You can only query future_fe_index information on cells "
               "that are locally owned."));

  return dealii::internal::DoFCellAccessorImplementation::Implementation::
    future_fe_index(*this);
}



template <typename DoFHandlerType, bool level_dof_access>
inline void
DoFCellAccessor<DoFHandlerType, level_dof_access>::set_future_fe_index(
  const unsigned int i) const
{
  Assert(
    (dynamic_cast<const dealii::DoFHandler<DoFHandlerType::dimension,
                                           DoFHandlerType::space_dimension> *>(
       this->dof_handler) != nullptr) ||
      (this->has_children() == false),
    ExcMessage("You can not set the future_fe_index on a cell that has "
               "children because no degrees of freedom will be assigned "
               "to this cell."));

  Assert(
    (dynamic_cast<const dealii::DoFHandler<DoFHandlerType::dimension,
                                           DoFHandlerType::space_dimension> *>(
       this->dof_handler) != nullptr) ||
      this->is_locally_owned(),
    ExcMessage("You can only set future_fe_index information on cells "
               "that are locally owned."));

  dealii::internal::DoFCellAccessorImplementation::Implementation::
    set_future_fe_index(*this, i);
}



template <typename DoFHandlerType, bool level_dof_access>
inline bool
DoFCellAccessor<DoFHandlerType, level_dof_access>::future_fe_index_set() const
{
  Assert(
    (dynamic_cast<const dealii::DoFHandler<DoFHandlerType::dimension,
                                           DoFHandlerType::space_dimension> *>(
       this->dof_handler) != nullptr) ||
      (this->has_children() == false),
    ExcMessage("You can not ask for the future_fe_index on a cell that has "
               "children because no degrees of freedom are assigned "
               "to this cell and, consequently, no finite element "
               "is associated with it."));
  Assert(
    (dynamic_cast<const dealii::DoFHandler<DoFHandlerType::dimension,
                                           DoFHandlerType::space_dimension> *>(
       this->dof_handler) != nullptr) ||
      (this->is_locally_owned()),
    ExcMessage("You can only query future_fe_index information on cells "
               "that are locally owned."));

  return dealii::internal::DoFCellAccessorImplementation::Implementation::
    future_fe_index_set(*this);
}



template <typename DoFHandlerType, bool level_dof_access>
inline void
DoFCellAccessor<DoFHandlerType, level_dof_access>::clear_future_fe_index() const
{
  Assert(
    (dynamic_cast<const dealii::DoFHandler<DoFHandlerType::dimension,
                                           DoFHandlerType::space_dimension> *>(
       this->dof_handler) != nullptr) ||
      (this->has_children() == false),
    ExcMessage("You can not ask for the future_fe_index on a cell that has "
               "children because no degrees of freedom are assigned "
               "to this cell and, consequently, no finite element "
               "is associated with it."));
  Assert(
    (dynamic_cast<const dealii::DoFHandler<DoFHandlerType::dimension,
                                           DoFHandlerType::space_dimension> *>(
       this->dof_handler) != nullptr) ||
      (this->is_locally_owned()),
    ExcMessage("You can only query future_fe_index information on cells "
               "that are locally owned."));

  dealii::internal::DoFCellAccessorImplementation::Implementation::
    clear_future_fe_index(*this);
}



template <typename DoFHandlerType, bool level_dof_access>
template <typename number, typename OutputVector>
inline void
DoFCellAccessor<DoFHandlerType, level_dof_access>::distribute_local_to_global(
  const Vector<number> &local_source,
  OutputVector &        global_destination) const
{
  this->distribute_local_to_global(local_source.begin(),
                                   local_source.end(),
                                   global_destination);
}



template <typename DoFHandlerType, bool level_dof_access>
template <typename ForwardIterator, typename OutputVector>
inline void
DoFCellAccessor<DoFHandlerType, level_dof_access>::distribute_local_to_global(
  ForwardIterator local_source_begin,
  ForwardIterator local_source_end,
  OutputVector &  global_destination) const
{
  Assert(this->dof_handler != nullptr,
         (typename std::decay<decltype(*this)>::type::ExcInvalidObject()));
  Assert(static_cast<unsigned int>(local_source_end - local_source_begin) ==
           this->get_fe().dofs_per_cell,
         (typename std::decay<decltype(*this)>::type::ExcVectorDoesNotMatch()));
  Assert(this->dof_handler->n_dofs() == global_destination.size(),
         (typename std::decay<decltype(*this)>::type::ExcVectorDoesNotMatch()));

  Assert(!this->has_children(), ExcMessage("Cell must be active"));

  const unsigned int n_dofs = local_source_end - local_source_begin;

  const types::global_dof_index *dofs =
    this->dof_handler->levels[this->level()]->get_cell_cache_start(
      this->present_index, n_dofs);

  // distribute cell vector
  global_destination.add(n_dofs, dofs, local_source_begin);
}



template <typename DoFHandlerType, bool level_dof_access>
template <typename ForwardIterator, typename OutputVector>
inline void
DoFCellAccessor<DoFHandlerType, level_dof_access>::distribute_local_to_global(
  const AffineConstraints<typename OutputVector::value_type> &constraints,
  ForwardIterator local_source_begin,
  ForwardIterator local_source_end,
  OutputVector &  global_destination) const
{
  Assert(this->dof_handler != nullptr,
         (typename std::decay<decltype(*this)>::type::ExcInvalidObject()));
  Assert(local_source_end - local_source_begin == this->get_fe().dofs_per_cell,
         (typename std::decay<decltype(*this)>::type::ExcVectorDoesNotMatch()));
  Assert(this->dof_handler->n_dofs() == global_destination.size(),
         (typename std::decay<decltype(*this)>::type::ExcVectorDoesNotMatch()));

  Assert(!this->has_children(), ExcMessage("Cell must be active."));

  const unsigned int n_dofs = local_source_end - local_source_begin;

  const types::global_dof_index *dofs =
    this->dof_handler->levels[this->level()]->get_cell_cache_start(
      this->present_index, n_dofs);

  // distribute cell vector
  constraints.distribute_local_to_global(local_source_begin,
                                         local_source_end,
                                         dofs,
                                         global_destination);
}



template <typename DoFHandlerType, bool level_dof_access>
template <typename number, typename OutputMatrix>
inline void
DoFCellAccessor<DoFHandlerType, level_dof_access>::distribute_local_to_global(
  const FullMatrix<number> &local_source,
  OutputMatrix &            global_destination) const
{
  Assert(this->dof_handler != nullptr,
         (typename std::decay<decltype(*this)>::type::ExcInvalidObject()));
  Assert(local_source.m() == this->get_fe().dofs_per_cell,
         (typename std::decay<decltype(*this)>::type::ExcMatrixDoesNotMatch()));
  Assert(local_source.n() == this->get_fe().dofs_per_cell,
         (typename std::decay<decltype(*this)>::type::ExcMatrixDoesNotMatch()));
  Assert(this->dof_handler->n_dofs() == global_destination.m(),
         (typename std::decay<decltype(*this)>::type::ExcMatrixDoesNotMatch()));
  Assert(this->dof_handler->n_dofs() == global_destination.n(),
         (typename std::decay<decltype(*this)>::type::ExcMatrixDoesNotMatch()));

  Assert(!this->has_children(), ExcMessage("Cell must be active."));

  const unsigned int n_dofs = local_source.m();

  const types::global_dof_index *dofs =
    this->dof_handler->levels[this->level()]->get_cell_cache_start(
      this->present_index, n_dofs);

  // distribute cell matrix
  for (unsigned int i = 0; i < n_dofs; ++i)
    global_destination.add(dofs[i], n_dofs, dofs, &local_source(i, 0));
}



template <typename DoFHandlerType, bool level_dof_access>
template <typename number, typename OutputMatrix, typename OutputVector>
inline void
DoFCellAccessor<DoFHandlerType, level_dof_access>::distribute_local_to_global(
  const FullMatrix<number> &local_matrix,
  const Vector<number> &    local_vector,
  OutputMatrix &            global_matrix,
  OutputVector &            global_vector) const
{
  Assert(this->dof_handler != nullptr,
         (typename std::decay<decltype(*this)>::type::ExcInvalidObject()));
  Assert(local_matrix.m() == this->get_fe().dofs_per_cell,
         (typename std::decay<decltype(*this)>::type::ExcMatrixDoesNotMatch()));
  Assert(local_matrix.n() == this->get_fe().dofs_per_cell,
         (typename std::decay<decltype(*this)>::type::ExcVectorDoesNotMatch()));
  Assert(this->dof_handler->n_dofs() == global_matrix.m(),
         (typename std::decay<decltype(*this)>::type::ExcMatrixDoesNotMatch()));
  Assert(this->dof_handler->n_dofs() == global_matrix.n(),
         (typename std::decay<decltype(*this)>::type::ExcMatrixDoesNotMatch()));
  Assert(local_vector.size() == this->get_fe().dofs_per_cell,
         (typename std::decay<decltype(*this)>::type::ExcVectorDoesNotMatch()));
  Assert(this->dof_handler->n_dofs() == global_vector.size(),
         (typename std::decay<decltype(*this)>::type::ExcVectorDoesNotMatch()));

  Assert(!this->has_children(), ExcMessage("Cell must be active."));

  const unsigned int             n_dofs = this->get_fe().dofs_per_cell;
  const types::global_dof_index *dofs =
    this->dof_handler->levels[this->level()]->get_cell_cache_start(
      this->present_index, n_dofs);

  // distribute cell matrices
  for (unsigned int i = 0; i < n_dofs; ++i)
    {
      global_matrix.add(dofs[i], n_dofs, dofs, &local_matrix(i, 0));
      global_vector(dofs[i]) += local_vector(i);
    }
}



DEAL_II_NAMESPACE_CLOSE

#endif
