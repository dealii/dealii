// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2014 by the deal.II authors
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

#ifndef __deal2__dof_accessor_templates_h
#define __deal2__dof_accessor_templates_h


#include <deal.II/base/config.h>
#include <deal.II/base/types.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_levels.h>
#include <deal.II/dofs/dof_faces.h>
#include <deal.II/hp/dof_level.h>
#include <deal.II/hp/dof_faces.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_iterator.templates.h>

#include <vector>
#include <limits>

DEAL_II_NAMESPACE_OPEN


/*------------------------- Functions: DoFAccessor ---------------------------*/


template <int structdim, class DH, bool level_dof_access>
inline
DoFAccessor<structdim,DH,level_dof_access>::DoFAccessor ()
{
  Assert (false, ExcInvalidObject());
}



template <int structdim, class DH, bool level_dof_access>
inline
DoFAccessor<structdim,DH,level_dof_access>::DoFAccessor (
  const Triangulation<DH::dimension,DH::space_dimension> *tria,
  const int                 level,
  const int                 index,
  const DH                 *dof_handler)
  :
  dealii::internal::DoFAccessor::Inheritance<structdim,DH::dimension,
  DH::space_dimension>::BaseClass (tria,
                                   level,
                                   index),
  dof_handler(const_cast<DH *>(dof_handler))
{}



template <int structdim, class DH, bool level_dof_access>
template <int structdim2, int dim2, int spacedim2>
DoFAccessor<structdim,DH,level_dof_access>::DoFAccessor (const InvalidAccessor<structdim2,dim2,spacedim2> &)
{
  Assert (false, ExcInvalidObject());
}



template <int structdim, class DH, bool level_dof_access>
template <int dim2, class DH2, bool level_dof_access2>
inline
DoFAccessor<structdim,DH,level_dof_access>::DoFAccessor (const DoFAccessor<dim2, DH2, level_dof_access2> &other)
  : BaseClass(other),
    dof_handler(0)
{
  Assert (false, ExcMessage("You are trying to assign iterators that are incompatible. "
                            "Reasons for incompatibility are that they point to different "
                            "types of DoFHandlers (e.g., dealii::DoFHandler and "
                            "dealii::hp::DoFHandler) or that the refer to objects of "
                            "different dimensionality (e.g., assigning a line iterator "
                            "to a quad iterator)."));
}



template <int structdim, class DH, bool level_dof_access>
template <bool level_dof_access2>
inline
DoFAccessor<structdim,DH,level_dof_access>::DoFAccessor (const DoFAccessor<structdim, DH, level_dof_access2> &other)
  : BaseClass(other),
    dof_handler(const_cast<DH *>(other.dof_handler))
{
}



template <int structdim, class DH, bool level_dof_access>
inline
void
DoFAccessor<structdim,DH,level_dof_access>::set_dof_handler (DH *dh)
{
  Assert (dh != 0, ExcInvalidObject());
  this->dof_handler = dh;
}



template <int structdim, class DH, bool level_dof_access>
inline
const DH &
DoFAccessor<structdim,DH,level_dof_access>::get_dof_handler () const
{
  return *this->dof_handler;
}



template <int structdim, class DH, bool level_dof_access>
inline
void
DoFAccessor<structdim,DH,level_dof_access>::copy_from (
  const TriaAccessorBase<structdim, DH::dimension, DH::space_dimension> &da)
{
  Assert (this->dof_handler != 0, ExcInvalidObject());
  BaseClass::copy_from(da);
}



template <int structdim, class DH, bool level_dof_access>
template <bool level_dof_access2>
inline
void
DoFAccessor<structdim,DH,level_dof_access>::copy_from (const DoFAccessor<structdim,DH,level_dof_access2> &a)
{
  BaseClass::copy_from (a);
  set_dof_handler (a.dof_handler);
}



template <int structdim, class DH, bool level_dof_access>
template <int dim2, class DH2, bool level_dof_access2>
inline
bool
DoFAccessor<structdim,DH,level_dof_access>::operator == (const DoFAccessor<dim2,DH2,level_dof_access2> &a) const
{
  Assert (structdim == dim2, ExcCantCompareIterators());
  Assert (this->dof_handler == a.dof_handler, ExcCantCompareIterators());
  return (BaseClass::operator == (a));
}



template <int structdim, class DH, bool level_dof_access>
template <int dim2, class DH2, bool level_dof_access2>
inline
bool
DoFAccessor<structdim,DH,level_dof_access>::operator != (const DoFAccessor<dim2,DH2,level_dof_access2> &a) const
{
  Assert (structdim == dim2, ExcCantCompareIterators());
  Assert (this->dof_handler == a.dof_handler, ExcCantCompareIterators());
  return (BaseClass::operator != (a));
}



template <int structdim, class DH, bool level_dof_access>
inline
TriaIterator<DoFAccessor<structdim,DH,level_dof_access> >
DoFAccessor<structdim,DH,level_dof_access>::child (const unsigned int i) const
{
  Assert (static_cast<unsigned int>(this->level()) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));

  TriaIterator<TriaAccessor<structdim,DH::dimension,DH::space_dimension> >
  t = TriaAccessor<structdim,DH::dimension,DH::space_dimension>::child(i);

  TriaIterator<DoFAccessor<structdim,DH,level_dof_access> > q (*t, this->dof_handler);
  return q;
}


namespace internal
{
  namespace DoFAccessor
  {
    /**
     * A class like the one with same
     * name in tria.cc. See there for
     * more information.
     */
    struct Implementation
    {
      /**
       * Implementations of the
       * get_dof_index/set_dof_index functions.
       */
      template <int spacedim>
      static
      types::global_dof_index
      get_dof_index (const dealii::DoFHandler<1,spacedim>   &dof_handler,
                     const unsigned int                      obj_level,
                     const unsigned int                      obj_index,
                     const unsigned int                      fe_index,
                     const unsigned int                      local_index,
                     dealii::internal::int2type<1>)
      {
        return dof_handler.levels[obj_level]->dof_object.
               get_dof_index (dof_handler,
                              obj_index,
                              fe_index,
                              local_index);
      }


      template <int spacedim>
      static
      void
      set_dof_index (const dealii::DoFHandler<1,spacedim>   &dof_handler,
                     const unsigned int              obj_level,
                     const unsigned int              obj_index,
                     const unsigned int              fe_index,
                     const unsigned int              local_index,
                     dealii::internal::int2type<1>,
                     const types::global_dof_index              global_index)
      {
        dof_handler.levels[obj_level]->dof_object.
        set_dof_index (dof_handler,
                       obj_index,
                       fe_index,
                       local_index,
                       global_index);
      }


      template <int spacedim>
      static
      types::global_dof_index
      get_dof_index (const dealii::DoFHandler<2,spacedim>   &dof_handler,
                     const unsigned int              obj_level,
                     const unsigned int              obj_index,
                     const unsigned int              fe_index,
                     const unsigned int              local_index,
                     dealii::internal::int2type<1>)
      {
        // faces have no levels
        Assert (obj_level == 0, ExcInternalError());
        return dof_handler.faces->lines.
               get_dof_index (dof_handler,
                              obj_index,
                              fe_index,
                              local_index);
      }


      template <int spacedim>
      static
      void
      set_dof_index (const dealii::DoFHandler<2,spacedim>   &dof_handler,
                     const unsigned int              obj_level,
                     const unsigned int              obj_index,
                     const unsigned int              fe_index,
                     const unsigned int              local_index,
                     dealii::internal::int2type<1>,
                     const types::global_dof_index              global_index)
      {
        // faces have no levels
        Assert (obj_level == 0, ExcInternalError());
        dof_handler.faces->lines.
        set_dof_index (dof_handler,
                       obj_index,
                       fe_index,
                       local_index,
                       global_index);
      }


      template <int spacedim>
      static
      types::global_dof_index
      get_dof_index (const dealii::DoFHandler<2,spacedim>   &dof_handler,
                     const unsigned int              obj_level,
                     const unsigned int              obj_index,
                     const unsigned int              fe_index,
                     const unsigned int              local_index,
                     dealii::internal::int2type<2>)
      {
        return dof_handler.levels[obj_level]->dof_object.
               get_dof_index (dof_handler,
                              obj_index,
                              fe_index,
                              local_index);
      }


      template <int spacedim>
      static
      void
      set_dof_index (const dealii::DoFHandler<2,spacedim>   &dof_handler,
                     const unsigned int              obj_level,
                     const unsigned int              obj_index,
                     const unsigned int              fe_index,
                     const unsigned int              local_index,
                     dealii::internal::int2type<2>,
                     const types::global_dof_index              global_index)
      {
        dof_handler.levels[obj_level]->dof_object.
        set_dof_index (dof_handler,
                       obj_index,
                       fe_index,
                       local_index,
                       global_index);
      }


      template <int spacedim>
      static
      types::global_dof_index
      get_dof_index (const dealii::DoFHandler<3,spacedim>   &dof_handler,
                     const unsigned int              obj_level,
                     const unsigned int              obj_index,
                     const unsigned int              fe_index,
                     const unsigned int              local_index,
                     dealii::internal::int2type<1>)
      {
        // faces have no levels
        Assert (obj_level == 0, ExcInternalError());
        return dof_handler.faces->lines.
               get_dof_index (dof_handler,
                              obj_index,
                              fe_index,
                              local_index);
      }


      template <int spacedim>
      static
      void
      set_dof_index (const dealii::DoFHandler<3,spacedim>   &dof_handler,
                     const unsigned int              obj_level,
                     const unsigned int              obj_index,
                     const unsigned int              fe_index,
                     const unsigned int              local_index,
                     dealii::internal::int2type<1>,
                     const types::global_dof_index              global_index)
      {
        // faces have no levels
        Assert (obj_level == 0, ExcInternalError());
        dof_handler.faces->lines.
        set_dof_index (dof_handler,
                       obj_index,
                       fe_index,
                       local_index,
                       global_index);
      }



      template <int spacedim>
      static
      types::global_dof_index
      get_dof_index (const dealii::DoFHandler<3,spacedim>   &dof_handler,
                     const unsigned int              obj_level,
                     const unsigned int              obj_index,
                     const unsigned int              fe_index,
                     const unsigned int              local_index,
                     dealii::internal::int2type<2>)
      {
        // faces have no levels
        Assert (obj_level == 0, ExcInternalError());
        return dof_handler.faces->quads.
               get_dof_index (dof_handler,
                              obj_index,
                              fe_index,
                              local_index);
      }


      template <int spacedim>
      static
      void
      set_dof_index (const dealii::DoFHandler<3,spacedim>   &dof_handler,
                     const unsigned int              obj_level,
                     const unsigned int              obj_index,
                     const unsigned int              fe_index,
                     const unsigned int              local_index,
                     dealii::internal::int2type<2>,
                     const types::global_dof_index              global_index)
      {
        // faces have no levels
        Assert (obj_level == 0, ExcInternalError());
        dof_handler.faces->quads.
        set_dof_index (dof_handler,
                       obj_index,
                       fe_index,
                       local_index,
                       global_index);
      }



      template <int spacedim>
      static
      types::global_dof_index
      get_dof_index (const dealii::DoFHandler<3,spacedim>   &dof_handler,
                     const unsigned int              obj_level,
                     const unsigned int              obj_index,
                     const unsigned int              fe_index,
                     const unsigned int              local_index,
                     dealii::internal::int2type<3>)
      {
        return dof_handler.levels[obj_level]->dof_object.
               get_dof_index (dof_handler,
                              obj_index,
                              fe_index,
                              local_index);
      }


      template <int spacedim>
      static
      void
      set_dof_index (const dealii::DoFHandler<3,spacedim>   &dof_handler,
                     const unsigned int              obj_level,
                     const unsigned int              obj_index,
                     const unsigned int              fe_index,
                     const unsigned int              local_index,
                     dealii::internal::int2type<3>,
                     const types::global_dof_index              global_index)
      {
        dof_handler.levels[obj_level]->dof_object.
        set_dof_index (dof_handler,
                       obj_index,
                       fe_index,
                       local_index,
                       global_index);
      }


      template <int spacedim>
      static
      types::global_dof_index
      get_dof_index (const dealii::hp::DoFHandler<1,spacedim> &dof_handler,
                     const unsigned int       obj_level,
                     const unsigned int       obj_index,
                     const unsigned int       fe_index,
                     const unsigned int       local_index,
                     const dealii::internal::int2type<1> &)
      {
        return dof_handler.levels[obj_level]->
               get_dof_index (obj_index,
                              fe_index,
                              local_index);
      }


      template <int spacedim>
      static
      void
      set_dof_index (const dealii::hp::DoFHandler<1,spacedim> &dof_handler,
                     const unsigned int       obj_level,
                     const unsigned int       obj_index,
                     const unsigned int       fe_index,
                     const unsigned int       local_index,
                     const dealii::internal::int2type<1> &,
                     const types::global_dof_index       global_index)
      {
        dof_handler.levels[obj_level]->
        set_dof_index (obj_index,
                       fe_index,
                       local_index,
                       global_index);
      }


      template <int spacedim>
      static
      types::global_dof_index
      get_dof_index (const dealii::hp::DoFHandler<2,spacedim> &dof_handler,
                     const unsigned int       obj_level,
                     const unsigned int       obj_index,
                     const unsigned int       fe_index,
                     const unsigned int       local_index,
                     const dealii::internal::int2type<1> &)
      {
        return dof_handler.faces->lines.
               get_dof_index (dof_handler,
                              obj_index,
                              fe_index,
                              local_index,
                              obj_level);
      }


      template <int spacedim>
      static
      void
      set_dof_index (const dealii::hp::DoFHandler<2,spacedim> &dof_handler,
                     const unsigned int       obj_level,
                     const unsigned int       obj_index,
                     const unsigned int       fe_index,
                     const unsigned int       local_index,
                     const dealii::internal::int2type<1> &,
                     const types::global_dof_index       global_index)
      {
        dof_handler.faces->lines.
        set_dof_index (dof_handler,
                       obj_index,
                       fe_index,
                       local_index,
                       global_index,
                       obj_level);
      }


      template <int spacedim>
      static
      types::global_dof_index
      get_dof_index (const dealii::hp::DoFHandler<2,spacedim> &dof_handler,
                     const unsigned int       obj_level,
                     const unsigned int       obj_index,
                     const unsigned int       fe_index,
                     const unsigned int       local_index,
                     const dealii::internal::int2type<2> &)
      {
        return dof_handler.levels[obj_level]->
               get_dof_index (obj_index,
                              fe_index,
                              local_index);
      }


      template <int spacedim>
      static
      void
      set_dof_index (const dealii::hp::DoFHandler<2,spacedim> &dof_handler,
                     const unsigned int       obj_level,
                     const unsigned int       obj_index,
                     const unsigned int       fe_index,
                     const unsigned int       local_index,
                     const dealii::internal::int2type<2> &,
                     const types::global_dof_index       global_index)
      {
        dof_handler.levels[obj_level]->
        set_dof_index (obj_index,
                       fe_index,
                       local_index,
                       global_index);
      }


      template <int spacedim>
      static
      types::global_dof_index
      get_dof_index (const dealii::hp::DoFHandler<3,spacedim> &dof_handler,
                     const unsigned int       obj_level,
                     const unsigned int       obj_index,
                     const unsigned int       fe_index,
                     const unsigned int       local_index,
                     const dealii::internal::int2type<1> &)
      {
        return dof_handler.faces->lines.
               get_dof_index (dof_handler,
                              obj_index,
                              fe_index,
                              local_index,
                              obj_level);
      }


      template <int spacedim>
      static
      void
      set_dof_index (const dealii::hp::DoFHandler<3,spacedim> &dof_handler,
                     const unsigned int       obj_level,
                     const unsigned int       obj_index,
                     const unsigned int       fe_index,
                     const unsigned int       local_index,
                     const dealii::internal::int2type<1> &,
                     const types::global_dof_index       global_index)
      {
        dof_handler.faces->lines.
        set_dof_index (dof_handler,
                       obj_index,
                       fe_index,
                       local_index,
                       global_index,
                       obj_level);
      }


      template <int spacedim>
      static
      types::global_dof_index
      get_dof_index (const dealii::hp::DoFHandler<3,spacedim> &dof_handler,
                     const unsigned int       obj_level,
                     const unsigned int       obj_index,
                     const unsigned int       fe_index,
                     const unsigned int       local_index,
                     const dealii::internal::int2type<2> &)
      {
        return dof_handler.faces->quads.
               get_dof_index (dof_handler,
                              obj_index,
                              fe_index,
                              local_index,
                              obj_level);
      }


      template <int spacedim>
      static
      void
      set_dof_index (const dealii::hp::DoFHandler<3,spacedim> &dof_handler,
                     const unsigned int       obj_level,
                     const unsigned int       obj_index,
                     const unsigned int       fe_index,
                     const unsigned int       local_index,
                     const dealii::internal::int2type<2> &,
                     const types::global_dof_index       global_index)
      {
        dof_handler.faces->quads.
        set_dof_index (dof_handler,
                       obj_index,
                       fe_index,
                       local_index,
                       global_index,
                       obj_level);
      }


      template <int spacedim>
      static
      types::global_dof_index
      get_dof_index (const dealii::hp::DoFHandler<3,spacedim> &dof_handler,
                     const unsigned int       obj_level,
                     const unsigned int       obj_index,
                     const unsigned int       fe_index,
                     const unsigned int       local_index,
                     const dealii::internal::int2type<3> &)
      {
        return dof_handler.levels[obj_level]->
               get_dof_index (obj_index,
                              fe_index,
                              local_index);
      }


      template <int spacedim>
      static
      void
      set_dof_index (const dealii::hp::DoFHandler<3,spacedim> &dof_handler,
                     const unsigned int       obj_level,
                     const unsigned int       obj_index,
                     const unsigned int       fe_index,
                     const unsigned int       local_index,
                     const dealii::internal::int2type<3> &,
                     const types::global_dof_index       global_index)
      {
        dof_handler.levels[obj_level]->
        set_dof_index (obj_index,
                       fe_index,
                       local_index,
                       global_index);
      }


      template <int structdim, int dim, int spacedim>
      static
      bool
      fe_index_is_active (const dealii::DoFHandler<dim,spacedim> &,
                          const unsigned int,
                          const unsigned int,
                          const unsigned int fe_index,
                          const dealii::internal::int2type<structdim> &)
      {
        return (fe_index == 0);
      }



      template <int structdim, int dim, int spacedim>
      static
      unsigned int
      n_active_fe_indices (const dealii::DoFHandler<dim,spacedim> &dof_handler,
                           const unsigned int obj_level,
                           const unsigned int obj_index,
                           const dealii::internal::int2type<structdim> &)
      {
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
        Assert ((dim==structdim
                 ?
                 typename
                 internal::Triangulation::Iterators<dim,spacedim>::
                 raw_cell_iterator (&dof_handler.get_tria(),
                                    obj_level,
                                    obj_index)->used()
                 :
                 (structdim==1
                  ?
                  typename
                  internal::Triangulation::Iterators<dim,spacedim>::
                  raw_line_iterator (&dof_handler.get_tria(),
                                     obj_level,
                                     obj_index)->used()
                  :
                  true))
                == true,
                ExcMessage ("This cell is not active and therefore can't be "
                            "queried for its active FE indices"));
        return 1;
      }



      template <int structdim, int dim, int spacedim>
      static
      unsigned int
      nth_active_fe_index (const dealii::DoFHandler<dim,spacedim> &dof_handler,
                           const unsigned int obj_level,
                           const unsigned int obj_index,
                           const unsigned int n,
                           const dealii::internal::int2type<structdim> &)
      {
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
        Assert ((dim==structdim
                 ?
                 typename
                 internal::Triangulation::Iterators<dim,spacedim>::
                 raw_cell_iterator (&dof_handler.get_tria(),
                                    obj_level,
                                    obj_index)->used()
                 :
                 (structdim==1
                  ?
                  typename
                  internal::Triangulation::Iterators<dim,spacedim>::
                  raw_line_iterator (&dof_handler.get_tria(),
                                     obj_level,
                                     obj_index)->used()
                  :
                  true))
                == true,
                ExcMessage ("This cell is not active and therefore can't be "
                            "queried for its active FE indices"));
        Assert (n == 0, ExcIndexRange (n, 0, 1));

        return dealii::DoFHandler<dim,spacedim>::default_fe_index;
      }


      template <int spacedim>
      static
      bool
      fe_index_is_active (const dealii::hp::DoFHandler<1,spacedim> &dof_handler,
                          const unsigned int obj_level,
                          const unsigned int obj_index,
                          const unsigned int fe_index,
                          const dealii::internal::int2type<1> &)
      {
        return dof_handler.levels[obj_level]->fe_index_is_active(obj_index,
                                                                 fe_index);
      }


      template <int spacedim>
      static
      unsigned int
      n_active_fe_indices (const dealii::hp::DoFHandler<1,spacedim> &,
                           const unsigned int obj_level,
                           const unsigned int obj_index,
                           const dealii::internal::int2type<1> &)
      {
        // on a cell, the number of active elements is one
        return 1;
      }



      template <int spacedim>
      static
      unsigned int
      nth_active_fe_index (const dealii::hp::DoFHandler<1,spacedim> &dof_handler,
                           const unsigned int obj_level,
                           const unsigned int obj_index,
                           const unsigned int n,
                           const dealii::internal::int2type<1> &)
      {
        Assert (n==0, ExcMessage("On cells, there can only be one active FE index"));
        return dof_handler.levels[obj_level]->active_fe_index (obj_index);
      }


      template <int spacedim>
      static
      bool
      fe_index_is_active (const dealii::hp::DoFHandler<2,spacedim> &dof_handler,
                          const unsigned int obj_level,
                          const unsigned int obj_index,
                          const unsigned int fe_index,
                          const dealii::internal::int2type<1> &)
      {
        return dof_handler.faces->lines.fe_index_is_active(dof_handler,
                                                           obj_index,
                                                           fe_index,
                                                           obj_level);
      }


      template <int spacedim>
      static
      unsigned int
      n_active_fe_indices (const dealii::hp::DoFHandler<2,spacedim> &dof_handler,
                           const unsigned int ,
                           const unsigned int obj_index,
                           const dealii::internal::int2type<1> &)
      {
        return dof_handler.faces->lines.n_active_fe_indices (dof_handler,
                                                             obj_index);
      }


      template <int spacedim>
      static
      unsigned int
      nth_active_fe_index (const dealii::hp::DoFHandler<2,spacedim> &dof_handler,
                           const unsigned int obj_level,
                           const unsigned int obj_index,
                           const unsigned int n,
                           const dealii::internal::int2type<1> &)
      {
        return dof_handler.faces->lines.nth_active_fe_index (dof_handler,
                                                             obj_level,
                                                             obj_index,
                                                             n);
      }



      template <int spacedim>
      static
      bool
      fe_index_is_active (const dealii::hp::DoFHandler<2,spacedim> &dof_handler,
                          const unsigned int obj_level,
                          const unsigned int obj_index,
                          const unsigned int fe_index,
                          const dealii::internal::int2type<2> &)
      {
        return dof_handler.levels[obj_level]->fe_index_is_active(obj_index,
                                                                 fe_index);
      }


      template <int spacedim>
      static
      unsigned int
      n_active_fe_indices (const dealii::hp::DoFHandler<2,spacedim> &,
                           const unsigned int obj_level,
                           const unsigned int obj_index,
                           const dealii::internal::int2type<2> &)
      {
        // on a cell, the number of active elements is one
        return 1;
      }



      template <int spacedim>
      static
      unsigned int
      nth_active_fe_index (const dealii::hp::DoFHandler<2,spacedim> &dof_handler,
                           const unsigned int obj_level,
                           const unsigned int obj_index,
                           const unsigned int n,
                           const dealii::internal::int2type<2> &)
      {
        Assert (n==0, ExcMessage("On cells, there can only be one active FE index"));
        return dof_handler.levels[obj_level]->active_fe_index (obj_index);
      }



      template <int spacedim>
      static
      bool
      fe_index_is_active (const dealii::hp::DoFHandler<3,spacedim> &dof_handler,
                          const unsigned int obj_level,
                          const unsigned int obj_index,
                          const unsigned int fe_index,
                          const dealii::internal::int2type<1> &)
      {
        return dof_handler.faces->lines.fe_index_is_active(dof_handler,
                                                           obj_index,
                                                           fe_index,
                                                           obj_level);
      }


      template <int spacedim>
      static
      unsigned int
      n_active_fe_indices (const dealii::hp::DoFHandler<3,spacedim> &dof_handler,
                           const unsigned int,
                           const unsigned int obj_index,
                           const dealii::internal::int2type<1> &)
      {
        return dof_handler.faces->lines.n_active_fe_indices (dof_handler,
                                                             obj_index);
      }



      template <int spacedim>
      static
      unsigned int
      nth_active_fe_index (const dealii::hp::DoFHandler<3,spacedim> &dof_handler,
                           const unsigned int obj_level,
                           const unsigned int obj_index,
                           const unsigned int n,
                           const dealii::internal::int2type<1> &)
      {
        return dof_handler.faces->lines.nth_active_fe_index (dof_handler,
                                                             obj_level,
                                                             obj_index,
                                                             n);
      }



      template <int spacedim>
      static
      bool
      fe_index_is_active (const dealii::hp::DoFHandler<3,spacedim> &dof_handler,
                          const unsigned int obj_level,
                          const unsigned int obj_index,
                          const unsigned int fe_index,
                          const dealii::internal::int2type<2> &)
      {
        return dof_handler.faces->quads.fe_index_is_active(dof_handler,
                                                           obj_index,
                                                           fe_index,
                                                           obj_level);
      }

      template <int spacedim>
      static
      bool
      fe_index_is_active (const dealii::hp::DoFHandler<3,spacedim> &dof_handler,
                          const unsigned int obj_level,
                          const unsigned int obj_index,
                          const unsigned int fe_index,
                          const dealii::internal::int2type<3> &)
      {
        return dof_handler.levels[obj_level]->fe_index_is_active(obj_index,
                                                                 fe_index);
      }


      template <int spacedim>
      static
      unsigned int
      n_active_fe_indices (const dealii::hp::DoFHandler<3,spacedim> &dof_handler,
                           const unsigned int ,
                           const unsigned int obj_index,
                           const dealii::internal::int2type<2> &)
      {
        return dof_handler.faces->quads.n_active_fe_indices (dof_handler,
                                                             obj_index);
      }



      template <int spacedim>
      static
      unsigned int
      nth_active_fe_index (const dealii::hp::DoFHandler<3,spacedim> &dof_handler,
                           const unsigned int obj_level,
                           const unsigned int obj_index,
                           const unsigned int n,
                           const dealii::internal::int2type<2> &)
      {
        return dof_handler.faces->quads.nth_active_fe_index (dof_handler,
                                                             obj_level,
                                                             obj_index,
                                                             n);
      }



      template <int spacedim>
      static
      unsigned int
      n_active_fe_indices (const dealii::hp::DoFHandler<3,spacedim> &,
                           const unsigned int obj_level,
                           const unsigned int obj_index,
                           const dealii::internal::int2type<3> &)
      {
        // on a cell, the number of active elements is one
        return 1;
      }



      template <int spacedim>
      static
      unsigned int
      nth_active_fe_index (const dealii::hp::DoFHandler<3,spacedim> &dof_handler,
                           const unsigned int obj_level,
                           const unsigned int obj_index,
                           const unsigned int n,
                           const dealii::internal::int2type<3> &)
      {
        Assert (n==0, ExcMessage("On cells, there can only be one active FE index"));
        return dof_handler.levels[obj_level]->active_fe_index (obj_index);
      }

      /**
       * Set the @p local_index-th
       * degree of freedom
       * corresponding to the finite
       * element specified by @p
       * fe_index on the vertex with
       * global number @p
       * vertex_index to @p
       * global_index.
       */
      template <int dim, int spacedim>
      static
      void
      set_vertex_dof_index (dealii::DoFHandler<dim,spacedim> &dof_handler,
                            const unsigned int vertex_index,
                            const unsigned int fe_index,
                            const unsigned int local_index,
                            const types::global_dof_index global_index)
      {
        Assert ((fe_index == dealii::DoFHandler<dim,spacedim>::default_fe_index),
                ExcMessage ("Only the default FE index is allowed for non-hp DoFHandler objects"));
        Assert (dof_handler.selected_fe != 0,
                ExcMessage ("No finite element collection is associated with "
                            "this DoFHandler"));
        Assert (local_index < dof_handler.selected_fe->dofs_per_vertex,
                ExcIndexRange(local_index, 0,
                              dof_handler.selected_fe->dofs_per_vertex));

        dof_handler.vertex_dofs[vertex_index *
                                dof_handler.selected_fe->dofs_per_vertex
                                + local_index]
          = global_index;
      }


      template <int dim, int spacedim>
      static
      void
      set_vertex_dof_index (dealii::hp::DoFHandler<dim,spacedim> &dof_handler,
                            const unsigned int vertex_index,
                            const unsigned int fe_index,
                            const unsigned int local_index,
                            const types::global_dof_index global_index)
      {
        Assert ( (fe_index != dealii::hp::DoFHandler<dim,spacedim>::default_fe_index),
                 ExcMessage ("You need to specify a FE index when working "
                             "with hp DoFHandlers"));
        Assert (dof_handler.finite_elements != 0,
                ExcMessage ("No finite element collection is associated with "
                            "this DoFHandler"));
        Assert (local_index < (*dof_handler.finite_elements)[fe_index].dofs_per_vertex,
                ExcIndexRange(local_index, 0,
                              (*dof_handler.finite_elements)[fe_index].dofs_per_vertex));
        Assert (fe_index < dof_handler.finite_elements->size(),
                ExcInternalError());
        Assert (dof_handler.vertex_dofs_offsets[vertex_index] !=
                numbers::invalid_dof_index,
                ExcMessage ("This vertex is unused and has no DoFs associated with it"));

        // hop along the list of index
        // sets until we find the one
        // with the correct fe_index, and
        // then poke into that
        // part. trigger an exception if
        // we can't find a set for this
        // particular fe_index
        const types::global_dof_index starting_offset = dof_handler.vertex_dofs_offsets[vertex_index];
        types::global_dof_index *pointer = &dof_handler.vertex_dofs[starting_offset];
        while (true)
          {
            Assert (pointer <= &dof_handler.vertex_dofs.back(), ExcInternalError());

            // a fe index is always small
            Assert((*pointer)<std::numeric_limits<unsigned int>::max(), ExcInternalError());
            const types::global_dof_index this_fe_index = *pointer;

            Assert (this_fe_index != numbers::invalid_dof_index,
                    ExcInternalError());
            Assert (this_fe_index < dof_handler.finite_elements->size(),
                    ExcInternalError());

            if (this_fe_index == fe_index)
              {
                *(pointer + 1 + local_index) = global_index;
                return;
              }
            else
              pointer += static_cast<types::global_dof_index>(
                           (*dof_handler.finite_elements)[this_fe_index].dofs_per_vertex + 1);
          }
      }


      /**
       * Get the @p local_index-th
       * degree of freedom
       * corresponding to the finite
       * element specified by @p
       * fe_index on the vertex with
       * global number @p
       * vertex_index to @p
       * global_index.
       */

      template <int dim, int spacedim>
      static
      types::global_dof_index
      get_vertex_dof_index (const dealii::DoFHandler<dim,spacedim> &dof_handler,
                            const unsigned int vertex_index,
                            const unsigned int fe_index,
                            const unsigned int local_index)
      {
        Assert ((fe_index == dealii::DoFHandler<dim,spacedim>::default_fe_index),
                ExcMessage ("Only the default FE index is allowed for non-hp DoFHandler objects"));
        Assert (dof_handler.selected_fe != 0,
                ExcMessage ("No finite element collection is associated with "
                            "this DoFHandler"));
        Assert (local_index < dof_handler.selected_fe->dofs_per_vertex,
                ExcIndexRange(local_index, 0,
                              dof_handler.selected_fe->dofs_per_vertex));

        return
          dof_handler.vertex_dofs[vertex_index *
                                  dof_handler.selected_fe->dofs_per_vertex
                                  + local_index];
      }


      template<int dim, int spacedim>
      static
      types::global_dof_index
      get_vertex_dof_index (const dealii::hp::DoFHandler<dim,spacedim> &dof_handler,
                            const unsigned int vertex_index,
                            const unsigned int fe_index,
                            const unsigned int local_index)
      {
        Assert ( (fe_index != dealii::hp::DoFHandler<dim,spacedim>::default_fe_index),
                 ExcMessage ("You need to specify a FE index when working "
                             "with hp DoFHandlers"));
        Assert (dof_handler.finite_elements != 0,
                ExcMessage ("No finite element collection is associated with "
                            "this DoFHandler"));
        Assert (local_index < (*dof_handler.finite_elements)[fe_index].dofs_per_vertex,
                ExcIndexRange(local_index, 0,
                              (*dof_handler.finite_elements)[fe_index].dofs_per_vertex));
        Assert (vertex_index < dof_handler.vertex_dofs_offsets.size(),
                ExcIndexRange (vertex_index, 0,
                               dof_handler.vertex_dofs_offsets.size()));
        Assert (dof_handler.vertex_dofs_offsets[vertex_index] !=
                numbers::invalid_dof_index,
                ExcMessage ("This vertex is unused and has no DoFs associated with it"));

        // hop along the list of index
        // sets until we find the one
        // with the correct fe_index, and
        // then poke into that
        // part. trigger an exception if
        // we can't find a set for this
        // particular fe_index
        const types::global_dof_index starting_offset = dof_handler.vertex_dofs_offsets[vertex_index];
        const types::global_dof_index *pointer = &dof_handler.vertex_dofs[starting_offset];
        while (true)
          {
            Assert (pointer <= &dof_handler.vertex_dofs.back(), ExcInternalError());

            Assert((*pointer)<std::numeric_limits<types::global_dof_index>::max(), ExcInternalError());
            const types::global_dof_index this_fe_index = *pointer;

            Assert (this_fe_index != numbers::invalid_dof_index,
                    ExcInternalError());
            Assert (this_fe_index < dof_handler.finite_elements->size(),
                    ExcInternalError());

            if (this_fe_index == fe_index)
              return *(pointer + 1 + local_index);
            else
              pointer += static_cast<types::global_dof_index>(
                           (*dof_handler.finite_elements)[this_fe_index].dofs_per_vertex + 1);
          }
      }


      /**
       * Return the number of
       * different finite elements
       * that are active on a given
       * vertex.
       */
      template<int dim, int spacedim>
      static
      unsigned int
      n_active_vertex_fe_indices (const dealii::hp::DoFHandler<dim,spacedim> &dof_handler,
                                  const unsigned int vertex_index)
      {
        Assert (dof_handler.finite_elements != 0,
                ExcMessage ("No finite element collection is associated with "
                            "this DoFHandler"));

        // if this vertex is unused, return 0
        if (dof_handler.vertex_dofs_offsets[vertex_index] == numbers::invalid_dof_index)
          return 0;

        // hop along the list of index
        // sets and count the number of
        // hops
        const types::global_dof_index starting_offset = dof_handler.vertex_dofs_offsets[vertex_index];
        const types::global_dof_index *pointer = &dof_handler.vertex_dofs[starting_offset];

        Assert (*pointer != numbers::invalid_dof_index,
                ExcInternalError());

        unsigned int counter = 0;
        while (true)
          {
            Assert (pointer <= &dof_handler.vertex_dofs.back(), ExcInternalError());

            const types::global_dof_index this_fe_index = *pointer;

            if (this_fe_index == numbers::invalid_dof_index)
              return counter;
            else
              {
                pointer += static_cast<types::global_dof_index>(
                             (*dof_handler.finite_elements)[this_fe_index].dofs_per_vertex + 1);
                ++counter;
              }
          }
      }



      /**
       * Return the fe index of the
       * n-th finite element active
       * on a given vertex.
       */
      template<int dim, int spacedim>
      static
      unsigned int
      nth_active_vertex_fe_index (const dealii::hp::DoFHandler<dim,spacedim> &dof_handler,
                                  const unsigned int vertex_index,
                                  const unsigned int n)
      {
        Assert (dof_handler.finite_elements != 0,
                ExcMessage ("No finite element collection is associated with "
                            "this DoFHandler"));
        Assert (n < n_active_vertex_fe_indices(dof_handler, vertex_index),
                ExcIndexRange (n, 0, n_active_vertex_fe_indices(dof_handler,
                                                                vertex_index)));
        // make sure we don't ask on
        // unused vertices
        Assert (dof_handler.vertex_dofs_offsets[vertex_index] !=
                numbers::invalid_dof_index,
                ExcInternalError());

        // hop along the list of index
        // sets and count the number of
        // hops
        const types::global_dof_index starting_offset = dof_handler.vertex_dofs_offsets[vertex_index];
        const types::global_dof_index *pointer = &dof_handler.vertex_dofs[starting_offset];

        Assert (*pointer != numbers::invalid_dof_index,
                ExcInternalError());

        unsigned int counter = 0;
        while (true)
          {
            Assert (pointer <= &dof_handler.vertex_dofs.back(), ExcInternalError());

            Assert((*pointer)<std::numeric_limits<unsigned int>::max(), ExcInternalError());
            const types::global_dof_index this_fe_index = *pointer;

            Assert (this_fe_index < dof_handler.finite_elements->size(),
                    ExcInternalError());

            if (counter == n)
              return this_fe_index;

            Assert (this_fe_index != numbers::invalid_dof_index,
                    ExcInternalError());

            pointer += static_cast<types::global_dof_index>(
                         (*dof_handler.finite_elements)[this_fe_index].dofs_per_vertex + 1);
            ++counter;
          }
      }



      /**
       * Return whether a particular
       * finite element index is
       * active on the specified
       * vertex.
       */
      template<int dim, int spacedim>
      static
      bool
      fe_is_active_on_vertex (const dealii::hp::DoFHandler<dim,spacedim> &dof_handler,
                              const unsigned int vertex_index,
                              const unsigned int fe_index)
      {
        Assert ( (fe_index != dealii::hp::DoFHandler<dim,spacedim>::default_fe_index),
                 ExcMessage ("You need to specify a FE index when working "
                             "with hp DoFHandlers"));
        Assert (dof_handler.finite_elements != 0,
                ExcMessage ("No finite element collection is associated with "
                            "this DoFHandler"));
        Assert (fe_index < dof_handler.finite_elements->size(),
                ExcInternalError());

        // make sure we don't ask on
        // unused vertices
        Assert (dof_handler.vertex_dofs_offsets[vertex_index] !=
                numbers::invalid_dof_index,
                ExcInternalError());

        // hop along the list of index
        // sets and see whether we find
        // the given index
        const types::global_dof_index starting_offset = dof_handler.vertex_dofs_offsets[vertex_index];
        const types::global_dof_index *pointer = &dof_handler.vertex_dofs[starting_offset];

        Assert (*pointer != numbers::invalid_dof_index,
                ExcInternalError());

        while (true)
          {
            Assert (pointer <= &dof_handler.vertex_dofs.back(), ExcInternalError());

            Assert((*pointer)<std::numeric_limits<types::global_dof_index>::max(), ExcInternalError());
            const types::global_dof_index this_fe_index = *pointer;

            Assert (this_fe_index < dof_handler.finite_elements->size(),
                    ExcInternalError());

            if (this_fe_index == numbers::invalid_dof_index)
              return false;
            else if (this_fe_index == fe_index)
              return true;
            else
              pointer += (*dof_handler.finite_elements)[this_fe_index].dofs_per_vertex + 1;
          }
      }

      template<class DH, bool level_dof_access>
      static
      void set_mg_dof_indices (const dealii::DoFAccessor<1,DH,level_dof_access> &,
                               const int,
                               const std::vector<types::global_dof_index> &,
                               const unsigned int)
      {
        AssertThrow (false, ExcNotImplemented ()); //TODO[TH]: implement
      }



      template<class DH, bool level_dof_access>
      static
      void set_mg_dof_indices (dealii::DoFAccessor<2, DH,level_dof_access> &accessor,
                               const int level,
                               const std::vector<types::global_dof_index> &dof_indices,
                               const unsigned int fe_index)
      {
        const FiniteElement<DH::dimension, DH::space_dimension> &fe = accessor.get_dof_handler ().get_fe ()[fe_index];
        std::vector<types::global_dof_index>::const_iterator next = dof_indices.begin ();

        for (unsigned int vertex = 0; vertex < GeometryInfo<2>::vertices_per_cell; ++vertex)
          for (unsigned int dof = 0; dof < fe.dofs_per_vertex; ++dof)
            accessor.set_mg_vertex_dof_index(level, vertex, dof, *next++, fe_index);

        for (unsigned int line = 0; line < GeometryInfo<2>::lines_per_cell; ++line)
          for (unsigned int dof = 0; dof < fe.dofs_per_line; ++dof)
            accessor.line(line)->set_mg_dof_index(level, dof, *next++);

        for (unsigned int dof = 0; dof < fe.dofs_per_quad; ++dof)
          accessor.set_mg_dof_index(level, dof, *next++);

        Assert (next == dof_indices.end (), ExcInternalError ());
      }



      template<class DH, bool level_dof_access>
      static
      void set_mg_dof_indices (const dealii::DoFAccessor<3, DH,level_dof_access> &,
                               const int,
                               const std::vector<types::global_dof_index> &,
                               const unsigned int)
      {
        AssertThrow (false, ExcNotImplemented ()); //TODO[TH]: implement
      }

    };
  }
}



template <int dim, class DH, bool level_dof_access>
inline
types::global_dof_index
DoFAccessor<dim,DH,level_dof_access>::dof_index (const unsigned int i,
                                                 const unsigned int fe_index) const
{
  // access the respective DoF
  return dealii::internal::DoFAccessor::Implementation::get_dof_index (*this->dof_handler,
         this->level(),
         this->present_index,
         fe_index,
         i,
         dealii::internal::int2type<dim>());
}


template<int structdim, class DH, bool level_dof_access>
inline
types::global_dof_index
DoFAccessor<structdim, DH,level_dof_access>::mg_dof_index (const int level,
                                                           const unsigned int i) const
{
  return this->dof_handler->template get_dof_index<structdim> (level, this->present_index, 0, i);
}


template <int dim, class DH, bool level_dof_access>
inline
void
DoFAccessor<dim,DH,level_dof_access>::set_dof_index (const unsigned int i,
                                                     const types::global_dof_index index,
                                                     const unsigned int fe_index) const
{
  // access the respective DoF
  dealii::internal::DoFAccessor::Implementation::set_dof_index (*this->dof_handler,
      this->level(),
      this->present_index,
      fe_index,
      i,
      dealii::internal::int2type<dim>(),
      index);
}



template <int dim, class DH, bool level_dof_access>
inline
unsigned int
DoFAccessor<dim,DH,level_dof_access>::n_active_fe_indices () const
{
  // access the respective DoF
  return
    dealii::internal::DoFAccessor::Implementation::
    n_active_fe_indices (*this->dof_handler,
                         this->level(),
                         this->present_index,
                         dealii::internal::int2type<dim>());
}



template <int dim, class DH, bool level_dof_access>
inline
unsigned int
DoFAccessor<dim,DH,level_dof_access>::nth_active_fe_index (const unsigned int n) const
{
  // access the respective DoF
  return
    dealii::internal::DoFAccessor::Implementation::
    nth_active_fe_index (*this->dof_handler,
                         this->level(),
                         this->present_index,
                         n,
                         dealii::internal::int2type<dim>());
}



template <int dim, class DH, bool level_dof_access>
inline
bool
DoFAccessor<dim,DH,level_dof_access>::fe_index_is_active (const unsigned int fe_index) const
{
  // access the respective DoF
  return
    dealii::internal::DoFAccessor::Implementation::
    fe_index_is_active (*this->dof_handler,
                        this->level(),
                        this->present_index,
                        fe_index,
                        dealii::internal::int2type<dim>());
}



template <int structdim, class DH, bool level_dof_access>
inline
types::global_dof_index
DoFAccessor<structdim, DH,level_dof_access>::vertex_dof_index (const unsigned int vertex,
    const unsigned int i,
    const unsigned int fe_index) const
{
  return
    dealii::internal::DoFAccessor::Implementation::get_vertex_dof_index
    (*this->dof_handler,
     this->vertex_index(vertex),
     fe_index,
     i);
}


template<int structdim, class DH, bool level_dof_access>
inline
types::global_dof_index
DoFAccessor<structdim, DH,level_dof_access>::mg_vertex_dof_index (const int level,
    const unsigned int vertex,
    const unsigned int i,
    const unsigned int fe_index) const
{
  Assert (this->dof_handler != 0, ExcInvalidObject ());
  Assert (&this->dof_handler->get_fe () != 0, ExcInvalidObject ());
  Assert (vertex < GeometryInfo<structdim>::vertices_per_cell, ExcIndexRange (vertex, 0, GeometryInfo<structdim>::vertices_per_cell));
  Assert (i < this->dof_handler->get_fe ()[fe_index].dofs_per_vertex, ExcIndexRange (i, 0, this->dof_handler->get_fe ()[fe_index].dofs_per_vertex));
  return this->dof_handler->mg_vertex_dofs[this->vertex_index (vertex)].get_index (level, i);
}


template <int structdim, class DH, bool level_dof_access>
inline
void
DoFAccessor<structdim, DH,level_dof_access>::set_vertex_dof_index (const unsigned int vertex,
    const unsigned int i,
    const types::global_dof_index index,
    const unsigned int fe_index) const
{
  dealii::internal::DoFAccessor::Implementation::set_vertex_dof_index
  (*this->dof_handler,
   this->vertex_index(vertex),
   fe_index,
   i,
   index);
}


template<int structdim, class DH, bool level_dof_access>
inline
void
DoFAccessor<structdim, DH,level_dof_access>::set_mg_vertex_dof_index (const int level,
    const unsigned int vertex,
    const unsigned int i,
    const types::global_dof_index index,
    const unsigned int fe_index) const
{
  Assert (this->dof_handler != 0, ExcInvalidObject ());
  Assert (&this->dof_handler->get_fe () != 0, ExcInvalidObject ());
  Assert (vertex < GeometryInfo<structdim>::vertices_per_cell, ExcIndexRange (vertex, 0, GeometryInfo<structdim>::vertices_per_cell));
  Assert (i < this->dof_handler->get_fe ()[fe_index].dofs_per_vertex, ExcIndexRange (i, 0, this->dof_handler->get_fe ()[fe_index].dofs_per_vertex));
  this->dof_handler->mg_vertex_dofs[this->vertex_index (vertex)].set_index (level, i, index);
}


template<int structdim, class DH, bool level_dof_access>
inline
void
DoFAccessor<structdim, DH,level_dof_access>::set_mg_dof_index (const int level,
    const unsigned int i,
    const types::global_dof_index index) const
{
  this->dof_handler->template set_dof_index<structdim> (level, this->present_index, 0, i, index);
}


namespace internal
{
  namespace DoFAccessor
  {
    template <int dim, int spacedim>
    inline
    const FiniteElement<dim,spacedim> &
    get_fe (const FiniteElement<dim,spacedim> &fe,
            const unsigned int)
    {
      return fe;
    }



    template <int dim, int spacedim>
    inline
    const FiniteElement<dim,spacedim> &
    get_fe (const dealii::hp::FECollection<dim,spacedim> &fe,
            const unsigned int                            index)
    {
      return fe[index];
    }
  }
}


template <int dim, class DH, bool level_dof_access>
inline
const FiniteElement<DH::dimension,DH::space_dimension> &
DoFAccessor<dim,DH,level_dof_access>::get_fe (const unsigned int fe_index) const
{
  Assert (fe_index_is_active (fe_index) == true,
          ExcMessage ("This function can only be called for active fe indices"));

  return dealii::internal::DoFAccessor::get_fe (this->dof_handler->get_fe(), fe_index);
}



namespace internal
{
  namespace DoFAccessor
  {
    template <class DH, bool level_dof_access>
    void get_dof_indices (const dealii::DoFAccessor<1,DH,level_dof_access>   &accessor,
                          std::vector<types::global_dof_index> &dof_indices,
                          const unsigned int         fe_index)
    {
      const unsigned int dofs_per_vertex = accessor.get_fe(fe_index).dofs_per_vertex,
                         dofs_per_line   = accessor.get_fe(fe_index).dofs_per_line;
      std::vector<types::global_dof_index>::iterator next = dof_indices.begin();
      for (unsigned int vertex=0; vertex<2; ++vertex)
        for (unsigned int d=0; d<dofs_per_vertex; ++d)
          *next++ = accessor.vertex_dof_index(vertex,d,fe_index);
      for (unsigned int d=0; d<dofs_per_line; ++d)
        *next++ = accessor.dof_index(d,fe_index);
    }



    template <class DH, bool level_dof_access>
    void get_dof_indices (const dealii::DoFAccessor<2,DH,level_dof_access>   &accessor,
                          std::vector<types::global_dof_index> &dof_indices,
                          const unsigned int         fe_index)
    {
      const unsigned int dofs_per_vertex = accessor.get_fe(fe_index).dofs_per_vertex,
                         dofs_per_line   = accessor.get_fe(fe_index).dofs_per_line,
                         dofs_per_quad   = accessor.get_fe(fe_index).dofs_per_quad;
      std::vector<types::global_dof_index>::iterator next = dof_indices.begin();
      for (unsigned int vertex=0; vertex<4; ++vertex)
        for (unsigned int d=0; d<dofs_per_vertex; ++d)
          *next++ = accessor.vertex_dof_index(vertex,d,fe_index);
      // now copy dof numbers from the line. for
      // lines with the wrong orientation (which
      // might occur in 3d), we have already made
      // sure that we're ok by picking the correct
      // vertices (this happens automatically in
      // the vertex() function). however, if the
      // line is in wrong orientation, we look at
      // it in flipped orientation and we will have
      // to adjust the shape function indices that
      // we see to correspond to the correct
      // (face-local) ordering.
      for (unsigned int line=0; line<4; ++line)
        for (unsigned int d=0; d<dofs_per_line; ++d)
          *next++ = accessor.line(line)->dof_index(accessor.get_fe(fe_index).
                                                   adjust_line_dof_index_for_line_orientation(d,
                                                       accessor.line_orientation(line)),
                                                   fe_index);
      for (unsigned int d=0; d<dofs_per_quad; ++d)
        *next++ = accessor.dof_index(d,fe_index);
    }



    template <class DH, bool level_dof_access>
    void get_dof_indices (const dealii::DoFAccessor<3,DH,level_dof_access>   &accessor,
                          std::vector<types::global_dof_index> &dof_indices,
                          const unsigned int         fe_index)
    {
      const unsigned int dofs_per_vertex = accessor.get_fe(fe_index).dofs_per_vertex,
                         dofs_per_line   = accessor.get_fe(fe_index).dofs_per_line,
                         dofs_per_quad   = accessor.get_fe(fe_index).dofs_per_quad,
                         dofs_per_hex    = accessor.get_fe(fe_index).dofs_per_hex;
      std::vector<types::global_dof_index>::iterator next = dof_indices.begin();
      for (unsigned int vertex=0; vertex<8; ++vertex)
        for (unsigned int d=0; d<dofs_per_vertex; ++d)
          *next++ = accessor.vertex_dof_index(vertex,d,fe_index);
      // now copy dof numbers from the line. for
      // lines with the wrong orientation, we have
      // already made sure that we're ok by picking
      // the correct vertices (this happens
      // automatically in the vertex()
      // function). however, if the line is in
      // wrong orientation, we look at it in
      // flipped orientation and we will have to
      // adjust the shape function indices that we
      // see to correspond to the correct
      // (cell-local) ordering.
      for (unsigned int line=0; line<12; ++line)
        for (unsigned int d=0; d<dofs_per_line; ++d)
          *next++ = accessor.line(line)->dof_index(accessor.get_fe(fe_index).
                                                   adjust_line_dof_index_for_line_orientation(d,
                                                       accessor.line_orientation(line)),fe_index);
      // now copy dof numbers from the face. for
      // faces with the wrong orientation, we
      // have already made sure that we're ok by
      // picking the correct lines and vertices
      // (this happens automatically in the
      // line() and vertex() functions). however,
      // if the face is in wrong orientation, we
      // look at it in flipped orientation and we
      // will have to adjust the shape function
      // indices that we see to correspond to the
      // correct (cell-local) ordering. The same
      // applies, if the face_rotation or
      // face_orientation is non-standard
      for (unsigned int quad=0; quad<6; ++quad)
        for (unsigned int d=0; d<dofs_per_quad; ++d)
          *next++ = accessor.quad(quad)->dof_index(accessor.get_fe(fe_index).
                                                   adjust_quad_dof_index_for_face_orientation(d,
                                                       accessor.face_orientation(quad),
                                                       accessor.face_flip(quad),
                                                       accessor.face_rotation(quad)),
                                                   fe_index);
      for (unsigned int d=0; d<dofs_per_hex; ++d)
        *next++ = accessor.dof_index(d,fe_index);
    }



    template<class DH, bool level_dof_access>
    void get_mg_dof_indices (const dealii::DoFAccessor<1, DH,level_dof_access> &accessor,
                             const int level,
                             std::vector<types::global_dof_index> &dof_indices,
                             const unsigned int fe_index)
    {
      const DH &handler = accessor.get_dof_handler();
      Assert(handler.n_dofs(level) != numbers::invalid_dof_index,
             ExcNotInitialized());

      const FiniteElement<DH::dimension, DH::space_dimension> &fe
        = handler.get_fe ()[fe_index];
      std::vector<types::global_dof_index>::iterator next = dof_indices.begin ();

      for (unsigned int vertex = 0; vertex < GeometryInfo<1>::vertices_per_cell; ++vertex)
        for (unsigned int dof = 0; dof < fe.dofs_per_vertex; ++dof)
          *next++ = accessor.mg_vertex_dof_index (level, vertex, dof);

      for (unsigned int dof = 0; dof < fe.dofs_per_line; ++dof)
        *next++ = accessor.mg_dof_index (level, dof);

      Assert (next == dof_indices.end (), ExcInternalError ());
    }



    template<class DH, bool level_dof_access>
    void get_mg_dof_indices (const dealii::DoFAccessor<2, DH,level_dof_access> &accessor,
                             const int level,
                             std::vector<types::global_dof_index> &dof_indices,
                             const unsigned int fe_index)
    {
      const DH &handler = accessor.get_dof_handler();
      Assert(handler.n_dofs(level) != numbers::invalid_dof_index,
             ExcNotInitialized());

      const FiniteElement<DH::dimension, DH::space_dimension> &fe = handler.get_fe ()[fe_index];
      std::vector<types::global_dof_index>::iterator next = dof_indices.begin ();

      for (unsigned int vertex = 0; vertex < GeometryInfo<2>::vertices_per_cell; ++vertex)
        for (unsigned int dof = 0; dof < fe.dofs_per_vertex; ++dof)
          *next++ = accessor.mg_vertex_dof_index (level, vertex, dof);

      for (unsigned int line = 0; line < GeometryInfo<2>::lines_per_cell; ++line)
        for (unsigned int dof = 0; dof < fe.dofs_per_line; ++dof)
          *next++ = accessor.line (line)->mg_dof_index (level, dof);

      for (unsigned int dof = 0; dof < fe.dofs_per_quad; ++dof)
        *next++ = accessor.mg_dof_index (level, dof);

      Assert (next == dof_indices.end (), ExcInternalError ());
    }



    template<class DH, bool level_dof_access>
    void get_mg_dof_indices (const dealii::DoFAccessor<3, DH,level_dof_access> &accessor,
                             const int level,
                             std::vector<types::global_dof_index> &dof_indices,
                             const unsigned int fe_index)
    {
      const DH &handler = accessor.get_dof_handler();
      Assert(handler.n_dofs(level) != numbers::invalid_dof_index,
             ExcNotInitialized());

      const FiniteElement<DH::dimension, DH::space_dimension> &fe = handler.get_fe ()[fe_index];
      std::vector<types::global_dof_index>::iterator next = dof_indices.begin ();

      for (unsigned int vertex = 0; vertex < GeometryInfo<3>::vertices_per_cell; ++vertex)
        for (unsigned int dof = 0; dof < fe.dofs_per_vertex; ++dof)
          *next++ = accessor.mg_vertex_dof_index (level, vertex, dof);

      for (unsigned int line = 0; line < GeometryInfo<3>::lines_per_cell; ++line)
        for (unsigned int dof = 0; dof < fe.dofs_per_line; ++dof)
          *next++ = accessor.line (line)->mg_dof_index (level, dof);

      for (unsigned int quad = 0; quad < GeometryInfo<3>::quads_per_cell; ++quad)
        for (unsigned int dof = 0; dof < fe.dofs_per_quad; ++dof)
          *next++ = accessor.quad (quad)->mg_dof_index (level, dof);

      for (unsigned int dof = 0; dof < fe.dofs_per_hex; ++dof)
        *next++ = accessor.mg_dof_index (level, dof);

      Assert (next == dof_indices.end (), ExcInternalError ());
    }


  }
}


template <int structdim, class DH, bool level_dof_access>
inline
void
DoFAccessor<structdim,DH,level_dof_access>::get_dof_indices (std::vector<types::global_dof_index> &dof_indices,
    const unsigned int         fe_index) const
{
  Assert (this->dof_handler != 0, ExcNotInitialized());
  Assert (&this->dof_handler->get_fe() != 0, ExcMessage ("DoFHandler not initialized"));
  Assert (static_cast<unsigned int>(this->level()) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));

  switch (structdim)
    {
    case 1:
      Assert (dof_indices.size() ==
              (2*this->dof_handler->get_fe()[fe_index].dofs_per_vertex +
               this->dof_handler->get_fe()[fe_index].dofs_per_line),
              ExcVectorDoesNotMatch());
      break;
    case 2:
      Assert (dof_indices.size() ==
              (4*this->dof_handler->get_fe()[fe_index].dofs_per_vertex +
               4*this->dof_handler->get_fe()[fe_index].dofs_per_line +
               this->dof_handler->get_fe()[fe_index].dofs_per_quad),
              ExcVectorDoesNotMatch());
      break;
    case 3:
      Assert (dof_indices.size() ==
              (8*this->dof_handler->get_fe()[fe_index].dofs_per_vertex +
               12*this->dof_handler->get_fe()[fe_index].dofs_per_line +
               6*this->dof_handler->get_fe()[fe_index].dofs_per_quad +
               this->dof_handler->get_fe()[fe_index].dofs_per_hex),
              ExcVectorDoesNotMatch());
      break;
    default:
      Assert (false, ExcNotImplemented());
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
  Assert (this->fe_index_is_active (fe_index)
          ||
          (this->dof_handler->get_fe()[fe_index].dofs_per_cell ==
           GeometryInfo<structdim>::vertices_per_cell *
           this->dof_handler->get_fe()[fe_index].dofs_per_vertex),
          ExcInternalError());

  // now do the actual work
  dealii::internal::DoFAccessor::get_dof_indices (*this, dof_indices, fe_index);
}



template<int structdim, class DH, bool level_dof_access>
inline
void DoFAccessor<structdim, DH,level_dof_access>::get_mg_dof_indices (const int level,
    std::vector<types::global_dof_index> &dof_indices,
    const unsigned int fe_index) const
{
  Assert (this->dof_handler != 0, ExcInvalidObject ());
  Assert (&this->dof_handler->get_fe () != 0, ExcInvalidObject ());

  switch (structdim)
    {
    case 1:
    {
      Assert (dof_indices.size () ==
              2 * this->dof_handler->get_fe ()[fe_index].dofs_per_vertex +
              this->dof_handler->get_fe ()[fe_index].dofs_per_line,
              ExcVectorDoesNotMatch ());
      break;
    }

    case 2:
    {
      Assert (dof_indices.size () ==
              4 * (this->dof_handler->get_fe ()[fe_index].dofs_per_vertex +
                   this->dof_handler->get_fe ()[fe_index].dofs_per_line) +
              this->dof_handler->get_fe ()[fe_index].dofs_per_quad,
              ExcVectorDoesNotMatch ());
      break;
    }

    case 3:
    {
      Assert (dof_indices.size () ==
              8 * this->dof_handler->get_fe ()[fe_index].dofs_per_vertex +
              12 * this->dof_handler->get_fe ()[fe_index].dofs_per_line +
              6 * this->dof_handler->get_fe ()[fe_index].dofs_per_quad +
              this->dof_handler->get_fe ()[fe_index].dofs_per_hex,
              ExcVectorDoesNotMatch ());
      break;
    }

    default:
      Assert (false, ExcNotImplemented ());
    }

  internal::DoFAccessor::get_mg_dof_indices (*this,
                                             level,
                                             dof_indices,
                                             fe_index);
}


template<int structdim, class DH, bool level_dof_access>
inline
void DoFAccessor<structdim, DH,level_dof_access>::set_mg_dof_indices (const int level,
    const std::vector<types::global_dof_index> &dof_indices,
    const unsigned int fe_index)
{
  Assert (this->dof_handler != 0, ExcInvalidObject ());
  Assert (&this->dof_handler->get_fe () != 0, ExcInvalidObject ());

  switch (structdim)
    {
    case 1:
    {
      Assert (dof_indices.size () ==
              2 * this->dof_handler->get_fe ()[fe_index].dofs_per_vertex +
              this->dof_handler->get_fe ()[fe_index].dofs_per_line,
              ExcVectorDoesNotMatch ());
      break;
    }

    case 2:
    {
      Assert (dof_indices.size () ==
              4 * (this->dof_handler->get_fe ()[fe_index].dofs_per_vertex +
                   this->dof_handler->get_fe ()[fe_index].dofs_per_line) +
              this->dof_handler->get_fe ()[fe_index].dofs_per_quad,
              ExcVectorDoesNotMatch ());
      break;
    }

    case 3:
    {
      Assert (dof_indices.size () ==
              8 * this->dof_handler->get_fe ()[fe_index].dofs_per_vertex +
              12 * this->dof_handler->get_fe ()[fe_index].dofs_per_line +
              6 * this->dof_handler->get_fe ()[fe_index].dofs_per_quad +
              this->dof_handler->get_fe ()[fe_index].dofs_per_hex,
              ExcVectorDoesNotMatch ());
      break;
    }

    default:
      Assert (false, ExcNotImplemented ());
    }

  internal::DoFAccessor::Implementation::set_mg_dof_indices (*this,
                                                             level,
                                                             dof_indices,
                                                             fe_index);
}


namespace internal
{
  namespace DoFAccessor
  {
    template <bool level_dof_access, class DH>
    inline
    typename dealii::internal::DoFHandler::Iterators<DH, level_dof_access>::quad_iterator
    get_quad(const dealii::Triangulation<DH::dimension, DH::space_dimension> *tria,
             unsigned int index,
             DH *dof_handler)
    {
    }


    template<bool level_dof_access>
    inline
    typename dealii::internal::DoFHandler::Iterators<dealii::DoFHandler<2,2>, level_dof_access>::quad_iterator
    get_quad(const dealii::Triangulation<2,2> *,
             unsigned int,
             dealii::DoFHandler<2,2> *)
    {
      Assert(false, ExcNotImplemented());
      return typename dealii::internal::DoFHandler::Iterators<dealii::DoFHandler<2,2>, level_dof_access>::line_iterator();
    }

    template<bool level_dof_access>
    inline
    typename dealii::internal::DoFHandler::Iterators<dealii::DoFHandler<2,3>, level_dof_access>::quad_iterator
    get_quad(const dealii::Triangulation<2,3> *,
             unsigned int,
             dealii::DoFHandler<2,3> *)
    {
      Assert(false, ExcNotImplemented());
      return typename dealii::internal::DoFHandler::Iterators<dealii::DoFHandler<2,3>, level_dof_access>::line_iterator();
    }

    template<bool level_dof_access>
    inline
    typename dealii::internal::DoFHandler::Iterators<dealii::hp::DoFHandler<2,2>, level_dof_access>::quad_iterator
    get_quad(const dealii::Triangulation<2,2> *,
             unsigned int,
             dealii::hp::DoFHandler<2,2> *)
    {
      Assert(false, ExcNotImplemented());
      return typename dealii::internal::DoFHandler::Iterators<dealii::hp::DoFHandler<2,2>, level_dof_access>::line_iterator();
    }

    template<bool level_dof_access>
    inline
    typename dealii::internal::DoFHandler::Iterators<dealii::hp::DoFHandler<2,3>, level_dof_access>::quad_iterator
    get_quad(const dealii::Triangulation<2,3> *,
             unsigned int,
             dealii::hp::DoFHandler<2,3> *)
    {
      Assert(false, ExcNotImplemented());
      return typename dealii::internal::DoFHandler::Iterators<dealii::hp::DoFHandler<2,3>, level_dof_access>::line_iterator();
    }
  }
}


template <int structdim, class DH, bool level_dof_access>
inline
typename dealii::internal::DoFHandler::Iterators<DH,level_dof_access>::line_iterator
DoFAccessor<structdim,DH,level_dof_access>::line (const unsigned int i) const
{
  // if we are asking for a particular line and this object refers to
  // a line, then the only valid index is i==0 and we should return
  // *this
  if (structdim == 1)
    {
      Assert (i==0, ExcMessage ("You can only ask for line zero if the "
                                "current object is a line itself."));
      return
        typename dealii::internal::DoFHandler::Iterators<DH,level_dof_access>::cell_iterator
        (&this->get_triangulation(),
         this->level(),
         this->index(),
         &this->get_dof_handler());
    }

  // otherwise we need to be in structdim>=2
  Assert (structdim > 1, ExcImpossibleInDim(structdim));
  Assert (DH::dimension > 1, ExcImpossibleInDim(DH::dimension));

  // checking of 'i' happens in line_index(i)
  return typename dealii::internal::DoFHandler::Iterators<DH,level_dof_access>::line_iterator
         (this->tria,
          0,  // only sub-objects are allowed, which have no level
          this->line_index(i),
          this->dof_handler);
}


template <int structdim, class DH, bool level_dof_access>
inline
typename dealii::internal::DoFHandler::Iterators<DH,level_dof_access>::quad_iterator
DoFAccessor<structdim,DH,level_dof_access>::quad (const unsigned int i) const
{
  // if we are asking for a
  // particular quad and this object
  // refers to a quad, then the only
  // valid index is i==0 and we
  // should return *this
  if (structdim == 2)
    {
      Assert (i==0, ExcMessage ("You can only ask for quad zero if the "
                                "current object is a quad itself."));
      return
        typename dealii::internal::DoFHandler::Iterators<DH>::cell_iterator
        (&this->get_triangulation(),
         this->level(),
         this->index(),
         &this->get_dof_handler());
    }

  // otherwise we need to be in structdim>=3
  Assert (structdim > 2, ExcImpossibleInDim(structdim));
  Assert (DH::dimension > 2, ExcImpossibleInDim(DH::dimension));

  // checking of 'i' happens in quad_index(i)
  return typename dealii::internal::DoFHandler::Iterators<DH,level_dof_access>::quad_iterator
         (this->tria,
          0,  // only sub-objects are allowed, which have no level
          this->quad_index(i),
          this->dof_handler);
}


/*------------------------- Functions: DoFAccessor<0,1,spacedim> ---------------------------*/


template <template <int, int> class DH, int spacedim, bool level_dof_access>
inline
DoFAccessor<0,DH<1,spacedim>, level_dof_access>::DoFAccessor ()
{
  Assert (false, ExcInvalidObject());
}



template <template <int, int> class DH, int spacedim, bool level_dof_access>
inline
DoFAccessor<0,DH<1,spacedim>, level_dof_access>::
DoFAccessor (const Triangulation<1,spacedim> *tria,
             const typename TriaAccessor<0,1,spacedim>::VertexKind vertex_kind,
             const unsigned int    vertex_index,
             const DH<1,spacedim> *dof_handler)
  :
  BaseClass (tria,
             vertex_kind,
             vertex_index),
  dof_handler(const_cast<DH<1,spacedim>*>(dof_handler))
{}



template <template <int, int> class DH, int spacedim, bool level_dof_access>
inline
DoFAccessor<0,DH<1,spacedim>, level_dof_access>::
DoFAccessor (const Triangulation<1,spacedim> *,
             const int                 ,
             const int                 ,
             const DH<1,spacedim> *)
  :
  dof_handler(0)
{
  Assert (false,
          ExcMessage ("This constructor can not be called for face iterators in 1d."));
}



template <template <int, int> class DH, int spacedim, bool level_dof_access>
template <int structdim2, int dim2, int spacedim2>
DoFAccessor<0,DH<1,spacedim>, level_dof_access>::DoFAccessor (const InvalidAccessor<structdim2,dim2,spacedim2> &)
{
  Assert (false, ExcInvalidObject());
}



template <template <int, int> class DH, int spacedim, bool level_dof_access>
template <int dim2, class DH2, bool level_dof_access2>
inline
DoFAccessor<0,DH<1,spacedim>, level_dof_access>::DoFAccessor (const DoFAccessor<dim2, DH2, level_dof_access2> &)
{
  Assert (false, ExcInvalidObject());
}



template <template <int, int> class DH, int spacedim, bool level_dof_access>
inline
void
DoFAccessor<0,DH<1,spacedim>, level_dof_access>::set_dof_handler (DH<1,spacedim> *dh)
{
  Assert (dh != 0, ExcInvalidObject());
  this->dof_handler = dh;
}



template <template <int, int> class DH, int spacedim, bool level_dof_access>
inline
const DH<1,spacedim> &
DoFAccessor<0,DH<1,spacedim>, level_dof_access>::get_dof_handler () const
{
  return *this->dof_handler;
}



template <template <int, int> class DH, int spacedim, bool level_dof_access>
inline
void
DoFAccessor<0,DH<1,spacedim>, level_dof_access>::get_dof_indices (
  std::vector<types::global_dof_index> &dof_indices,
  const unsigned int fe_index) const
{
  for (unsigned int i=0; i<dof_indices.size(); ++i)
    dof_indices[i]
      = dealii::internal::DoFAccessor::Implementation::get_vertex_dof_index (
          *dof_handler,
          this->global_vertex_index,
          fe_index,
          i);
}



template <template <int, int> class DH, int spacedim, bool level_dof_access>
inline
types::global_dof_index
DoFAccessor<0,DH<1,spacedim>, level_dof_access>::
vertex_dof_index (const unsigned int vertex,
                  const unsigned int i,
                  const unsigned int fe_index) const
{
  Assert (vertex == 0, ExcIndexRange (vertex, 0, 1));
  return dealii::internal::DoFAccessor::Implementation::get_vertex_dof_index (
           *dof_handler,
           this->global_vertex_index,
           fe_index,
           i);
}



template <template <int, int> class DH, int spacedim, bool level_dof_access>
inline
void
DoFAccessor<0,DH<1,spacedim>, level_dof_access>::copy_from (const TriaAccessorBase<0,1,spacedim> &da)
{
  Assert (this->dof_handler != 0, ExcInvalidObject());
  BaseClass::copy_from(da);
}



template <template <int, int> class DH, int spacedim, bool level_dof_access>
template <bool level_dof_access2>
inline
void
DoFAccessor<0,DH<1,spacedim>, level_dof_access>::copy_from (const DoFAccessor<0,DH<1,spacedim>, level_dof_access2> &a)
{
  BaseClass::copy_from (a);
  set_dof_handler (a.dof_handler);
}



template <template <int, int> class DH, int spacedim, bool level_dof_access>
template <int dim2, class DH2, bool level_dof_access2>
inline
bool
DoFAccessor<0,DH<1,spacedim>, level_dof_access>::operator == (const DoFAccessor<dim2,DH2,level_dof_access2> &a) const
{
  Assert (dim2 == 0, ExcCantCompareIterators());
  Assert (this->dof_handler == a.dof_handler, ExcCantCompareIterators());
  return (BaseClass::operator == (a));
}



template <template <int, int> class DH, int spacedim, bool level_dof_access>
template <int dim2, class DH2, bool level_dof_access2>
inline
bool
DoFAccessor<0,DH<1,spacedim>, level_dof_access>::operator != (const DoFAccessor<dim2,DH2,level_dof_access2> &a) const
{
  Assert (dim2 == 0, ExcCantCompareIterators());
  Assert (this->dof_handler == a.dof_handler, ExcCantCompareIterators());
  return (BaseClass::operator != (a));
}



/*------------------------- Functions: DoFCellAccessor -----------------------*/


namespace internal
{
  namespace DoFCellAccessor
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
       * Implement the updating of the
       * cache. Currently not
       * implemented for hp::DoFHandler
       * objects.
       */
      template <int spacedim, bool level_dof_access>
      static
      void
      update_cell_dof_indices_cache (const DoFCellAccessor<DoFHandler<1,spacedim>, level_dof_access> &accessor)
      {
        // check as in documentation that
        // cell is either active, or dofs
        // are only in vertices. otherwise
        // simply don't update the cache at
        // all. the get_dof_indices
        // function will then make sure we
        // don't access the invalid data
        if (accessor.has_children()
            &&
            (accessor.get_fe().dofs_per_cell !=
             accessor.get_fe().dofs_per_vertex * GeometryInfo<1>::vertices_per_cell))
          return;

        const unsigned int dofs_per_vertex = accessor.get_fe().dofs_per_vertex,
                           dofs_per_line   = accessor.get_fe().dofs_per_line,
                           dofs_per_cell   = accessor.get_fe().dofs_per_cell;

        // make sure the cache is at least
        // as big as we need it when
        // writing to the last element of
        // this cell
        Assert (accessor.present_index * dofs_per_cell + dofs_per_cell
                <=
                accessor.dof_handler->levels[accessor.present_level]
                ->cell_dof_indices_cache.size(),
                ExcInternalError());

        std::vector<types::global_dof_index>::iterator next
          = (accessor.dof_handler->levels[accessor.present_level]
             ->cell_dof_indices_cache.begin() + accessor.present_index * dofs_per_cell);

        for (unsigned int vertex=0; vertex<2; ++vertex)
          for (unsigned int d=0; d<dofs_per_vertex; ++d)
            *next++ = accessor.vertex_dof_index(vertex,d);
        for (unsigned int d=0; d<dofs_per_line; ++d)
          *next++ = accessor.dof_index(d);
      }



      template <int spacedim, bool level_dof_access>
      static
      void
      update_cell_dof_indices_cache (const DoFCellAccessor<DoFHandler<2,spacedim>, level_dof_access> &accessor)
      {
        // check as in documentation that
        // cell is either active, or dofs
        // are only in vertices. otherwise
        // simply don't update the cache at
        // all. the get_dof_indices
        // function will then make sure we
        // don't access the invalid data
        if (accessor.has_children()
            &&
            (accessor.get_fe().dofs_per_cell !=
             accessor.get_fe().dofs_per_vertex * GeometryInfo<2>::vertices_per_cell))
          return;

        const unsigned int dofs_per_vertex = accessor.get_fe().dofs_per_vertex,
                           dofs_per_line   = accessor.get_fe().dofs_per_line,
                           dofs_per_quad   = accessor.get_fe().dofs_per_quad,
                           dofs_per_cell   = accessor.get_fe().dofs_per_cell;

        // make sure the cache is at least
        // as big as we need it when
        // writing to the last element of
        // this cell
        Assert (accessor.present_index * dofs_per_cell + dofs_per_cell
                <=
                accessor.dof_handler->levels[accessor.present_level]
                ->cell_dof_indices_cache.size(),
                ExcInternalError());

        std::vector<types::global_dof_index>::iterator next
          = (accessor.dof_handler->levels[accessor.present_level]
             ->cell_dof_indices_cache.begin() + accessor.present_index * dofs_per_cell);

        for (unsigned int vertex=0; vertex<4; ++vertex)
          for (unsigned int d=0; d<dofs_per_vertex; ++d)
            *next++ = accessor.vertex_dof_index(vertex,d);
        for (unsigned int line=0; line<4; ++line)
          for (unsigned int d=0; d<dofs_per_line; ++d)
            *next++ = accessor.line(line)->dof_index(d);
        for (unsigned int d=0; d<dofs_per_quad; ++d)
          *next++ = accessor.dof_index(d);
      }


      template <int spacedim, bool level_dof_access>
      static
      void
      update_cell_dof_indices_cache (const DoFCellAccessor<DoFHandler<3,spacedim>, level_dof_access> &accessor)
      {
        // check as in documentation that
        // cell is either active, or dofs
        // are only in vertices. otherwise
        // simply don't update the cache at
        // all. the get_dof_indices
        // function will then make sure we
        // don't access the invalid data
        if (accessor.has_children()
            &&
            (accessor.get_fe().dofs_per_cell !=
             accessor.get_fe().dofs_per_vertex * GeometryInfo<3>::vertices_per_cell))
          return;

        const unsigned int dofs_per_vertex = accessor.get_fe().dofs_per_vertex,
                           dofs_per_line   = accessor.get_fe().dofs_per_line,
                           dofs_per_quad   = accessor.get_fe().dofs_per_quad,
                           dofs_per_hex    = accessor.get_fe().dofs_per_hex,
                           dofs_per_cell   = accessor.get_fe().dofs_per_cell;

        // make sure the cache is at least
        // as big as we need it when
        // writing to the last element of
        // this cell
        Assert (accessor.present_index * dofs_per_cell + dofs_per_cell
                <=
                accessor.dof_handler->levels[accessor.present_level]
                ->cell_dof_indices_cache.size(),
                ExcInternalError());

        std::vector<types::global_dof_index>::iterator next
          = (accessor.dof_handler->levels[accessor.present_level]
             ->cell_dof_indices_cache.begin() + accessor.present_index * dofs_per_cell);

        for (unsigned int vertex=0; vertex<8; ++vertex)
          for (unsigned int d=0; d<dofs_per_vertex; ++d)
            *next++ = accessor.vertex_dof_index(vertex,d);
        // now copy dof numbers from the line. for
        // lines with the wrong orientation, we have
        // already made sure that we're ok by picking
        // the correct vertices (this happens
        // automatically in the vertex()
        // function). however, if the line is in
        // wrong orientation, we look at it in
        // flipped orientation and we will have to
        // adjust the shape function indices that we
        // see to correspond to the correct
        // (cell-local) ordering.
        for (unsigned int line=0; line<12; ++line)
          for (unsigned int d=0; d<dofs_per_line; ++d)
            *next++ = accessor.line(line)->dof_index(accessor.dof_handler->get_fe().
                                                     adjust_line_dof_index_for_line_orientation(d,
                                                         accessor.line_orientation(line)));
        // now copy dof numbers from the face. for
        // faces with the wrong orientation, we
        // have already made sure that we're ok by
        // picking the correct lines and vertices
        // (this happens automatically in the
        // line() and vertex() functions). however,
        // if the face is in wrong orientation, we
        // look at it in flipped orientation and we
        // will have to adjust the shape function
        // indices that we see to correspond to the
        // correct (cell-local) ordering. The same
        // applies, if the face_rotation or
        // face_orientation is non-standard
        for (unsigned int quad=0; quad<6; ++quad)
          for (unsigned int d=0; d<dofs_per_quad; ++d)
            *next++ = accessor.quad(quad)->dof_index(accessor.dof_handler->get_fe().
                                                     adjust_quad_dof_index_for_face_orientation(d,
                                                         accessor.face_orientation(quad),
                                                         accessor.face_flip(quad),
                                                         accessor.face_rotation(quad)));
        for (unsigned int d=0; d<dofs_per_hex; ++d)
          *next++ = accessor.dof_index(d);
      }


      // implementation for the case of
      // hp::DoFHandler objects. it's
      // not implemented there, for no
      // space dimension
      template <int dim, int spacedim, bool level_dof_access>
      static
      void
      update_cell_dof_indices_cache (const DoFCellAccessor<dealii::hp::DoFHandler<dim,spacedim>, level_dof_access> &accessor)
      {
        // caches are only for cells with DoFs, i.e., for active ones
        if (accessor.has_children())
          return;

        const unsigned int dofs_per_cell   = accessor.get_fe().dofs_per_cell;

        // make sure the cache is at least
        // as big as we need it when
        // writing to the last element of
        // this cell
        Assert (static_cast<unsigned int>(accessor.present_index)
                <
                accessor.dof_handler->levels[accessor.present_level]
                ->cell_cache_offsets.size(),
                ExcInternalError());
        Assert (accessor.dof_handler->levels[accessor.present_level]
                ->cell_cache_offsets[accessor.present_index]
                <=
                accessor.dof_handler->levels[accessor.present_level]
                ->cell_dof_indices_cache.size(),
                ExcInternalError());

        std::vector<types::global_dof_index> dof_indices (dofs_per_cell);
        static_cast<const dealii::DoFAccessor<dim,dealii::hp::DoFHandler<dim,spacedim>,level_dof_access> &>
        (accessor).get_dof_indices (dof_indices, accessor.active_fe_index());

        types::global_dof_index *next_dof_index
          = &accessor.dof_handler->levels[accessor.present_level]
            ->cell_dof_indices_cache[accessor.dof_handler->levels[accessor.present_level]
                                     ->cell_cache_offsets[accessor.present_index]];
        for (unsigned int i=0; i<dofs_per_cell; ++i, ++next_dof_index)
          *next_dof_index = dof_indices[i];
      }



      /**
       * Implement setting dof
       * indices on a
       * cell. Currently not
       * implemented for
       * hp::DoFHandler objects.
       */
      template <int spacedim, bool level_dof_access>
      static
      void
      set_dof_indices (DoFCellAccessor<DoFHandler<1,spacedim>, level_dof_access> &accessor,
                       const std::vector<types::global_dof_index>          &local_dof_indices)
      {
        Assert (accessor.has_children() == false,
                ExcInternalError());

        const unsigned int dofs_per_vertex = accessor.get_fe().dofs_per_vertex,
                           dofs_per_line   = accessor.get_fe().dofs_per_line,
                           dofs_per_cell   = accessor.get_fe().dofs_per_cell;

        Assert (local_dof_indices.size() == dofs_per_cell,
                ExcInternalError());

        unsigned int index = 0;

        for (unsigned int vertex=0; vertex<2; ++vertex)
          for (unsigned int d=0; d<dofs_per_vertex; ++d, ++index)
            accessor.set_vertex_dof_index(vertex,d,
                                          local_dof_indices[index]);

        for (unsigned int d=0; d<dofs_per_line; ++d, ++index)
          accessor.set_dof_index(d, local_dof_indices[index]);

        Assert (index == dofs_per_cell,
                ExcInternalError());
      }



      template <int spacedim, bool level_dof_access>
      static
      void
      set_dof_indices (DoFCellAccessor<DoFHandler<2,spacedim>, level_dof_access> &accessor,
                       const std::vector<types::global_dof_index>          &local_dof_indices)
      {
        Assert (accessor.has_children() == false,
                ExcInternalError());

        const unsigned int dofs_per_vertex = accessor.get_fe().dofs_per_vertex,
                           dofs_per_line   = accessor.get_fe().dofs_per_line,
                           dofs_per_quad   = accessor.get_fe().dofs_per_quad,
                           dofs_per_cell   = accessor.get_fe().dofs_per_cell;

        Assert (local_dof_indices.size() == dofs_per_cell,
                ExcInternalError());

        unsigned int index = 0;

        for (unsigned int vertex=0; vertex<4; ++vertex)
          for (unsigned int d=0; d<dofs_per_vertex; ++d, ++index)
            accessor.set_vertex_dof_index(vertex,d,
                                          local_dof_indices[index]);
        for (unsigned int line=0; line<4; ++line)
          for (unsigned int d=0; d<dofs_per_line; ++d, ++index)
            accessor.line(line)->set_dof_index(d, local_dof_indices[index]);

        for (unsigned int d=0; d<dofs_per_quad; ++d, ++index)
          accessor.set_dof_index(d, local_dof_indices[index]);

        Assert (index == dofs_per_cell,
                ExcInternalError());
      }



      template <int spacedim, bool level_dof_access>
      static
      void
      set_dof_indices (DoFCellAccessor<DoFHandler<3,spacedim>, level_dof_access> &accessor,
                       const std::vector<types::global_dof_index>          &local_dof_indices)
      {
        Assert (accessor.has_children() == false,
                ExcInternalError());

        const unsigned int dofs_per_vertex = accessor.get_fe().dofs_per_vertex,
                           dofs_per_line   = accessor.get_fe().dofs_per_line,
                           dofs_per_quad   = accessor.get_fe().dofs_per_quad,
                           dofs_per_hex    = accessor.get_fe().dofs_per_hex,
                           dofs_per_cell   = accessor.get_fe().dofs_per_cell;

        Assert (local_dof_indices.size() == dofs_per_cell,
                ExcInternalError());

        unsigned int index = 0;

        for (unsigned int vertex=0; vertex<8; ++vertex)
          for (unsigned int d=0; d<dofs_per_vertex; ++d, ++index)
            accessor.set_vertex_dof_index(vertex,d,
                                          local_dof_indices[index]);
        // now copy dof numbers into the line. for
        // lines with the wrong orientation, we have
        // already made sure that we're ok by picking
        // the correct vertices (this happens
        // automatically in the vertex()
        // function). however, if the line is in
        // wrong orientation, we look at it in
        // flipped orientation and we will have to
        // adjust the shape function indices that we
        // see to correspond to the correct
        // (cell-local) ordering.
        for (unsigned int line=0; line<12; ++line)
          for (unsigned int d=0; d<dofs_per_line; ++d, ++index)
            accessor.line(line)->set_dof_index(accessor.dof_handler->get_fe().
                                               adjust_line_dof_index_for_line_orientation(d,
                                                   accessor.line_orientation(line)),
                                               local_dof_indices[index]);
        // now copy dof numbers into the face. for
        // faces with the wrong orientation, we
        // have already made sure that we're ok by
        // picking the correct lines and vertices
        // (this happens automatically in the
        // line() and vertex() functions). however,
        // if the face is in wrong orientation, we
        // look at it in flipped orientation and we
        // will have to adjust the shape function
        // indices that we see to correspond to the
        // correct (cell-local) ordering. The same
        // applies, if the face_rotation or
        // face_orientation is non-standard
        for (unsigned int quad=0; quad<6; ++quad)
          for (unsigned int d=0; d<dofs_per_quad; ++d, ++index)
            accessor.quad(quad)->set_dof_index(accessor.dof_handler->get_fe().
                                               adjust_quad_dof_index_for_face_orientation(d,
                                                   accessor.face_orientation(quad),
                                                   accessor.face_flip(quad),
                                                   accessor.face_rotation(quad)),
                                               local_dof_indices[index]);
        for (unsigned int d=0; d<dofs_per_hex; ++d, ++index)
          accessor.set_dof_index(d, local_dof_indices[index]);

        Assert (index == dofs_per_cell,
                ExcInternalError());
      }


      // implementation for the case of
      // hp::DoFHandler objects. it's
      // not implemented there, for no
      // space dimension
      template <int dim, int spacedim, bool level_dof_access>
      static
      void
      set_dof_indices (const DoFCellAccessor<dealii::hp::DoFHandler<dim,spacedim>, level_dof_access> &,
                       const std::vector<types::global_dof_index> &)
      {
        Assert (false, ExcNotImplemented());
      }



      /**
       * Do what the active_fe_index
       * function in the parent class
       * is supposed to do.
       */
      template <int dim, int spacedim, bool level_dof_access>
      static
      unsigned int
      active_fe_index (const DoFCellAccessor<DoFHandler<dim,spacedim>, level_dof_access> &)
      {
        // ::DoFHandler only supports a
        // single active fe with index
        // zero
        return 0;
      }



      template <int dim, int spacedim, bool level_dof_access>
      static
      unsigned int
      active_fe_index (const DoFCellAccessor<dealii::hp::DoFHandler<dim,spacedim>, level_dof_access> &accessor)
      {
        Assert (static_cast<unsigned int>(accessor.level()) < accessor.dof_handler->levels.size(),
                ExcMessage ("DoFHandler not initialized"));

        return accessor.dof_handler->levels[accessor.level()]
               ->active_fe_index(accessor.present_index);
      }



      /**
       * Do what the
       * set_active_fe_index function
       * in the parent class is
       * supposed to do.
       */
      template <int dim, int spacedim, bool level_dof_access>
      static
      void
      set_active_fe_index (const DoFCellAccessor<DoFHandler<dim,spacedim>, level_dof_access> &,
                           const unsigned int                                i)
      {
        // ::DoFHandler only supports a
        // single active fe with index
        // zero
        typedef dealii::DoFAccessor<dim,DoFHandler<dim,spacedim>, level_dof_access> BaseClass;
        Assert (i == 0, typename BaseClass::ExcInvalidObject());
      }



      template <int dim, int spacedim, bool level_dof_access>
      static
      void
      set_active_fe_index (DoFCellAccessor<dealii::hp::DoFHandler<dim,spacedim>, level_dof_access> &accessor,
                           const unsigned int                                      i)
      {
        typedef dealii::DoFAccessor<dim,DoFHandler<dim,spacedim>, level_dof_access> BaseClass;
        Assert (accessor.dof_handler != 0,
                typename BaseClass::ExcInvalidObject());
        Assert (static_cast<unsigned int>(accessor.level()) <
                accessor.dof_handler->levels.size(),
                ExcMessage ("DoFHandler not initialized"));

        accessor.dof_handler->levels[accessor.level()]
        ->set_active_fe_index (accessor.present_index, i);
      }



      template <int dim, int spacedim, bool level_dof_access, typename ForwardIterator, class OutputVector>
      static
      void
      distribute_local_to_global (const DoFCellAccessor<dealii::DoFHandler<dim,spacedim>, level_dof_access> &accessor,
                                  ForwardIterator local_source_begin,
                                  ForwardIterator local_source_end,
                                  OutputVector   &global_destination)
      {
        typedef dealii::DoFAccessor<dim,DoFHandler<dim,spacedim>, level_dof_access> BaseClass;
        Assert (accessor.dof_handler != 0,
                typename BaseClass::ExcInvalidObject());
        Assert (&accessor.get_fe() != 0,
                typename BaseClass::ExcInvalidObject());
        Assert (static_cast<unsigned int>(local_source_end-local_source_begin)
                ==
                accessor.get_fe().dofs_per_cell,
                typename BaseClass::ExcVectorDoesNotMatch());
        Assert (accessor.dof_handler->n_dofs() == global_destination.size(),
                typename BaseClass::ExcVectorDoesNotMatch());

        Assert (!accessor.has_children(),
                ExcMessage ("Cell must be active"));

        const unsigned int n_dofs = local_source_end - local_source_begin;

        types::global_dof_index *dofs = &accessor.dof_handler->levels[accessor.level()]
                                        ->cell_dof_indices_cache[accessor.present_index * n_dofs];

        // distribute cell vector
        global_destination.add(n_dofs, dofs, local_source_begin);
      }



      template <int dim, int spacedim, bool level_dof_access, typename ForwardIterator, class OutputVector>
      static
      void
      distribute_local_to_global (const DoFCellAccessor<dealii::hp::DoFHandler<dim,spacedim>, level_dof_access> &accessor,
                                  ForwardIterator local_source_begin,
                                  ForwardIterator local_source_end,
                                  OutputVector   &global_destination)
      {
        typedef dealii::DoFAccessor<dim,DoFHandler<dim,spacedim>, level_dof_access> BaseClass;
        Assert (accessor.dof_handler != 0,
                typename BaseClass::ExcInvalidObject());
        Assert (&accessor.get_fe() != 0,
                typename BaseClass::ExcInvalidObject());
        Assert (local_source_end-local_source_begin == accessor.get_fe().dofs_per_cell,
                typename BaseClass::ExcVectorDoesNotMatch());
        Assert (accessor.dof_handler->n_dofs() == global_destination.size(),
                typename BaseClass::ExcVectorDoesNotMatch());

        const unsigned int n_dofs = local_source_end - local_source_begin;

//TODO[WB/MK]: This function could me made more efficient because it allocates memory, which could be avoided by passing in another argument as a scratch array. This should be fixed eventually. another option would be to let the surrounding class have a (static, mutable) scratch array that is thread-local

        // get indices of dofs
        std::vector<types::global_dof_index> dofs (n_dofs);
        accessor.get_dof_indices (dofs);

        // distribute cell vector
        global_destination.add (n_dofs, dofs.begin(), local_source_begin);
      }



      template <int dim, int spacedim, bool level_dof_access, typename ForwardIterator, class OutputVector>
      static
      void
      distribute_local_to_global (const DoFCellAccessor<dealii::DoFHandler<dim,spacedim>, level_dof_access> &accessor,
                                  const ConstraintMatrix &constraints,
                                  ForwardIterator         local_source_begin,
                                  ForwardIterator         local_source_end,
                                  OutputVector           &global_destination)
      {
        typedef dealii::DoFAccessor<dim,DoFHandler<dim,spacedim>, level_dof_access> BaseClass;
        Assert (accessor.dof_handler != 0,
                typename BaseClass::ExcInvalidObject());
        Assert (&accessor.get_fe() != 0,
                typename BaseClass::ExcInvalidObject());
        Assert (local_source_end-local_source_begin == accessor.get_fe().dofs_per_cell,
                typename BaseClass::ExcVectorDoesNotMatch());
        Assert (accessor.dof_handler->n_dofs() == global_destination.size(),
                typename BaseClass::ExcVectorDoesNotMatch());

        Assert (!accessor.has_children(),
                ExcMessage ("Cell must be active."));

        const unsigned int n_dofs = local_source_end - local_source_begin;

        types::global_dof_index *dofs = &accessor.dof_handler->levels[accessor.level()]
                                        ->cell_dof_indices_cache[accessor.present_index * n_dofs];

        // distribute cell vector
        constraints.distribute_local_to_global (local_source_begin, local_source_end,
                                                dofs, global_destination);
      }



      template <int dim, int spacedim, bool level_dof_access, typename ForwardIterator, class OutputVector>
      static
      void
      distribute_local_to_global (const DoFCellAccessor<dealii::hp::DoFHandler<dim,spacedim>, level_dof_access> &accessor,
                                  const ConstraintMatrix &constraints,
                                  ForwardIterator         local_source_begin,
                                  ForwardIterator         local_source_end,
                                  OutputVector           &global_destination)
      {
        typedef dealii::DoFAccessor<dim,DoFHandler<dim,spacedim>, level_dof_access> BaseClass;
        Assert (accessor.dof_handler != 0,
                typename BaseClass::ExcInvalidObject());
        Assert (&accessor.get_fe() != 0,
                typename BaseClass::ExcInvalidObject());
        Assert (local_source_end-local_source_begin == accessor.get_fe().dofs_per_cell,
                typename BaseClass::ExcVectorDoesNotMatch());
        Assert (accessor.dof_handler->n_dofs() == global_destination.size(),
                typename BaseClass::ExcVectorDoesNotMatch());

        const unsigned int n_dofs = local_source_end - local_source_begin;

//TODO[WB/MK]: This function could me made more efficient because it allocates memory, which could be avoided by passing in another argument as a scratch array. This should be fixed eventually

        // get indices of dofs
        std::vector<types::global_dof_index> dofs (n_dofs);
        accessor.get_dof_indices (dofs);

        // distribute cell vector
        constraints.distribute_local_to_global (local_source_begin, local_source_end,
                                                dofs.begin(), global_destination);
      }



      template <int dim, int spacedim, bool level_dof_access, typename number, class OutputMatrix>
      static
      void
      distribute_local_to_global (const DoFCellAccessor<dealii::DoFHandler<dim,spacedim>, level_dof_access> &accessor,
                                  const dealii::FullMatrix<number> &local_source,
                                  OutputMatrix                     &global_destination)
      {
        typedef dealii::DoFAccessor<dim,DoFHandler<dim,spacedim>, level_dof_access> BaseClass;
        Assert (accessor.dof_handler != 0,
                typename BaseClass::ExcInvalidObject());
        Assert (&accessor.get_fe() != 0,
                typename BaseClass::ExcInvalidObject());
        Assert (local_source.m() == accessor.get_fe().dofs_per_cell,
                typename BaseClass::ExcMatrixDoesNotMatch());
        Assert (local_source.n() == accessor.get_fe().dofs_per_cell,
                typename BaseClass::ExcMatrixDoesNotMatch());
        Assert (accessor.dof_handler->n_dofs() == global_destination.m(),
                typename BaseClass::ExcMatrixDoesNotMatch());
        Assert (accessor.dof_handler->n_dofs() == global_destination.n(),
                typename BaseClass::ExcMatrixDoesNotMatch());

        Assert (!accessor.has_children(),
                ExcMessage ("Cell must be active."));

        const unsigned int n_dofs = local_source.m();

        types::global_dof_index *dofs = &accessor.dof_handler->levels[accessor.level()]
                                        ->cell_dof_indices_cache[accessor.present_index * n_dofs];

        // distribute cell matrix
        for (unsigned int i=0; i<n_dofs; ++i)
          global_destination.add(dofs[i], n_dofs, dofs,
                                 &local_source(i,0));
      }



      template <int dim, int spacedim, bool level_dof_access, typename number, class OutputMatrix>
      static
      void
      distribute_local_to_global (const DoFCellAccessor<dealii::hp::DoFHandler<dim,spacedim>, level_dof_access> &accessor,
                                  const dealii::FullMatrix<number> &local_source,
                                  OutputMatrix                     &global_destination)
      {
        typedef dealii::DoFAccessor<dim,DoFHandler<dim,spacedim>, level_dof_access> BaseClass;
        Assert (accessor.dof_handler != 0,
                typename BaseClass::ExcInvalidObject());
        Assert (&accessor.get_fe() != 0,
                typename BaseClass::ExcInvalidObject());
        Assert (local_source.m() == accessor.get_fe().dofs_per_cell,
                typename BaseClass::ExcMatrixDoesNotMatch());
        Assert (local_source.n() == accessor.get_fe().dofs_per_cell,
                typename BaseClass::ExcVectorDoesNotMatch());
        Assert (accessor.dof_handler->n_dofs() == global_destination.m(),
                typename BaseClass::ExcMatrixDoesNotMatch());
        Assert (accessor.dof_handler->n_dofs() == global_destination.n(),
                typename BaseClass::ExcMatrixDoesNotMatch());

        const unsigned int n_dofs = local_source.size();

//TODO[WB/MK]: This function could me made more efficient because it allocates memory, which could be avoided by passing in another argument as a scratch array.

        // get indices of dofs
        std::vector<types::global_dof_index> dofs (n_dofs);
        accessor.get_dof_indices (dofs);

        // distribute cell matrix
        global_destination.add(dofs,local_source);
      }



      template <int dim, int spacedim, bool level_dof_access, typename number,
                class OutputMatrix, typename OutputVector>
      static
      void
      distribute_local_to_global (const DoFCellAccessor<dealii::DoFHandler<dim,spacedim>, level_dof_access> &accessor,
                                  const dealii::FullMatrix<number> &local_matrix,
                                  const dealii::Vector<number>     &local_vector,
                                  OutputMatrix                     &global_matrix,
                                  OutputVector                     &global_vector)
      {
        typedef dealii::DoFAccessor<dim,DoFHandler<dim,spacedim>, level_dof_access> BaseClass;
        Assert (accessor.dof_handler != 0,
                typename BaseClass::ExcInvalidObject());
        Assert (&accessor.get_fe() != 0,
                typename BaseClass::ExcInvalidObject());
        Assert (local_matrix.m() == accessor.get_fe().dofs_per_cell,
                typename BaseClass::ExcMatrixDoesNotMatch());
        Assert (local_matrix.n() == accessor.get_fe().dofs_per_cell,
                typename BaseClass::ExcVectorDoesNotMatch());
        Assert (accessor.dof_handler->n_dofs() == global_matrix.m(),
                typename BaseClass::ExcMatrixDoesNotMatch());
        Assert (accessor.dof_handler->n_dofs() == global_matrix.n(),
                typename BaseClass::ExcMatrixDoesNotMatch());
        Assert (local_vector.size() == accessor.get_fe().dofs_per_cell,
                typename BaseClass::ExcVectorDoesNotMatch());
        Assert (accessor.dof_handler->n_dofs() == global_vector.size(),
                typename BaseClass::ExcVectorDoesNotMatch());

        Assert (!accessor.has_children(),
                ExcMessage ("Cell must be active."));

        const unsigned int n_dofs = accessor.get_fe().dofs_per_cell;
        types::global_dof_index *dofs = &accessor.dof_handler->levels[accessor.level()]
                                        ->cell_dof_indices_cache[accessor.present_index *n_dofs];

        // distribute cell matrices
        for (unsigned int i=0; i<n_dofs; ++i)
          {
            global_matrix.add(dofs[i], n_dofs, dofs, &local_matrix(i,0));
            global_vector(dofs[i]) += local_vector(i);
          }
      }



      template <int dim, int spacedim, bool level_dof_access, typename number,
                class OutputMatrix, typename OutputVector>
      static
      void
      distribute_local_to_global (const DoFCellAccessor<dealii::hp::DoFHandler<dim,spacedim>, level_dof_access> &accessor,
                                  const dealii::FullMatrix<number> &local_matrix,
                                  const dealii::Vector<number>     &local_vector,
                                  OutputMatrix                     &global_matrix,
                                  OutputVector                     &global_vector)
      {
        typedef dealii::DoFAccessor<dim,DoFHandler<dim,spacedim>, level_dof_access> BaseClass;
        Assert (accessor.dof_handler != 0,
                typename BaseClass::ExcInvalidObject());
        Assert (&accessor.get_fe() != 0,
                typename BaseClass::ExcInvalidObject());
        Assert (local_matrix.m() == accessor.get_fe().dofs_per_cell,
                typename BaseClass::ExcMatrixDoesNotMatch());
        Assert (local_matrix.n() == accessor.get_fe().dofs_per_cell,
                typename BaseClass::ExcVectorDoesNotMatch());
        Assert (accessor.dof_handler->n_dofs() == global_matrix.m(),
                typename BaseClass::ExcMatrixDoesNotMatch());
        Assert (accessor.dof_handler->n_dofs() == global_matrix.n(),
                typename BaseClass::ExcMatrixDoesNotMatch());
        Assert (local_vector.size() == accessor.get_fe().dofs_per_cell,
                typename BaseClass::ExcVectorDoesNotMatch());
        Assert (accessor.dof_handler->n_dofs() == global_vector.size(),
                typename BaseClass::ExcVectorDoesNotMatch());

        const unsigned int n_dofs = local_matrix.size();

//TODO[WB/MK]: This function could me made more efficient because it
//allocates memory, which could be avoided by passing in another
//argument as a scratch array. Comment(GK) Do not bother and leva this
//to ConstraintMatrix or MeshWorker::Assembler

        // get indices of dofs
        std::vector<types::global_dof_index> dofs (n_dofs);
        accessor.get_dof_indices (dofs);

        // distribute cell matrix and vector
        global_matrix.add(dofs,local_matrix);
        global_vector.add(dofs,local_vector);
      }
    };
  }
}


template <class DH, bool level_dof_access>
inline
DoFCellAccessor<DH,level_dof_access>::DoFCellAccessor (
  const Triangulation<DH::dimension,DH::space_dimension> *tria,
  const int                 level,
  const int                 index,
  const AccessorData       *local_data)
  :
  DoFAccessor<DH::dimension,DH,level_dof_access> (tria,level,index, local_data)
{}


template <class DH, bool level_dof_access>
template <int structdim2, int dim2, int spacedim2>
inline
DoFCellAccessor<DH,level_dof_access>::DoFCellAccessor (const InvalidAccessor<structdim2,dim2,spacedim2> &)
{
  Assert (false, typename BaseClass::ExcInvalidObject());
}



template <class DH, bool level_dof_access>
template <int dim2, class DH2, bool level_dof_access2>
inline
DoFCellAccessor<DH,level_dof_access>::DoFCellAccessor (const DoFAccessor<dim2,DH2,level_dof_access2> &other)
  :
  BaseClass(other)
{}


template <class DH, bool level_dof_access>
inline
TriaIterator<DoFCellAccessor<DH,level_dof_access> >
DoFCellAccessor<DH,level_dof_access>::neighbor (const unsigned int i) const
{
  TriaIterator<DoFCellAccessor<DH,level_dof_access> >
  q (this->tria,
     this->neighbor_level (i),
     this->neighbor_index (i),
     this->dof_handler);

#ifdef DEBUG
  if (q.state() != IteratorState::past_the_end)
    Assert (q->used(), TriaAccessorExceptions::ExcUnusedCellAsNeighbor());
#endif
  return q;
}


template <class DH, bool level_dof_access>
inline
TriaIterator<DoFCellAccessor<DH,level_dof_access> >
DoFCellAccessor<DH,level_dof_access>::child (const unsigned int i) const
{
  TriaIterator<DoFCellAccessor<DH,level_dof_access> >
  q (this->tria,
     this->level()+1,
     this->child_index (i),
     this->dof_handler);

#ifdef DEBUG
  if (q.state() != IteratorState::past_the_end)
    Assert (q->used(), TriaAccessorExceptions::ExcUnusedCellAsChild());
#endif
  return q;
}


template <class DH, bool level_dof_access>
inline
TriaIterator<DoFCellAccessor<DH,level_dof_access> >
DoFCellAccessor<DH,level_dof_access>::parent () const
{
  TriaIterator<DoFCellAccessor<DH,level_dof_access> >
  q (this->tria,
     this->level() - 1,
     this->parent_index (),
     this->dof_handler);

  return q;
}


namespace internal
{
  namespace DoFCellAccessor
  {
    template <class DH, bool level_dof_access>
    inline
    TriaIterator<dealii::DoFAccessor<DH::dimension-1,DH,level_dof_access> >
    get_face (const dealii::DoFCellAccessor<DH,level_dof_access> &cell,
              const unsigned int i,
              const dealii::internal::int2type<1>)
    {
      dealii::DoFAccessor<0, DH,level_dof_access>
      a (&cell.get_triangulation(),
         ((i == 0) && cell.at_boundary(0)
          ?
          dealii::TriaAccessor<0, 1, DH::space_dimension>::left_vertex
          :
          ((i == 1) && cell.at_boundary(1)
           ?
           dealii::TriaAccessor<0, 1, DH::space_dimension>::right_vertex
           :
           dealii::TriaAccessor<0, 1, DH::space_dimension>::interior_vertex)),
         cell.vertex_index(i),
         &cell.get_dof_handler());
      return dealii::TriaIterator<dealii::DoFAccessor<0,DH,level_dof_access> > (a);
    }


    template <class DH, bool level_dof_access>
    inline
    TriaIterator<dealii::DoFAccessor<DH::dimension-1,DH,level_dof_access> >
    get_face (const dealii::DoFCellAccessor<DH,level_dof_access> &cell,
              const unsigned int i,
              const dealii::internal::int2type<2>)
    {
      return cell.line(i);
    }


    template <class DH, bool level_dof_access>
    inline
    TriaIterator<dealii::DoFAccessor<DH::dimension-1,DH,level_dof_access> >
    get_face (const dealii::DoFCellAccessor<DH,level_dof_access> &cell,
              const unsigned int i,
              const dealii::internal::int2type<3>)
    {
      return cell.quad(i);
    }
  }
}


template <class DH, bool level_dof_access>
inline
typename DoFCellAccessor<DH,level_dof_access>::face_iterator
DoFCellAccessor<DH,level_dof_access>::face (const unsigned int i) const
{
  Assert (i<GeometryInfo<dim>::faces_per_cell, ExcIndexRange (i, 0, GeometryInfo<dim>::faces_per_cell));
  Assert (static_cast<unsigned int>(this->level()) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));

  const unsigned int dim = DH::dimension;
  return dealii::internal::DoFCellAccessor::get_face (*this, i, dealii::internal::int2type<dim>());
}



template <class DH, bool level_dof_access>
inline
void
DoFCellAccessor<DH,level_dof_access>::get_dof_indices (std::vector<types::global_dof_index> &dof_indices) const
{
  Assert (this->active(), ExcMessage ("get_dof_indices() only works on active cells."));
  Assert (this->is_artificial() == false,
          ExcMessage ("Can't ask for DoF indices on artificial cells."));
  AssertDimension (dof_indices.size(), this->get_fe().dofs_per_cell);

  const types::global_dof_index *cache
    = this->dof_handler->levels[this->present_level]
      ->get_cell_cache_start (this->present_index, this->get_fe().dofs_per_cell);
  for (unsigned int i=0; i<this->get_fe().dofs_per_cell; ++i, ++cache)
    dof_indices[i] = *cache;
}



template<class DH, bool level_dof_access>
inline
void DoFCellAccessor<DH,level_dof_access>::get_mg_dof_indices (std::vector<types::global_dof_index> &dof_indices) const
{
  DoFAccessor<dim, DH,level_dof_access>::get_mg_dof_indices (this->level (), dof_indices);
}



template<class DH, bool level_dof_access>
inline
void DoFCellAccessor<DH,level_dof_access>::set_mg_dof_indices (const std::vector<types::global_dof_index> &dof_indices)
{
  DoFAccessor<dim, DH,level_dof_access>::set_mg_dof_indices (this->level (), dof_indices);
}



template<class DH, bool level_dof_access>
inline
void DoFCellAccessor<DH,level_dof_access>::get_active_or_mg_dof_indices (std::vector<types::global_dof_index> &dof_indices) const
{
  if (level_dof_access)
    get_mg_dof_indices (dof_indices);
  else
    get_dof_indices (dof_indices);
}



template <class DH, bool level_dof_access>
template <class InputVector, typename number>
inline
void
DoFCellAccessor<DH,level_dof_access>::get_dof_values (const InputVector &values,
                                                      Vector<number>    &local_values) const
{
  get_dof_values (values, local_values.begin(), local_values.end());
}



template <class DH, bool level_dof_access>
template <class InputVector, typename ForwardIterator>
inline
void
DoFCellAccessor<DH,level_dof_access>::get_dof_values (const InputVector &values,
                                                      ForwardIterator    local_values_begin,
                                                      ForwardIterator    local_values_end) const
{
  Assert (this->is_artificial() == false,
          ExcMessage ("Can't ask for DoF indices on artificial cells."));
  Assert (!this->has_children(),
          ExcMessage ("Cell must be active."));

  Assert (static_cast<unsigned int>(local_values_end-local_values_begin)
          == this->get_fe().dofs_per_cell,
          typename DoFCellAccessor::ExcVectorDoesNotMatch());
  Assert (values.size() == this->get_dof_handler().n_dofs(),
          typename DoFCellAccessor::ExcVectorDoesNotMatch());

  const types::global_dof_index *cache
    = this->dof_handler->levels[this->present_level]
      ->get_cell_cache_start (this->present_index, this->get_fe().dofs_per_cell);

  values.extract_subvector_to (cache,
                               cache + this->get_fe().dofs_per_cell,
                               local_values_begin);
}



template <class DH, bool level_dof_access>
template <class InputVector, typename ForwardIterator>
inline
void
DoFCellAccessor<DH,level_dof_access>::get_dof_values (const ConstraintMatrix &constraints,
                                                      const InputVector      &values,
                                                      ForwardIterator         local_values_begin,
                                                      ForwardIterator         local_values_end) const
{
  Assert (this->is_artificial() == false,
          ExcMessage ("Can't ask for DoF indices on artificial cells."));
  Assert (!this->has_children(),
          ExcMessage ("Cell must be active."));

  Assert (static_cast<unsigned int>(local_values_end-local_values_begin)
          == this->get_fe().dofs_per_cell,
          typename DoFCellAccessor::ExcVectorDoesNotMatch());
  Assert (values.size() == this->get_dof_handler().n_dofs(),
          typename DoFCellAccessor::ExcVectorDoesNotMatch());


  const types::global_dof_index *cache
    = this->dof_handler->levels[this->present_level]
      ->get_cell_cache_start (this->present_index, this->get_fe().dofs_per_cell);

  constraints.get_dof_values(values, *cache, local_values_begin,
                             local_values_end);
}



template <class DH, bool level_dof_access>
template <class OutputVector, typename number>
inline
void
DoFCellAccessor<DH,level_dof_access>::set_dof_values (const Vector<number> &local_values,
                                                      OutputVector         &values) const
{
  Assert (this->is_artificial() == false,
          ExcMessage ("Can't ask for DoF indices on artificial cells."));
  Assert (!this->has_children(),
          ExcMessage ("Cell must be active."));

  Assert (static_cast<unsigned int>(local_values.size())
          == this->get_fe().dofs_per_cell,
          typename DoFCellAccessor::ExcVectorDoesNotMatch());
  Assert (values.size() == this->get_dof_handler().n_dofs(),
          typename DoFCellAccessor::ExcVectorDoesNotMatch());


  const types::global_dof_index *cache
    = this->dof_handler->levels[this->present_level]
      ->get_cell_cache_start (this->present_index, this->get_fe().dofs_per_cell);

  for (unsigned int i=0; i<this->get_fe().dofs_per_cell; ++i, ++cache)
    values(*cache) = local_values(i);
}



template <class DH, bool level_dof_access>
inline
const FiniteElement<DH::dimension,DH::space_dimension> &
DoFCellAccessor<DH,level_dof_access>::get_fe () const
{
  Assert ((dynamic_cast<const dealii::DoFHandler<DH::dimension,DH::space_dimension>*>
           (this->dof_handler) != 0)
          ||
          (this->has_children() == false),
          ExcMessage ("In hp::DoFHandler objects, finite elements are only associated "
                      "with active cells. Consequently, you can not ask for the "
                      "active finite element on cells with children."));
  return dealii::internal::DoFAccessor::get_fe (this->dof_handler->get_fe(), active_fe_index());
}



template <class DH, bool level_dof_access>
inline
unsigned int
DoFCellAccessor<DH,level_dof_access>::active_fe_index () const
{
  Assert ((dynamic_cast<const dealii::DoFHandler<DH::dimension,DH::space_dimension>*>
           (this->dof_handler) != 0)
          ||
          (this->has_children() == false),
          ExcMessage ("You can not ask for the active_fe_index on a cell that has "
                      "children because no degrees of freedom are assigned "
                      "to this cell and, consequently, no finite element "
                      "is associated with it."));
  return dealii::internal::DoFCellAccessor::Implementation::active_fe_index (*this);
}



template <class DH, bool level_dof_access>
inline
void
DoFCellAccessor<DH,level_dof_access>::set_active_fe_index (const unsigned int i)
{
  Assert ((dynamic_cast<const dealii::DoFHandler<DH::dimension,DH::space_dimension>*>
           (this->dof_handler) != 0)
          ||
          (this->has_children() == false),
          ExcMessage ("You can not set the active_fe_index on a cell that has "
                      "children because no degrees of freedom will be assigned "
                      "to this cell."));
  dealii::internal::DoFCellAccessor::Implementation::set_active_fe_index (*this, i);
}




template <class DH, bool level_dof_access>
template <typename number, typename OutputVector>
inline
void
DoFCellAccessor<DH,level_dof_access>::distribute_local_to_global (
  const Vector<number> &local_source,
  OutputVector         &global_destination) const
{
  dealii::internal::DoFCellAccessor::Implementation::
  distribute_local_to_global (*this, local_source.begin(),
                              local_source.end(), global_destination);
}



template <class DH, bool level_dof_access>
template <typename ForwardIterator, typename OutputVector>
inline
void
DoFCellAccessor<DH,level_dof_access>::distribute_local_to_global (
  ForwardIterator  local_source_begin,
  ForwardIterator  local_source_end,
  OutputVector    &global_destination) const
{
  dealii::internal::DoFCellAccessor::Implementation::
  distribute_local_to_global (*this, local_source_begin,
                              local_source_end, global_destination);
}



template <class DH, bool level_dof_access>
template <typename ForwardIterator, typename OutputVector>
inline
void
DoFCellAccessor<DH,level_dof_access>::distribute_local_to_global (
  const ConstraintMatrix &constraints,
  ForwardIterator         local_source_begin,
  ForwardIterator         local_source_end,
  OutputVector           &global_destination) const
{
  dealii::internal::DoFCellAccessor::Implementation::
  distribute_local_to_global (*this, constraints, local_source_begin,
                              local_source_end, global_destination);
}



template <class DH, bool level_dof_access>
template <typename number, typename OutputMatrix>
inline
void
DoFCellAccessor<DH,level_dof_access>::distribute_local_to_global (
  const FullMatrix<number> &local_source,
  OutputMatrix             &global_destination) const
{
  dealii::internal::DoFCellAccessor::Implementation::
  distribute_local_to_global (*this,local_source,global_destination);
}



template <class DH, bool level_dof_access>
template <typename number, typename OutputMatrix, typename OutputVector>
inline
void
DoFCellAccessor<DH,level_dof_access>::distribute_local_to_global (
  const FullMatrix<number> &local_matrix,
  const Vector<number>     &local_vector,
  OutputMatrix             &global_matrix,
  OutputVector             &global_vector) const
{
  dealii::internal::DoFCellAccessor::Implementation::
  distribute_local_to_global (*this,local_matrix,local_vector,
                              global_matrix,global_vector);
}


DEAL_II_NAMESPACE_CLOSE

#endif
