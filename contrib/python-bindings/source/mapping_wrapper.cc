// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2020 by the deal.II authors
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

#include <cell_accessor_wrapper.h>
#include <mapping_wrapper.h>
#include <point_wrapper.h>
#include <triangulation_wrapper.h>

DEAL_II_NAMESPACE_OPEN

namespace python
{
  namespace internal
  {
    template <int dim, int spacedim>
    PointWrapper
    transform_unit_to_real_cell(void *mapping_ptr,
                                void *cell_accessor_ptr,
                                void *point_ptr)
    {
      const MappingQGeneric<dim, spacedim> *mapping =
        static_cast<const MappingQGeneric<dim, spacedim> *>(mapping_ptr);

      const CellAccessor<dim, spacedim> *cell_accessor =
        static_cast<const CellAccessor<dim, spacedim> *>(cell_accessor_ptr);

      const Point<dim> *point = static_cast<const Point<dim> *>(point_ptr);

      typename Triangulation<dim, spacedim>::active_cell_iterator cell(
        &cell_accessor->get_triangulation(),
        cell_accessor->level(),
        cell_accessor->index());

      Point<spacedim> p_real =
        mapping->transform_unit_to_real_cell(cell, *point);

      boost::python::list coord_list;
      for (int i = 0; i < spacedim; ++i)
        coord_list.append(p_real[i]);

      return PointWrapper(coord_list);
    }



    template <int dim, int spacedim>
    PointWrapper
    transform_real_to_unit_cell(void *mapping_ptr,
                                void *cell_accessor_ptr,
                                void *point_ptr)
    {
      const MappingQGeneric<dim, spacedim> *mapping =
        static_cast<const MappingQGeneric<dim, spacedim> *>(mapping_ptr);

      const CellAccessor<dim, spacedim> *cell_accessor =
        static_cast<const CellAccessor<dim, spacedim> *>(cell_accessor_ptr);

      const Point<spacedim> *point =
        static_cast<const Point<spacedim> *>(point_ptr);

      typename Triangulation<dim, spacedim>::active_cell_iterator cell(
        &cell_accessor->get_triangulation(),
        cell_accessor->level(),
        cell_accessor->index());

      Point<dim> p_real = mapping->transform_real_to_unit_cell(cell, *point);

      boost::python::list coord_list;
      for (int i = 0; i < dim; ++i)
        coord_list.append(p_real[i]);

      return PointWrapper(coord_list);
    }



    template <int dim, int spacedim>
    PointWrapper
    project_real_point_to_unit_point_on_face(void *mapping_ptr,
                                             void *cell_accessor_ptr,
                                             const unsigned int face_no,
                                             void *             point_ptr)
    {
      const MappingQGeneric<dim, spacedim> *mapping =
        static_cast<const MappingQGeneric<dim, spacedim> *>(mapping_ptr);

      const CellAccessor<dim, spacedim> *cell_accessor =
        static_cast<const CellAccessor<dim, spacedim> *>(cell_accessor_ptr);

      const Point<spacedim> *point =
        static_cast<const Point<spacedim> *>(point_ptr);

      typename Triangulation<dim, spacedim>::active_cell_iterator cell(
        &cell_accessor->get_triangulation(),
        cell_accessor->level(),
        cell_accessor->index());

      Point<dim - 1> p_face =
        mapping->project_real_point_to_unit_point_on_face(cell,
                                                          face_no,
                                                          *point);

      boost::python::list coord_list;
      for (int i = 0; i < dim - 1; ++i)
        coord_list.append(p_face[i]);
      coord_list.append(0.);

      return PointWrapper(coord_list);
    }

  } // namespace internal



  MappingQGenericWrapper::MappingQGenericWrapper()
    : dim(-1)
    , spacedim(-1)
    , degree(-1)
    , mapping_ptr(nullptr)
  {}



  MappingQGenericWrapper::MappingQGenericWrapper(const int dim,
                                                 const int spacedim,
                                                 const int degree)
    : dim(dim)
    , spacedim(spacedim)
    , degree(degree)
  {
    if ((dim == 2) && (spacedim == 2))
      {
        mapping_ptr = new MappingQGeneric<2, 2>(degree);
      }
    else if ((dim == 2) && (spacedim == 3))
      {
        mapping_ptr = new MappingQGeneric<2, 3>(degree);
      }
    else if ((dim == 3) && (spacedim == 3))
      {
        mapping_ptr = new MappingQGeneric<3, 3>(degree);
      }
    else
      AssertThrow(false, ExcMessage("Wrong dim-spacedim combination."));
  }



  MappingQGenericWrapper::MappingQGenericWrapper(
    const MappingQGenericWrapper &other)
  {
    dim      = other.dim;
    spacedim = other.spacedim;
    degree   = other.degree;

    AssertThrow(other.mapping_ptr != nullptr,
                ExcMessage("Underlying mapping does not exist."));

    if ((dim == 2) && (spacedim == 2))
      {
        mapping_ptr = new MappingQGeneric<2, 2>(other.degree);
      }
    else if ((dim == 2) && (spacedim == 3))
      {
        mapping_ptr = new MappingQGeneric<2, 3>(other.degree);
      }
    else if ((dim == 3) && (spacedim == 3))
      {
        mapping_ptr = new MappingQGeneric<3, 3>(other.degree);
      }
    else
      AssertThrow(false, ExcMessage("Wrong dim-spacedim combination."));
  }



  MappingQGenericWrapper::~MappingQGenericWrapper()
  {
    if (dim != -1)
      {
        if ((dim == 2) && (spacedim == 2))
          {
            // We cannot call delete on a void pointer so cast the void pointer
            // back first.
            MappingQGeneric<2, 2> *tmp =
              static_cast<MappingQGeneric<2, 2> *>(mapping_ptr);
            delete tmp;
          }
        else if ((dim == 2) && (spacedim == 3))
          {
            MappingQGeneric<2, 3> *tmp =
              static_cast<MappingQGeneric<2, 3> *>(mapping_ptr);
            delete tmp;
          }
        else
          {
            MappingQGeneric<3, 3> *tmp =
              static_cast<MappingQGeneric<3, 3> *>(mapping_ptr);
            delete tmp;
          }

        dim         = -1;
        spacedim    = -1;
        degree      = -1;
        mapping_ptr = nullptr;
      }
  }



  PointWrapper
  MappingQGenericWrapper::transform_unit_to_real_cell(CellAccessorWrapper &cell,
                                                      PointWrapper &       p)
  {
    AssertThrow(
      dim == p.get_dim(),
      ExcMessage(
        "Dimension of the point is not equal to the dimension of the mapping."));

    if ((dim == 2) && (spacedim == 2))
      return internal::transform_unit_to_real_cell<2, 2>(mapping_ptr,
                                                         cell.cell_accessor,
                                                         p.point);
    else if ((dim == 2) && (spacedim == 3))
      return internal::transform_unit_to_real_cell<2, 3>(mapping_ptr,
                                                         cell.cell_accessor,
                                                         p.point);
    else
      return internal::transform_unit_to_real_cell<3, 3>(mapping_ptr,
                                                         cell.cell_accessor,
                                                         p.point);
  }



  PointWrapper
  MappingQGenericWrapper::transform_real_to_unit_cell(CellAccessorWrapper &cell,
                                                      PointWrapper &       p)
  {
    AssertThrow(
      spacedim == p.get_dim(),
      ExcMessage(
        "Dimension of the point is not equal to the space dimension of the mapping."));

    if ((dim == 2) && (spacedim == 2))
      return internal::transform_real_to_unit_cell<2, 2>(mapping_ptr,
                                                         cell.cell_accessor,
                                                         p.point);
    else if ((dim == 2) && (spacedim == 3))
      return internal::transform_real_to_unit_cell<2, 3>(mapping_ptr,
                                                         cell.cell_accessor,
                                                         p.point);
    else
      return internal::transform_real_to_unit_cell<3, 3>(mapping_ptr,
                                                         cell.cell_accessor,
                                                         p.point);
  }



  PointWrapper
  MappingQGenericWrapper::project_real_point_to_unit_point_on_face(
    CellAccessorWrapper &cell,
    const unsigned int   face_no,
    PointWrapper &       p)
  {
    AssertThrow(
      spacedim == p.get_dim(),
      ExcMessage(
        "Dimension of the point is not equal to the space dimension of the mapping."));

    if ((dim == 2) && (spacedim == 2))
      return internal::project_real_point_to_unit_point_on_face<2, 2>(
        mapping_ptr, cell.cell_accessor, face_no, p.point);
    else if ((dim == 2) && (spacedim == 3))
      return internal::project_real_point_to_unit_point_on_face<2, 3>(
        mapping_ptr, cell.cell_accessor, face_no, p.point);
    else
      return internal::project_real_point_to_unit_point_on_face<3, 3>(
        mapping_ptr, cell.cell_accessor, face_no, p.point);
  }



  void *
  MappingQGenericWrapper::get_mapping() const
  {
    return mapping_ptr;
  }

} // namespace python

DEAL_II_NAMESPACE_CLOSE
