// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2021 by the deal.II authors
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


#include <deal.II/boost_adaptors/bounding_box.h>

#include <deal.II/fe/mapping.h>

#include <deal.II/grid/tria.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#ifdef DEAL_II_BOOST_HAS_BROKEN_HEADER_DEPRECATIONS
#  define BOOST_ALLOW_DEPRECATED_HEADERS
#endif
#include <boost/geometry.hpp>
#ifdef DEAL_II_BOOST_HAS_BROKEN_HEADER_DEPRECATIONS
#  undef BOOST_ALLOW_DEPRECATED_HEADERS
#endif
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

DEAL_II_NAMESPACE_OPEN


template <int dim, int spacedim>
boost::container::small_vector<Point<spacedim>,
                               GeometryInfo<dim>::vertices_per_cell>
Mapping<dim, spacedim>::get_vertices(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell) const
{
  boost::container::small_vector<Point<spacedim>,
                                 GeometryInfo<dim>::vertices_per_cell>
    vertices;
  for (const unsigned int i : cell->vertex_indices())
    vertices.push_back(cell->vertex(i));

  return vertices;
}



template <int dim, int spacedim>
Point<spacedim>
Mapping<dim, spacedim>::get_center(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const bool map_center_of_reference_cell) const
{
  if (map_center_of_reference_cell)
    {
      Point<dim> reference_center;
      for (unsigned int d = 0; d < dim; ++d)
        reference_center[d] = .5;
      return transform_unit_to_real_cell(cell, reference_center);
    }
  else
    {
      const auto      vertices = get_vertices(cell);
      Point<spacedim> center;
      for (const auto &v : vertices)
        center += v;
      return center / GeometryInfo<dim>::vertices_per_cell;
    }
}



template <int dim, int spacedim>
BoundingBox<spacedim>
Mapping<dim, spacedim>::get_bounding_box(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell) const
{
  if (preserves_vertex_locations())
    return cell->bounding_box();
  else
    return BoundingBox<spacedim>(get_vertices(cell));
}



template <int dim, int spacedim>
void
Mapping<dim, spacedim>::fill_fe_immersed_surface_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &,
  const NonMatching::ImmersedSurfaceQuadrature<dim> &,
  const typename Mapping<dim, spacedim>::InternalDataBase &,
  dealii::internal::FEValuesImplementation::MappingRelatedData<dim, spacedim> &)
  const
{
  AssertThrow(false, ExcNotImplemented());
}



template <int dim, int spacedim>
void
Mapping<dim, spacedim>::transform_points_real_to_unit_cell(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const ArrayView<const Point<spacedim>> &                    real_points,
  const ArrayView<Point<dim>> &                               unit_points) const
{
  AssertDimension(real_points.size(), unit_points.size());
  for (unsigned int i = 0; i < real_points.size(); ++i)
    {
      try
        {
          unit_points[i] = transform_real_to_unit_cell(cell, real_points[i]);
        }
      catch (typename Mapping<dim>::ExcTransformationFailed &)
        {
          unit_points[i]    = Point<dim>();
          unit_points[i][0] = std::numeric_limits<double>::infinity();
        }
    }
}



template <int dim, int spacedim>
Point<dim - 1>
Mapping<dim, spacedim>::project_real_point_to_unit_point_on_face(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const unsigned int                                          face_no,
  const Point<spacedim> &                                     p) const
{
  // The function doesn't make physical sense for dim=1
  Assert(dim > 1, ExcNotImplemented());
  // Not implemented for higher dimensions
  Assert(dim <= 3, ExcNotImplemented());

  Point<dim> unit_cell_pt = transform_real_to_unit_cell(cell, p);

  const unsigned int unit_normal_direction =
    GeometryInfo<dim>::unit_normal_direction[face_no];

  if (dim == 2)
    {
      if (unit_normal_direction == 0)
        return Point<dim - 1>{unit_cell_pt(1)};
      else if (unit_normal_direction == 1)
        return Point<dim - 1>{unit_cell_pt(0)};
    }
  else if (dim == 3)
    {
      if (unit_normal_direction == 0)
        return Point<dim - 1>{unit_cell_pt(1), unit_cell_pt(2)};
      else if (unit_normal_direction == 1)
        return Point<dim - 1>{unit_cell_pt(0), unit_cell_pt(2)};
      else if (unit_normal_direction == 2)
        return Point<dim - 1>{unit_cell_pt(0), unit_cell_pt(1)};
    }

  // We should never get here
  Assert(false, ExcInternalError());
  return {};
}



#ifndef DOXYGEN
template <int dim, int spacedim>
void
Mapping<dim, spacedim>::fill_fe_face_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const unsigned int                                          face_no,
  const hp::QCollection<dim - 1> &                            quadrature,
  const typename Mapping<dim, spacedim>::InternalDataBase &   internal_data,
  dealii::internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
    &output_data) const
{
  // base class version, implement overridden function in derived classes
  AssertDimension(quadrature.size(), 1);
  fill_fe_face_values(cell, face_no, quadrature[0], internal_data, output_data);
}



template <int dim, int spacedim>
void
Mapping<dim, spacedim>::fill_fe_face_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const unsigned int                                          face_no,
  const Quadrature<dim - 1> &                                 quadrature,
  const typename Mapping<dim, spacedim>::InternalDataBase &   internal_data,
  dealii::internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
    &output_data) const
{
  Assert(false,
         ExcMessage("Use of a deprecated interface, please implement "
                    "fill_fe_face_values taking a hp::QCollection argument"));
  (void)cell;
  (void)face_no;
  (void)quadrature;
  (void)internal_data;
  (void)output_data;
}



template <int dim, int spacedim>
std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
Mapping<dim, spacedim>::get_face_data(
  const UpdateFlags               update_flags,
  const hp::QCollection<dim - 1> &quadrature) const
{
  // base class version, implement overridden function in derived classes
  return get_face_data(update_flags, quadrature[0]);
}



template <int dim, int spacedim>
std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
Mapping<dim, spacedim>::get_face_data(
  const UpdateFlags          update_flags,
  const Quadrature<dim - 1> &quadrature) const
{
  Assert(false,
         ExcMessage("Use of a deprecated interface, please implement "
                    "fill_fe_face_values taking a hp::QCollection argument"));
  (void)update_flags;
  (void)quadrature;

  return std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>();
}
#endif

/* ---------------------------- InternalDataBase --------------------------- */


template <int dim, int spacedim>
Mapping<dim, spacedim>::InternalDataBase::InternalDataBase()
  : update_each(update_default)
{}



template <int dim, int spacedim>
std::size_t
Mapping<dim, spacedim>::InternalDataBase::memory_consumption() const
{
  return sizeof(*this);
}


/* ------------------------------ Global functions ------------------------- */

template <int dim, int spacedim>
const Mapping<dim, spacedim> &
get_default_linear_mapping(const Triangulation<dim, spacedim> &triangulation)
{
  const auto &reference_cells = triangulation.get_reference_cells();
  Assert(reference_cells.size() == 1,
         ExcMessage(
           "This function can only work for triangulations that "
           "use only a single cell type -- for example, only triangles "
           "or only quadrilaterals. For mixed meshes, there is no "
           "single linear mapping object that can be used for all "
           "cells of the triangulation. The triangulation you are "
           "passing to this function uses multiple cell types."));

  return reference_cells.front()
    .template get_default_linear_mapping<dim, spacedim>();
}



// explicit instantiations
#include "mapping.inst"


DEAL_II_NAMESPACE_CLOSE
