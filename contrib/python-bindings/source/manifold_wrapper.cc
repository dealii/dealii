// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2017 by the deal.II authors
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

#include <deal.II/grid/manifold_lib.h>

#include <manifold_wrapper.h>
#include <point_wrapper.h>

DEAL_II_NAMESPACE_OPEN

namespace python
{
  namespace internal
  {
    template <int dim, int spacedim>
    Manifold<dim, spacedim> *
    create_spherical_manifold(const PointWrapper &center)
    {
      const Point<spacedim> *point =
        static_cast<const Point<spacedim> *>(center.get_point());
      return new SphericalManifold<dim, spacedim>(*point);
    }



    template <int dim, int spacedim>
    Manifold<dim, spacedim> *
    create_polar_manifold(const PointWrapper &center)
    {
      const Point<spacedim> *point =
        static_cast<const Point<spacedim> *>(center.get_point());
      return new PolarManifold<dim, spacedim>(*point);
    }



    template <int dim, int spacedim>
    Manifold<dim, spacedim> *
    create_cylindrical_manifold(const int axis, const double tolerance)
    {
      return new CylindricalManifold<dim, spacedim>(axis, tolerance);
    }



    template <int dim, int spacedim>
    Manifold<dim, spacedim> *
    clone(void *manifold_ptr)
    {
      const Manifold<dim, spacedim> *manifold =
        static_cast<const Manifold<dim, spacedim> *>(manifold_ptr);

      if (const SphericalManifold<dim, spacedim> *d =
            dynamic_cast<const SphericalManifold<dim, spacedim> *>(manifold))
        {
          return new SphericalManifold<dim, spacedim>(*d);
        }
      else if (const PolarManifold<dim, spacedim> *d =
                 dynamic_cast<const PolarManifold<dim, spacedim> *>(manifold))
        {
          return new PolarManifold<dim, spacedim>(*d);
        }
      else if (const CylindricalManifold<dim, spacedim> *d =
                 dynamic_cast<const CylindricalManifold<dim, spacedim> *>(
                   manifold))
        {
          return new CylindricalManifold<dim, spacedim>(*d);
        }
      else
        ExcMessage("Unsupported manifold type in clone.");

      return nullptr;
    }
  } // namespace internal



  ManifoldWrapper::ManifoldWrapper(const int dim, const int spacedim)
    : dim(dim)
    , spacedim(spacedim)
  {
    AssertThrow(((dim == 2) && (spacedim == 2)) ||
                  ((dim == 2) && (spacedim == 3)) ||
                  ((dim == 3) && (spacedim == 3)),
                ExcMessage("Wrong dim-spacedim combination."));
  }



  ManifoldWrapper::ManifoldWrapper(const ManifoldWrapper &other)
  {
    dim      = other.dim;
    spacedim = other.spacedim;

    AssertThrow(other.manifold_ptr != nullptr,
                ExcMessage("Underlying manifold does not exist."));

    if ((dim == 2) && (spacedim == 2))
      {
        manifold_ptr = internal::clone<2, 2>(other.manifold_ptr);
      }
    else if ((dim == 2) && (spacedim == 3))
      {
        manifold_ptr = internal::clone<2, 3>(other.manifold_ptr);
      }
    else if ((dim == 3) && (spacedim == 3))
      {
        manifold_ptr = internal::clone<3, 3>(other.manifold_ptr);
      }
    else
      AssertThrow(false, ExcMessage("Wrong dim-spacedim combination."));
  }



  ManifoldWrapper::~ManifoldWrapper()
  {
    if (dim != -1)
      {
        if ((dim == 2) && (spacedim == 2))
          {
            // We cannot call delete on a void pointer so cast the void pointer
            // back first.
            Manifold<2, 2> *tmp = static_cast<Manifold<2, 2> *>(manifold_ptr);
            delete tmp;
          }
        else if ((dim == 2) && (spacedim == 3))
          {
            Manifold<2, 3> *tmp = static_cast<Manifold<2, 3> *>(manifold_ptr);
            delete tmp;
          }
        else
          {
            Manifold<3, 3> *tmp = static_cast<Manifold<3, 3> *>(manifold_ptr);
            delete tmp;
          }

        dim          = -1;
        spacedim     = -1;
        manifold_ptr = nullptr;
      }
  }



  void
  ManifoldWrapper::create_spherical(const PointWrapper center)
  {
    if ((dim == 2) && (spacedim == 2))
      {
        manifold_ptr = internal::create_spherical_manifold<2, 2>(center);
      }
    else if ((dim == 2) && (spacedim == 3))
      {
        manifold_ptr = internal::create_spherical_manifold<2, 3>(center);
      }
    else if ((dim == 3) && (spacedim == 3))
      {
        manifold_ptr = internal::create_spherical_manifold<3, 3>(center);
      }
    else
      AssertThrow(false, ExcMessage("Wrong dim-spacedim combination."));
  }



  void
  ManifoldWrapper::create_polar(const PointWrapper center)
  {
    if ((dim == 2) && (spacedim == 2))
      {
        manifold_ptr = internal::create_polar_manifold<2, 2>(center);
      }
    else if ((dim == 2) && (spacedim == 3))
      {
        manifold_ptr = internal::create_polar_manifold<2, 3>(center);
      }
    else if ((dim == 3) && (spacedim == 3))
      {
        manifold_ptr = internal::create_polar_manifold<3, 3>(center);
      }
    else
      AssertThrow(false, ExcMessage("Wrong dim-spacedim combination."));
  }



  void
  ManifoldWrapper::create_cylindrical(const int axis, const double tolerance)
  {
    if ((dim == 2) && (spacedim == 3))
      {
        manifold_ptr =
          internal::create_cylindrical_manifold<2, 3>(axis, tolerance);
      }
    else if ((dim == 3) && (spacedim == 3))
      {
        manifold_ptr =
          internal::create_cylindrical_manifold<3, 3>(axis, tolerance);
      }
    else
      AssertThrow(false, ExcMessage("Wrong dim-spacedim combination."));
  }



  void *
  ManifoldWrapper::get_manifold() const
  {
    return manifold_ptr;
  }



  int
  ManifoldWrapper::get_dim() const
  {
    return dim;
  }



  int
  ManifoldWrapper::get_spacedim() const
  {
    return spacedim;
  }


} // namespace python

DEAL_II_NAMESPACE_CLOSE
