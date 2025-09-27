// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/grid/manifold_lib.h>

#include <manifold_wrapper.h>
#include <point_wrapper.h>

#include <memory>

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
    create_cylindrical_manifold(const boost::python::list &direction_list,
                                const boost::python::list &axial_point_list)
    {
      Tensor<1, spacedim> direction;
      for (int d = 0; d < spacedim; ++d)
        direction[d] = boost::python::extract<double>(direction_list[d]);

      Point<spacedim> axial_point;
      for (int d = 0; d < spacedim; ++d)
        axial_point[d] = boost::python::extract<double>(axial_point_list[d]);

      return new CylindricalManifold<dim, spacedim>(direction, axial_point);
    }



    template <int dim, int spacedim>
    Manifold<dim, spacedim> *
    create_function_manifold(const std::string &push_forward,
                             const std::string &pull_back)
    {
      return new FunctionManifold<dim, spacedim>(push_forward, pull_back);
    }



    template <int dim, int spacedim>
    Manifold<dim, spacedim> *
    create_function_manifold(boost::python::object &push_forward,
                             boost::python::object &pull_back)
    {
      return new FunctionManifold<dim, spacedim>(
        std::make_unique<FunctionWrapper<dim>>(push_forward, spacedim),
        std::make_unique<FunctionWrapper<spacedim>>(pull_back, dim));
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
      else if (const FunctionManifold<dim, spacedim> *d =
                 dynamic_cast<const FunctionManifold<dim, spacedim> *>(
                   manifold))
        {
          return new FunctionManifold<dim, spacedim>(*d);
        }
      else
        ExcMessage("Unsupported manifold type in clone.");

      return nullptr;
    }
  } // namespace internal



  ManifoldWrapper::ManifoldWrapper(const int dim, const int spacedim)
    : dim(dim)
    , spacedim(spacedim)
    , manifold_ptr(nullptr)
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
      AssertThrow(false,
                  ExcMessage(
                    "Given dim-spacedim combination is not implemented."));
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
      AssertThrow(false,
                  ExcMessage(
                    "Given dim-spacedim combination is not implemented."));
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
      AssertThrow(false,
                  ExcMessage(
                    "Given dim-spacedim combination is not implemented."));
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
      AssertThrow(false,
                  ExcMessage(
                    "Given dim-spacedim combination is not implemented."));
  }



  void
  ManifoldWrapper::create_cylindrical(const boost::python::list &direction,
                                      const boost::python::list &axial_point)
  {
    if ((dim == 2) && (spacedim == 3))
      {
        manifold_ptr =
          internal::create_cylindrical_manifold<2, 3>(direction, axial_point);
      }
    else if ((dim == 3) && (spacedim == 3))
      {
        manifold_ptr =
          internal::create_cylindrical_manifold<3, 3>(direction, axial_point);
      }
    else
      AssertThrow(false,
                  ExcMessage(
                    "Given dim-spacedim combination is not implemented."));
  }



  void
  ManifoldWrapper::create_function_string(const std::string &push_forward,
                                          const std::string &pull_back)
  {
    if ((dim == 2) && (spacedim == 2))
      {
        manifold_ptr =
          internal::create_function_manifold<2, 2>(push_forward, pull_back);
      }
    else if ((dim == 2) && (spacedim == 3))
      {
        manifold_ptr =
          internal::create_function_manifold<2, 3>(push_forward, pull_back);
      }
    else if ((dim == 3) && (spacedim == 3))
      {
        manifold_ptr =
          internal::create_function_manifold<3, 3>(push_forward, pull_back);
      }
    else
      AssertThrow(false,
                  ExcMessage(
                    "Given dim-spacedim combination is not implemented."));
  }



  void
  ManifoldWrapper::create_function(boost::python::object &push_forward,
                                   boost::python::object &pull_back)
  {
    if ((dim == 2) && (spacedim == 2))
      {
        manifold_ptr =
          internal::create_function_manifold<2, 2>(push_forward, pull_back);
      }
    else if ((dim == 2) && (spacedim == 3))
      {
        manifold_ptr =
          internal::create_function_manifold<2, 3>(push_forward, pull_back);
      }
    else if ((dim == 3) && (spacedim == 3))
      {
        manifold_ptr =
          internal::create_function_manifold<3, 3>(push_forward, pull_back);
      }
    else
      AssertThrow(false,
                  ExcMessage(
                    "Given dim-spacedim combination is not implemented."));
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
