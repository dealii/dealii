// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2018 by the deal.II authors
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

#include <deal.II/base/exceptions.h>
#include <deal.II/base/geometric_utilities.h>

DEAL_II_NAMESPACE_OPEN

namespace GeometricUtilities
{
  namespace Coordinates
  {
    DeclException1(NegativeRadius,
                   double,
                   << "The radius <" << arg1 << "> can not be negative.");

    DeclException1(SphericalAzimuth,
                   double,
                   << "The azimuth angle <" << arg1 << "> is not in [0,2Pi).");

    DeclException1(SphericalPolar,
                   double,
                   << "The polar angle <" << arg1 << "> is not in [0,Pi].");

    template <int dim>
    std::array<double, dim>
    to_spherical(const Point<dim>& position)
    {
      std::array<double, dim> scoord;

      // radius
      scoord[0] = position.norm();
      // azimuth angle \theta:
      scoord[1] = std::atan2(position(1), position(0));
      // correct to [0,2*pi)
      if(scoord[1] < 0.0)
        scoord[1] += 2.0 * numbers::PI;

      // polar angle \phi:
      if(dim == 3)
        {
          // acos returns the angle in the range [0,\pi]
          if(scoord[0] > std::numeric_limits<double>::min())
            scoord[2] = std::acos(position(2) / scoord[0]);
          else
            scoord[2] = 0.0;
        }
      return scoord;
    }

    template <std::size_t dim>
    Point<dim>
    from_spherical(const std::array<double, dim>& scoord)
    {
      Point<dim> ccoord;

      Assert(scoord[0] >= 0., NegativeRadius(scoord[0]));

      Assert(scoord[1] >= 0. && scoord[1] < 2. * numbers::PI,
             SphericalAzimuth(scoord[1]));

      switch(dim)
        {
          case 2:
            {
              ccoord[0] = scoord[0] * std::cos(scoord[1]);
              ccoord[1] = scoord[0] * std::sin(scoord[1]);
              break;
            }
          case 3:
            {
              Assert(scoord[2] >= 0. && scoord[2] <= numbers::PI,
                     SphericalPolar(scoord[2]));

              ccoord[0] = scoord[0] * std::sin(scoord[2]) * std::cos(scoord[1]);
              ccoord[1] = scoord[0] * std::sin(scoord[2]) * std::sin(scoord[1]);
              ccoord[2] = scoord[0] * std::cos(scoord[2]);
              break;
            }
          default:
            Assert(false, ExcNotImplemented());
            break;
        }

      return ccoord;
    }

    // explicit instantiations
#include "geometric_utilities.inst"

  } // namespace Coordinates
} // namespace GeometricUtilities

DEAL_II_NAMESPACE_CLOSE
