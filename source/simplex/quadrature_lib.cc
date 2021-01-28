// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/polynomial.h>

#include <deal.II/simplex/quadrature_lib.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>


DEAL_II_NAMESPACE_OPEN

namespace Simplex
{
  template <int dim>
  QGauss<dim>::QGauss(const unsigned int n_points_1D)
    : QSimplex<dim>(Quadrature<dim>())
  {
    // fill quadrature points and quadrature weights
    if (dim == 1)
      {
        const dealii::QGauss<dim> quad(n_points_1D);

        this->quadrature_points = quad.get_points();
        this->weights           = quad.get_weights();
      }
    else if (dim == 2)
      {
        if (n_points_1D == 1)
          {
            const double p = 1.0 / 3.0;
            this->quadrature_points.emplace_back(p, p);
            this->weights.emplace_back(0.5);
          }
        else if (n_points_1D == 2)
          {
            const double Q23 = 2.0 / 3.0;
            const double Q16 = 1.0 / 6.0;

            this->quadrature_points.emplace_back(Q23, Q16);
            this->quadrature_points.emplace_back(Q16, Q23);
            this->quadrature_points.emplace_back(Q16, Q16);
            this->weights.emplace_back(Q16);
            this->weights.emplace_back(Q16);
            this->weights.emplace_back(Q16);
          }
        else if (n_points_1D == 3)
          {
            const double q12 = 0.5;

            // clang-format off
            this->quadrature_points.emplace_back(0.3333333333330, 0.3333333333330);
            this->quadrature_points.emplace_back(0.7974269853530, 0.1012865073230);
            this->quadrature_points.emplace_back(0.1012865073230, 0.7974269853530);
            this->quadrature_points.emplace_back(0.1012865073230, 0.1012865073230);
            this->quadrature_points.emplace_back(0.0597158717898, 0.4701420641050);
            this->quadrature_points.emplace_back(0.4701420641050, 0.0597158717898);
            this->quadrature_points.emplace_back(0.4701420641050, 0.4701420641050);
            // clang-format on

            this->weights.emplace_back(q12 * 0.225);
            this->weights.emplace_back(q12 * 0.125939180545);
            this->weights.emplace_back(q12 * 0.125939180545);
            this->weights.emplace_back(q12 * 0.125939180545);
            this->weights.emplace_back(q12 * 0.132394152789);
            this->weights.emplace_back(q12 * 0.132394152789);
            this->weights.emplace_back(q12 * 0.132394152789);
          }
        else if (n_points_1D == 4)
          {
            const double q12 = 0.5;

            // clang-format off
            this->quadrature_points.emplace_back(0.06493051315916486, 0.06493051315916486);
            this->quadrature_points.emplace_back(0.06493051315916486, 0.8701389736816703);
            this->quadrature_points.emplace_back(0.8701389736816703, 0.06493051315916486);
            this->quadrature_points.emplace_back(0.2845755842491703, 0.517039939069323);
            this->quadrature_points.emplace_back(0.3135591843849315, 0.04386347179237249);
            this->quadrature_points.emplace_back(0.19838447668150672, 0.2845755842491703);
            this->quadrature_points.emplace_back(0.642577343822696, 0.3135591843849315);
            this->quadrature_points.emplace_back(0.517039939069323, 0.19838447668150672);
            this->quadrature_points.emplace_back(0.04386347179237249, 0.642577343822696);
            this->quadrature_points.emplace_back(0.517039939069323, 0.2845755842491703);
            this->quadrature_points.emplace_back(0.04386347179237249, 0.3135591843849315);
            this->quadrature_points.emplace_back(0.19838447668150672, 0.517039939069323);
            this->quadrature_points.emplace_back(0.642577343822696, 0.04386347179237249);
            this->quadrature_points.emplace_back(0.2845755842491703, 0.19838447668150672);
            this->quadrature_points.emplace_back(0.3135591843849315, 0.642577343822696);
            // clang-format on

            this->weights.emplace_back(q12 * 0.053077801790232484);
            this->weights.emplace_back(q12 * 0.053077801790232484);
            this->weights.emplace_back(q12 * 0.053077801790232484);
            this->weights.emplace_back(q12 * 0.07085308369213365);
            this->weights.emplace_back(q12 * 0.069274682079417);
            this->weights.emplace_back(q12 * 0.07085308369213365);
            this->weights.emplace_back(q12 * 0.069274682079417);
            this->weights.emplace_back(q12 * 0.07085308369213365);
            this->weights.emplace_back(q12 * 0.069274682079417);
            this->weights.emplace_back(q12 * 0.07085308369213365);
            this->weights.emplace_back(q12 * 0.069274682079417);
            this->weights.emplace_back(q12 * 0.07085308369213365);
            this->weights.emplace_back(q12 * 0.069274682079417);
            this->weights.emplace_back(q12 * 0.07085308369213365);
            this->weights.emplace_back(q12 * 0.069274682079417);
          }
      }
    else if (dim == 3)
      {
        if (n_points_1D == 1)
          {
            const double Q14 = 1.0 / 4.0;
            const double Q16 = 1.0 / 6.0;

            this->quadrature_points.emplace_back(Q14, Q14, Q14);
            this->weights.emplace_back(Q16);
          }
        else if (n_points_1D == 2)
          {
            const double Q124 = 1.0 / 6.0 / 4.0;

            const double palpha = (5.0 + 3.0 * sqrt(5.0)) / 20.0;
            const double pbeta  = (5.0 - sqrt(5.0)) / 20.0;
            this->quadrature_points.emplace_back(pbeta, pbeta, pbeta);
            this->quadrature_points.emplace_back(palpha, pbeta, pbeta);
            this->quadrature_points.emplace_back(pbeta, palpha, pbeta);
            this->quadrature_points.emplace_back(pbeta, pbeta, palpha);
            this->weights.emplace_back(Q124);
            this->weights.emplace_back(Q124);
            this->weights.emplace_back(Q124);
            this->weights.emplace_back(Q124);
          }
        else if (n_points_1D == 3)
          {
            const double Q16 = 1.0 / 6.0;

            // clang-format off
            this->quadrature_points.emplace_back(0.5684305841968444, 0.1438564719343852, 0.1438564719343852);
            this->quadrature_points.emplace_back(0.1438564719343852, 0.1438564719343852, 0.1438564719343852);
            this->quadrature_points.emplace_back(0.1438564719343852, 0.1438564719343852, 0.5684305841968444);
            this->quadrature_points.emplace_back(0.1438564719343852, 0.5684305841968444, 0.1438564719343852);
            this->quadrature_points.emplace_back(0.0000000000000000, 0.5000000000000000, 0.5000000000000000);
            this->quadrature_points.emplace_back(0.5000000000000000, 0.0000000000000000, 0.5000000000000000);
            this->quadrature_points.emplace_back(0.5000000000000000, 0.5000000000000000, 0.0000000000000000);
            this->quadrature_points.emplace_back(0.5000000000000000, 0.0000000000000000, 0.0000000000000000);
            this->quadrature_points.emplace_back(0.0000000000000000, 0.5000000000000000, 0.0000000000000000);
            this->quadrature_points.emplace_back(0.0000000000000000, 0.0000000000000000, 0.5000000000000000);
            // clang-format on

            this->weights.emplace_back(0.2177650698804054 * Q16);
            this->weights.emplace_back(0.2177650698804054 * Q16);
            this->weights.emplace_back(0.2177650698804054 * Q16);
            this->weights.emplace_back(0.2177650698804054 * Q16);
            this->weights.emplace_back(0.0214899534130631 * Q16);
            this->weights.emplace_back(0.0214899534130631 * Q16);
            this->weights.emplace_back(0.0214899534130631 * Q16);
            this->weights.emplace_back(0.0214899534130631 * Q16);
            this->weights.emplace_back(0.0214899534130631 * Q16);
            this->weights.emplace_back(0.0214899534130631 * Q16);
          }
        else if (n_points_1D == 4)
          {
            const double Q16 = 1.0 / 6.0;

            // clang-format off
            this->quadrature_points.emplace_back(0.25, 0.25, 0.25);
            this->quadrature_points.emplace_back(0.3157011497782028, 0.3157011497782028, 0.3157011497782028);
            this->quadrature_points.emplace_back(0.3157011497782028, 0.3157011497782028, 0.05289655066539156);
            this->quadrature_points.emplace_back(0.3157011497782028, 0.05289655066539156, 0.3157011497782028);
            this->quadrature_points.emplace_back(0.05289655066539156, 0.3157011497782028, 0.3157011497782028);
            this->quadrature_points.emplace_back(0.05048982259839635, 0.05048982259839635, 0.44951017740160365);
            this->quadrature_points.emplace_back(0.05048982259839635, 0.44951017740160365, 0.05048982259839635);
            this->quadrature_points.emplace_back(0.44951017740160365, 0.05048982259839635, 0.05048982259839635);
            this->quadrature_points.emplace_back(0.05048982259839635, 0.44951017740160365, 0.44951017740160365);
            this->quadrature_points.emplace_back(0.44951017740160365, 0.05048982259839635, 0.44951017740160365);
            this->quadrature_points.emplace_back(0.44951017740160365, 0.44951017740160365, 0.05048982259839635);
            this->quadrature_points.emplace_back(0.188833831026001, 0.188833831026001, 0.575171637587);
            this->quadrature_points.emplace_back(0.02126547254148314, 0.02126547254148314, 0.8108302410985486);
            this->quadrature_points.emplace_back(0.188833831026001, 0.575171637587, 0.188833831026001);
            this->quadrature_points.emplace_back(0.02126547254148314, 0.8108302410985486, 0.02126547254148314);
            this->quadrature_points.emplace_back(0.575171637587, 0.188833831026001, 0.188833831026001);
            this->quadrature_points.emplace_back(0.8108302410985486, 0.02126547254148314, 0.02126547254148314);
            this->quadrature_points.emplace_back(0.188833831026001, 0.575171637587, 0.047160700360998065);
            this->quadrature_points.emplace_back(0.02126547254148314, 0.8108302410985486, 0.1466388138184851);
            this->quadrature_points.emplace_back(0.575171637587, 0.188833831026001, 0.047160700360998065);
            this->quadrature_points.emplace_back(0.8108302410985486, 0.02126547254148314, 0.1466388138184851);
            this->quadrature_points.emplace_back(0.575171637587, 0.047160700360998065, 0.188833831026001);
            this->quadrature_points.emplace_back(0.8108302410985486, 0.1466388138184851, 0.02126547254148314);
            this->quadrature_points.emplace_back(0.188833831026001, 0.188833831026001, 0.047160700360998065);
            this->quadrature_points.emplace_back(0.02126547254148314, 0.02126547254148314, 0.1466388138184851);
            this->quadrature_points.emplace_back(0.188833831026001, 0.047160700360998065, 0.188833831026001);
            this->quadrature_points.emplace_back(0.02126547254148314, 0.1466388138184851, 0.02126547254148314);
            this->quadrature_points.emplace_back(0.047160700360998065, 0.188833831026001, 0.188833831026001);
            this->quadrature_points.emplace_back(0.1466388138184851, 0.02126547254148314, 0.02126547254148314);
            this->quadrature_points.emplace_back(0.188833831026001, 0.047160700360998065, 0.575171637587);
            this->quadrature_points.emplace_back(0.02126547254148314, 0.1466388138184851, 0.8108302410985486);
            this->quadrature_points.emplace_back(0.047160700360998065, 0.188833831026001, 0.575171637587);
            this->quadrature_points.emplace_back(0.1466388138184851, 0.02126547254148314, 0.8108302410985486);
            this->quadrature_points.emplace_back(0.047160700360998065, 0.575171637587, 0.188833831026001);
            this->quadrature_points.emplace_back(0.1466388138184851, 0.8108302410985486, 0.02126547254148314);
            // clang-format on

            this->weights.emplace_back(Q16 * 0.09548528946413085);
            this->weights.emplace_back(Q16 * 0.04232958120996703);
            this->weights.emplace_back(Q16 * 0.04232958120996703);
            this->weights.emplace_back(Q16 * 0.04232958120996703);
            this->weights.emplace_back(Q16 * 0.04232958120996703);
            this->weights.emplace_back(Q16 * 0.03189692783285758);
            this->weights.emplace_back(Q16 * 0.03189692783285758);
            this->weights.emplace_back(Q16 * 0.03189692783285758);
            this->weights.emplace_back(Q16 * 0.03189692783285758);
            this->weights.emplace_back(Q16 * 0.03189692783285758);
            this->weights.emplace_back(Q16 * 0.03189692783285758);
            this->weights.emplace_back(Q16 * 0.03720713072833462);
            this->weights.emplace_back(Q16 * 0.008110770829903342);
            this->weights.emplace_back(Q16 * 0.03720713072833462);
            this->weights.emplace_back(Q16 * 0.008110770829903342);
            this->weights.emplace_back(Q16 * 0.03720713072833462);
            this->weights.emplace_back(Q16 * 0.008110770829903342);
            this->weights.emplace_back(Q16 * 0.03720713072833462);
            this->weights.emplace_back(Q16 * 0.008110770829903342);
            this->weights.emplace_back(Q16 * 0.03720713072833462);
            this->weights.emplace_back(Q16 * 0.008110770829903342);
            this->weights.emplace_back(Q16 * 0.03720713072833462);
            this->weights.emplace_back(Q16 * 0.008110770829903342);
            this->weights.emplace_back(Q16 * 0.03720713072833462);
            this->weights.emplace_back(Q16 * 0.008110770829903342);
            this->weights.emplace_back(Q16 * 0.03720713072833462);
            this->weights.emplace_back(Q16 * 0.008110770829903342);
            this->weights.emplace_back(Q16 * 0.03720713072833462);
            this->weights.emplace_back(Q16 * 0.008110770829903342);
            this->weights.emplace_back(Q16 * 0.03720713072833462);
            this->weights.emplace_back(Q16 * 0.008110770829903342);
            this->weights.emplace_back(Q16 * 0.03720713072833462);
            this->weights.emplace_back(Q16 * 0.008110770829903342);
            this->weights.emplace_back(Q16 * 0.03720713072833462);
            this->weights.emplace_back(Q16 * 0.008110770829903342);
          }
      }

    AssertDimension(this->quadrature_points.size(), this->weights.size());
    Assert(this->quadrature_points.size() > 0,
           ExcNotImplemented(
             "Simplex::QGauss is currently only implemented for "
             "n_points_1D = 1, 2, 3, and 4 while you are asking for "
             "n_points_1D = " +
             Utilities::to_string(n_points_1D)));
  }



  template <int dim>
  QGaussWedge<dim>::QGaussWedge(const unsigned int n_points)
    : Quadrature<dim>()
  {
    AssertDimension(dim, 3);

    Simplex::QGauss<2> quad_tri(n_points);
    QGauss<1>          quad_line(n_points);

    for (unsigned int i = 0; i < quad_line.size(); ++i)
      for (unsigned int j = 0; j < quad_tri.size(); ++j)
        {
          this->quadrature_points.emplace_back(quad_tri.point(j)[0],
                                               quad_tri.point(j)[1],
                                               quad_line.point(i)[0]);
          this->weights.emplace_back(quad_tri.weight(j) * quad_line.weight(i));
        }

    AssertDimension(this->quadrature_points.size(), this->weights.size());
    Assert(this->quadrature_points.size() > 0,
           ExcMessage("No valid quadrature points!"));
  }



  template <int dim>
  QGaussPyramid<dim>::QGaussPyramid(const unsigned int n_points_1D)
    : Quadrature<dim>()
  {
    AssertDimension(dim, 3);

    if (n_points_1D == 1)
      {
        const double Q14 = 1.0 / 4.0;
        const double Q43 = 4.0 / 3.0;

        this->quadrature_points.emplace_back(0, 0, Q14);
        this->weights.emplace_back(Q43);
      }
    else if (n_points_1D == 2)
      {
        // clang-format off
        this->quadrature_points.emplace_back(-0.26318405556971, -0.26318405556971, 0.54415184401122);
        this->quadrature_points.emplace_back(-0.50661630334979, -0.50661630334979, 0.12251482265544);
        this->quadrature_points.emplace_back(-0.26318405556971, +0.26318405556971, 0.54415184401122);
        this->quadrature_points.emplace_back(-0.50661630334979, +0.50661630334979, 0.12251482265544);
        this->quadrature_points.emplace_back(+0.26318405556971, -0.26318405556971, 0.54415184401122);
        this->quadrature_points.emplace_back(+0.50661630334979, -0.50661630334979, 0.12251482265544);
        this->quadrature_points.emplace_back(+0.26318405556971, +0.26318405556971, 0.54415184401122);
        this->quadrature_points.emplace_back(+0.50661630334979, +0.50661630334979, 0.12251482265544);
        // clang-format on

        this->weights.emplace_back(0.10078588207983);
        this->weights.emplace_back(0.23254745125351);
        this->weights.emplace_back(0.10078588207983);
        this->weights.emplace_back(0.23254745125351);
        this->weights.emplace_back(0.10078588207983);
        this->weights.emplace_back(0.23254745125351);
        this->weights.emplace_back(0.10078588207983);
        this->weights.emplace_back(0.23254745125351);
      }

    AssertDimension(this->quadrature_points.size(), this->weights.size());
    Assert(this->quadrature_points.size() > 0,
           ExcMessage("No valid quadrature points!"));
  }

} // namespace Simplex


template class Simplex::QGauss<1>;
template class Simplex::QGauss<2>;
template class Simplex::QGauss<3>;
template class Simplex::QGaussWedge<1>;
template class Simplex::QGaussWedge<2>;
template class Simplex::QGaussWedge<3>;
template class Simplex::QGaussPyramid<1>;
template class Simplex::QGaussPyramid<2>;
template class Simplex::QGaussPyramid<3>;

DEAL_II_NAMESPACE_CLOSE
