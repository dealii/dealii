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
            Quadrature<dim>::operator=(QWitherdenVincent<dim>(n_points_1D));
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
            Quadrature<dim>::operator=(QWitherdenVincent<dim>(n_points_1D));
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

  namespace
  {
    template <std::size_t b_dim>
    std::vector<std::array<double, b_dim>>
    all_permutations(const std::array<double, b_dim> &b_point)
    {
      std::vector<std::array<double, b_dim>> output;

      // We want all possible permutations of the barycentric coordinates.
      // The easiest way to get all of them is to sort the input first and
      // then use next_permutation to cycle through them all.
      std::array<double, b_dim> temp = b_point;
      std::sort(temp.begin(), temp.end());
      do
        {
          output.push_back(temp);
        }
      while (std::next_permutation(temp.begin(), temp.end()));

      return output;
    }
  } // namespace



  template <int dim>
  QWitherdenVincent<dim>::QWitherdenVincent(const unsigned int n_points_1D)
    : QSimplex<dim>(Quadrature<dim>())
  {
    Assert(1 <= dim && dim <= 3, ExcNotImplemented());
    // Just use Gauss in 1D: this is a high-order open rule so this is a
    // reasonable equivalent for generic programming.
    if (dim == 1)
      {
        Quadrature<dim>::operator=(dealii::QGauss<dim>(n_points_1D));
        return;
      }

    std::array<double, dim + 1> centroid;
    std::fill(centroid.begin(), centroid.end(), 1.0 / (dim + 1.0));
    std::vector<std::vector<std::array<double, dim + 1>>> b_point_permutations;
    std::vector<double>                                   b_weights;

    // We can simplify the implementation of these quadrature rules
    // by quite a bit by exploiting symmetry - we do essentially the
    // same thing for each barycentric coordinate, so we can express
    // our quadrature rule as permutations of barycentric points
    // instead of writing things out explicitly.

    // Apply a Barycentric permutation where one point is different.
    auto process_point_1 = [&](const double a, const double w) {
      const double                b = 1.0 - dim * a;
      std::array<double, dim + 1> b_point;
      std::fill(b_point.begin(), b_point.begin() + dim, a);
      b_point[dim] = b;

      b_weights.push_back(w);
      b_point_permutations.push_back(all_permutations(b_point));
    };

    // Apply a Barycentric permutation where two points (in 3D) are different.
    auto process_point_2 = [&](const double a, const double w) {
      Assert(dim == 3, ExcInternalError());
      const double                b = (1.0 - 2.0 * a) / 2.0;
      std::array<double, dim + 1> b_point;
      std::fill(b_point.begin(), b_point.begin() + dim - 1, a);
      b_point[dim - 1] = b;
      b_point[dim]     = b;

      b_weights.push_back(w);
      b_point_permutations.push_back(all_permutations(b_point));
    };

    // Apply a Barycentric permutation where three (or four) points
    // are different (since there are two inputs).
    auto process_point_3 = [&](const double a, const double b, const double w) {
      const double                c = 1.0 - (dim - 1.0) * a - b;
      std::array<double, dim + 1> b_point;
      std::fill(b_point.begin(), b_point.begin() + dim - 1, a);
      b_point[dim - 1] = b;
      b_point[dim]     = c;

      b_weights.push_back(w);
      b_point_permutations.push_back(all_permutations(b_point));
    };

    if (n_points_1D == 1)
      {
        b_point_permutations.push_back({centroid});
        b_weights.push_back(1.0);
      }
    else if (n_points_1D == 2)
      {
        // This is WV-4 in 2D and WV-3 in 3D
        if (dim == 2)
          {
            process_point_1(9.1576213509770743e-02, 1.0995174365532187e-01);
            process_point_1(4.4594849091596489e-01, 2.2338158967801147e-01);
          }
        else if (dim == 3)
          {
            process_point_1(3.281633025163817e-01, 1.362178425370874e-01);
            process_point_1(1.080472498984286e-01, 1.137821574629126e-01);
          }
      }
    else if (n_points_1D == 3)
      {
        // This is the WV-5 rule in both 2D and 3D
        if (dim == 2)
          {
            b_weights.push_back(0.225);
            b_point_permutations.push_back({centroid});

            process_point_1(1.0128650732345634e-01, 1.2593918054482714e-01);
            process_point_1(4.7014206410511511e-01, 1.3239415278850619e-01);
          }
        else if (dim == 3)
          {
            process_point_1(3.108859192633006e-01, 1.126879257180159e-01);
            process_point_1(9.273525031089125e-02, 7.349304311636196e-02);

            process_point_2(4.550370412564964e-02, 4.254602077708147e-02);
          }
      }
    else if (n_points_1D == 4)
      {
        // This is the WV-7 rule in both 2D and 3D
        if (dim == 2)
          {
            process_point_1(3.3730648554587850e-02, 1.6545050110792131e-02);
            process_point_1(4.7430969250471822e-01, 7.7086646185986069e-02);
            process_point_1(2.4157738259540357e-01, 1.2794417123015558e-01);
            process_point_3(4.7036644652595216e-02,
                            1.9868331479735168e-01,
                            5.5878732903199779e-02);
          }
        else if (dim == 3)
          {
            b_point_permutations.push_back({centroid});
            b_weights.push_back(9.548528946413085e-02);

            process_point_1(3.157011497782028e-01, 4.232958120996703e-02);
            process_point_2(5.048982259839635e-02, 3.189692783285758e-02);

            process_point_3(1.888338310260010e-01,
                            5.751716375870000e-01,
                            3.720713072833462e-02);
            process_point_3(2.126547254148314e-02,
                            8.108302410985486e-01,
                            8.110770829903342e-03);
          }
      }
    else if (n_points_1D == 5)
      {
        // This is the WV-9 rule in both 2D and 3D
        if (dim == 2)
          {
            b_point_permutations.push_back({centroid});
            b_weights.push_back(9.7135796282798836e-02);

            process_point_1(4.4729513394452691e-02, 2.5577675658698031e-02);
            process_point_1(4.8968251919873762e-01, 3.1334700227139071e-02);
            process_point_1(4.3708959149293664e-01, 7.7827541004774278e-02);
            process_point_1(1.8820353561903275e-01, 7.9647738927210249e-02);

            process_point_3(3.6838412054736258e-02,
                            2.2196298916076568e-01,
                            4.3283539377289376e-02);
          }
        else if (dim == 3)
          {
            b_point_permutations.push_back({centroid});
            b_weights.push_back(5.801054891248025e-02);

            process_point_1(6.198169755222693e-10, 6.431928175925639e-05);
            process_point_1(1.607745353952616e-01, 2.317333846242546e-02);
            process_point_1(3.222765218214210e-01, 2.956291233542929e-02);
            process_point_1(4.510891834541358e-02, 8.063979979616182e-03);

            process_point_2(1.122965460043761e-01, 3.813408010370246e-02);

            process_point_3(4.588714487524592e-01,
                            2.554579233041310e-03,
                            8.384422198298552e-03);
            process_point_3(3.377587068533860e-02,
                            7.183503264420745e-01,
                            1.023455935274533e-02);
            process_point_3(1.836413698099279e-01,
                            3.441591057817528e-02,
                            2.052491596798814e-02);
          }
      }
    else if (n_points_1D == 6)
      {
        // There is no WV-11 rule in 3D yet
        if (dim == 2)
          {
            b_point_permutations.push_back({centroid});
            b_weights.push_back(8.5761179732224219e-02);

            process_point_1(2.8485417614371900e-02, 1.0431870512894697e-02);
            process_point_1(4.9589190096589092e-01, 1.6606273054585369e-02);
            process_point_1(1.0263548271224643e-01, 3.8630759237019321e-02);
            process_point_1(4.3846592676435220e-01, 6.7316154079468296e-02);
            process_point_1(2.1021995670317828e-01, 7.0515684111716576e-02);

            process_point_3(7.3254276860644785e-03,
                            1.4932478865208237e-01,
                            1.0290289572953278e-02);
            process_point_3(4.6010500165429957e-02,
                            2.8958112563770588e-01,
                            4.0332476640500554e-02);
          }
        else if (dim == 3)
          {
            Assert(false, ExcNotImplemented());
          }
      }
    else
      {
        Assert(false, ExcNotImplemented());
      }

    Assert(b_point_permutations.size() == b_weights.size(), ExcInternalError());
    for (unsigned int permutation_n = 0; permutation_n < b_weights.size();
         ++permutation_n)
      {
        for (const std::array<double, dim + 1> &b_point :
             b_point_permutations[permutation_n])
          {
            const double volume = (dim == 2 ? 1.0 / 2.0 : 1.0 / 6.0);
            this->weights.emplace_back(volume * b_weights[permutation_n]);
            Point<dim> c_point;
            std::copy(b_point.begin(),
                      b_point.begin() + dim,
                      c_point.begin_raw());
            this->quadrature_points.emplace_back(c_point);
          }
      }
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

template class Simplex::QWitherdenVincent<1>;
template class Simplex::QWitherdenVincent<2>;
template class Simplex::QWitherdenVincent<3>;

DEAL_II_NAMESPACE_CLOSE
