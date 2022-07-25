// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2021 by the deal.II authors
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

#include <deal.II/base/config.h>

#include <deal.II/base/polynomials_barycentric.h>
#include <deal.II/base/qprojector.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p_bubbles.h>
#include <deal.II/fe/fe_tools.h>

DEAL_II_NAMESPACE_OPEN

namespace FE_P_BubblesImplementation
{
  template <int dim>
  std::vector<unsigned int>
  get_dpo_vector(const unsigned int degree)
  {
    std::vector<unsigned int> dpo(dim + 1);
    if (degree == 0)
      {
        dpo[dim] = 1; // single interior dof
      }
    else
      {
        Assert(degree == 1 || degree == 2, ExcNotImplemented());
        dpo[0] = 1; // vertex dofs

        if (degree == 2)
          {
            dpo[1] = 1; // line dofs

            if (dim > 1)
              dpo[dim] = 1; // the internal bubble function
            if (dim == 3)
              dpo[dim - 1] = 1; // face bubble functions
          }
      }

    return dpo;
  }



  template <int dim>
  std::vector<Point<dim>>
  unit_support_points(const unsigned int degree)
  {
    Assert(degree < 3, ExcNotImplemented());
    // Start with the points used by FE_SimplexP, and then add bubbles.
    FE_SimplexP<dim>        fe_p(degree);
    std::vector<Point<dim>> points = fe_p.get_unit_support_points();

    Point<dim> centroid;
    std::fill(centroid.begin_raw(), centroid.end_raw(), 1.0 / double(dim + 1));

    switch (dim)
      {
        case 1:
          // nothing more to do
          return points;
        case 2:
          {
            if (degree == 2)
              points.push_back(centroid);
            return points;
          }
        case 3:
          {
            if (degree == 2)
              {
                const double q13 = 1.0 / 3.0;
                points.emplace_back(q13, q13, 0.0);
                points.emplace_back(q13, 0.0, q13);
                points.emplace_back(0.0, q13, q13);
                points.emplace_back(q13, q13, q13);
                points.push_back(centroid);
              }
            return points;
          }
        default:
          Assert(false, ExcNotImplemented());
      }
    return points;
  }



  template <>
  std::vector<Point<0>>
  unit_support_points<0>(const unsigned int /*degree*/)
  {
    std::vector<Point<0>> points;
    points.emplace_back();
    return points;
  }



  template <int dim>
  BarycentricPolynomials<dim>
  get_basis(const unsigned int degree)
  {
    Point<dim> centroid;
    std::fill(centroid.begin_raw(), centroid.end_raw(), 1.0 / double(dim + 1));

    auto M = [](const unsigned int d) {
      return BarycentricPolynomial<dim, double>::monomial(d);
    };

    switch (degree)
      {
        // we don't need to add bubbles to P0 or P1
        case 0:
        case 1:
          return BarycentricPolynomials<dim>::get_fe_p_basis(degree);
        case 2:
          {
            const auto fe_p =
              BarycentricPolynomials<dim>::get_fe_p_basis(degree);
            // no further work is needed in 1D
            if (dim == 1)
              return fe_p;

            // in 2D and 3D we add a centroid bubble function
            auto c_bubble = BarycentricPolynomial<dim>() + 1;
            for (unsigned int d = 0; d < dim + 1; ++d)
              c_bubble = c_bubble * M(d);
            c_bubble = c_bubble / c_bubble.value(centroid);

            std::vector<BarycentricPolynomial<dim>> bubble_functions;
            if (dim == 2)
              {
                bubble_functions.push_back(c_bubble);
              }
            else if (dim == 3)
              {
                // need 'face bubble' functions in addition to the centroid.
                // Furthermore we need to subtract them off from the other
                // functions so that we end up with an interpolatory basis
                auto b0 = 27 * M(0) * M(1) * M(2);
                bubble_functions.push_back(b0 - b0.value(centroid) * c_bubble);
                auto b1 = 27 * M(0) * M(1) * M(3);
                bubble_functions.push_back(b1 - b1.value(centroid) * c_bubble);
                auto b2 = 27 * M(0) * M(2) * M(3);
                bubble_functions.push_back(b2 - b2.value(centroid) * c_bubble);
                auto b3 = 27 * M(1) * M(2) * M(3);
                bubble_functions.push_back(b3 - b3.value(centroid) * c_bubble);

                bubble_functions.push_back(c_bubble);
              }

            // Extract out the support points for the extra bubble (both
            // volume and face) functions:
            const std::vector<Point<dim>> support_points =
              unit_support_points<dim>(degree);
            const std::vector<Point<dim>> bubble_support_points(
              support_points.begin() + fe_p.n(), support_points.end());
            Assert(bubble_support_points.size() == bubble_functions.size(),
                   ExcInternalError());
            const unsigned int n_bubbles = bubble_support_points.size();

            // Assemble the final basis:
            std::vector<BarycentricPolynomial<dim>> lump_polys;
            for (unsigned int i = 0; i < fe_p.n(); ++i)
              {
                BarycentricPolynomial<dim> p = fe_p[i];

                for (unsigned int j = 0; j < n_bubbles; ++j)
                  {
                    p = p -
                        p.value(bubble_support_points[j]) * bubble_functions[j];
                  }

                lump_polys.push_back(p);
              }

            for (auto &p : bubble_functions)
              lump_polys.push_back(std::move(p));

              // Sanity check:
#ifdef DEBUG
            BarycentricPolynomial<dim> unity;
            for (const auto &p : lump_polys)
              unity = unity + p;

            Point<dim> test;
            for (unsigned int d = 0; d < dim; ++d)
              test[d] = 2.0;
            Assert(std::abs(unity.value(test) - 1.0) < 1e-10,
                   ExcInternalError());
#endif

            return BarycentricPolynomials<dim>(lump_polys);
          }
        default:
          Assert(degree < 3, ExcNotImplemented());
      }

    Assert(degree < 3, ExcNotImplemented());
    // bogus return to placate compilers
    return BarycentricPolynomials<dim>::get_fe_p_basis(degree);
  }



  template <int dim>
  FiniteElementData<dim>
  get_fe_data(const unsigned int degree)
  {
    // It's not efficient, but delegate computation of the degree of the
    // finite element (which is different from the input argument) to the
    // basis.
    const auto polys = get_basis<dim>(degree);
    return FiniteElementData<dim>(get_dpo_vector<dim>(degree),
                                  ReferenceCells::get_simplex<dim>(),
                                  1, // n_components
                                  polys.degree(),
                                  FiniteElementData<dim>::H1);
  }
} // namespace FE_P_BubblesImplementation



template <int dim, int spacedim>
FE_SimplexP_Bubbles<dim, spacedim>::FE_SimplexP_Bubbles(
  const unsigned int degree)
  : FE_SimplexPoly<dim, spacedim>(
      FE_P_BubblesImplementation::get_basis<dim>(degree),
      FE_P_BubblesImplementation::get_fe_data<dim>(degree),
      FE_P_BubblesImplementation::unit_support_points<dim>(degree),
      {FE_P_BubblesImplementation::unit_support_points<dim - 1>(degree)},
      // Interface constraints are not yet implemented
      FullMatrix<double>())
  , approximation_degree(degree)
{}



template <int dim, int spacedim>
std::string
FE_SimplexP_Bubbles<dim, spacedim>::get_name() const
{
  return "FE_SimplexP_Bubbles<" + Utilities::dim_string(dim, spacedim) + ">" +
         "(" + std::to_string(approximation_degree) + ")";
}



template <int dim, int spacedim>
std::unique_ptr<FiniteElement<dim, spacedim>>
FE_SimplexP_Bubbles<dim, spacedim>::clone() const
{
  return std::make_unique<FE_SimplexP_Bubbles<dim, spacedim>>(*this);
}

// explicit instantiations
#include "fe_simplex_p_bubbles.inst"

DEAL_II_NAMESPACE_CLOSE
