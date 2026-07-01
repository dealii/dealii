// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check that the hp identity functions for mixed meshes by checking if the
// expected number of entries are returned


#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_wedge_p.h>

#include <deal.II/hp/fe_collection.h>

#include "../tests.h"


template <int dim>
void
test_hp_vertex_dof_identities(const FiniteElement<dim> &fe,
                              const FiniteElement<dim> &fe_other)
{
  const auto result = fe.hp_vertex_dof_identities(fe_other);
  // expect one entry for all cases
  if (result.size() == 1)
    deallog << "hp vertex dof identities between " << fe.get_name() << " and "
            << fe_other.get_name() << " ok " << result.size() << std::endl;
  else
    {
      deallog << "hp vertex dof identities between " << fe.get_name() << " and "
              << fe_other.get_name() << " expected size 1 but got "
              << result.size() << std::endl;

      DEAL_II_ASSERT_UNREACHABLE();
    }
}

template <int dim>
void
test_hp_line_dof_identities(const FiniteElement<dim> &fe,
                            const FiniteElement<dim> &fe_other)
{
  const auto results = fe.hp_line_dof_identities(fe_other);

  unsigned int n_expected = 0;
  // with the same degree the DoFs could agree but do not have to
  if (fe.degree == fe_other.degree)
    {
      // for degree 1 nothing to do
      if (fe.degree > 1)
        {
          // check if points are equidistant
          bool fe_is_equidistant, fe_other_is_equidistant;

          const double scale =
            (fe.unit_support_point(1) - fe.unit_support_point(0))[0];
          const Point<dim> first_line_point =
            fe.unit_support_point(fe.get_first_line_index());
          const double distance =
            fe.unit_support_point(0).distance(first_line_point);
          if (std::fabs(distance - (scale / fe.degree)) < 1e-12)
            fe_is_equidistant = true;
          else
            fe_is_equidistant = false;

          const double     scale_other = (fe_other.unit_support_point(1) -
                                      fe_other.unit_support_point(0))[0];
          const Point<dim> first_line_point_other =
            fe_other.unit_support_point(fe_other.get_first_line_index());
          const double distance_other =
            fe_other.unit_support_point(0).distance(first_line_point_other);
          if (std::fabs(distance_other - (scale_other / fe_other.degree)) <
              1e-12)
            fe_other_is_equidistant = true;
          else
            fe_other_is_equidistant = false;

          // if both are equidistant or both are not equidistant then we have
          // degree - 1 DoFs on each line
          if (fe_is_equidistant == fe_other_is_equidistant)
            n_expected = fe.degree - 1;
          // if they do not agree the midpoints could coincide if they have one
          else if (fe.degree % 2 == 0 && fe_other.degree % 2 == 0)
            {
              n_expected = 1;
            }
        }
    }
  // if they have different degrees the midpoints could coincide if they have
  // one
  else if (fe.degree % 2 == 0 && fe_other.degree % 2 == 0)
    {
      // in this case they have a shared midpoint between them
      n_expected = 1;
    }

  if (results.size() == n_expected)
    deallog << "hp line dof identities between " << fe.get_name() << " and "
            << fe_other.get_name() << " ok " << results.size() << std::endl;
  else
    {
      deallog << "hp line dof identities between " << fe.get_name() << " and "
              << fe_other.get_name() << " expected size " << n_expected
              << " but got " << results.size() << std::endl;

      DEAL_II_ASSERT_UNREACHABLE();
    }
}



template <int dim>
void
test_hp_quad_dof_identities(const FiniteElement<dim> &fe,
                            const FiniteElement<dim> &fe_other)
{
  deallog << "hp quad dof identities between " << fe.get_name() << " and "
          << fe_other.get_name() << " ";

  const auto reference_cell       = fe.reference_cell();
  const auto reference_cell_other = fe_other.reference_cell();

  // go over all faces and see if they are quads or tris
  // if they are the same then check if there are identities
  for (const auto f : reference_cell.face_indices())
    for (const auto f_other : reference_cell_other.face_indices())
      if (reference_cell.face_reference_cell(f) ==
          reference_cell_other.face_reference_cell(f_other))
        {
          const auto results = fe.hp_quad_dof_identities(fe_other, f, f_other);

          unsigned int n_expected = 0;
          // they could have the same support points if they have the same
          // degree
          if (fe.degree == fe_other.degree)
            {
              // nothing to do for degree 1
              if (fe.degree > 1)
                {
                  // check if points are equidistant
                  bool fe_is_equidistant, fe_other_is_equidistant;

                  const double scale =
                    (fe.unit_support_point(1) - fe.unit_support_point(0))[0];
                  const Point<dim> first_line_point =
                    fe.unit_support_point(fe.get_first_line_index());
                  const double distance =
                    fe.unit_support_point(0).distance(first_line_point);
                  if (std::fabs(distance - (scale / fe.degree)) < 1e-12)
                    fe_is_equidistant = true;
                  else
                    fe_is_equidistant = false;

                  const double scale_other =
                    (fe_other.unit_support_point(1) -
                     fe_other.unit_support_point(0))[0];
                  const Point<dim> first_line_point_other =
                    fe_other.unit_support_point(
                      fe_other.get_first_line_index());
                  const double distance_other =
                    fe_other.unit_support_point(0).distance(
                      first_line_point_other);
                  if (std::fabs(distance_other -
                                (scale_other / fe_other.degree)) < 1e-12)
                    fe_other_is_equidistant = true;
                  else
                    fe_other_is_equidistant = false;

                  // if both are equidistant or both are not equidistant then
                  // they should be the same
                  if (fe_is_equidistant == fe_other_is_equidistant)
                    n_expected = fe.n_dofs_per_quad(f);
                  // if they do not agree the midpoints could coincide if they
                  // have one
                  else if (fe.degree % 2 == 0 && fe_other.degree % 2 == 0)
                    {
                      n_expected = 1;
                    }
                }
            }
          // both can have the same midpoint
          else if (reference_cell.face_reference_cell(f).is_hyper_cube())
            if (fe.degree % 2 == 0 && fe_other.degree % 2 == 0)
              {
                // only midpoint is identical
                n_expected = 1;
              }

          if (results.size() == n_expected)
            deallog << "ok " << results.size() << ", ";
          else
            {
              deallog << " on face " << f << " and " << f_other
                      << " expected size is " << n_expected << " but got "
                      << results.size() << std::endl;

              DEAL_II_ASSERT_UNREACHABLE();
            }
        }
  deallog << std::endl;
}

template <int dim>
void
test()
{
  hp::FECollection<dim> fes;
  for (unsigned int i = 1; i <= 3; ++i)
    fes.push_back(FE_SimplexP<dim>(i));
  for (unsigned int i = 1; i <= 4; ++i)
    fes.push_back(FE_Q<dim>(i));
  for (unsigned int i = 1; i <= 4; ++i)
    fes.push_back(FE_Q<dim>(QIterated<1>(QTrapezoid<1>(), i)));
  if constexpr (dim == 3)
    {
      for (unsigned int i = 1; i <= 2; ++i)
        fes.push_back(FE_WedgeP<dim>(i));
      for (unsigned int i = 1; i <= 3; ++i)
        fes.push_back(FE_PyramidP<dim>(i));
    }

  for (unsigned int i = 0; i < fes.size(); ++i)
    for (unsigned int j = 0; j < fes.size(); ++j)
      test_hp_vertex_dof_identities(fes[i], fes[j]);
  deallog << std::endl;

  for (unsigned int i = 0; i < fes.size(); ++i)
    for (unsigned int j = 0; j < fes.size(); ++j)
      test_hp_line_dof_identities(fes[i], fes[j]);
  deallog << std::endl;

  if constexpr (dim == 3)
    for (unsigned int i = 0; i < fes.size(); ++i)
      for (unsigned int j = 0; j < fes.size(); ++j)
        test_hp_quad_dof_identities(fes[i], fes[j]);
  deallog << std::endl;
}



int
main()
{
  initlog();

  deallog.push("2d");
  test<2>();
  deallog.pop();

  deallog.push("3d");
  test<3>();
  deallog.pop();
}
