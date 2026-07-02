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


// verify symmetric consistency of hp identity functions on mixed meshes
// fe_1.hp_vertex_dof_identities(fe_2) should give the same result as
// fe_2.hp_vertex_dof_identities(fe_1)

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_wedge_p.h>

#include <deal.II/hp/fe_collection.h>

#include "../tests.h"


bool
check_vectors_for_consistency(
  const std::vector<std::pair<unsigned int, unsigned int>> &result_fe,
  const std::vector<std::pair<unsigned int, unsigned int>> &result_fe_other)
{
  bool is_consistent = true;

  // check if they have the same length
  if (result_fe.size() != result_fe_other.size())
    is_consistent = false;

  // check if they have the same entries
  if (is_consistent)
    for (unsigned int i = 0; i < result_fe.size(); ++i)
      {
        if (result_fe[i].first == result_fe_other[i].second &&
            result_fe[i].second == result_fe_other[i].first)
          {
            // consistent result
          }
        else
          is_consistent = false;
      }
  return is_consistent;
}

template <int dim>
void
test_hp_vertex_dof_identities(const FiniteElement<dim> &fe,
                              const FiniteElement<dim> &fe_other)
{
  const auto result_fe       = fe.hp_vertex_dof_identities(fe_other);
  const auto result_fe_other = fe_other.hp_vertex_dof_identities(fe);

  const bool is_consistent =
    check_vectors_for_consistency(result_fe, result_fe_other);

  if (is_consistent)
    deallog << "hp vertex identities between " << fe.get_name() << " and "
            << fe_other.get_name() << " are consistent" << std::endl;
  else
    {
      deallog << "hp vertex identities between " << fe.get_name() << " and "
              << fe_other.get_name() << " are not consistent!" << std::endl;
      for (auto &r : result_fe)
        deallog << r.first << " " << r.second << "  ";
      deallog << std::endl;
      for (auto &r : result_fe_other)
        deallog << r.first << " " << r.second << "  ";
      deallog << std::endl;

      DEAL_II_ASSERT_UNREACHABLE();
    }
}

template <int dim>
void
test_hp_line_dof_identities(const FiniteElement<dim> &fe,
                            const FiniteElement<dim> &fe_other)
{
  const auto results_fe    = fe.hp_line_dof_identities(fe_other);
  const auto results_other = fe_other.hp_line_dof_identities(fe);

  const bool is_consistent =
    check_vectors_for_consistency(results_fe, results_other);

  if (is_consistent)
    deallog << "hp line dof identities between " << fe.get_name() << " and "
            << fe_other.get_name() << " are consistent" << std::endl;
  else
    {
      deallog << "hp line dof identities between " << fe.get_name() << " and "
              << fe_other.get_name() << " are not consistent!" << std::endl;
      for (auto &r : results_fe)
        deallog << r.first << " " << r.second << "  ";
      deallog << std::endl;
      for (auto &r : results_other)
        deallog << r.first << " " << r.second << "  ";
      deallog << std::endl;

      DEAL_II_ASSERT_UNREACHABLE();
    }
}

template <int dim>
void
test_hp_quad_dof_identities(const FiniteElement<dim> &fe,
                            const FiniteElement<dim> &fe_other)
{
  std::vector<std::pair<unsigned int, unsigned int>> results_fe;
  std::vector<std::pair<unsigned int, unsigned int>> results_fe_other;

  const auto reference_cell       = fe.reference_cell();
  const auto reference_cell_other = fe_other.reference_cell();
  // go over all faces and see if they are quads or tris
  // if they are the same then check if there are identities
  for (const auto f : reference_cell.face_indices())
    for (const auto f_other : reference_cell_other.face_indices())
      if (reference_cell.face_reference_cell(f) ==
          reference_cell_other.face_reference_cell(f_other))
        {
          const auto identities =
            fe.hp_quad_dof_identities(fe_other, f, f_other);
          const auto identities_other =
            fe_other.hp_quad_dof_identities(fe, f_other, f);

          for (const auto &r : identities)
            results_fe.emplace_back(r);

          for (const auto &r : identities_other)
            results_fe_other.emplace_back(r);
        }

  const bool is_consistent =
    check_vectors_for_consistency(results_fe, results_fe_other);

  if (is_consistent)
    deallog << "hp quad dof identities between " << fe.get_name() << " and "
            << fe_other.get_name() << " are consistent" << std::endl;
  else
    {
      deallog << "hp quad dof identities between " << fe.get_name() << " and "
              << fe_other.get_name() << " are not consistent!" << std::endl;
      for (auto &r : results_fe)
        deallog << r.first << " " << r.second << "  ";
      deallog << std::endl;
      for (auto &r : results_fe_other)
        deallog << r.first << " " << r.second << "  ";
      deallog << std::endl;

      DEAL_II_ASSERT_UNREACHABLE();
    }
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
    for (unsigned int j = i; j < fes.size(); ++j)
      test_hp_vertex_dof_identities(fes[i], fes[j]);
  deallog << std::endl;

  for (unsigned int i = 0; i < fes.size(); ++i)
    for (unsigned int j = i; j < fes.size(); ++j)
      test_hp_line_dof_identities(fes[i], fes[j]);
  deallog << std::endl;

  if constexpr (dim == 3)
    for (unsigned int i = 0; i < fes.size(); ++i)
      for (unsigned int j = i; j < fes.size(); ++j)
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
