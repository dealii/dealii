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



// verify hp identity functions for mixed meshes
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
  deallog << "hp vertex dof identities between " << fe.get_name() << " and "
          << fe_other.get_name() << ": " << result[0].first << " "
          << result[0].second << std::endl;
}

template <int dim>
void
test_hp_line_dof_identities(const FiniteElement<dim> &fe,
                            const FiniteElement<dim> &fe_other)
{
  const auto results = fe.hp_line_dof_identities(fe_other);

  deallog << "hp line dof identities between " << fe.get_name() << " and "
          << fe_other.get_name() << ": ";
  for (const auto &result : results)
    deallog << result.first << " " << result.second << ", ";
  deallog << std::endl;
}

template <int dim>
void
test_hp_quad_dof_identities(const FiniteElement<dim> &fe,
                            const FiniteElement<dim> &fe_other)
{
  std::vector<std::pair<unsigned int, unsigned int>> results;

  if (((dynamic_cast<const FE_SimplexP<dim, dim> *>(&fe)) &&
       (dynamic_cast<const FE_Q<dim, dim> *>(&fe_other))) ||
      ((dynamic_cast<const FE_SimplexP<dim, dim> *>(&fe_other)) &&
       (dynamic_cast<const FE_Q<dim, dim> *>(&fe))))
    return;

  if ((dynamic_cast<const FE_WedgeP<dim, dim> *>(&fe)))
    {
      unsigned int face_no = 0;
      if ((dynamic_cast<const FE_Q<dim, dim> *>(&fe_other)))
        face_no = 2;
      results = fe.hp_quad_dof_identities(fe_other, face_no);

      if ((dynamic_cast<const FE_WedgeP<dim, dim> *>(&fe_other)) ||
          (dynamic_cast<const FE_PyramidP<dim, dim> *>(&fe_other)))
        {
          const auto additional_results =
            fe.hp_quad_dof_identities(fe_other, 2);
          for (const auto &r : additional_results)
            results.emplace_back(r);
        }
    }
  else if ((dynamic_cast<const FE_PyramidP<dim, dim> *>(&fe)))
    {
      unsigned int face_no = 0;
      if ((dynamic_cast<const FE_SimplexP<dim, dim> *>(&fe_other)))
        face_no = 1;
      results = fe.hp_quad_dof_identities(fe_other, face_no);

      if ((dynamic_cast<const FE_WedgeP<dim, dim> *>(&fe_other)) ||
          (dynamic_cast<const FE_PyramidP<dim, dim> *>(&fe_other)))
        {
          const auto additional_results =
            fe.hp_quad_dof_identities(fe_other, 1);
          for (const auto &r : additional_results)
            results.emplace_back(r);
        }
    }
  else
    {
      results = fe.hp_quad_dof_identities(fe_other);
    }

  deallog << "hp quad dof identities between " << fe.get_name() << " and "
          << fe_other.get_name() << ": ";
  for (const auto &result : results)
    deallog << result.first << " " << result.second << ", ";
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
      for (unsigned int i = 1; i <= 1; ++i)
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
