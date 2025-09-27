// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2002 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Compute all refinement and restriction matrices for Tetrahedral elements.
// Each element has 8 children and there are a total of 3 different
// isotropic refinement choices for Tets. Also check for the DG case.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <string>
#include <vector>

#include "../tests.h"

#define PRECISION 6



template <int dim>
void
test_prolongation(int degree)
{
  FE_SimplexP<dim, dim> fe(degree);

  const unsigned int                           nc = 8;
  std::vector<std::vector<FullMatrix<double>>> matrices(
    static_cast<unsigned int>(IsotropicRefinementChoice::cut_tet_49),
    std::vector<FullMatrix<double>>(
      GeometryInfo<dim>::n_children(RefinementCase<dim>::isotropic_refinement),
      FullMatrix<double>(fe.n_dofs_per_cell(), fe.n_dofs_per_cell())));

  for (unsigned int refinement_direction =
         static_cast<unsigned int>(IsotropicRefinementChoice::cut_tet_68);
       refinement_direction <=
       static_cast<unsigned int>(IsotropicRefinementChoice::cut_tet_49);
       refinement_direction++)
    {
      deallog << "Refinement Direction: " << refinement_direction << std::endl;
      for (unsigned int i = 0; i < nc; ++i)
        {
          matrices[refinement_direction - 1][i] = fe.get_prolongation_matrix(
            i, RefinementCase<dim>(refinement_direction));

          deallog << "Child: " << i << std::endl;
          for (unsigned int j = 0; j < fe.n_dofs_per_cell(); ++j)
            {
              for (unsigned int k = 0; k < fe.n_dofs_per_cell(); ++k)
                deallog << matrices[refinement_direction - 1][i][j][k] << " ";
              deallog << std::endl;
            }
          deallog << std::endl;
        }
    }
}

template <int dim>
void
test_restriction(int degree)
{
  FE_SimplexP<dim, dim> fe(degree);

  const unsigned int                           nc = 8;
  std::vector<std::vector<FullMatrix<double>>> matrices(
    static_cast<unsigned int>(IsotropicRefinementChoice::cut_tet_49),
    std::vector<FullMatrix<double>>(
      GeometryInfo<dim>::n_children(RefinementCase<dim>::isotropic_refinement),
      FullMatrix<double>(fe.n_dofs_per_cell(), fe.n_dofs_per_cell())));

  for (unsigned int refinement_direction =
         static_cast<unsigned int>(IsotropicRefinementChoice::cut_tet_68);
       refinement_direction <=
       static_cast<unsigned int>(IsotropicRefinementChoice::cut_tet_49);
       refinement_direction++)
    {
      deallog << "Refinement Direction: " << refinement_direction << std::endl;
      for (unsigned int i = 0; i < nc; ++i)
        {
          matrices[refinement_direction - 1][i] = fe.get_restriction_matrix(
            i, RefinementCase<dim>(refinement_direction));

          deallog << "Child: " << i << std::endl;
          for (unsigned int j = 0; j < fe.n_dofs_per_cell(); ++j)
            {
              for (unsigned int k = 0; k < fe.n_dofs_per_cell(); ++k)
                deallog << matrices[refinement_direction - 1][i][j][k] << " ";
              deallog << std::endl;
            }
          deallog << std::endl;
        }
    }
}

template <int dim>
void
test_restriction_additive(int degree)
{
  FE_SimplexDGP<dim, dim> fe(degree);

  const unsigned int                           nc = 8;
  std::vector<std::vector<FullMatrix<double>>> matrices(
    static_cast<unsigned int>(IsotropicRefinementChoice::cut_tet_49),
    std::vector<FullMatrix<double>>(
      GeometryInfo<dim>::n_children(RefinementCase<dim>::isotropic_refinement),
      FullMatrix<double>(fe.n_dofs_per_cell(), fe.n_dofs_per_cell())));

  for (unsigned int refinement_direction =
         static_cast<unsigned int>(IsotropicRefinementChoice::cut_tet_68);
       refinement_direction <=
       static_cast<unsigned int>(IsotropicRefinementChoice::cut_tet_49);
       refinement_direction++)
    {
      deallog << "Refinement Direction: " << refinement_direction << std::endl;
      for (unsigned int i = 0; i < nc; ++i)
        {
          matrices[refinement_direction - 1][i] = fe.get_restriction_matrix(
            i, RefinementCase<dim>(refinement_direction));

          deallog << "Child: " << i << std::endl;
          for (unsigned int j = 0; j < fe.n_dofs_per_cell(); ++j)
            {
              for (unsigned int k = 0; k < fe.n_dofs_per_cell(); ++k)
                deallog << matrices[refinement_direction - 1][i][j][k] << " ";
              deallog << std::endl;
            }
          deallog << std::endl;
        }
    }
}



int
main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(PRECISION);
  deallog << std::fixed;
  deallog.attach(logfile);

  for (int degree = 1; degree < 4; degree++)
    {
      test_prolongation<3>(degree);
      test_restriction<3>(degree);
      test_restriction_additive<3>(degree);
    }

  return 0;
}
