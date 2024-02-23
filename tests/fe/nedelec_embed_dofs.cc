// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Checks if get_embedding_dofs() correctly maps the local DoF indices
// of an FE_Nedelec element of sub_degree to an element of sup_degree,
// where sub_degree <= sup_degree

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_nedelec.h>

#include "interpolate_common.h"

template <int dim>
void
check_nedelec_embed(const unsigned int sub_degree,
                    const unsigned int sup_degree)
{
  Assert((sub_degree > 0) && (sub_degree <= sup_degree),
         ExcIndexRange(sub_degree, 1, sup_degree));

  const FE_Nedelec<dim> fe_sub(sub_degree - 1);
  const FE_Nedelec<dim> fe_sup(sup_degree - 1);

  deallog << fe_sub.get_name() << ' ' << fe_sup.get_name() << ' '
          << fe_sub.dofs_per_cell << ' ' << fe_sup.dofs_per_cell;

  // Set the quadrature sufficiently high for the enriched finite element space
  const QGauss<dim> quadrature(fe_sup.degree + 1);

  std::vector<double> dofs_sub(fe_sub.dofs_per_cell, 0.);
  std::vector<double> dofs_sup(fe_sup.dofs_per_cell, 0.);

  // Assign arbitrary values to dofs_sub
  for (unsigned int i = 0; i < dofs_sub.size(); ++i)
    dofs_sub[i] = std::sin(i + 0.5);


  // Map the DoFs from fe_sub to fe_sup
  const std::vector<unsigned int> embedding_dofs =
    fe_sup.get_embedding_dofs(fe_sub.degree);
  for (unsigned int i = 0; i < dofs_sub.size(); ++i)
    dofs_sup[embedding_dofs[i]] = dofs_sub[i];

  // Compare the values at each quadrature point
  double result = 0.;
  for (unsigned int k = 0; k < quadrature.size(); ++k)
    for (unsigned int comp = 0; comp < fe_sup.n_components(); ++comp)
      {
        double eval_sub = 0., eval_sup = 0.;

        for (unsigned int i = 0; i < dofs_sub.size(); ++i)
          eval_sub +=
            dofs_sub[i] *
            fe_sub.shape_value_component(i, quadrature.point(k), comp);

        for (unsigned int i = 0; i < dofs_sup.size(); ++i)
          eval_sup +=
            dofs_sup[i] *
            fe_sup.shape_value_component(i, quadrature.point(k), comp);

        const double diff = std::abs(eval_sub - eval_sup);
        result            = std::max(result, diff);
      }
  deallog << " max diff " << result << std::endl;
}

int
main()
{
  initlog();

  // For 2-D
  {
    for (unsigned int low_deg = 1; low_deg <= 5; ++low_deg)
      for (unsigned int high_deg = low_deg; high_deg <= 10; ++high_deg)
        check_nedelec_embed<2>(low_deg, high_deg);
  }
}
