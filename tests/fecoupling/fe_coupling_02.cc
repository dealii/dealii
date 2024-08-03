// ---------------------------------------------------------------------
//
// Copyright (C) 2023 - 2024 by the deal.II authors
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


// Test QuadratureCouplingType::tensor_product integration in FECouplingValues
// with DoFCouplingType::independent. Typical for BEM coupling.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_coupling_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>

#include <fstream>
#include <iomanip>
#include <iostream>

#include "../tests.h"

#include "../test_grids.h"

int
main()
{
  initlog(0);
  constexpr unsigned int dim = 1;

  Triangulation<dim> tria;
  TestGrids::hyper_line(tria, 2);

  DoFHandler<dim> dofh(tria);
  FE_Q<dim>       fe(1);
  dofh.distribute_dofs(fe);

  UpdateFlags update_flags =
    update_quadrature_points | update_values | update_JxW_values;

  FEValues<dim> fv1(fe, QGauss<dim>(fe.degree + 1), update_flags);
  FEValues<dim> fv2(fe, QGauss<dim>(fe.degree + 1), update_flags);


  std::vector<types::global_dof_index> v1(fe.dofs_per_cell);
  std::vector<types::global_dof_index> v2(fe.dofs_per_cell);

  const auto &cell1 = dofh.begin();
  fv1.reinit(cell1);
  cell1->get_dof_indices(v1);

  unsigned int boundary_face = numbers::invalid_unsigned_int;
  unsigned int internal_face = numbers::invalid_unsigned_int;

  for (const unsigned int f : cell1->face_indices())
    {
      if (cell1->at_boundary(f))
        {
          boundary_face = f;
        }
      else if (!cell1->at_boundary(f))
        {
          internal_face = f;
        }
      if (boundary_face != numbers::invalid_unsigned_int &&
          internal_face != numbers::invalid_unsigned_int)
        break;
    }

  const auto &cell2 = cell1->neighbor(internal_face);
  fv2.reinit(cell2);
  cell2->get_dof_indices(v2);

  // Now we have two cells, cell1 and cell2, and we want to couple them. We
  // integrate over the tensor product of the two cells:
  //
  // \int_T1 \int_T2 |x1-x2| v_i(x1) v_j(x2) dx1 dx2
  //

  FECouplingValues<dim> fcv(fv1,
                            fv2,
                            DoFCouplingType::independent,
                            QuadratureCouplingType::tensor_product);

  deallog << "dofs_per_cell = " << fe.dofs_per_cell << std::endl
          << "n_q_points_per_cell = " << fv2.n_quadrature_points << std::endl
          << "n_q_points = " << fcv.n_quadrature_points() << std::endl;

  FullMatrix<double> matrix(fcv.n_first_dofs(), fcv.n_second_dofs());

  const auto first  = fcv.get_first_extractor(FEValuesExtractors::Scalar(0));
  const auto second = fcv.get_second_extractor(FEValuesExtractors::Scalar(0));

  // We need to loop over all the coupling quadrature points
  for (const auto q : fcv.quadrature_point_indices())
    {
      const auto &[x, y] = fcv.quadrature_point(q);
      for (const auto i : fcv.first_dof_indices())
        {
          const auto &vi = fcv[first].value(i, q);
          for (const auto j : fcv.second_dof_indices())
            {
              const auto &vj = fcv[second].value(j, q);

              deallog << std::left                                           //
                      << "vi[" << i << "]: " << std::setw(10) << vi          //
                      << "vj[" << j << "]: " << std::setw(10) << vj          //
                      << "x[" << q << "]: " << std::setw(10) << x            //
                      << "y[" << q << "]: " << std::setw(10) << y            //
                      << "JxW[" << q << "]: " << std::setw(10) << fcv.JxW(q) //
                      << "distance = " << x.distance(y) << std::endl;
              matrix(i, j) += x.distance(y) * vi * vj * fcv.JxW(q);
            }
        }
    }
  // Now print the matrix:
  deallog << "matrix = " << std::endl;
  for (unsigned int i = 0; i < matrix.m(); ++i)
    {
      for (unsigned int j = 0; j < matrix.n(); ++j)
        {
          deallog << matrix(i, j) << " ";
        }
      deallog << std::endl;
    }
}
