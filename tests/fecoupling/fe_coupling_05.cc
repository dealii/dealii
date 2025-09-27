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


// Test QuadratureCouplingType::matching integration in FECouplingValues
// with DoFCouplingType::contiguous. Typical for bulk-surface coupling.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_coupling_values.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_interface_values.h>
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
  initlog();
  constexpr unsigned int dim   = 2;
  constexpr unsigned int codim = dim - 1;

  Triangulation<dim> tria;
  TestGrids::hyper_line(tria, 2);

  Triangulation<codim, dim> codim_tria;
  TestGrids::hyper_line(codim_tria, 2);

  DoFHandler<dim> dofh(tria);
  FE_Q<dim>       fe(1);
  dofh.distribute_dofs(fe);

  DoFHandler<codim, dim> codim_dofh(codim_tria);
  FE_Q<codim, dim>       codim_fe(1);
  codim_dofh.distribute_dofs(codim_fe);

  UpdateFlags update_flags =
    update_quadrature_points | update_values | update_JxW_values;

  QGauss<codim> quad(fe.degree + 1);

  FEFaceValues<dim>    fv1(fe, quad, update_flags);
  FEValues<codim, dim> fv2(codim_fe, quad, update_flags);

  std::vector<types::global_dof_index> v1(fe.dofs_per_cell);
  std::vector<types::global_dof_index> v2(codim_fe.dofs_per_cell);

  deallog << "Bulk dofs: " << fe.dofs_per_cell
          << ", Surface dofs: " << codim_fe.dofs_per_cell << std::endl;

  const auto &cell1 = dofh.begin();
  cell1->get_dof_indices(v1);
  // Coincides with the first cell of codim tria.
  const unsigned int face = 2;

  const auto &cell2 = codim_dofh.begin();
  cell2->get_dof_indices(v2);

  fv1.reinit(cell1, face);
  fv2.reinit(cell2);

  // Now we have two cells, cell1 and cell2, and we want to integrate on the
  // face.
  //
  // \int_f v_i(x) * v_j(x) dx
  //

  FECouplingValues<dim, codim, dim> fcv(fv1,
                                        fv2,
                                        DoFCouplingType::contiguous,
                                        QuadratureCouplingType::matching);

  FullMatrix<double> matrix(fcv.n_coupling_dofs(), fcv.n_coupling_dofs());
  const FEValuesExtractors::Scalar scalar(0);

  const auto bulk    = fcv.get_first_extractor(scalar);
  const auto surface = fcv.get_second_extractor(scalar);

  deallog << "Renumbering vectors: " << std::endl;

  for (const auto i : fcv.coupling_dof_indices())
    {
      auto id_first  = fcv.coupling_dof_to_dof_indices(i).first;
      auto id_second = fcv.coupling_dof_to_dof_indices(i).second;
      deallog << "i: " << i << " bulk: " << (int)id_first
              << " surface: " << (int)id_second << std::endl;
    }

  // We need to loop over all the coupling quadrature points
  for (const auto q : fcv.quadrature_point_indices())
    {
      const auto &[x, y] = fcv.quadrature_point(q);
      for (const auto i : fcv.coupling_dof_indices())
        {
          const auto &bulk_vi    = fcv[bulk].value(i, q);
          const auto &surface_vi = fcv[surface].value(i, q);

          for (const auto j : fcv.coupling_dof_indices())
            {
              const auto &bulk_vj    = fcv[bulk].value(j, q);
              const auto &surface_vj = fcv[surface].value(j, q);

              deallog << std::left                                           //
                      << "bvi[" << i << "]: " << std::setw(10) << bulk_vi    //
                      << "bvj[" << j << "]: " << std::setw(10) << bulk_vj    //
                      << "svi[" << i << "]: " << std::setw(10) << surface_vi //
                      << "svj[" << j << "]: " << std::setw(10) << surface_vj //
                      << "x[" << q << "]: " << std::setw(10) << x            //
                      << "  y[" << q << "]: " << std::setw(10) << y          //
                      << std::endl;
              matrix(i, j) += (bulk_vi * bulk_vj +             //
                               10 * bulk_vi * surface_vj +     //
                               100 * surface_vi * bulk_vj +    //
                               1000 * surface_vi * surface_vj) //
                              * fcv.JxW(q);
            }
        }
    }
  // Now print the matrix:
  deallog << "matrix = " << std::endl << std::left;
  for (unsigned int i = 0; i < matrix.m(); ++i)
    {
      for (unsigned int j = 0; j < matrix.n(); ++j)
        {
          deallog << std::setw(8) << matrix(i, j) << " ";
        }
      deallog << std::endl;
    }
}
