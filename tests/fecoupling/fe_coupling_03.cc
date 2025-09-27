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
// with DoFCouplingType::contiguous. Typical for DG interfaces.

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
  constexpr unsigned int dim = 2;

  Triangulation<dim> tria;
  TestGrids::hyper_line(tria, 2);

  DoFHandler<dim> dofh(tria);
  FE_DGQ<dim>     fe(1);
  dofh.distribute_dofs(fe);

  UpdateFlags update_flags =
    update_quadrature_points | update_values | update_JxW_values;

  QGauss<dim - 1> quad(fe.degree + 1);

  FEFaceValues<dim>      fv1(fe, quad, update_flags);
  FEFaceValues<dim>      fv2(fe, quad, update_flags);
  FEInterfaceValues<dim> fiv(fe, quad, update_flags);

  std::vector<types::global_dof_index> v1(fe.dofs_per_cell);
  std::vector<types::global_dof_index> v2(fe.dofs_per_cell);

  const auto &cell1 = dofh.begin();
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
  fv1.reinit(cell1, internal_face);

  const auto &cell2 = cell1->neighbor(internal_face);
  cell2->get_dof_indices(v2);
  fv2.reinit(cell2, cell1->neighbor_of_neighbor(internal_face));

  // Now reinit the interface values
  fiv.reinit(cell1,
             internal_face,
             numbers::invalid_unsigned_int,
             cell2,
             cell1->neighbor_of_neighbor(internal_face),
             numbers::invalid_unsigned_int);

  // Now we have two cells, cell1 and cell2, and we want to integrate on the
  // face.
  //
  // \int_f [v_i](x1) * [v_j](x2)) dx1 dx2
  //

  FECouplingValues<dim> fcv(fv1,
                            fv2,
                            DoFCouplingType::contiguous,
                            QuadratureCouplingType::matching);

  deallog << "dofs_per_cell = " << fe.dofs_per_cell << std::endl
          << "n_q_points_per_cell = " << fv2.n_quadrature_points << std::endl
          << "n_q_points = " << fcv.n_quadrature_points() << std::endl
          << "n_interface_dofs = " << fiv.n_current_interface_dofs()
          << std::endl
          << "n_coupling_dofs = " << fcv.n_coupling_dofs() << std::endl;

  FullMatrix<double> matrix(fcv.n_coupling_dofs(), fcv.n_coupling_dofs());
  const FEValuesExtractors::Scalar scalar(0);

  const auto first  = fcv.get_first_extractor(scalar);
  const auto second = fcv.get_second_extractor(scalar);

  // We need to loop over all the coupling quadrature points
  for (const auto q : fcv.quadrature_point_indices())
    {
      const auto &[x, y] = fcv.quadrature_point(q);
      for (const auto i : fcv.coupling_dof_indices())
        {
          const auto &lvi       = fcv[first].value(i, q);
          const auto &rvi       = fcv[second].value(i, q);
          const auto  jump_i    = lvi - rvi;
          const auto &jump_ifiv = fiv[scalar].jump_in_values(i, q);


          for (const auto j : fcv.coupling_dof_indices())
            {
              const auto &lvj = fcv[first].value(j, q);
              const auto &rvj = fcv[second].value(j, q);

              const auto  jump_j    = lvj - rvj;
              const auto &jump_jfiv = fiv[scalar].jump_in_values(j, q);

              deallog << std::left << "fcv -- "                           //
                      << "[vi][" << i << "]: " << std::setw(10) << jump_i //
                      << "[vj][" << j << "]: " << std::setw(10) << jump_j //
                      << "x[" << q << "]: " << std::setw(10) << x         //
                      << "  y[" << q << "]: " << std::setw(10) << y       //
                      << "JxW[" << q << "]: " << std::setw(10) << fcv.JxW(q)
                      << std::endl;
              deallog << std::left << "fiv -- "                              //
                      << "[vi][" << i << "]: " << std::setw(10) << jump_ifiv //
                      << "[vj][" << j << "]: " << std::setw(10) << jump_jfiv //
                      << std::endl;
              matrix(i, j) += jump_i * jump_j * fcv.JxW(q);
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
