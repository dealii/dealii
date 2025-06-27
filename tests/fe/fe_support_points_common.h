// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/logstream.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_bdm.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_q_bubbles.h>
#include <deal.II/fe/fe_q_dg0.h>
#include <deal.II/fe/fe_q_hierarchical.h>
#include <deal.II/fe/fe_q_iso_q1.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/reference_cell.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "../tests.h"



template <int dim>
inline void
check_support(const FiniteElement<dim> &finel, const char *name)
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr, 0., 1.);
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(finel);

  deallog << name << '<' << dim << '>' << " cell support points" << std::endl;

  const std::vector<Point<dim>> &cell_points = finel.get_unit_support_points();

  for (unsigned int k = 0; k < cell_points.size(); ++k)
    deallog << std::setprecision(3) << cell_points[k] << std::endl;

  const std::vector<Point<dim - 1>> &face_points =
    finel.get_unit_face_support_points();
  const std::vector<double> dummy_weights(face_points.size());

  Quadrature<dim - 1> q(face_points, dummy_weights);

  for (const unsigned int i : GeometryInfo<dim>::face_indices())
    {
      const auto qp = QProjector<dim>::project_to_face(
        ReferenceCells::get_hypercube<dim>(),
        q,
        i,
        numbers::default_geometric_orientation);
      deallog << name << '<' << dim << '>' << " face " << i << " support points"
              << std::endl;

      for (unsigned int k = 0; k < face_points.size(); ++k)
        deallog << std::setprecision(3) << qp.point(k) << std::endl;
    }

  deallog << name << '<' << dim << '>' << " cell generalized support points"
          << std::endl;

  const std::vector<Point<dim>> &cell_g_points =
    finel.get_generalized_support_points();

  for (unsigned int k = 0; k < cell_g_points.size(); ++k)
    deallog << std::setprecision(3) << cell_g_points[k] << std::endl;
}


#define CHECK_ALL(EL, deg, dim)  \
  {                              \
    FE_##EL<dim> EL(deg);        \
    check_support(EL, #EL #deg); \
  }
#define CHECK_SYS1(sub1, N1, dim) \
  {                               \
    FESystem<dim> q(sub1, N1);    \
    check_support(q, #sub1 #N1);  \
  }
#define CHECK_SYS2(sub1, N1, sub2, N2, dim) \
  {                                         \
    FESystem<dim> q(sub1, N1, sub2, N2);    \
    check_support(q, #sub1 #N1 #sub2 #N2);  \
  }
#define CHECK_SYS3(sub1, N1, sub2, N2, sub3, N3, dim) \
  {                                                   \
    FESystem<dim> q(sub1, N1, sub2, N2, sub3, N3);    \
    check_support(q, #sub1 #N1 #sub2 #N2 #sub3 #N3);  \
  }
