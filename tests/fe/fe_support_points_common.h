// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#include "../tests.h"
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_q_hierarchical.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q_dg0.h>
#include <deal.II/fe/fe_q_iso_q1.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_bdm.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>



template <int dim>
inline void
check_support (const FiniteElement<dim> &finel, const char *name)
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr, 0., 1.);
  DoFHandler<dim> dof (tr);
  dof.distribute_dofs (finel);

  deallog << name << '<' << dim << '>' << " cell support points" << std::endl;

  const std::vector<Point<dim> > &cell_points = finel.get_unit_support_points ();

  for (unsigned int k=0; k<cell_points.size(); ++k)
    deallog << std::setprecision(3) << cell_points[k] << std::endl;

  const std::vector<Point<dim-1> > &face_points = finel.get_unit_face_support_points ();
  const std::vector<double> dummy_weights (face_points.size());

  Quadrature<dim-1> q(face_points, dummy_weights);

  for (unsigned int i=0; i<GeometryInfo<dim>::faces_per_cell; ++i)
    {
      std::vector<Point<dim> > q_points (q.get_points().size());
      QProjector<dim>::project_to_face (q, i, q_points);
      Quadrature<dim> qp(q_points);
      deallog << name << '<' << dim << '>' << " face " << i
              << " support points" << std::endl;

      for (unsigned int k=0; k<face_points.size(); ++k)
        deallog << std::setprecision(3) << qp.point(k)
                << std::endl;
    }
}


#define CHECK_ALL(EL,deg,dim) { FE_ ## EL<dim> EL(deg); check_support(EL, #EL #deg); }
#define CHECK_SYS1(sub1,N1,dim) { FESystem<dim> q(sub1, N1); check_support(q, #sub1 #N1); }
#define CHECK_SYS2(sub1,N1,sub2,N2,dim) { FESystem<dim> q(sub1, N1, sub2, N2); \
    check_support(q, #sub1 #N1 #sub2 #N2); }
#define CHECK_SYS3(sub1,N1,sub2,N2,sub3,N3,dim) { FESystem<dim> q(sub1, N1, sub2, N2, sub3, N3); \
    check_support(q, #sub1 #N1 #sub2 #N2 #sub3 #N3); }


