// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2013 by the deal.II authors
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


// Test projection onto lines

#include "../tests.h"
#include <iomanip>
#include <fstream>

#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/geometry_info.h>
#include <cmath>

template<int dim>
void check_line(Quadrature<1> &quadrature)
{
  Point<dim> p1;
  Point<dim> p2;
  p1(0) = 1.;
  p2(0) = 7.;
  if (dim>1)
    {
      p1(1) = 3;
      p2(1) = -5.;
    }
  if (dim>2)
    {
      p1(2) = 0;
      p2(2) = 10.;
    }
  Quadrature<dim> q = QProjector<dim>::project_to_line(quadrature, p1, p2);
  double s = 0.;

  for (unsigned int k=0; k<q.size(); ++k)
    {
      deallog << k << '\t' << q.point(k) << std::endl;
      s += q.weight(k);
    }
  deallog << "length: " << s << std::endl;
}

template<int dim>
void check_face(Quadrature<1> &q1)
{
  deallog << "Checking dim " << dim
          << " 1d-points " << q1.size()
          << std::endl;

  Quadrature<dim-1> subquadrature(q1);
  for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
    {
      deallog << "Face " << f
              << std::endl;

      Quadrature<dim> quadrature
        = QProjector<dim>::project_to_face(subquadrature, f);
      for (unsigned int k=0; k<quadrature.size(); ++k)
        deallog << quadrature.point(k) << std::endl;
    }

  for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
    for (unsigned int s=0; s<GeometryInfo<dim>::max_children_per_face; ++s)
      {
        deallog << "Face " << f << " subface " << s
                << std::endl;

        Quadrature<dim> quadrature
          = QProjector<dim>::project_to_face(subquadrature, f);
        for (unsigned int k=0; k<quadrature.size(); ++k)
          deallog << quadrature.point(k) << std::endl;
      }
}

template<int dim>
void check_faces (Quadrature<1> &q1)
{
  const unsigned int nq = q1.size();

  deallog << "Checking dim " << dim
          << " 1d-points " << nq
          << std::endl;


  Quadrature<dim-1> subquadrature(q1);
  const unsigned int nqs = subquadrature.size();

  Quadrature<dim> faces = QProjector<dim>::project_to_all_faces(subquadrature);

  for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
    {

      deallog << "Face " << f
              << " orientation false"
              << std::endl;

      unsigned int
      offset = QProjector<dim>::DataSetDescriptor::face(f, false, false, false, nqs);

      for (unsigned int k=0; k<nqs; ++k)
        deallog << faces.point(offset+k) << std::endl;

      deallog << "Face " << f
              << " orientation true"
              << std::endl;

      offset = QProjector<dim>::DataSetDescriptor::face(f, true, false, false, nqs);

      for (unsigned int k=0; k<nqs; ++k)
        deallog << faces.point(offset+k) << std::endl;
    }

  /*
  Quadrature<dim> subs = QProjector<dim>::project_to_all_subfaces(subquadrature);


  for (unsigned int f=0;f<GeometryInfo<dim>::faces_per_cell;++f)
  for (unsigned int s=0;s<GeometryInfo<dim>::max_children_per_face;++s)
  {

  deallog << "Face " << f << " subface " << s
  << " orientation false"
  << std::endl;

  unsigned int
  offset = QProjector<dim>::DataSetDescriptor::subface(f, s, false, nqs);

  for (unsigned int k=0;k<nqs;++k)
  deallog << faces.point(offset+k) << std::endl;

  deallog << "Face " << f << " subface " << s
  << " orientation true"
  << std::endl;

  offset = QProjector<dim>::DataSetDescriptor::subface(f, s, true, nqs);

  for (unsigned int k=0;k<nqs;++k)
  deallog << faces.point(offset+k) << std::endl;
  }
  */
}


void check(Quadrature<1> &q)
{
  deallog << std::endl;
  deallog.push("line");
  check_line<1> (q);
  check_line<2> (q);
  check_line<3> (q);
  deallog.pop();
  deallog.push("face");
  check_face<1>(q);
  check_face<2>(q);
  check_face<3>(q);
  deallog.pop();
  deallog.push("all");
  check_faces<2>(q);
  check_faces<3>(q);
  deallog.pop();
}

int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog << std::setprecision(2);

  deallog.threshold_double(1.e-10);

  Quadrature<1> none(0);
  check(none);

  QGauss<1> midpoint(1);
  check(midpoint);

  QTrapez<1> trapez;
  check(trapez);

  QSimpson<1> simpson;
  check(simpson);

  QMilne<1> milne;
  check(milne);
}
