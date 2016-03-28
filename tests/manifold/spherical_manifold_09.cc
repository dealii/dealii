// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
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

// Check SphericalManifold on faces.

#include "../tests.h"

#include <deal.II/base/utilities.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/mapping_q_generic.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/manifold_lib.h>

// Check that the spherical manifold finds the right intermediate
// points both using face quadratures, as well as using face
// quadratures projected to cells.

#include <numeric>

template<int dim, int spacedim>
void test()
{
  deallog << "dim=" << dim << ", spacedim=" << spacedim << std::endl;

  Triangulation<dim, spacedim>   triangulation;

  GridGenerator::hyper_cube (triangulation, 2.0, 4.0);

  typename Triangulation<dim,spacedim>::active_cell_iterator
  cell = triangulation.begin_active();

  // Center and radius of the Ball
  Point<spacedim> center = cell->center();
  double radius = center.distance(cell->vertex(0));

  static  const SphericalManifold<dim,spacedim> manifold(cell->center());

  triangulation.set_all_manifold_ids(0);
  triangulation.set_manifold (0, manifold);

  const QGauss<dim-1> face_quad(3);

  const FE_Q<dim-1> face_fe_q(1);

  // Compute the points with the faces
  deallog << "Face quadratures ("
          << face_quad.size() * GeometryInfo<dim>::faces_per_cell
          << ")" << std::endl;

  for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
    {
      std::vector<Point<spacedim> > vertices;
      std::vector<double> weights(GeometryInfo<dim>::vertices_per_face);

      for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_face; ++i)
        vertices.push_back(cell->face(f)->vertex(i));

      for (unsigned int i=0; i<face_quad.size(); ++i)
        {
          for (unsigned int v=0; v<weights.size(); ++v)
            weights[v] = face_fe_q.shape_value(v, face_quad.point(i));

          deallog << manifold.get_new_point(Quadrature<spacedim>(vertices, weights))
                  << std::endl;
        }
    }

  // Project the face quads
  Quadrature<dim> quad = QProjector<dim>::project_to_all_faces(face_quad);

  deallog << "Face quadratures projected on cell "
          << "(" << quad.size() << ")" << std::endl;

  std::vector<Point<spacedim> > vertices;
  for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
    vertices.push_back(cell->vertex(i));

  // Used to compute weights.
  FE_Q<dim> fe_q(1);

  std::vector<double> weights(GeometryInfo<dim>::vertices_per_cell);

  for (unsigned int i=0; i<quad.size(); ++i)
    {
      for (unsigned int v=0; v<weights.size(); ++v)
        weights[v] = fe_q.shape_value(v, quad.point(i));

      deallog << manifold.get_new_point(Quadrature<spacedim>(vertices, weights))
              << std::endl;
    }

  // for(unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f) {
  //   std::vector<Point<spacedim> > vertices;
  //   std::vector<double> weights(GeometryInfo<dim>::vertices_per_face);

  //   for(unsigned int i=0; i<GeometryInfo<dim>::vertices_per_face; ++i)
  //     vertices.push_back(cell->face(f)->vertex(i));

  //   for (unsigned int i=0; i<face_quad.size(); ++i)
  //   {
  //     for(unsigned int v=0; v<weights.size(); ++v)
  //  weights[v] = face_fe_q.shape_value(v, face_quad.point(i));

  //     deallog << manifold.get_new_point(Quadrature<spacedim>(vertices, weights))
  //        << std::endl;
  //   }
  // }
}


int
main()
{
  std::ofstream logfile ("output");
  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);

  test<2,2>();
  test<2,3>();

  // stest<3,3>();

  return 0;
}



