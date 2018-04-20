// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2018 by the deal.II authors
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

// like grid_transform, but use a spatially variable coefficient


#include "../tests.h"
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/fe/mapping_q.h>



template <int dim>
class Coefficient : public Function<dim>
{
public:
  virtual double value (const Point<dim> &p,
                        const unsigned int) const
  {
    return (p[0]>0 ? 10 : 1);
  }
};


int main ()
{
  const unsigned int dim=2;
  Point<dim> origin;
  SphericalManifold<dim> boundary(origin);
  MappingQ<dim> mapping(2);
  Triangulation<dim> tria;
  const double inner_radius=1.;
  const double outer_radius=5.;
  GridGenerator::hyper_shell(tria, origin, inner_radius, outer_radius, 8);
  tria.set_all_manifold_ids(numbers::flat_manifold_id);
  GridTools::copy_boundary_to_manifold_id(tria);
  tria.set_manifold(0, boundary);
  tria.refine_global(2);

  // build up a map of vertex indices
  // of boundary vertices to the new
  // boundary points
  std::map<unsigned int,Point<dim> > new_points;

  // new center and new radius
  // of the inner circle.
  const Point<dim> n_center(0,-1);
  const double n_radius=0.5;

  Triangulation<dim>::cell_iterator cell=tria.begin_active(),
                                    endc=tria.end();
  Triangulation<dim>::face_iterator face;
  for (; cell!=endc; ++cell)
    {
      if (cell->at_boundary())
        for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
          {
            face=cell->face(face_no);
            if (face->at_boundary())
              for (unsigned int vertex_no=0;
                   vertex_no<GeometryInfo<dim>::vertices_per_face; ++vertex_no)
                {
                  const Point<dim> &v=face->vertex(vertex_no);
                  if (std::fabs(std::sqrt(v.square())-outer_radius)<1e-12)
                    {
                      // leave the
                      // point, where
                      // they are.
                      new_points.insert(std::pair<types::global_dof_index, Point<dim> > (
                                          face->vertex_index(vertex_no), v));
                    }
                  else if (std::fabs(std::sqrt(v.square())-inner_radius)<1e-12)
                    {
                      // move the
                      // center of
                      // the inner
                      // circle to
                      // (-1,0) and
                      // take half
                      // the radius
                      // of the
                      // circle.
                      new_points.insert(std::pair<types::global_dof_index, Point<dim> > (
                                          face->vertex_index(vertex_no), n_radius/inner_radius*v+n_center));
                      face->set_manifold_id(1);
                    }
                  else
                    Assert(false, ExcInternalError());
                }
          }
    }

  Coefficient<dim> c;
  GridTools::laplace_transform (new_points, tria, &c, true);
  SphericalManifold<dim> inner_ball(n_center);
  tria.set_manifold(1, inner_ball);

  GridOut grid_out;
  std::ofstream eps_stream2("output");
  grid_out.write_eps(tria, eps_stream2, &mapping);

  tria.clear();

  return 0;
}
