//-------------------------------------------------------------------
//    Copyright (C) 2017 - 2018 by the deal.II authors.
//
//    This file is subject to LGPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-------------------------------------------------------------------


// Test that transfinite interpolation manifold returns valid results on the
// same geometry as transfinite_manifold_05 (ball inside square, curved
// spherical surface). In an initial version, the line search in Newton for
// the transfinite interpolation would eagerly search too far outside the
// valid chart domain, leading to failures in the spherical manifold.

#include "../tests.h"
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>


int main ()
{
  initlog();

  const int dim = 3;
  Triangulation<dim> tria1, tria2, tria;
  GridGenerator::hyper_shell(tria1, Point<dim>(), 0.4, std::sqrt(dim), 6);
  GridGenerator::hyper_ball(tria2, Point<dim>(), 0.4);
  GridGenerator::merge_triangulations(tria1, tria2, tria);
  tria.set_all_manifold_ids(0);
  for (typename Triangulation<dim>::cell_iterator cell = tria.begin();
       cell != tria.end(); ++cell)
    {
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
        {
          bool face_at_sphere_boundary = true;
          for (unsigned int v=0; v<GeometryInfo<dim-1>::vertices_per_cell; ++v)
            if (std::abs(cell->face(f)->vertex(v).norm()-0.4) > 1e-12)
              face_at_sphere_boundary = false;
          if (face_at_sphere_boundary)
            cell->face(f)->set_all_manifold_ids(1);
        }
    }
  static const SphericalManifold<dim> spherical_manifold;
  tria.set_manifold(1, spherical_manifold);
  static TransfiniteInterpolationManifold<dim> transfinite;
  transfinite.initialize(tria);
  tria.set_manifold(0, transfinite);

  tria.refine_global(1);

  deallog.precision(10);
  deallog << "Cell centers" << std::endl;
  for (auto cell : tria.cell_iterators())
    deallog << cell->id() << " has center " << cell->center(/*respect_manifold*/true) << std::endl;

  return 0;
}
