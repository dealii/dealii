//-------------------------------------------------------------------
//    Copyright (C) 2016 by the deal.II authors.
//
//    This file is subject to LGPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-------------------------------------------------------------------

// Test that a Manifold with projection only works

#include "../tests.h"
#include <fstream>
#include <deal.II/base/logstream.h>


// all include files you need here
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_out.h>

template <int dim, int spacedim>
class MyManifold : public Manifold<dim, spacedim>
{
public:
  Point<spacedim>
  project_to_manifold (const std::vector<Point<spacedim> > &vertices,
                       const Point<spacedim> &candidate) const
  {
    // Shift the y coordinate to 4*x*(1-x)
    Point<spacedim> p = candidate;
    if (spacedim > 1)
      p[1] = 4*p[0]*(1-p[0]);
    return p;
  }
};

// Helper function
template <int dim, int spacedim>
void test(unsigned int ref=1)
{
  deallog << "Testing dim=" << dim
          << ", spacedim="<< spacedim << std::endl;
  MyManifold<dim,spacedim> manifold;

  // Test that we get the right thing with two simple points.
  Point<spacedim> p0, p1;
  p1[0] = 1.0;

  Point<spacedim> p2 = manifold.get_intermediate_point(p0,p1,.5);
  deallog << "p2(0.5): " << p2 << std::endl;

  p2 = manifold.get_intermediate_point(p0,p1,.2);
  deallog << "p2(0.2): " << p2 << std::endl;
}

int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);

  test<1,1>();
  test<1,2>();
  test<2,2>();
  test<2,3>();
  test<3,3>();

  return 0;
}

