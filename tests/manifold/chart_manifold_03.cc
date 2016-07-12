//-------------------------------------------------------------------
//    Copyright (C) 2016 by the deal.II authors.
//
//    This file is subject to LGPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-------------------------------------------------------------------


// Test periodicity of FlatManifold, where the
// flat manifold is implemented as a ChartManifold with identity
// pull-back and push-forward

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
class MyFlatManifold : public ChartManifold<dim,spacedim,spacedim>
{
public:
  MyFlatManifold (const Tensor<1,spacedim> &periodicity)
    :
    ChartManifold<dim,spacedim,spacedim> (periodicity)
  {}


  virtual
  Point<spacedim>
  pull_back(const Point<spacedim> &space_point) const
  {
    return space_point;
  }


  virtual
  Point<spacedim>
  push_forward(const Point<spacedim> &chart_point) const
  {
    return chart_point;
  }

  virtual
  DerivativeForm<1,spacedim,spacedim>
  push_forward_gradient(const Point<spacedim> &chart_point) const
  {
    DerivativeForm<1,spacedim,spacedim> x;
    for (unsigned int d=0; d<spacedim; ++d)
      x[d][d] = 1;
    return x;
  }
};



// Helper function
template <int dim, int spacedim>
void test(unsigned int ref=1)
{
  deallog << "Testing dim=" << dim
          << ", spacedim="<< spacedim << std::endl;

  Tensor<1,spacedim> periodicity;
  periodicity[0] = 5.0;

  MyFlatManifold<dim,spacedim> manifold(periodicity);

  Quadrature<spacedim> quad;
  std::vector<std::vector<Point<spacedim> > >ps(10,std::vector<Point<spacedim> >(2));
  Point<spacedim> middle;
  std::vector<double > ws(2, 0.5);

  // Case 1: both points are close to left boundary of periodicity
  ps[0][0][0] = 1;
  ps[0][1][0] = 2;
  // Case 2: same, with different order
  ps[1][0][0] = 2;
  ps[1][1][0] = 1;

  // Case 3: one is close to left, one to right
  ps[2][0][0] = 1;
  ps[2][1][0] = 4;
  // Case 4: same, opposite order
  ps[3][0][0] = 4;
  ps[3][1][0] = 1;

  // Case 5: both close to right
  ps[4][0][0] = 3;
  ps[4][1][0] = 4;
  // Case 6: same, opposite order
  ps[5][0][0] = 4;
  ps[5][1][0] = 3;

  // Case 7: both close to middle
  ps[6][0][0] = 2;
  ps[6][1][0] = 3;
  // Case 8: same, opposite order
  ps[7][0][0] = 3;
  ps[7][1][0] = 2;

  // Case 9: Corner cases
  ps[8][0][0] = -1e-10;
  ps[8][1][0] = 5+1e-10;
  // Case 10: same, opposite order
  ps[9][0][0] = 5+1e-10;
  ps[9][1][0] = -1e-10;

  for (unsigned int i=0; i<ps.size(); ++i)
    {
      quad = Quadrature<spacedim>(ps[i],ws);
      middle = manifold.get_new_point(quad);
      deallog << "P0: " << ps[i][0] << " , P1: " << ps[i][1] << " , Middle: " << middle << std::endl;
    }

}

int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.threshold_double(1.e-8);

  test<1,1>();
  test<1,2>();
  test<2,2>();
  test<2,3>();
  test<3,3>();

  return 0;
}

