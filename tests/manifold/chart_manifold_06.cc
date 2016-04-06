//-------------------------------------------------------------------
//    Copyright (C) 2016 by the deal.II authors.
//
//    This file is subject to LGPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-------------------------------------------------------------------


// Test direction vector of flat manifold with periodicity, where the
// flat manifold is implemented as a ChartManifold with identity
// pull-back and push-forward

#include "../tests.h"
#include <fstream>
#include <deal.II/base/logstream.h>
#include <deal.II/grid/manifold.h>



template <int dim, int spacedim>
class MyFlatManifold : public ChartManifold<dim,spacedim,spacedim>
{
public:
  MyFlatManifold (const Tensor<1,spacedim> &periodicity)
    :
    ChartManifold<dim,spacedim,spacedim>(periodicity)
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
void test()
{
  deallog << "Testing dim=" << dim
          << ", spacedim="<< spacedim << std::endl;

  // make the domain periodic in the first direction with periodicity 1.1
  Tensor<1,spacedim> periodicity;
  periodicity[0] = 1.1;
  MyFlatManifold<dim,spacedim> manifold(periodicity);

  Point<spacedim> x1, x2;
  for (unsigned int d=0; d<spacedim; ++d)
    {
      x1[d] = 0.1;
      x2[d] = 0.9;
    }

  // get the connecting vector between these two points. because we
  // have no periodicity, this should simply be the vector with
  // components all equal to 0.8 except for the first, which ought to be -0.3
  deallog << manifold.get_tangent_vector (x1, x2) << std::endl;

  // then also test the opposite direction
  deallog << manifold.get_tangent_vector (x2, x1) << std::endl;
}

int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);

  test<1,1>();
  test<1,2>();
  test<1,3>();
  test<2,2>();
  test<2,3>();
  test<3,3>();

  return 0;
}

