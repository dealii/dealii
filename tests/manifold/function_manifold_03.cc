//----------------------------  function_manifold_chart ---------------------------
//    Copyright (C) 2011 - 2015 by the mathLab team.
//
//    This file is subject to LGPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------- function_manifold_chart ---------------------------


// Test a simple parabolic manifold, including gradients and tangent vector

#include "../tests.h"
#include <fstream>
#include <deal.II/base/logstream.h>


// all include files you need here
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_out.h>

// Helper function
template <int dim, int spacedim>
void test(unsigned int ref=1)
{
  deallog << "Testing dim " << dim
          << ", spacedim " << spacedim << std::endl;

  // Here the only allowed axis is z. In cylinder the default is x.
  std::string push_forward_expression;
  std::string pull_back_expression;

  switch (spacedim)
    {
    case 2:
      push_forward_expression = "x; x^2";
      pull_back_expression = "x";
      break;
    case 3:
      push_forward_expression = "x; x^2; 0";
      pull_back_expression = "x";
      break;
    default:
      Assert(false, ExcInternalError());
    }

  FunctionManifold<dim,spacedim,1> manifold(push_forward_expression,
                                            pull_back_expression);

  // Two points and two weights
  std::vector<Point<spacedim> > p(2);
  p[1][0] = 1.0;
  p[1][1] = 1.0;
  std::vector<double> w(2);

  unsigned int n_intermediates = 16;

  deallog << "P0: " << p[0]
          << ", P1: " << p[1] << std::endl;

  for (unsigned int i=0; i<n_intermediates+1; ++i)
    {
      w[0] = 1.0-(double)i/((double)n_intermediates);
      w[1] = 1.0 - w[0];

      Point<spacedim> ip = manifold.get_new_point(Quadrature<spacedim>(p, w));
      Tensor<1,spacedim> t1 = manifold.get_tangent_vector(ip, p[0]);
      Tensor<1,spacedim> t2 = manifold.get_tangent_vector(ip, p[1]);

      deallog << "P: " << ip
              << ", T(P, P0): " << t1
              << ", T(P, P1): " << t2 << std::endl;

    }
}

int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);


  test<2,2>();
  test<2,3>();
  test<3,3>();

  return 0;
}

