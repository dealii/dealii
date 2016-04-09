//----------------------------  function_manifold_chart ---------------------------
//    Copyright (C) 2011 - 2015 by the mathLab team.
//
//    This file is subject to LGPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------- composition_manifold ---------------------------


// Stress periodicity in CompositionManifold. Compose SphericalManifold with
// the identity, and make sure periodicity is respected.

#include "../tests.h"
#include <fstream>
#include <deal.II/base/logstream.h>


// all include files you need here
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/composition_manifold.h>


int main ()
{
  initlog();
  std::ostream &out = deallog.get_file_stream();

  const int dim=2, spacedim=2;

  Point<spacedim> center;

  SphericalManifold<dim,spacedim> S(center);
  FunctionManifold<dim,spacedim,spacedim> F("x;y", "x;y");

  CompositionManifold<dim> manifold(S,F);

  // Chart points.
  Point<spacedim> cp[2];

  // Radius
  cp[0][0] = 1.0;
  cp[1][0] = 1.0;

  // Last point
  cp[0][1] = -numbers::PI/4;
  cp[1][1] =  numbers::PI/4;

  // Spacedim points
  std::vector<Point<spacedim> > sp(2);

  // Weights
  std::vector<double> w(2);

  sp[0] = manifold.push_forward(cp[0]);
  sp[1] = manifold.push_forward(cp[1]);

  unsigned int n_intermediates = 32;

  out << "set terminal aqua " << 0 << std::endl
      << "set size ratio -1" << std::endl
      << "plot '-' with vectors " << std::endl;

  for (unsigned int v=0; v<sp.size(); ++v)
    out << center << " "
        << sp[v] << std::endl;


  for (unsigned int i=0; i<n_intermediates+1; ++i)
    {
      w[0] = 1.0-(double)i/((double)n_intermediates);
      w[1] = 1.0 - w[0];

      Point<spacedim> ip = manifold.get_new_point(Quadrature<spacedim>(sp, w));
      Tensor<1,spacedim> t1 = manifold.get_tangent_vector(ip, sp[0]);
      Tensor<1,spacedim> t2 = manifold.get_tangent_vector(ip, sp[1]);

      out << ip << " "
          << t2 << std::endl;
    }

  out << "e" << std::endl;

  out << "set terminal aqua " << 1 << std::endl
      << "set size ratio -1" << std::endl
      << "plot '-' w lp " << std::endl;

  for (unsigned int i=0; i<n_intermediates+1; ++i)
    {
      w[0] = 1.0-(double)i/((double)n_intermediates);
      w[1] = 1.0 - w[0];

      Point<spacedim> ip = manifold.
                           pull_back(manifold.get_new_point(Quadrature<spacedim>(sp, w)));

      ip[0] = w[1];

      out << ip << std::endl;
    }
  out << "e" << std::endl;

  return 0;
}
