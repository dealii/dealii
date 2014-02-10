//----------------------------  manifold_id_01.cc  ---------------------------
//    Copyright (C) 2011, 2013 by the mathLab team.
//
//    This file is subject to LGPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  flat_manifold_01.cc  ---------------------------


// Test periodicity conditions using FlatManifold

#include "../tests.h"
#include <fstream>
#include <base/logstream.h>


// all include files you need here
#include <deal.II/grid/manifold.h>
#include <deal.II/grid/manifold_lib.h>

template <int spacedim> 
void test() {
  Point<spacedim> period;
  for(unsigned int i=0; i<spacedim; ++i)
    period[i] = 1.0;
  FlatManifold<spacedim> manifold(period);
  deallog << "Testing spacedim = " << spacedim << std::endl;
  deallog << "Period: " << period << std::endl;
  std::vector<Point<spacedim> > p(2);
  std::vector<double> w(2, 0.5);
  double eps = .1;
  p[0] = period*eps*2;
  p[1] = period*(1.0-eps);
  Point<spacedim> middle = manifold.get_new_point(p, w);
  deallog << "P(0): " << p[0] << std::endl;
  deallog << "P(1): " << p[1] << std::endl;
  deallog << "Middle: " << middle << std::endl;
}

int main () 
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-5);
  
  test<1>();
  test<2>();
  test<3>();

  return 0;
}
                  
