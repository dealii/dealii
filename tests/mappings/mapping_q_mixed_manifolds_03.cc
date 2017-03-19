// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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


// Prints the mapped points (quadrature points) of the mesh in
// mapping_q_mixed_manifolds_01, once using smoothing and once without.

#include "../tests.h"
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/fe/mapping_q_generic.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/grid/manifold_lib.h>
#include <vector>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>



template <int dim>
void test()
{
  Point<dim> direction;
  direction[dim-1] = 1.;

  std_cxx11::shared_ptr<Manifold<dim> > cylinder_manifold
  (dim == 2 ? static_cast<Manifold<dim>*>(new SphericalManifold<dim>(Point<dim>())) :
   static_cast<Manifold<dim>*>(new CylindricalManifold<dim>(direction, Point<dim>())));

  // create cube and shift it to position 0.7
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, -0.5, 0.5);
  Tensor<1,dim> shift;
  shift[0] = 1.;
  GridTools::shift(shift, tria);
  tria.begin()->face(0)->set_all_manifold_ids(1);
  tria.set_manifold(1, *cylinder_manifold);

  FE_Nothing<dim> fe;
  for (unsigned int degree = 6; degree < 7; ++degree)
    {
      QIterated<dim> quad(QTrapez<1>(), degree);
      {
        MappingQGeneric<dim> mapping(degree, true);
        FEValues<dim> fe_values(mapping, fe, quad, update_quadrature_points);
        fe_values.reinit(tria.begin());
        deallog << "Points with smoothing in " << dim << "D" << std::endl;
        for (unsigned int q=0; q<quad.size(); ++q)
          deallog << fe_values.quadrature_point(q) << std::endl;
      }
      {
        MappingQGeneric<dim> mapping(degree, false);
        FEValues<dim> fe_values(mapping, fe, quad, update_quadrature_points);
        fe_values.reinit(tria.begin());
        deallog << "Points without smoothing in " << dim << "D" << std::endl;
        for (unsigned int q=0; q<quad.size(); ++q)
          deallog << fe_values.quadrature_point(q) << std::endl;
      }
    }
}


int main()
{
  std::ofstream logfile ("output");
  deallog.attach(logfile);
  deallog.threshold_double(1.e-14);

  test<2>();
  test<3>();
}
