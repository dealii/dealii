// ---------------------------------------------------------------------
//
// Copyright (C) 2013 by the deal.II authors
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


#include "../tests.h"
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/fe_values.h>
#include <vector>
#include <fstream>
#include <string>


template<int dim>
inline void
plot_derivatives(Mapping<dim> &mapping,
                 FiniteElement<dim> &finel,
                 const char *name)
{
  deallog.push (name);

  Triangulation<dim> tr;
  DoFHandler<dim> dof(tr);
  GridGenerator::hyper_cube(tr);
  typename DoFHandler<dim>::cell_iterator c = dof.begin();
  dof.distribute_dofs(finel);

  QGauss<dim-1> q(1);
//  QIterated<dim> q(q_trapez, div);
  FEFaceValues<dim> fe(mapping, finel, q, UpdateFlags(update_gradients
                                                      | update_second_derivatives));
  for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
    {
      fe.reinit(c, face);

      for (unsigned int k=0; k<q.size(); ++k)
        {
          deallog << "Face " << face
                  << " Point " << q.point(k) << std::endl;
          for (unsigned int i=0; i<finel.dofs_per_cell; ++i)
            {
              deallog << "\tGrad " << fe.shape_grad(i,k);
              deallog << "\t2nd " << fe.shape_hessian(i,k);
              deallog << std::endl;
            }
        }
    }
  deallog.pop ();
}



template<int dim>
void plot_FE_Q_shape_functions()
{
  MappingQ1<dim> m;
//  FE_Q<dim> q1(1);
//  plot_derivatives(m, q1, "Q1");
//  plot_face_shape_functions(m, q1, "Q1");
  FE_Q<dim> q2(2);
  plot_derivatives(m, q2, "Q2");
  FE_Q<dim> q3(3);
  plot_derivatives(m, q3, "Q3");
  FE_Q<dim> q4(4);
  plot_derivatives(m, q4, "Q4");
}


template<int dim>
void plot_FE_DGQ_shape_functions()
{
  MappingQ1<dim> m;
  FE_DGQ<dim> q1(1);
  plot_derivatives(m, q1, "DGQ1");
  FE_DGQ<dim> q2(2);
  plot_derivatives(m, q2, "DGQ2");
  FE_DGQ<dim> q3(3);
  plot_derivatives(m, q3, "DGQ3");
  FE_DGQ<dim> q4(4);
  plot_derivatives(m, q4, "DGQ4");
}


int
main()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision(2);
  deallog << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console(0);

  deallog.push ("2d");
  plot_FE_Q_shape_functions<2>();
//  plot_FE_DGQ_shape_functions<2>();
  deallog.pop ();

  deallog.push ("3d");
//  plot_FE_Q_shape_functions<3>();
  deallog.pop ();
  return 0;
}



