// ---------------------------------------------------------------------
//
// Copyright (C) 2013 - 2015 by the deal.II authors
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
  GridGenerator::hyper_cube(tr, 2., 5.);
  typename DoFHandler<dim>::cell_iterator c = dof.begin();
  dof.distribute_dofs(finel);

  const unsigned int div = 1;

  QTrapez<dim> q;
//  QIterated<dim> q(q_trapez, div);
  FEValues<dim> fe(mapping, finel, q, UpdateFlags(update_gradients
                                                  | update_second_derivatives));
  fe.reinit(c);

  unsigned int k=0;
  for (unsigned int mz=0; mz<=((dim>2) ? div : 0) ; ++mz)
    {
      for (unsigned int my=0; my<=((dim>1) ? div : 0) ; ++my)
        {
          for (unsigned int mx=0; mx<=div; ++mx)
            {
              deallog << q.point(k) << std::endl;

              for (unsigned int i=0; i<finel.dofs_per_cell; ++i)
                {
                  deallog << "\tGrad " << fe.shape_grad(i,k);
                  deallog << "\t2nd " << fe.shape_hessian(i,k);
                  deallog << std::endl;
                }
              k++;
            }
        }
    }
  deallog.pop ();
}



template<int dim>
void plot_FE_Q_shape_functions()
{
  MappingQGeneric<dim> m(1);
  FE_Q<dim> q1(1);
  plot_derivatives(m, q1, "Q1");
//  plot_face_shape_functions(m, q1, "Q1");
  FE_Q<dim> q2(2);
  plot_derivatives(m, q2, "Q2");
//  plot_face_shape_functions(m, q2, "Q2");
  FE_Q<dim> q3(QIterated<1>(QTrapez<1>(),3));
  plot_derivatives(m, q3, "Q3");
//  plot_face_shape_functions(m, q3, "Q3");
  FE_Q<dim> q4(QIterated<1>(QTrapez<1>(),4));
  plot_derivatives(m, q4, "Q4");
//  plot_face_shape_functions(m, q4, "Q4");
//    FE_Q<dim> q5(5);
//    plot_derivatives(m, q5, "Q5");
//    FE_Q<dim> q6(6);
//    plot_derivatives(m, q6, "Q6");
//    FE_Q<dim> q7(7);
//    plot_derivatives(m, q7, "Q7");
//    FE_Q<dim> q8(8);
//    plot_derivatives(m, q8, "Q8");
//    FE_Q<dim> q9(9);
//    plot_derivatives(m, q9, "Q9");
//    FE_Q<dim> q10(10);
//    plot_derivatives(m, q10, "Q10");
}


template<int dim>
void plot_FE_DGQ_shape_functions()
{
  MappingQGeneric<dim> m(1);
  FE_DGQ<dim> q1(1);
  plot_derivatives(m, q1, "DGQ1");
//  plot_face_shape_functions(m, q1, "DGQ1");
  FE_DGQ<dim> q2(2);
  plot_derivatives(m, q2, "DGQ2");
//  plot_face_shape_functions(m, q2, "DGQ2");
  FE_DGQArbitraryNodes<dim> q3(QIterated<1>(QTrapez<1>(),3));
  plot_derivatives(m, q3, "DGQ3");
//  plot_face_shape_functions(m, q3, "DGQ3");
  FE_DGQArbitraryNodes<dim> q4(QIterated<1>(QTrapez<1>(),4));
  plot_derivatives(m, q4, "DGQ4");
//  plot_face_shape_functions(m, q4, "DGQ4");
//    FE_DGQ<dim> q5(5);
//    plot_derivatives(m, q5, "DGQ5");
//    FE_DGQ<dim> q6(6);
//    plot_derivatives(m, q6, "DGQ6");
//    FE_DGQ<dim> q7(7);
//    plot_derivatives(m, q7, "DGQ7");
//    FE_DGQ<dim> q8(8);
//    plot_derivatives(m, q8, "DGQ8");
//    FE_DGQ<dim> q9(9);
//    plot_derivatives(m, q9, "DGQ9");
//    FE_DGQ<dim> q10(10);
//    plot_derivatives(m, q10, "DGQ10");
}


int
main()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision(8);
  deallog << std::fixed;
  deallog.attach(logfile);

  deallog.push ("1d");
  plot_FE_Q_shape_functions<1>();
  deallog.pop ();

  deallog.push ("2d");
  plot_FE_Q_shape_functions<2>();
  plot_FE_DGQ_shape_functions<2>();
  deallog.pop ();

  deallog.push ("3d");
//  plot_FE_Q_shape_functions<3>();
  deallog.pop ();


  // FESystem test.
  MappingQGeneric<2> m(1);
  FESystem<2> q2_q3(FE_Q<2>(2), 1,
                    FE_Q<2>(QIterated<1>(QTrapez<1>(),3)), 1);
//  plot_derivatives(m, q2_q3, "Q2_Q3");

  return 0;
}
