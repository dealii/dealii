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
#include <deal.II/fe/mapping_cartesian.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/fe_values.h>
#include <vector>
#include <fstream>
#include <string>


template <int dim>
void performance (Triangulation<dim> &tr,
                  const Mapping<dim> &mapping,
                  const FiniteElement<dim> &fe,
                  const Quadrature<dim> &quadrature,
                  UpdateFlags flags)
{
  deallog << "Create dofs" << std::endl;
  DoFHandler<dim> dof (tr);
  dof.distribute_dofs (fe);

  deallog << "Create FEValues" << std::endl;

  FEValues<dim> val (mapping, fe, quadrature, flags);

  DoFHandler<dim>::active_cell_iterator cell = dof.begin_active();
  DoFHandler<dim>::active_cell_iterator end = dof.end();

  deallog << "Loop" << std::endl;

  for (; cell != end ; ++cell)
    val.reinit(cell);

  deallog << "End" << std::endl;
}

template <int dim>
void loop (std::vector<FiniteElement<dim> *> elements,
           const Mapping<dim> &mapping,
           Triangulation<dim> &tr)
{
  QGauss<dim> gauss (4);

  typename std::vector<FiniteElement<dim> *>::iterator elementp = elements.begin ();
  typename std::vector<FiniteElement<dim> *>::iterator end = elements.end ();

  for (; elementp != end; ++elementp)
    {
      const FiniteElement<dim> &element = **elementp;

      char dofs[20];
      std::ostrstream ost (dofs, 19);
      ost << element.n_dofs_per_cell() << std::ends;

      deallog.push(dofs);

      deallog.push("points");
      performance (tr, mapping, element, gauss, update_q_points);
      deallog.pop();

      deallog.push("values");
      performance (tr, mapping, element, gauss, update_values);
      deallog.pop();

      deallog.push("grads-");
      performance (tr, mapping, element, gauss, update_gradients);
      deallog.pop();

      deallog.push("2nd---");
      performance (tr, mapping, element, gauss, update_second_derivatives);
      deallog.pop();

      deallog.push("matrix");
      performance (tr, mapping, element, gauss, update_q_points
                   | update_JxW_values
                   | update_values
                   | update_gradients);
      deallog.pop();

      deallog.push("all---");
      performance (tr, mapping, element, gauss, update_q_points
                   | update_JxW_values
                   | update_values
                   | update_gradients
                   | update_second_derivatives);
      deallog.pop();
      deallog.pop();
    }
}

int main ()
{
  std::ofstream of ("performance.log");
  deallog.attach (of);
  deallog.log_execution_time(true);
  deallog.log_time_differences(true);
  Triangulation<2> tr;
  GridGenerator::hyper_ball (tr);
  tr.refine_global (8);

  MappingCartesian<2> cartesian;
  MappingQ1<2> mapping;
  MappingQ<2> mappingq1(1);
  MappingQ<2> mappingq2(2);
  std::vector<FiniteElement<2>*> el2d;
  el2d.push_back (new FE_Q<2> (1));
  el2d.push_back (new FE_Q<2> (2));
  el2d.push_back (new FE_Q<2> (3));
  el2d.push_back (new FE_Q<2> (4));

  loop (el2d, cartesian, tr);
  loop (el2d, mapping, tr);
  loop (el2d, mappingq1, tr);
  loop (el2d, mappingq2, tr);
}
