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


// Test the different FEValuesBase::get_function_values

#include "../tests.h"
#include <fstream>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>

#include <deal.II/lac/vector.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_tools.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/grid/grid_generator.h>


// Call this function with a system consisting of several copies of
// the SAME element
template<int dim>
void vector_values(const FiniteElement<dim> &fe)
{
  Assert(fe.n_base_elements() == 1, ExcNotImplemented());
  deallog.push(fe.get_name());

  QTrapez<dim> quadrature;
  std::vector<unsigned int> renumbering(fe.dofs_per_cell);
  std::vector<std::vector<unsigned int> > component_start;
  FETools::compute_component_wise(fe, renumbering, component_start);

  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);
  tr.refine_global(1);

  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);
  DoFRenumbering::component_wise(dof);

  Vector<float> v(dof.n_dofs());
  for (unsigned int i=0; i<v.size(); ++i)
    v(i) = i;

  FEValues<dim> feval(fe, quadrature, update_values);
  std::vector<Vector<double> > local(quadrature.size(),
                                     Vector<double>(fe.n_components()));

  typename DoFHandler<dim>::active_cell_iterator cell = dof.begin_active();
  const typename DoFHandler<dim>::active_cell_iterator end = dof.end();

  unsigned int cell_no = 0;
  while (cell != end)
    {
      deallog << "Cell " << cell_no++ << std::endl;
      feval.reinit(cell);
      feval.get_function_values(v, local);
      for (unsigned int c=0; c<fe.n_components(); ++c)
        {
          deallog << "Component " << c;
          for (unsigned int k=0; k<quadrature.size(); ++k)
            deallog << '\t' << (int) local[k](c);
          deallog << std::endl;
        }
      ++cell;
    }
  deallog.pop();
}

template<int dim>
void test_vectors()
{
  FE_Q<dim> q1(1);

  FESystem<dim> q1_3(q1,3);
  vector_values(q1_3);
}


int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test_vectors<2>();

  return 0;
}
