// ---------------------------------------------------------------------
//
// Copyright (C) 2014 by the Alexander Grayver & deal.II authors
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

// The test checks that project_boundary_values_curl_conforming
// works correctly for high-order FE_Nedelec elements used via
// FESystem. This requires the produced constraints to be the same
// for FE_Nedelec and FESystem(FE_Nedelec, 1).

#include "../tests.h"
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/numerics/vector_tools.h>

std::ofstream logfile ("output");

template <int dim>
class BoundaryFunction: public Function<dim>
{
public:
  BoundaryFunction ();
  virtual void vector_value (const Point<dim> &p, Vector<double> &values) const;
};

template <int dim>
BoundaryFunction<dim>::BoundaryFunction (): Function<dim> (dim)
{
}

template <int dim>
void BoundaryFunction<dim>::vector_value (const Point<dim> &,
                                          Vector<double> &values) const
{
  for (unsigned int d = 0; d < dim; ++d)
    values (d) = d + 1.0;
}

template <int dim>
void test_boundary_values (const FiniteElement<dim> &fe, ConstraintMatrix &constraints)
{
  Triangulation<dim> triangulation;
  GridGenerator::subdivided_hyper_cube (triangulation, 2);
  DoFHandler<dim> dof_handler (triangulation);
  dof_handler.distribute_dofs (fe);
  BoundaryFunction<dim> boundary_function;
  constraints.clear ();
  VectorTools::project_boundary_values_curl_conforming (dof_handler, 0, boundary_function, 0, constraints);
  constraints.close ();
}

template <int dim>
void test(unsigned order)
{
  deallog << "dim:" << dim << " order:" << order << "\t";

  ConstraintMatrix constraints_fe, constraints_fes;

  {
    FE_Nedelec<3> fe (order);
    test_boundary_values (fe, constraints_fe);
  }

  {
    FESystem<3> fe(FE_Nedelec<3>(order),1);
    test_boundary_values (fe, constraints_fes);
  }

  if(constraints_fes.n_constraints() == constraints_fe.n_constraints())
  {
    const IndexSet& lines = constraints_fes.get_local_lines ();

    for(unsigned i = 0; i < lines.n_elements(); ++i)
    {
      if(!constraints_fe.is_constrained(lines.nth_index_in_set(i)))
      {
        deallog << "Failed" << std::endl;
        return;
      }

      const std::vector<std::pair<types::global_dof_index,double> >& c1
              = *constraints_fes.get_constraint_entries(lines.nth_index_in_set(i));
      const std::vector<std::pair<types::global_dof_index,double> >& c2
              = *constraints_fe.get_constraint_entries(lines.nth_index_in_set(i));

      for(size_t j = 0; j < c1.size(); ++j)
        if((c1[j].first != c2[j].first) || (fabs(c1[j].second - c2[j].second) > 1e-14))
        {
          deallog << "Failed" << std::endl;
          return;
        }
    }
  }
  else
  {
    deallog << "Failed" << std::endl;
    return;
  }

  deallog << "OK" << std::endl;
}

int main ()
{
  deallog << std::setprecision (2);
  deallog.attach (logfile);
  deallog.depth_console (0);
  deallog.threshold_double (1e-12);

  test<2>(0);
  test<2>(1);
  test<2>(2);

  test<3>(0);
  test<3>(1);
  //test<3>(2);
}
