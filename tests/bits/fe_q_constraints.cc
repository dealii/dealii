// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
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



// This file implements a simple test procedure for the algorithm
// which constructs the constraint matrices for FE_Q<3> elements
// of arbitrary order. After interpolating a polynomial with the
// degree of the ansatzspace, the hanging node constraints are
// distributed to the solution vector. This procedure should not
// change the solution vector. Hence the difference between the
// exact function and the FE function is measured before and
// after the distribution of the hanging node constraints.

#include "../tests.h"
#include <deal.II/fe/fe_q.h>

#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/polynomial.h>

#include <deal.II/numerics/vector_tools.h>

#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>

#include <fstream>
#include <iomanip>

template <int dim>
class TestFunction : public Function<dim>
{
private:
  std::vector<Polynomials::Polynomial<double> > base;
  const unsigned int p_order;

public:
  TestFunction (const unsigned int p_order);

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component = 0) const;
};


template <int dim>
TestFunction<dim>::TestFunction (const unsigned int p_order) :
  p_order (p_order)
{
  std::vector<double> coeff(p_order);

  for (unsigned int d = 0; d < dim; ++d)
    {
      for (unsigned int c = 0; c < p_order; ++c)
        coeff[c] = (double) Testing::rand() / (double) RAND_MAX;
      base.push_back (Polynomials::Polynomial<double> (coeff));
    }
}


template <int dim>
double TestFunction<dim>::value (const Point<dim>   &p,
                                 const unsigned int  /*component*/) const
{
  double val = base[0].value(p(0));
  for (unsigned int i = 1; i < dim; ++i)
    val *= base[i].value(p(i));
  return val;
}


template <int dim>
class TestFEQConstraints
{
public:
  TestFEQConstraints (unsigned int p_order, unsigned int refinements);
  void run ();

private:
  const unsigned int refinements;
  const unsigned int p_order;

  void refine_grid_random ();
  void make_grid_and_dofs ();

  void test ();

  Triangulation<dim>     triangulation;
  FE_Q<dim>              fe;
  DoFHandler<dim>        dof_handler;
  ConstraintMatrix     hanging_node_constraints;
};

template <int dim>
TestFEQConstraints<dim>::TestFEQConstraints (unsigned int p_order,
                                             unsigned int refinements) :
  refinements (refinements),
  p_order (p_order),
  fe (p_order),
  dof_handler (triangulation)
{
}


// Actually this function creates with pseudo-random numbers errors,
// which are then used to refine the grid. This should be a tough test for
// the hanging nodes.
template <int dim>
void TestFEQConstraints<dim>::refine_grid_random ()
{
  const unsigned int n_cells = triangulation.n_active_cells();
  Vector<double> estimated_error_per_cell (n_cells);

  for (unsigned int i = 0; i < n_cells; ++i)
    estimated_error_per_cell(i) = (double) Testing::rand() / (double) RAND_MAX;

  GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                   estimated_error_per_cell,
                                                   0.3, 0.03);
  triangulation.execute_coarsening_and_refinement ();
}


template <int dim>
void TestFEQConstraints<dim>::make_grid_and_dofs ()
{
  GridGenerator::hyper_cube (triangulation, -1, 1);
  triangulation.refine_global (2);

  for (unsigned int i = 0; i < refinements; ++i)
    refine_grid_random ();

  dof_handler.distribute_dofs (fe);

  deallog << "---------------------------------------------------------"
          << std::endl;
  deallog << "P-Order: " << p_order
          << "  Number of degrees of freedom: "
          << dof_handler.n_dofs()
          << std::endl;

  hanging_node_constraints.clear ();
  DoFTools::make_hanging_node_constraints (dof_handler,
                                           hanging_node_constraints);
  hanging_node_constraints.close ();
}


template <int dim>
void TestFEQConstraints<dim>::test ()
{
  TestFunction<dim> test_function (p_order);
  double l2test,
         l2norm1,
         l2norm2,
         l2error;

  Vector<double> solution;
  solution.reinit (dof_handler.n_dofs ());
  hanging_node_constraints.distribute (solution);

  QGauss<dim> quadrature (p_order + 1);

  Vector<double> norm_per_cell (triangulation.n_active_cells ());
  VectorTools::interpolate (dof_handler, test_function, solution);

  // First error check. Simply the interpolation error of the used FE-Space
  // on the given triangulation.
  VectorTools::integrate_difference (dof_handler, solution,
                                     ZeroFunction<dim>(1),
                                     norm_per_cell,
                                     quadrature,
                                     VectorTools::L2_norm);
  l2test = norm_per_cell.l2_norm ();
  deallog << "L2-Norm of test function "
          << l2test << std::endl;

  // First error check. Simply the interpolation error of the used FE-Space
  // on the given triangulation.
  VectorTools::integrate_difference (dof_handler, solution,
                                     test_function,
                                     norm_per_cell,
                                     quadrature,
                                     VectorTools::L2_norm);
  l2norm1 = norm_per_cell.l2_norm () / l2test;

  // Second error check. Interpolation error, after redistribution of the
  // values onto the DoFs on the hanging nodes.
  hanging_node_constraints.distribute (solution);
  VectorTools::integrate_difference (dof_handler, solution,
                                     test_function,
                                     norm_per_cell,
                                     quadrature,
                                     VectorTools::L2_norm);
  l2norm2 = norm_per_cell.l2_norm () / l2test;
  l2error = fabs (l2norm1 - l2norm2);
  deallog << "Normed L2-Error 1: " << l2norm1
          << "  Normed L2-Error 2: " << l2norm2 << std::endl
          << "Normed L2-Diff: " << l2error << "   ";
  if (l2error < 1.0e-16)
    deallog << "OK" << std::endl;
  else
    deallog << "Error !" << std::endl;
}


template <int dim>
void TestFEQConstraints<dim>::run ()
{
  make_grid_and_dofs ();
  test ();
}


int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  unsigned int ref_level[] = {4,3,2};

  // Initialise the random generator with an
  // arbitrary number. This should ensure
  // reproducible results.
  Testing::srand(4375384);

  // With these parameters, the elements are
  // tested up to a polynomial degree of
  // 6. The triangulation is refined 2 times
  // uniformly and several times (pseudo-)
  // randomly.
  TestFEQConstraints<3> *test;
  for (unsigned int p = 1; p < 4; ++p)
    {
      test = new TestFEQConstraints<3> (p, ref_level[p-1]);
      test->run ();
      delete test;
    }

  return 0;
}
