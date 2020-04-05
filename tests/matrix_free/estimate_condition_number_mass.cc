// ---------------------------------------------------------------------
//
// Copyright (C) 2013 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



// this function computes condition number estimates for the mass matrix at
// different polynomial degrees. The mesh uses a hypercube mesh with no
// hanging nodes and no other constraints

#include <deal.II/base/function_lib.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_q_hierarchical.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


void
output_double_number(double input, const std::string &text)
{
  deallog << text << input << std::endl;
}

template <int dim, int fe_degree, typename Number>
void
mass_operator(const MatrixFree<dim, Number> &              data,
              Vector<Number> &                             dst,
              const Vector<Number> &                       src,
              const std::pair<unsigned int, unsigned int> &cell_range)
{
  FEEvaluation<dim, fe_degree, fe_degree + 1, 1, Number> fe_eval(data);
  const unsigned int n_q_points = fe_eval.n_q_points;

  for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      fe_eval.reinit(cell);
      fe_eval.read_dof_values(src);
      fe_eval.evaluate(true, false, false);
      for (unsigned int q = 0; q < n_q_points; ++q)
        {
          fe_eval.submit_value(fe_eval.get_value(q), q);
        }
      fe_eval.integrate(true, false);
      fe_eval.distribute_local_to_global(dst);
    }
}



template <int dim, int fe_degree, typename Number>
class MatrixFreeTest
{
public:
  typedef VectorizedArray<Number> vector_t;

  MatrixFreeTest(const MatrixFree<dim, Number> &data_in)
    : data(data_in){};

  void
  vmult(Vector<Number> &dst, const Vector<Number> &src) const
  {
    dst = 0;
    const std::function<void(const MatrixFree<dim, Number> &,
                             Vector<Number> &,
                             const Vector<Number> &,
                             const std::pair<unsigned int, unsigned int> &)>
      wrap = mass_operator<dim, fe_degree, Number>;
    data.cell_loop(wrap, dst, src);
  };

private:
  const MatrixFree<dim, Number> &data;
};



template <int dim, int fe_degree>
void
test(const FiniteElement<dim> &fe, const unsigned int n_iterations)
{
  typedef double number;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);

  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);
  AffineConstraints<double> constraints;
  constraints.close();

  deallog << "Testing " << dof.get_fe().get_name() << std::endl;

  MatrixFree<dim, number> mf_data;
  {
    const QGauss<1>                                  quad(fe_degree + 1);
    typename MatrixFree<dim, number>::AdditionalData data;
    data.tasks_block_size = 8 / VectorizedArray<number>::size();
    mf_data.reinit(dof, constraints, quad, data);
  }

  MatrixFreeTest<dim, fe_degree, number> mf(mf_data);
  Vector<number>                         in(dof.n_dofs()), out(dof.n_dofs());

  VectorTools::create_right_hand_side(dof,
                                      QGauss<dim>(fe_degree + 1),
                                      Functions::CosineFunction<dim>(),
                                      in);

  // prescribe number of iterations in CG instead of tolerance. This is needed
  // to get the same test output on different platforms where roundoff errors
  // accumulate differently. Beware of this strange solver setting when seeing
  // "failure" in the output
  SolverControl control(n_iterations, 0);
  SolverCG<>    solver(control);
  solver.connect_condition_number_slot(
    std::bind(output_double_number,
              std::placeholders::_1,
              "Condition number estimate: "));
  try
    {
      solver.solve(mf, out, in, PreconditionIdentity());
    }
  catch (...)
    {}
  deallog << std::endl;
}


int
main()
{
  initlog();

  deallog << std::setprecision(2);

  {
    // iterations taken from results at tolerance 1e-9
    deallog.push("2d");
    test<2, 1>(FE_Q<2>(1), 15);
    test<2, 1>(FE_DGQ<2>(1), 3);
    test<2, 2>(FE_Q<2>(2), 34);
    test<2, 2>(FE_DGQ<2>(2), 6);
    test<2, 4>(FE_Q<2>(4), 56);
    test<2, 4>(FE_Q<2>(QIterated<1>(QTrapez<1>(), 4)), 75);
    test<2, 4>(FE_DGQ<2>(4), 18);
    test<2, 4>(FE_Q_Hierarchical<2>(4), 736);
    test<2, 6>(FE_Q<2>(6), 77);
    test<2, 6>(FE_Q<2>(QIterated<1>(QTrapez<1>(), 6)), 161);
    test<2, 6>(FE_DGQ<2>(6), 32);
    test<2, 6>(FE_DGQArbitraryNodes<2>(QIterated<1>(QTrapez<1>(), 6)), 38);
    test<2, 10>(FE_Q<2>(10), 108);
    test<2, 10>(FE_Q<2>(QIterated<1>(QTrapez<1>(), 10)), 782);
    test<2, 10>(FE_DGQ<2>(10), 70);
    test<2, 10>(FE_DGQArbitraryNodes<2>(QIterated<1>(QTrapez<1>(), 10)), 118);
    test<2, 16>(FE_Q<2>(16), 156);
    test<2, 25>(FE_Q<2>(25), 200);
    deallog.pop();
    deallog.push("3d");
    test<3, 1>(FE_Q<3>(1), 37);
    test<3, 2>(FE_Q<3>(2), 80);
    test<3, 5>(FE_Q<3>(5), 187);
    test<3, 5>(FE_Q<3>(QIterated<1>(QTrapez<1>(), 5)), 494);
    deallog.pop();
  }
}
