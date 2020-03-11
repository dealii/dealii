// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2018-2014 by the deal.II authors
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



// tests matrix-free face evaluation in a very simple form

#include <deal.II/base/function.h>

#include <deal.II/fe/fe_dgq.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



template <int dim, int fe_degree, typename number>
class MatrixFreeTest
{
public:
  MatrixFreeTest(const MatrixFree<dim, number> &data)
    : data(data)
  {}

  void
  action(Vector<number> &src) const
  {
    data.loop(&MatrixFreeTest::dummy,
              &MatrixFreeTest::dummy,
              &MatrixFreeTest::local_apply_boundary_face,
              this,
              src,
              src);
  }

private:
  void
  dummy(const MatrixFree<dim, number> &,
        Vector<number> &,
        const Vector<number> &,
        const std::pair<unsigned int, unsigned int> &) const
  {}
  void
  local_apply_boundary_face(
    const MatrixFree<dim, number> &data,
    Vector<number> &,
    const Vector<number> &                       src,
    const std::pair<unsigned int, unsigned int> &face_range) const
  {
    FEFaceEvaluation<dim, fe_degree, fe_degree + 1, 1, number> fe_eval(data,
                                                                       true);
    for (unsigned int face = face_range.first; face < face_range.second; face++)
      {
        fe_eval.reinit(face);
        fe_eval.read_dof_values(src);
        fe_eval.evaluate(true, true);

        // Only one vectorization component is filled in this case because we
        // only have one cell (otherwise the output will not be stable among
        // systems...)
        deallog << "Face " << face << ": ";
        for (unsigned int q = 0; q < fe_eval.n_q_points; ++q)
          {
            deallog << fe_eval.get_value(q)[0] << " "
                    << fe_eval.get_normal_derivative(q)[0] << "   ";
          }
        deallog << std::endl;
      }
  }

  const MatrixFree<dim, number> &data;
};



template <int dim>
class BoundaryFunction : public Function<dim>
{
public:
  virtual double
  value(const Point<dim> &p, const unsigned int) const
  {
    return p[0] + 0.2 * p[1];
  }
};



template <int dim>
Point<dim>
grid_transform(const Point<dim> &in)
{
  Point<dim> out;
  for (unsigned int d = 0; d < dim; ++d)
    out[d] = in[d] + 0.75 * d * in[0];
  return out;
}



template <int dim, int fe_degree>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, 0, 1);
  GridTools::transform(&grid_transform<dim>, tria);

  FE_DGQ<dim>     fe(fe_degree);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);
  AffineConstraints<double> constraints;
  constraints.close();

  MatrixFree<dim, double> mf_data;
  {
    const QGauss<1>                                  quad(fe_degree + 1);
    typename MatrixFree<dim, double>::AdditionalData data;
    data.tasks_parallel_scheme = MatrixFree<dim, double>::AdditionalData::none;
    data.mapping_update_flags_inner_faces =
      (update_gradients | update_JxW_values);
    data.mapping_update_flags_boundary_faces =
      (update_gradients | update_JxW_values);

    mf_data.reinit(dof, constraints, quad, data);
  }
  MatrixFreeTest<dim, fe_degree, double> mf(mf_data);
  Vector<double>                         in(dof.n_dofs());
  VectorTools::interpolate(dof, BoundaryFunction<dim>(), in);
  mf.action(in);
}


int
main()
{
  initlog();

  deallog << std::setprecision(3);

  {
    deallog.push("2d");
    test<2, 1>();
    test<2, 2>();
    deallog.pop();
    deallog.push("3d");
    test<3, 1>();
    deallog.pop();
  }
}
