// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// tests matrix-free face evaluation with the option to compress
// This is the same as face_value_optimization but uses fe_degree=-1

#include <deal.II/base/function.h>

#include <deal.II/fe/fe_dgq.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include "../tests.h"


template <int dim, int fe_degree, typename number>
class MatrixFreeTest
{
public:
  MatrixFreeTest(const MatrixFree<dim, number> &data)
    : data(data)
    , error(2)
  {}

  ~MatrixFreeTest()
  {
    deallog << "Error between variants: " << error[0] / error[1] << std::endl;
  }

  void
  check_error(Vector<number> &src) const
  {
    data.loop(&MatrixFreeTest::local_apply,
              &MatrixFreeTest::local_apply_face,
              &MatrixFreeTest::local_apply_face,
              this,
              src,
              src);
  }

private:
  void
  local_apply(const MatrixFree<dim, number> &,
              Vector<number> &,
              const Vector<number> &,
              const std::pair<unsigned int, unsigned int> &) const
  {}

  void
  local_apply_face(
    const MatrixFree<dim, number> &data,
    Vector<number> &,
    const Vector<number>                        &src,
    const std::pair<unsigned int, unsigned int> &face_range) const
  {
    FEFaceEvaluation<dim, -1, 0, 1, number> ref(data, true);
    FEFaceEvaluation<dim, -1, 0, 1, number> check(data, true);

    for (unsigned int face = face_range.first; face < face_range.second; ++face)
      {
        ref.reinit(face);
        check.reinit(face);

        ref.read_dof_values(src);
        ref.evaluate(EvaluationFlags::values);
        check.gather_evaluate(src, EvaluationFlags::values);

        for (unsigned int q = 0; q < ref.n_q_points; ++q)
          {
            VectorizedArray<number> diff =
              (ref.get_value(q) - check.get_value(q));
            for (unsigned int v = 0; v < VectorizedArray<number>::size(); ++v)
              {
                if (std::abs(diff[v]) > 1e-12)
                  {
                    deallog << "Error detected on face" << face << ", v=" << v
                            << '!' << std::endl;
                    deallog << "ref: ";
                    for (unsigned int i = 0; i < ref.dofs_per_cell; ++i)
                      deallog << ref.get_dof_value(i)[v] << ' ';
                    deallog << std::endl;
                    deallog << "done: " << check.get_value(q)[v]
                            << " instead of " << ref.get_value(q)[v]
                            << std::endl;

                    deallog
                      << data.get_face_info(face).cells_interior[v] << ' '
                      << (int)data.get_face_info(face).interior_face_no << ' '
                      << (int)data.get_face_info(face).face_orientation << ' '
                      << (int)data.get_face_info(face).subface_index
                      << std::endl;
                    deallog << std::endl;
                  }
                error[0] += std::abs(diff[v]);
                error[1] += std::abs(ref.get_value(q)[v]);
              }
          }
      }

    FEFaceEvaluation<dim, fe_degree, fe_degree + 1, 1, number> refr(data,
                                                                    false);
    FEFaceEvaluation<dim, fe_degree, fe_degree + 1, 1, number> checkr(data,
                                                                      false);

    for (unsigned int face = face_range.first;
         face < std::min(data.n_inner_face_batches(), face_range.second);
         face++)
      {
        refr.reinit(face);
        checkr.reinit(face);

        refr.read_dof_values(src);
        refr.evaluate(EvaluationFlags::values);
        checkr.gather_evaluate(src, EvaluationFlags::values);

        for (unsigned int q = 0; q < ref.n_q_points; ++q)
          {
            VectorizedArray<number> diff =
              (refr.get_value(q) - checkr.get_value(q));
            for (unsigned int v = 0; v < VectorizedArray<number>::size(); ++v)
              {
                if (std::abs(diff[v]) > 1e-12)
                  {
                    deallog << "Error detected on face" << face << ", v=" << v
                            << '!' << std::endl;
                    deallog << "ref: ";
                    for (unsigned int i = 0; i < ref.dofs_per_cell; ++i)
                      deallog << refr.get_dof_value(i)[v] << ' ';
                    deallog << std::endl;
                    deallog << "done: " << check.get_value(q)[v]
                            << " instead of " << ref.get_value(q)[v]
                            << std::endl;

                    deallog
                      << data.get_face_info(face).cells_exterior[v] << ' '
                      << (int)data.get_face_info(face).exterior_face_no << ' '
                      << (int)data.get_face_info(face).face_orientation << ' '
                      << (int)data.get_face_info(face).subface_index
                      << std::endl;
                    deallog << std::endl;
                  }
                error[0] += std::abs(diff[v]);
                error[1] += std::abs(refr.get_value(q)[v]);
              }
          }
      }
  }

  const MatrixFree<dim, number> &data;
  mutable std::vector<double>    error;
};



template <int dim, int fe_degree>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_ball(tria);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  typename Triangulation<dim>::active_cell_iterator cell, endc;
  cell = tria.begin_active();
  endc = tria.end();
  for (; cell != endc; ++cell)
    if (cell->center().norm() < 0.5)
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  tria.begin(tria.n_levels() - 1)->set_refine_flag();
  tria.last()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  tria.refine_global(1);
  cell = tria.begin_active();
  for (unsigned int i = 0; i < 10 - 3 * dim; ++i)
    {
      cell                 = tria.begin_active();
      endc                 = tria.end();
      unsigned int counter = 0;
      for (; cell != endc; ++cell, ++counter)
        if (counter % (7 - i) == 0)
          cell->set_refine_flag();
      tria.execute_coarsening_and_refinement();
    }

  FE_DGQ<dim>     fe(fe_degree);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);
  AffineConstraints<double> constraints;
  constraints.close();

  using number = double;

  MatrixFree<dim, number> mf_data;
  {
    const QGauss<1>                                  quad(fe_degree + 1);
    typename MatrixFree<dim, number>::AdditionalData data;
    data.tasks_parallel_scheme = MatrixFree<dim, number>::AdditionalData::none;
    data.tasks_block_size      = 3;
    data.mapping_update_flags_inner_faces =
      (update_gradients | update_JxW_values);
    data.mapping_update_flags_boundary_faces =
      (update_gradients | update_JxW_values);

    mf_data.reinit(MappingQ1<dim>{}, dof, constraints, quad, data);
  }

  MatrixFreeTest<dim, fe_degree, number> mf(mf_data);

  Vector<number> in(dof.n_dofs());

  for (unsigned int i = 0; i < dof.n_dofs(); ++i)
    {
      if (constraints.is_constrained(i))
        continue;
      const double entry = Testing::rand() / (double)RAND_MAX;
      in(i)              = entry;
    }

  mf.check_error(in);
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
    test<3, 2>();
    deallog.pop();
  }
}
