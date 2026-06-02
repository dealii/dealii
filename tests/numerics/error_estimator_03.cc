// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// Print some values for KellyErrorEstimator with a ComponentMask, AMR, and
// multiple input vectors to verify recent patches don't change the output.


#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/error_estimator.h>

#include "../tests.h"

template <int dim>
class TestFunction : public Function<dim>
{
public:
  TestFunction(const unsigned int n_components)
    : Function<dim>(n_components)
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const override
  {
    return (p[component % dim] + 1.0) * (p.norm() + 1.0) * (component + 2.0);
  }
};



template <int dim>
void
check()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, -1, 1, true);

  for (unsigned int i = 0; i < 4; ++i)
    {
      tria.begin_active()->set_refine_flag();
      tria.execute_coarsening_and_refinement();
    }

  DoFHandler<dim> dof_handler(tria);
  FESystem<dim>   fe(
    FE_Q<dim>(2), 1, FE_Q<dim>(1), 2, FE_Q<dim>(1), 2, FE_Q<dim>(1), 1);
  dof_handler.distribute_dofs(fe);

  deallog << "Number of components = " << fe.n_components() << std::endl;
  ComponentMask component_mask(fe.n_components(), false);
  component_mask.set(0, true);
  component_mask.set(3, true);
  component_mask.set(4, true);

  TestFunction<dim> function(fe.n_components());
  std::map<types::boundary_id, const Function<dim> *> neumann_bc;
  neumann_bc[0] = &function;

  std::vector<Vector<double>> dof_vectors;
  std::vector<Vector<float>>  error_vectors;
  for (unsigned int i = 0; i < 3; ++i)
    {
      dof_vectors.emplace_back(dof_handler.n_dofs());
      for (types::global_dof_index j = 0; j < dof_handler.n_dofs(); ++j)
        dof_vectors.back()[j] = random_value();

      error_vectors.emplace_back(tria.n_active_cells());
    }

  std::vector<const ReadVector<double> *> dof_vector_ptrs;
  std::vector<Vector<float> *>            error_vector_ptrs;
  for (unsigned int i = 0; i < 3; ++i)
    {
      dof_vector_ptrs.push_back(&dof_vectors[i]);
      error_vector_ptrs.push_back(&error_vectors[i]);
    }

  ArrayView<const ReadVector<double> *> dof_vector_view(dof_vector_ptrs);
  ArrayView<Vector<float> *>            error_view(error_vector_ptrs);
  KellyErrorEstimator<dim>::estimate(MappingQ<dim>(1),
                                     dof_handler,
                                     QGauss<dim - 1>(3),
                                     neumann_bc,
                                     dof_vector_view,
                                     error_view,
                                     component_mask);

  for (const auto &vec : error_vectors)
    {
      deallog << "estimated error:" << std::endl;
      vec.print(deallog.get_file_stream(), 8);
      deallog << std::endl;
    }

  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();
  deallog << std::setprecision(2) << std::fixed;

  deallog.push("1d");
  check<1>();
  deallog.pop();
  deallog.push("2d");
  check<2>();
  deallog.pop();
  deallog.push("3d");
  check<3>();
  deallog.pop();
}
