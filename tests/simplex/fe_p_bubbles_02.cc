// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Verify that FE_SimplexP_Bubbles can be used with a lumped mass matrix by
// computing a convergence rate.

#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_simplex_p_bubbles.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/reference_cell.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools_integrate_difference.h>
#include <deal.II/numerics/vector_tools_interpolate.h>

#include "../tests.h"


template <int dim, int spacedim = dim>
void
test_interpolate()
{
  deallog << std::endl << "test interpolation" << std::endl;
  Functions::CosineFunction<dim> func(1);

  const unsigned int                        n_refinements = 7 - dim;
  std::vector<Triangulation<dim, spacedim>> trias;
  for (unsigned int refinement_n = 0; refinement_n < n_refinements;
       ++refinement_n)
    {
      Triangulation<dim, spacedim> tria;
      Triangulation<dim, spacedim> tria_hypercube;
      GridGenerator::subdivided_hyper_cube(tria_hypercube,
                                           std::pow(2, refinement_n));
      if (dim == spacedim)
        GridTools::distort_random(0.2, tria_hypercube);
      GridGenerator::convert_hypercube_to_simplex_mesh(tria_hypercube, tria);
      trias.emplace_back(std::move(tria));
    }

  deallog << "dim = " << dim << std::endl;
  for (unsigned int degree = 1; degree < 3; ++degree)
    {
      deallog << "degree = " << degree << std::endl;
      double old_error = -1.0;
      for (unsigned int refinement_n = 0; refinement_n < n_refinements - degree;
           ++refinement_n)
        {
          const Triangulation<dim, spacedim> &tria = trias[refinement_n];
          deallog << "number of cells = " << tria.n_active_cells() << std::endl;

          FE_SimplexP_Bubbles<dim, spacedim> fe(degree);

          const ReferenceCell       type = fe.reference_cell();
          DoFHandler<dim, spacedim> dh(tria);
          dh.distribute_dofs(fe);
          deallog << "number of dofs = " << dh.n_dofs() << std::endl;

          const Mapping<dim, spacedim> &map =
            type.template get_default_linear_mapping<dim, spacedim>();

          Vector<double> solution(dh.n_dofs());
          VectorTools::interpolate(map, dh, func, solution);

          QGaussSimplex<dim> error_quad(4);
          Vector<double>     out_l2(tria.n_active_cells());
          VectorTools::integrate_difference(
            map, dh, solution, func, out_l2, error_quad, VectorTools::L2_norm);
          const double new_error =
            VectorTools::compute_global_error(tria,
                                              out_l2,
                                              VectorTools::L2_norm);
          deallog << "error = " << new_error << std::endl;
          if (old_error != -1.0)
            deallog << "ratio = " << old_error / new_error << std::endl;
          old_error = new_error;
        }
    }
}



template <int dim, int spacedim = dim>
void
test_lumped_project()
{
  deallog << std::endl << "test lumped project" << std::endl;
  const auto func = Functions::CosineFunction<dim>(1);

  const unsigned int                        n_refinements = 7 - dim;
  std::vector<Triangulation<dim, spacedim>> trias;
  for (unsigned int refinement_n = 0; refinement_n < n_refinements;
       ++refinement_n)
    {
      Triangulation<dim, spacedim> tria;
      Triangulation<dim, spacedim> tria_hypercube;
      GridGenerator::subdivided_hyper_cube(tria_hypercube,
                                           std::pow(2, refinement_n));
      if (dim == spacedim)
        GridTools::distort_random(0.2, tria_hypercube);
      GridGenerator::convert_hypercube_to_simplex_mesh(tria_hypercube, tria);
      trias.emplace_back(std::move(tria));
    }


  deallog << "dim = " << dim << std::endl;
  for (unsigned int degree = 1; degree < 3; ++degree)
    {
      deallog << "degree = " << degree << std::endl;
      double old_error = -1.0;
      for (unsigned int refinement_n = 0; refinement_n < n_refinements;
           ++refinement_n)
        {
          const Triangulation<dim, spacedim> &tria = trias[refinement_n];
          deallog << "number of cells = " << tria.n_active_cells() << std::endl;

          FE_SimplexP_Bubbles<dim, spacedim> fe(degree);

          const ReferenceCell       type = fe.reference_cell();
          DoFHandler<dim, spacedim> dh(tria);
          dh.distribute_dofs(fe);
          deallog << "number of dofs = " << dh.n_dofs() << std::endl;
          const Quadrature<dim> nodal_quad =
            FETools::compute_nodal_quadrature(fe);
          const Quadrature<dim> cell_quad = QGaussSimplex<dim>(
            std::max<unsigned int>(fe.tensor_degree() + 1, 2));

          Vector<double> lumped_mass(dh.n_dofs());
          Vector<double> consistent_rhs(dh.n_dofs());

          const Mapping<dim, spacedim> &map =
            type.template get_default_linear_mapping<dim, spacedim>();

          FEValues<dim> lumped_fev(map,
                                   fe,
                                   nodal_quad,
                                   update_values | update_JxW_values);
          FEValues<dim> consistent_fev(map,
                                       fe,
                                       cell_quad,
                                       update_quadrature_points |
                                         update_values | update_JxW_values);

          std::vector<types::global_dof_index> dofs(fe.dofs_per_cell);
          for (const auto &cell : dh.active_cell_iterators())
            {
              lumped_fev.reinit(cell);
              consistent_fev.reinit(cell);
              cell->get_dof_indices(dofs);

              for (unsigned int q = 0; q < nodal_quad.size(); ++q)
                {
                  for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
                    {
                      const double v = lumped_fev.shape_value(i, q);
                      Assert(std::abs(v - double(i == q)) < 1e-14,
                             ExcInternalError());
                      lumped_mass[dofs[i]] +=
                        lumped_fev.shape_value(i, q) * lumped_fev.JxW(q);
                    }
                }

              for (unsigned int q = 0; q < cell_quad.size(); ++q)
                for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
                  consistent_rhs[dofs[i]] +=
                    consistent_fev.shape_value(i, q) *
                    func.value(consistent_fev.quadrature_point(q)) *
                    consistent_fev.JxW(q);
            }

          Vector<double> solution(dh.n_dofs());
          for (std::size_t i = 0; i < solution.size(); ++i)
            solution[i] = consistent_rhs[i] / lumped_mass[i];

          QGaussSimplex<dim> error_quad(4);
          Vector<double>     out_l2(tria.n_active_cells());
          VectorTools::integrate_difference(
            map, dh, solution, func, out_l2, error_quad, VectorTools::L2_norm);

          const auto new_error =
            VectorTools::compute_global_error(tria,
                                              out_l2,
                                              VectorTools::L2_norm);
          deallog << "error = " << new_error << std::endl;
          if (old_error != -1.0)
            deallog << "ratio = " << old_error / new_error << std::endl;
          old_error = new_error;

#if 0
          static unsigned int counter = 0;
          for (unsigned int n_subdivisions = 1; n_subdivisions <= 2;
               ++n_subdivisions)
            {
              DataOut<dim> data_out;

              data_out.attach_dof_handler(dh);
              data_out.add_data_vector(solution, "solution",
                                       DataOut<dim>::type_dof_data);
              data_out.build_patches(map, n_subdivisions);
              std::ofstream output("test." + std::to_string(dim) + "." +
                                   std::to_string(counter++) + ".vtk");
              data_out.write_vtk(output);
            }
#endif
        }
    }
}



int
main()
{
  initlog();

  test_interpolate<2, 2>();
  test_interpolate<3, 3>();

  test_lumped_project<2, 2>();
  test_lumped_project<3, 3>();
}
