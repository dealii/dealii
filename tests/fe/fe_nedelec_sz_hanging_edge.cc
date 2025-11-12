/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2025 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------
 */

// @sect3{Include files}

#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nedelec_sz.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/vector_tools_boundary.h>

#include <cmath>
#include <fstream>
#include <iostream>

#include "../tests.h"

namespace HangingEdge
{

  using namespace dealii;

  template <int dim>
  class RHSValues : public Function<dim>
  {
  public:
    RHSValues()
      : Function<dim>(dim)
    {}

    void
    vector_value(const Point<dim> &p, Vector<double> &values) const override;

    void
    vector_value_list(const std::vector<Point<dim>> &points,
                      std::vector<Vector<double>>   &value_list) const override
    {
      Assert(value_list.size() == points.size(),
             ExcDimensionMismatch(value_list.size(), points.size()));

      for (unsigned int p = 0; p < points.size(); p++)
        {
          RHSValues<dim>::vector_value(points[p], value_list[p]);
        }
    }
  };

  // here we are using artificial boundary values,
  // the benefit is, that we know the exact solution for the
  // electric field E, i.e. this is also the exact solution
  // to the partial differential equation we aim to solve
  template <>
  void
  RHSValues<3>::vector_value(const Point<3> &p, Vector<double> &values) const
  {
    values(0) = 0;
    values(1) = sin(p[2]);
    values(2) = sin(p[0]);
  }



  template <int dim>
  class MaxwellProblem
  {
  public:
    MaxwellProblem(unsigned int poly_degree, unsigned int n_cells);

    void
    run();

  private:
    void
    make_grid();
    void
    setup_system();
    void
    assemble_system();
    void
    solve();
    void
    print();
    double
    output_error();
    void
    output_results();

    Triangulation<dim> triangulation;
    DoFHandler<dim>    dof_handler;
    FE_NedelecSZ<dim>  fe;

    AffineConstraints<double> constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;
    Vector<double>       solution, system_rhs;

    const unsigned int poly_degree;

    const unsigned int n_cells;
    const double       inner_radius;
    const double       outer_radius;
    const double       z_thickness;
  };



  template <int dim>
  MaxwellProblem<dim>::MaxwellProblem(const unsigned int poly_degree,
                                      const unsigned int n_cells)
    : dof_handler(triangulation)
    , fe(poly_degree)
    , poly_degree(poly_degree)
    , n_cells(n_cells)
    , inner_radius(2.0 / 3.0)
    , outer_radius(1.0)
    , z_thickness(0.5)
  {}



  template <int dim>
  void
  MaxwellProblem<dim>::make_grid()
  {
    // Vertices:
    std::vector<Point<3>> vertices;

    // center vertex
    vertices.push_back(Point<3>(0, 0, 0));

    // outer vertices
    for (unsigned int i = 0; i < 2 * n_cells; ++i)
      {
        double radius = (i % 2 == 0) ? outer_radius : inner_radius;
        vertices.push_back(
          Point<3>(radius * std::sin((numbers::PI * i) / (1.0 * n_cells)),
                   radius * std::cos((numbers::PI * i) / (1.0 * n_cells)),
                   0.0));
      }

    vertices.push_back(Point<3>(0, 0, z_thickness));
    for (unsigned int i = 0; i < 2 * n_cells; ++i)
      {
        double radius = (i % 2 == 0) ? outer_radius : inner_radius;
        vertices.push_back(
          Point<3>(radius * std::sin((numbers::PI * i) / (1.0 * n_cells)),
                   radius * std::cos((numbers::PI * i) / (1.0 * n_cells)),
                   z_thickness));
      }


    // CellData:
    std::vector<CellData<3>> cells(n_cells);
    for (unsigned int i = 0; i < n_cells; ++i)
      {
        // set the material id
        cells[i].material_id = 0;

        // set the vertices
        cells[i].vertices[0] = 0;
        cells[i].vertices[1] = 2 + (2 * i);
        cells[i].vertices[2] =
          ((2 * n_cells - 1) + (2 * i)) % (2 * n_cells) + 1;
        cells[i].vertices[3] = 1 + (2 * i);

        const int offset     = 2 * n_cells + 1;
        cells[i].vertices[4] = offset + 0;
        cells[i].vertices[5] = offset + 2 + (2 * i);
        cells[i].vertices[6] =
          offset + ((2 * n_cells - 1) + (2 * i)) % (2 * n_cells) + 1;
        cells[i].vertices[7] = offset + 1 + (2 * i);
      }

    triangulation.create_triangulation(vertices, cells, SubCellData());

    // Refine one cell except one
    unsigned int cell_counter = 0;
    for (auto &cell : triangulation.active_cell_iterators())
      {
        ++cell_counter;
        if (cell_counter >= n_cells)
          break;
        cell->set_refine_flag();
      }
    triangulation.execute_coarsening_and_refinement();
  }



  template <int dim>
  void
  MaxwellProblem<dim>::setup_system()
  {
    dof_handler.distribute_dofs(fe);

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());

    DoFTools::make_sparsity_pattern(dof_handler, dsp);

    constraints.clear();
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);

    VectorTools::project_boundary_values_curl_conforming_l2(
      dof_handler,
      0 /* vector component*/,
      RHSValues<dim>(),
      0 /* boundary id*/,
      constraints);
    constraints.close();

    constraints.condense(dsp);

    sparsity_pattern.copy_from(dsp);
    system_matrix.reinit(sparsity_pattern);
    system_rhs.reinit(dof_handler.n_dofs());
    solution.reinit(dof_handler.n_dofs());
  }

  template <int dim>
  void
  MaxwellProblem<dim>::assemble_system()
  {
    system_matrix = 0;
    system_rhs    = 0;

    const unsigned int curl_dim = (dim == 2) ? 1 : 3;

    // choose the quadrature formulas
    QGauss<dim> quadrature_formula(fe.degree + 2);

    // get the number of quadrature points and dofs
    const unsigned int n_q_points    = quadrature_formula.size(),
                       dofs_per_cell = fe.dofs_per_cell;

    // set update flags
    FEValues<dim> fe_values(fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

    // Extractors for the real part
    const FEValuesExtractors::Vector E_re(0);

    // create the local left hand side and right hand side
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    // loop over all cells
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        if (cell->is_locally_owned() == false)
          continue;

        // initialize values:
        cell_matrix = 0;
        cell_rhs    = 0;
        fe_values.reinit(cell);

        for (const unsigned int i : fe_values.dof_indices())
          {
            // only compute this once
            std::vector<Tensor<1, dim>>      phi_i(n_q_points);
            std::vector<Tensor<1, curl_dim>> curl_phi_i(n_q_points);
            for (unsigned int q_point = 0; q_point < n_q_points; q_point++)
              {
                phi_i[q_point]      = fe_values[E_re].value(i, q_point);
                curl_phi_i[q_point] = fe_values[E_re].curl(i, q_point);
              }

            // we use here, that the problem is symmetrical
            for (unsigned int j = i; j < dofs_per_cell; j++)
              {
                double mass_part = 0;
                double curl_part = 0;

                for (unsigned int q_point = 0; q_point < n_q_points; q_point++)
                  {
                    Tensor<1, dim> phi_j = fe_values[E_re].value(j, q_point);
                    Tensor<1, curl_dim> curl_phi_j =
                      fe_values[E_re].curl(j, q_point);

                    curl_part +=
                      curl_phi_i[q_point] * curl_phi_j * fe_values.JxW(q_point);

                    mass_part +=
                      phi_i[q_point] * phi_j * fe_values.JxW(q_point);
                  }

                double mass_term  = curl_part - mass_part;
                cell_matrix(i, j) = mass_term;
                cell_matrix(j, i) = mass_term;
              }
          }

        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(
          cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
      }
  }


  template <int dim>
  void
  MaxwellProblem<dim>::solve()
  {
    SparseDirectUMFPACK A_direct;
    A_direct.initialize(system_matrix);
    A_direct.vmult(solution, system_rhs);
    constraints.distribute(solution);
  }



  template <int dim>
  void
  MaxwellProblem<dim>::print()
  {
    // choose the quadrature formulas
    QGauss<dim>        quadrature_formula(fe.degree + 2);
    const unsigned int n_q_points = quadrature_formula.size();

    // set update flags
    FEValues<dim> fe_values(fe,
                            quadrature_formula,
                            update_values | update_quadrature_points |
                              update_JxW_values);

    // Extractors to real part
    const FEValuesExtractors::Vector E_re(0);

    for (auto &cell : dof_handler.active_cell_iterators())
      {
        fe_values.reinit(cell);

        std::vector<Tensor<1, dim>> E_value_re(n_q_points);
        fe_values[E_re].get_function_values(solution, E_value_re);

        for (unsigned int component = 0; component < dim; ++component)
          {
            double value = 0.0;
            for (unsigned int q_point = 0; q_point < n_q_points; q_point++)
              value += E_value_re[q_point][component] * fe_values.JxW(q_point);

            deallog << value << "\t";
          }

        deallog << std::endl;
      }
  }



  template <int dim>
  void
  MaxwellProblem<dim>::output_results()
  {
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, {"real_Ex", "real_Ey", "real_Ez"});
    data_out.build_patches();
    std::ofstream output("solution.vtk");
    data_out.write_vtk(output);
  }



  template <int dim>
  void
  MaxwellProblem<dim>::run()
  {
    deallog << "Testing for dim = " << dim
            << ", polynomial_degree p = " << poly_degree << std::endl;

    make_grid();
    setup_system();
    assemble_system();
    solve();
    print();
    // output_results();
  }
} // namespace HangingEdge

int
main()
{
  initlog();

  using namespace HangingEdge;

  const unsigned int poly_degree = 2;
  const unsigned int n_cells     = 7;

  MaxwellProblem<3> maxwell_3d(poly_degree, n_cells);
  maxwell_3d.run();

  return 0;
}
