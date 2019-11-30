/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 by the deal.II authors
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

 *
 * Author:
 */
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/meshworker/mesh_loop.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>


#include <fstream>
#include <iostream>
#include <cmath>


namespace StepBiharmonic
{
  using namespace dealii;


  namespace ExactSolution
  {
    using numbers::PI;

    /**
     * An exact solution of the form
     * $ u(x,y) = \sin(\pi x) \sin(\pi y) $.
     *
     * Note that this solution has zero boundary values for the *value*
     * of the solution, but not for its Laplacian. Consequently, the
     * boundary contribution to the penalty terms is not zero.
     */
    template <int dim>
    class Solution : public Function<dim>
    {
    public:
      static_assert(dim == 2, "Only dim==2 is implemented");

      virtual double value(const Point<dim> &p,
                           const unsigned int /*component*/ = 0) const override
      {
        return std::sin(PI * p[0]) * std::sin(PI * p[1]);
      }

      virtual Tensor<1, dim>
      gradient(const Point<dim> &p,
               const unsigned int /*component*/ = 0) const override
      {
        Tensor<1, dim> r;
        r[0] = PI * std::cos(PI * p[0]) * std::sin(PI * p[1]);
        r[1] = PI * std::cos(PI * p[1]) * std::sin(PI * p[0]);
        return r;
      }

      virtual void
      hessian_list(const std::vector<Point<dim>> &       points,
                   std::vector<SymmetricTensor<2, dim>> &hessians,
                   const unsigned int /*component*/ = 0) const override
      {
        for (unsigned i = 0; i < points.size(); ++i)
          {
            const double x = points[i][0];
            const double y = points[i][1];

            hessians[i][0][0] = -PI * PI * std::sin(PI * x) * std::sin(PI * y);
            hessians[i][0][1] = PI * PI * std::cos(PI * x) * std::cos(PI * y);
            hessians[i][1][1] = -PI * PI * std::sin(PI * x) * std::sin(PI * y);
          }
      }
    };


    /**
     * The corresponding right hand side.
     */
    template <int dim>
    class RightHandSide : public Function<dim>
    {
    public:
      static_assert(dim == 2, "Only dim==2 is implemented");

      virtual double value(const Point<dim> &p,
                           const unsigned int /*component*/ = 0) const override

      {
        return 4 * std::pow(PI, 4.0) * std::sin(PI * p[0]) *
               std::sin(PI * p[1]);
      }
    };
  } // namespace ExactSolution



  /*************************************************************/
  // @sect3{The main class}
  template <int dim>
  class BiharmonicProblem
  {
  public:
    BiharmonicProblem(const unsigned int fe_degree);

    void run();

  private:
    void make_grid();
    void setup_system();
    void assemble_system();
    void solve();
    void compute_errors();
    void output_results(const unsigned int iteration) const;

    Triangulation<dim>        triangulation;
    const MappingQ<dim>       mapping;
    const FE_Q<dim>           fe;
    DoFHandler<dim>           dof_handler;
    AffineConstraints<double> constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double> solution;
    Vector<double> system_rhs;
  };

  template <int dim>
  BiharmonicProblem<dim>::BiharmonicProblem(const unsigned int fe_degree)
    : mapping(1)
    , fe(fe_degree)
    , dof_handler(triangulation)
  {}



  template <int dim>
  void BiharmonicProblem<dim>::make_grid()
  {
    GridGenerator::hyper_cube(triangulation, 0., 1.);
    triangulation.refine_global(1);

    std::cout << "Number of active cells: " << triangulation.n_active_cells()
              << std::endl
              << "Total number of cells: " << triangulation.n_cells()
              << std::endl;
  }



  template <int dim>
  void BiharmonicProblem<dim>::setup_system()
  {
    dof_handler.distribute_dofs(fe);

    std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs()
              << std::endl;

    constraints.clear();
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);

    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             ExactSolution::Solution<dim>(),
                                             constraints);
    constraints.close();


    DynamicSparsityPattern c_sparsity(dof_handler.n_dofs());
    DoFTools::make_flux_sparsity_pattern(dof_handler,
                                         c_sparsity,
                                         constraints,
                                         true);
    sparsity_pattern.copy_from(c_sparsity);
    system_matrix.reinit(sparsity_pattern);

    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
  }



  template <int dim>
  struct ScratchData
  {
    ScratchData(const Mapping<dim> &      mapping,
                const FiniteElement<dim> &fe,
                const unsigned int        quadrature_degree,
                const UpdateFlags         update_flags = update_values |
                                                 update_gradients |
                                                 update_quadrature_points |
                                                 update_JxW_values,
                const UpdateFlags interface_update_flags =
                  update_values | update_gradients | update_quadrature_points |
                  update_JxW_values | update_normal_vectors)
      : fe_values(mapping, fe, QGauss<dim>(quadrature_degree), update_flags)
      , fe_interface_values(mapping,
                            fe,
                            QGauss<dim - 1>(quadrature_degree),
                            interface_update_flags)
    {}


    ScratchData(const ScratchData<dim> &scratch_data)
      : fe_values(scratch_data.fe_values.get_mapping(),
                  scratch_data.fe_values.get_fe(),
                  scratch_data.fe_values.get_quadrature(),
                  scratch_data.fe_values.get_update_flags())
      , fe_interface_values(scratch_data.fe_values.get_mapping(),
                            scratch_data.fe_values.get_fe(),
                            scratch_data.fe_interface_values.get_quadrature(),
                            scratch_data.fe_interface_values.get_update_flags())
    {}

    FEValues<dim>          fe_values;
    FEInterfaceValues<dim> fe_interface_values;
  };



  struct CopyDataFace
  {
    FullMatrix<double>                   cell_matrix;
    std::vector<types::global_dof_index> joint_dof_indices;
  };



  struct CopyData
  {
    FullMatrix<double>                   cell_matrix;
    Vector<double>                       cell_rhs;
    std::vector<types::global_dof_index> local_dof_indices;
    std::vector<CopyDataFace>            face_data;

    template <class Iterator>
    void reinit(const Iterator &cell, unsigned int dofs_per_cell)
    {
      cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
      cell_rhs.reinit(dofs_per_cell);

      local_dof_indices.resize(dofs_per_cell);
      cell->get_dof_indices(local_dof_indices);
    }
  };



  template <int dim>
  void BiharmonicProblem<dim>::assemble_system()
  {
    using Iterator = decltype(dof_handler.begin_active());
    const ExactSolution::RightHandSide<dim> right_hand_side;

    auto cell_worker = [&](const Iterator &  cell,
                           ScratchData<dim> &scratch_data,
                           CopyData &        copy_data) {
      const unsigned int n_dofs = scratch_data.fe_values.get_fe().dofs_per_cell;
      copy_data.reinit(cell, n_dofs);
      scratch_data.fe_values.reinit(cell);

      const auto &q_points = scratch_data.fe_values.get_quadrature_points();

      const FEValues<dim> &      fe_v = scratch_data.fe_values;
      const std::vector<double> &JxW  = fe_v.get_JxW_values();

      const double nu = 1.0;

      for (unsigned int point = 0; point < fe_v.n_quadrature_points; ++point)
        {
          for (unsigned int i = 0; i < n_dofs; ++i)
            {
              for (unsigned int j = 0; j < n_dofs; ++j)
                {
                  // \int_Z \nu \nabla^2 u \cdot \nabla^2 v \, dx.
                  copy_data.cell_matrix(i, j) +=
                    nu *
                    scalar_product(fe_v.shape_hessian(i, point),
                                   fe_v.shape_hessian(j, point)) *
                    JxW[point]; // dx
                }

              copy_data.cell_rhs(i) += fe_v.shape_value(i, point) *
                                       right_hand_side.value(q_points[point]) *
                                       JxW[point]; // dx
            }
        }
    };


    auto face_worker = [&](const Iterator &    cell,
                           const unsigned int &f,
                           const unsigned int &sf,
                           const Iterator &    ncell,
                           const unsigned int &nf,
                           const unsigned int &nsf,
                           ScratchData<dim> &  scratch_data,
                           CopyData &          copy_data) {
      FEInterfaceValues<dim> &fe_i = scratch_data.fe_interface_values;
      fe_i.reinit(cell, f, sf, ncell, nf, nsf);
      const auto &q_points = fe_i.get_quadrature_points();

      copy_data.face_data.emplace_back();
      CopyDataFace &copy_data_face = copy_data.face_data.back();

      const unsigned int n_dofs        = fe_i.n_current_interface_dofs();
      copy_data_face.joint_dof_indices = fe_i.get_interface_dof_indices();

      copy_data_face.cell_matrix.reinit(n_dofs, n_dofs);

      const std::vector<double> &        JxW     = fe_i.get_JxW_values();
      const std::vector<Tensor<1, dim>> &normals = fe_i.get_normal_vectors();

      // eta = 1/2 + 2C_2
      // gamma = eta/|e|

      double gamma = 1.0; // TODO:

      {
        int                degree = fe.tensor_degree();
        const unsigned int normal1 =
          GeometryInfo<dim>::unit_normal_direction[f];
        const unsigned int normal2 =
          GeometryInfo<dim>::unit_normal_direction[nf];
        const unsigned int deg1sq =
          degree * (degree + 1); //(deg1 == 0) ? 1 : deg1 * (deg1+1);
        const unsigned int deg2sq =
          degree * (degree + 1); //(deg2 == 0) ? 1 : deg2 * (deg2+1);

        double penalty1 = deg1sq / cell->extent_in_direction(normal1);
        double penalty2 = deg2sq / ncell->extent_in_direction(normal2);
        if (cell->has_children() ^ ncell->has_children())
          {
            penalty1 *= 8;
          }
        gamma = 0.5 * (penalty1 + penalty2);
      }


      for (unsigned int qpoint = 0; qpoint < q_points.size(); ++qpoint)
        {
          // \int_F -{grad^2 u n n } [grad v n]
          //   - {grad^2 v n n } [grad u n]
          //   +  gamma [grad u n ][grad v n]
          const auto &n = normals[qpoint];

          for (unsigned int i = 0; i < n_dofs; ++i)
            for (unsigned int j = 0; j < n_dofs; ++j)
              {
                copy_data_face.cell_matrix(i, j) +=
                  (-(fe_i.average_hessian(i, qpoint) * n *
                     n)                                    // - {grad^2 v n n }
                     * (fe_i.jump_gradient(j, qpoint) * n) // [grad u n]
                   - (fe_i.average_hessian(j, qpoint) * n *
                      n) // - {grad^2 u n n }
                       * (fe_i.jump_gradient(i, qpoint) * n) // [grad v n]
                   // gamma [grad u n ][grad v n]:
                   + gamma * (fe_i.jump_gradient(i, qpoint) * n) *
                       (fe_i.jump_gradient(j, qpoint) * n)) *
                  JxW[qpoint]; // dx
              }
        }
    };


    auto boundary_worker = [&](const Iterator &    cell,
                               const unsigned int &face_no,
                               ScratchData<dim> &  scratch_data,
                               CopyData &          copy_data) {
      // return;
      FEInterfaceValues<dim> &fe_i = scratch_data.fe_interface_values;
      fe_i.reinit(cell, face_no);
      const auto &q_points = fe_i.get_quadrature_points();

      copy_data.face_data.emplace_back();
      CopyDataFace &copy_data_face = copy_data.face_data.back();

      const unsigned int n_dofs        = fe_i.n_current_interface_dofs();
      copy_data_face.joint_dof_indices = fe_i.get_interface_dof_indices();

      copy_data_face.cell_matrix.reinit(n_dofs, n_dofs);

      const std::vector<double> &        JxW     = fe_i.get_JxW_values();
      const std::vector<Tensor<1, dim>> &normals = fe_i.get_normal_vectors();


      const ExactSolution::Solution<dim> exact_solution;
      std::vector<Tensor<1, dim>>        exact_gradients(q_points.size());
      exact_solution.gradient_list(q_points, exact_gradients);


      // eta = 1/2 + 2C_2
      // gamma = eta/|e|

      double gamma = 1.0;

      {
        int                degree = fe.tensor_degree();
        const unsigned int normal1 =
          GeometryInfo<dim>::unit_normal_direction[face_no];
        const unsigned int deg1sq =
          degree * (degree + 1); //(deg1 == 0) ? 1 : deg1 * (deg1+1);

        gamma = deg1sq / cell->extent_in_direction(normal1);
        //      gamma = 0.5*(penalty1 + penalty2);
      }

      for (unsigned int qpoint = 0; qpoint < q_points.size(); ++qpoint)
        {
          const auto &n = normals[qpoint];

          for (unsigned int i = 0; i < n_dofs; ++i)
            {
              for (unsigned int j = 0; j < n_dofs; ++j)
                copy_data_face.cell_matrix(i, j) +=
                  (-(fe_i.average_hessian(i, qpoint) * n *
                     n)                                    // - {grad^2 v n n }
                     * (fe_i.jump_gradient(j, qpoint) * n) // [grad u n]
                   //
                   - (fe_i.average_hessian(j, qpoint) * n *
                      n) // - {grad^2 u n n }
                       * (fe_i.jump_gradient(i, qpoint) * n) //  [grad v n]
                                                             //
                   + 2.0 * gamma *
                       (fe_i.jump_gradient(i, qpoint) * n) // 2 gamma [grad v n]
                       * (fe_i.jump_gradient(j, qpoint) * n) // [grad u n]
                   ) *
                  JxW[qpoint]; // dx

              copy_data.cell_rhs(i) +=
                (-(fe_i.average_hessian(i, qpoint) * n *
                   n) *                                    //  - {grad^2 v n n }
                   (exact_gradients[qpoint] * n)           // (grad u_exact n)
                 + 2.0 * gamma                             //
                     * (fe_i.jump_gradient(i, qpoint) * n) // [grad v n]
                     * (exact_gradients[qpoint] * n)       // (grad u_exact n)
                 ) *
                JxW[qpoint]; // dx
            }
        }
    };

    auto copier = [&](const CopyData &c) {
      constraints.distribute_local_to_global(c.cell_matrix,
                                             c.cell_rhs,
                                             c.local_dof_indices,
                                             system_matrix,
                                             system_rhs);

      for (auto &cdf : c.face_data)
        {
          constraints.distribute_local_to_global(cdf.cell_matrix,
                                                 cdf.joint_dof_indices,
                                                 system_matrix);
        }
    };

    const unsigned int n_gauss_points = dof_handler.get_fe().degree + 1;

    ScratchData<dim> scratch_data(mapping,
                                  fe,
                                  n_gauss_points,
                                  update_values | update_gradients |
                                    update_hessians | update_quadrature_points |
                                    update_JxW_values,
                                  update_values | update_gradients |
                                    update_hessians | update_quadrature_points |
                                    update_JxW_values | update_normal_vectors);
    CopyData         copy_data;
    MeshWorker::mesh_loop(dof_handler.begin_active(),
                          dof_handler.end(),
                          cell_worker,
                          copier,
                          scratch_data,
                          copy_data,
                          MeshWorker::assemble_own_cells |
                            MeshWorker::assemble_boundary_faces |
                            MeshWorker::assemble_own_interior_faces_once,
                          boundary_worker,
                          face_worker);
  }


  template <int dim>
  void BiharmonicProblem<dim>::solve()
  {
    std::cout << "   Solving system..." << std::endl;

    SparseDirectUMFPACK A_direct;
    A_direct.initialize(system_matrix);
    A_direct.vmult(solution, system_rhs);
    constraints.distribute(solution);
  }



  template <int dim>
  void BiharmonicProblem<dim>::compute_errors()
  {
    const unsigned int n_gauss_points =
      dof_handler.get_fe().tensor_degree() + 1;

    {
      Vector<float> norm_per_cell(triangulation.n_active_cells());
      VectorTools::integrate_difference(mapping,
                                        dof_handler,
                                        solution,
                                        ExactSolution::Solution<dim>(),
                                        norm_per_cell,
                                        QGauss<dim>(n_gauss_points + 1),
                                        VectorTools::L2_norm);
      const double error_norm =
        VectorTools::compute_global_error(triangulation,
                                          norm_per_cell,
                                          VectorTools::L2_norm);
      std::cout << "   Error in the L2 norm       :     " << error_norm
                << std::endl;
    }

    {
      Vector<float> norm_per_cell(triangulation.n_active_cells());
      VectorTools::integrate_difference(mapping,
                                        dof_handler,
                                        solution,
                                        ExactSolution::Solution<dim>(),
                                        norm_per_cell,
                                        QGauss<dim>(n_gauss_points + 1),
                                        VectorTools::H1_seminorm);
      const double error_norm =
        VectorTools::compute_global_error(triangulation,
                                          norm_per_cell,
                                          VectorTools::H1_seminorm);
      std::cout << "   Error in the H1 seminorm       : " << error_norm
                << std::endl;
    }

    // Now also compute the H2 seminorm error, integrating over the interiors
    // of the cells but not taking into account the interface jump terms.
    // This is *not* equivalent to the energy error for the problem.
    {
      const QGauss<dim>            quadrature_formula(fe.degree + 2);
      ExactSolution::Solution<dim> exact_solution;
      Vector<double> error_per_cell(triangulation.n_active_cells());

      FEValues<dim> fe_values(mapping,
                              fe,
                              quadrature_formula,
                              update_values | update_hessians |
                                update_quadrature_points | update_JxW_values);

      FEValuesExtractors::Scalar scalar(0);
      const unsigned int         n_q_points = quadrature_formula.size();

      std::vector<SymmetricTensor<2, dim>> exact_hessians(n_q_points);
      std::vector<Tensor<2, dim>>          hessians(n_q_points);
      for (auto cell : dof_handler.active_cell_iterators())
        {
          fe_values.reinit(cell);
          fe_values[scalar].get_function_hessians(solution, hessians);
          exact_solution.hessian_list(fe_values.get_quadrature_points(),
                                      exact_hessians);

          double diff = 0;
          for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
            {
              diff +=
                ((exact_hessians[q_point] - hessians[q_point]).norm_square() *
                 fe_values.JxW(q_point));
            }
          error_per_cell[cell->active_cell_index()] = std::sqrt(diff);
        }
      const double error_norm = error_per_cell.l2_norm();
      std::cout << "   Error in the broken H2 seminorm: " << error_norm
                << std::endl;
    }
  }


  template <int dim>
  void
  BiharmonicProblem<dim>::output_results(const unsigned int iteration) const
  {
    std::cout << "   Writing graphical output..." << std::endl;

    DataOut<dim> data_out;

    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "u");
    Vector<double>                     exact  = solution;
    unsigned int                       degree = fe.tensor_degree();
    const ExactSolution::Solution<dim> exact_solution;
    VectorTools::project(mapping,
                         dof_handler,
                         constraints,
                         QGauss<dim>(degree + 1),
                         exact_solution,
                         exact);
    data_out.add_data_vector(exact, "exact");

    data_out.build_patches();

    std::ofstream output_vtk(
      ("output_" + Utilities::int_to_string(iteration, 6) + ".vtk").c_str());
    data_out.write_vtk(output_vtk);
  }



  template <int dim>
  void BiharmonicProblem<dim>::run()
  {
    make_grid();

    const unsigned int n_cycles = 4;
    for (unsigned int cycle = 0; cycle < n_cycles; ++cycle)
      {
        std::cout << "Cycle: " << cycle << " of " << n_cycles << std::endl;



        triangulation.refine_global(1);
        setup_system();

        assemble_system();
        solve();

        output_results(cycle);

        compute_errors();
        std::cout << std::endl;
      }
  }
} // namespace StepBiharmonic



int main(int argc, char *argv[])
{
  try
    {
      using namespace dealii;
      using namespace StepBiharmonic;

      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);

      unsigned int degree = 2; // minimum degree 2

      // If provided on the command line, override the polynomial degree
      // by the one given there.
      if (argc > 1)
        degree = Utilities::string_to_int(argv[1]);

      BiharmonicProblem<2> my_bi(degree);
      my_bi.run();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
