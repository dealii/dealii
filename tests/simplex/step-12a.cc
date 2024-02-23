// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Step-12 with tetrahedron mesh.

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/precondition_block.h>
#include <deal.II/lac/solver_richardson.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/derivative_approximation.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

// Finally, the new include file for using the mesh_loop from the MeshWorker
// framework
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_simplex_p_bubbles.h>
#include <deal.II/fe/fe_wedge_p.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/precondition.h>

#include <deal.II/meshworker/mesh_loop.h>

#include <fstream>
#include <iostream>

// #define HEX

namespace Step12
{

  template <int dim>
  class BoundaryValues : public Function<dim>
  {
  public:
    BoundaryValues() = default;
    virtual void
    value_list(const std::vector<Point<dim>> &points,
               std::vector<double>           &values,
               const unsigned int             component = 0) const override;
  };

  template <int dim>
  void
  BoundaryValues<dim>::value_list(const std::vector<Point<dim>> &points,
                                  std::vector<double>           &values,
                                  const unsigned int component) const
  {
    (void)component;
    AssertIndexRange(component, 1);
    Assert(values.size() == points.size(),
           ExcDimensionMismatch(values.size(), points.size()));

    for (unsigned int i = 0; i < values.size(); ++i)
      {
        if (points[i][0] < 0.5)
          values[i] = 1.;
        else
          values[i] = 0.;
      }
  }


  template <int dim>
  Tensor<1, dim>
  beta(const Point<dim> &p)
  {
    Assert(dim >= 2, ExcNotImplemented());

    Point<dim> wind_field;
    wind_field[0] = -p[1];
    wind_field[1] = p[0];

    if (wind_field.norm() > 1e-6)
      wind_field /= wind_field.norm();

    return wind_field;
  }


  template <int dim>
  struct ScratchData
  {
    ScratchData(const Mapping<dim>        &mapping,
                const FiniteElement<dim>  &fe,
                const Quadrature<dim>     &quad,
                const Quadrature<dim - 1> &quad_face,
                const UpdateFlags          update_flags = update_values |
                                                 update_gradients |
                                                 update_quadrature_points |
                                                 update_JxW_values,
                const UpdateFlags interface_update_flags =
                  update_values | update_gradients | update_quadrature_points |
                  update_JxW_values | update_normal_vectors)
      : fe_values(mapping, fe, quad, update_flags)
      , fe_interface_values(mapping, fe, quad_face, interface_update_flags)
    {}


    ScratchData(const ScratchData<dim> &scratch_data)
      : fe_values(scratch_data.fe_values.get_mapping(),
                  scratch_data.fe_values.get_fe(),
                  scratch_data.fe_values.get_quadrature(),
                  scratch_data.fe_values.get_update_flags())
      , fe_interface_values(
          scratch_data.fe_values
            .get_mapping(), // TODO: implement for fe_interface_values
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
    void
    reinit(const Iterator &cell, unsigned int dofs_per_cell)
    {
      cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
      cell_rhs.reinit(dofs_per_cell);

      local_dof_indices.resize(dofs_per_cell);
      cell->get_dof_indices(local_dof_indices);
    }
  };


  template <int dim>
  class AdvectionProblem
  {
  public:
    AdvectionProblem();
    void
    run();

  private:
    void
    setup_system();
    void
    assemble_system();
    void
    solve();
    void
    refine_grid();
    void
    output_results(const unsigned int cycle) const;

    Triangulation<dim> triangulation;
#ifdef HEX
    const MappingQ1<dim> mapping;
#else
    FE_SimplexP<dim>     fe_mapping;
    const MappingFE<dim> mapping;
#endif

    // Furthermore we want to use DG elements.
#ifdef HEX
    FE_DGQ<dim> fe;
#else
    FE_SimplexDGP<dim>   fe;
#endif
    DoFHandler<dim> dof_handler;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double> solution;
    Vector<double> right_hand_side;
  };


  template <int dim>
  AdvectionProblem<dim>::AdvectionProblem()
#ifdef HEX
    : mapping()
    , fe(2)
#else
    : fe_mapping(1)
    , mapping(fe_mapping)
    , fe(2)
#endif
    , dof_handler(triangulation)
  {}


  template <int dim>
  void
  AdvectionProblem<dim>::setup_system()
  {
    dof_handler.distribute_dofs(fe);

    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_flux_sparsity_pattern(dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);

    system_matrix.reinit(sparsity_pattern);
    solution.reinit(dof_handler.n_dofs());
    right_hand_side.reinit(dof_handler.n_dofs());
  }

  template <int dim>
  void
  AdvectionProblem<dim>::assemble_system()
  {
    using Iterator = typename DoFHandler<dim>::active_cell_iterator;
    const BoundaryValues<dim> boundary_function;

    auto cell_worker = [&](const Iterator   &cell,
                           ScratchData<dim> &scratch_data,
                           CopyData         &copy_data) {
      const unsigned int n_dofs = scratch_data.fe_values.get_fe().dofs_per_cell;
      copy_data.reinit(cell, n_dofs);
      scratch_data.fe_values.reinit(cell);

      const auto &q_points = scratch_data.fe_values.get_quadrature_points();

      const FEValues<dim>       &fe_v = scratch_data.fe_values;
      const std::vector<double> &JxW  = fe_v.get_JxW_values();

      for (unsigned int point = 0; point < fe_v.n_quadrature_points; ++point)
        {
          auto beta_q = beta(q_points[point]);
          for (unsigned int i = 0; i < n_dofs; ++i)
            for (unsigned int j = 0; j < n_dofs; ++j)
              {
                copy_data.cell_matrix(i, j) +=
                  -beta_q                      // -\beta
                  * fe_v.shape_grad(i, point)  // \nabla \phi_i
                  * fe_v.shape_value(j, point) // \phi_j
                  * JxW[point];                // dx
              }
        }
    };

    auto boundary_worker = [&](const Iterator     &cell,
                               const unsigned int &face_no,
                               ScratchData<dim>   &scratch_data,
                               CopyData           &copy_data) {
      scratch_data.fe_interface_values.reinit(cell, face_no);
      const FEFaceValuesBase<dim> &fe_face =
        scratch_data.fe_interface_values.get_fe_face_values(0);

      const auto &q_points = fe_face.get_quadrature_points();

      const unsigned int n_facet_dofs = fe_face.get_fe().n_dofs_per_cell();
      const std::vector<double>         &JxW     = fe_face.get_JxW_values();
      const std::vector<Tensor<1, dim>> &normals = fe_face.get_normal_vectors();

      std::vector<double> g(q_points.size());
      boundary_function.value_list(q_points, g);

      for (unsigned int point = 0; point < q_points.size(); ++point)
        {
          const double beta_dot_n = beta(q_points[point]) * normals[point];

          if (beta_dot_n > 0)
            {
              for (unsigned int i = 0; i < n_facet_dofs; ++i)
                for (unsigned int j = 0; j < n_facet_dofs; ++j)
                  copy_data.cell_matrix(i, j) +=
                    fe_face.shape_value(i, point)   // \phi_i
                    * fe_face.shape_value(j, point) // \phi_j
                    * beta_dot_n                    // \beta . n
                    * JxW[point];                   // dx
            }
          else
            for (unsigned int i = 0; i < n_facet_dofs; ++i)
              copy_data.cell_rhs(i) += -fe_face.shape_value(i, point) // \phi_i
                                       * g[point]                     // g
                                       * beta_dot_n  // \beta . n
                                       * JxW[point]; // dx
        }
    };

    auto face_worker = [&](const Iterator     &cell,
                           const unsigned int &f,
                           const unsigned int &sf,
                           const Iterator     &ncell,
                           const unsigned int &nf,
                           const unsigned int &nsf,
                           ScratchData<dim>   &scratch_data,
                           CopyData           &copy_data) {
      FEInterfaceValues<dim> &fe_iv = scratch_data.fe_interface_values;
      fe_iv.reinit(cell, f, sf, ncell, nf, nsf);
      const auto &q_points = fe_iv.get_quadrature_points();

      copy_data.face_data.emplace_back();
      CopyDataFace &copy_data_face = copy_data.face_data.back();

      const unsigned int n_dofs        = fe_iv.n_current_interface_dofs();
      copy_data_face.joint_dof_indices = fe_iv.get_interface_dof_indices();

      copy_data_face.cell_matrix.reinit(n_dofs, n_dofs);

      const std::vector<double>         &JxW     = fe_iv.get_JxW_values();
      const std::vector<Tensor<1, dim>> &normals = fe_iv.get_normal_vectors();

      for (unsigned int qpoint = 0; qpoint < q_points.size(); ++qpoint)
        {
          const double beta_dot_n = beta(q_points[qpoint]) * normals[qpoint];
          for (unsigned int i = 0; i < n_dofs; ++i)
            for (unsigned int j = 0; j < n_dofs; ++j)
              copy_data_face.cell_matrix(i, j) +=
                fe_iv.jump_in_shape_values(i, qpoint) // [\phi_i]
                *
                fe_iv.shape_value((beta_dot_n > 0), j, qpoint) // phi_j^{upwind}
                * beta_dot_n                                   // (\beta . n)
                * JxW[qpoint];                                 // dx
        }
    };

    AffineConstraints<double> constraints;

    auto copier = [&](const CopyData &c) {
      constraints.distribute_local_to_global(c.cell_matrix,
                                             c.cell_rhs,
                                             c.local_dof_indices,
                                             system_matrix,
                                             right_hand_side);

      for (auto &cdf : c.face_data)
        {
          constraints.distribute_local_to_global(cdf.cell_matrix,
                                                 cdf.joint_dof_indices,
                                                 system_matrix);
        }
    };

    const unsigned int degree = dof_handler.get_fe().degree;

#ifdef HEX
    QGauss<dim> quad(degree + 1);

    QGauss<dim - 1> face_quad(degree + 1);
#else
    QGaussSimplex<dim> quad(degree + 1);

    QGaussSimplex<dim - 1> face_quad(degree + 1);
#endif

    ScratchData<dim> scratch_data(mapping, fe, quad, face_quad);
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
  void
  AdvectionProblem<dim>::solve()
  {
    SolverControl                    solver_control(1000, 1e-12);
    SolverRichardson<Vector<double>> solver(solver_control);

    PreconditionBlockSSOR<SparseMatrix<double>> preconditioner;

    preconditioner.initialize(system_matrix, fe.dofs_per_cell);

    solver.solve(system_matrix, solution, right_hand_side, preconditioner);

    deallog << "  Solver converged in " << solver_control.last_step()
            << " iterations." << std::endl;
  }


  template <int dim>
  void
  AdvectionProblem<dim>::refine_grid()
  {
    Vector<float> gradient_indicator(triangulation.n_active_cells());

    DerivativeApproximation::approximate_gradient(mapping,
                                                  dof_handler,
                                                  solution,
                                                  gradient_indicator);

    unsigned int cell_no = 0;
    for (const auto &cell : dof_handler.active_cell_iterators())
      gradient_indicator(cell_no++) *=
        std::pow(cell->diameter(), 1 + 1.0 * dim / 2);

    GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                                    gradient_indicator,
                                                    0.3,
                                                    0.1);

    triangulation.execute_coarsening_and_refinement();
  }


  template <int dim>
  void
  AdvectionProblem<dim>::output_results(const unsigned int cycle) const
  {
#if false
#  ifdef HEX
    const std::string filename =
      dim == 2 ? ("step12-quad-" + std::to_string(cycle) + ".vtk") :
                 ("step12-hex-" + std::to_string(cycle) + ".vtk");
#  else
    const std::string filename =
      dim == 2 ? ("step12-tri-" + std::to_string(cycle) + ".vtk") :
                 ("step12-tet-" + std::to_string(cycle) + ".vtk");
#  endif
    deallog << "  Writing solution to <" << filename << '>' << std::endl;
    std::ofstream output(filename);

    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "u", DataOut<dim>::type_dof_data);

    data_out.build_patches(mapping, 2);

    data_out.write_vtk(output);
#endif

    {
      Vector<float> values(triangulation.n_active_cells());
      VectorTools::integrate_difference(mapping,
                                        dof_handler,
                                        solution,
                                        Functions::ZeroFunction<dim>(),
                                        values,
                                        QGauss<dim>(fe.degree + 1),
                                        VectorTools::Linfty_norm);
      const double l_infty =
        VectorTools::compute_global_error(triangulation,
                                          values,
                                          VectorTools::Linfty_norm);
      deallog << "  L-infinity norm: " << l_infty << std::endl;
    }
  }


  template <int dim>
  void
  AdvectionProblem<dim>::run()
  {
    // #ifdef HEX
    //     for (unsigned int cycle = 0; cycle < 6; ++cycle)
    // #else
    for (unsigned int cycle = 0; cycle < 1; ++cycle)
      // #endif
      {
        deallog << "Cycle " << cycle << std::endl;

        if (cycle == 0)
          {
#ifdef HEX
            // GridGenerator::hyper_cube(triangulation);
            // triangulation.refine_global(3);
            GridGenerator::subdivided_hyper_cube(triangulation, 16);
#else
            GridGenerator::subdivided_hyper_cube_with_simplices(triangulation,
                                                                dim == 2 ? 32 :
                                                                           8);
#endif
          }
        else
          refine_grid();

        deallog << "  Number of active cells:       "
                << triangulation.n_active_cells() << std::endl;

        setup_system();

        deallog << "  Number of degrees of freedom: " << dof_handler.n_dofs()
                << std::endl;

        assemble_system();
        solve();

        output_results(cycle);
      }
  }
} // namespace Step12


int
main()
{
  initlog();

  {
    Step12::AdvectionProblem<2> dgmethod;
    dgmethod.run();
  }
  {
    Step12::AdvectionProblem<3> dgmethod;
    dgmethod.run();
  }

  return 0;
}
