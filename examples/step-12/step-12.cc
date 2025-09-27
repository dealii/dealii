/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2009 - 2024 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 *
 * Authors: Guido Kanschat, Texas A&M University, 2009
 *          Timo Heister, Clemson University, 2019
 */


// The first few files have already been covered in previous examples and will
// thus not be further commented on:
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/fe/mapping_q1.h>
// Here the discontinuous finite elements are defined. They are used in the same
// way as all other finite elements, though -- as you have seen in previous
// tutorial programs -- there isn't much user interaction with finite element
// classes at all: they are passed to <code>DoFHandler</code> and
// <code>FEValues</code> objects, and that is about it.
#include <deal.II/fe/fe_dgq.h>
// This header is needed for FEInterfaceValues to compute integrals on
// interfaces:
#include <deal.II/fe/fe_interface_values.h>
// We are going to use a standard solver, called Generalized minimal residual
// method (GMRES). It is an iterative solver which is applicable to arbitrary
// invertible matrices. This, in combination with a block SSOR preconditioner
// (defined in precondition_block.h), that uses the special block matrix
// structure of system matrices arising from DG discretizations.
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition_block.h>
// We are going to use gradients as refinement indicator.
#include <deal.II/numerics/derivative_approximation.h>

// Finally, the new include file for using the mesh_loop from the MeshWorker
// framework
#include <deal.II/meshworker/mesh_loop.h>

// Like in all programs, we finish this section by including the needed C++
// headers and declaring we want to use objects in the dealii namespace without
// prefix.
#include <iostream>
#include <fstream>


namespace Step12
{
  using namespace dealii;

  // @sect3{Equation data}
  //
  // First, we define a class describing the inhomogeneous boundary data. Since
  // only its values are used, we implement value_list(), but leave all other
  // functions of Function undefined.
  template <int dim>
  class BoundaryValues : public Function<dim>
  {
  public:
    BoundaryValues() = default;
    virtual void value_list(const std::vector<Point<dim>> &points,
                            std::vector<double>           &values,
                            const unsigned int component = 0) const override;
  };

  // Given the flow direction, the inflow boundary of the unit square $[0,1]^2$
  // are the right and the lower boundaries. We prescribe discontinuous boundary
  // values 1 and 0 on the x-axis and value 0 on the right boundary. The values
  // of this function on the outflow boundaries will not be used within the DG
  // scheme.
  template <int dim>
  void BoundaryValues<dim>::value_list(const std::vector<Point<dim>> &points,
                                       std::vector<double>           &values,
                                       const unsigned int component) const
  {
    (void)component;
    AssertIndexRange(component, 1);
    AssertDimension(values.size(), points.size());

    for (unsigned int i = 0; i < values.size(); ++i)
      {
        if (points[i][0] < 0.5)
          values[i] = 1.;
        else
          values[i] = 0.;
      }
  }


  // Finally, a function that computes and returns the wind field
  // $\beta=\beta(\mathbf x)$. As explained in the introduction, we will use a
  // rotational field around the origin in 2d. In 3d, we simply leave the
  // $z$-component unset (i.e., at zero), whereas the function can not be used
  // in 1d in its current implementation:
  template <int dim>
  Tensor<1, dim> beta(const Point<dim> &p)
  {
    Assert(dim >= 2, ExcNotImplemented());

    Tensor<1, dim> wind_field;
    wind_field[0] = -p[1];
    wind_field[1] = p[0];

    if (wind_field.norm() > 1e-10)
      wind_field /= wind_field.norm();

    return wind_field;
  }


  // @sect3{The ScratchData and CopyData classes}
  //
  // The following objects are the scratch and copy objects we use in the call
  // to MeshWorker::mesh_loop(). The new object is the FEInterfaceValues object,
  // that works similar to FEValues or FEFaceValues, except that it acts on
  // an interface between two cells and allows us to assemble the interface
  // terms in our weak form.

  template <int dim>
  struct ScratchData
  {
    ScratchData(const Mapping<dim>        &mapping,
                const FiniteElement<dim>  &fe,
                const Quadrature<dim>     &quadrature,
                const Quadrature<dim - 1> &quadrature_face,
                const UpdateFlags          update_flags = update_values |
                                                 update_gradients |
                                                 update_quadrature_points |
                                                 update_JxW_values,
                const UpdateFlags interface_update_flags =
                  update_values | update_gradients | update_quadrature_points |
                  update_JxW_values | update_normal_vectors)
      : fe_values(mapping, fe, quadrature, update_flags)
      , fe_interface_values(mapping,
                            fe,
                            quadrature_face,
                            interface_update_flags)
    {}


    ScratchData(const ScratchData<dim> &scratch_data)
      : fe_values(scratch_data.fe_values.get_mapping(),
                  scratch_data.fe_values.get_fe(),
                  scratch_data.fe_values.get_quadrature(),
                  scratch_data.fe_values.get_update_flags())
      , fe_interface_values(scratch_data.fe_interface_values.get_mapping(),
                            scratch_data.fe_interface_values.get_fe(),
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


  // @sect3{The AdvectionProblem class}
  //
  // After this preparations, we proceed with the main class of this program,
  // called AdvectionProblem.
  //
  // This should all be pretty familiar to you. Interesting details will only
  // come up in the implementation of the assemble function.
  template <int dim>
  class AdvectionProblem
  {
  public:
    AdvectionProblem();
    void run();

  private:
    void setup_system();
    void assemble_system();
    void solve();
    void refine_grid();
    void output_results(const unsigned int cycle) const;

    Triangulation<dim>   triangulation;
    const MappingQ1<dim> mapping;

    // Furthermore we want to use DG elements.
    const FE_DGQ<dim> fe;
    DoFHandler<dim>   dof_handler;

    const QGauss<dim>     quadrature;
    const QGauss<dim - 1> quadrature_face;

    // The next four members represent the linear system to be solved.
    // <code>system_matrix</code> and <code>right_hand_side</code> are generated
    // by <code>assemble_system()</code>, the <code>solution</code> is computed
    // in <code>solve()</code>. The <code>sparsity_pattern</code> is used to
    // determine the location of nonzero elements in <code>system_matrix</code>.
    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double> solution;
    Vector<double> right_hand_side;
  };


  // We start with the constructor. The 1 in the constructor call of
  // <code>fe</code> is the polynomial degree.
  template <int dim>
  AdvectionProblem<dim>::AdvectionProblem()
    : mapping()
    , fe(1)
    , dof_handler(triangulation)
    , quadrature(fe.tensor_degree() + 1)
    , quadrature_face(fe.tensor_degree() + 1)
  {}


  template <int dim>
  void AdvectionProblem<dim>::setup_system()
  {
    // In the function that sets up the usual finite element data structures, we
    // first need to distribute the DoFs.
    dof_handler.distribute_dofs(fe);

    // We start by generating the sparsity pattern. To this end, we first fill
    // an intermediate object of type DynamicSparsityPattern with the couplings
    // appearing in the system. After building the pattern, this object is
    // copied to <code>sparsity_pattern</code> and can be discarded.

    // To build the sparsity pattern for DG discretizations, we can call the
    // function analogue to DoFTools::make_sparsity_pattern, which is called
    // DoFTools::make_flux_sparsity_pattern:
    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_flux_sparsity_pattern(dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);

    // Finally, we set up the structure of all components of the linear system.
    system_matrix.reinit(sparsity_pattern);
    solution.reinit(dof_handler.n_dofs());
    right_hand_side.reinit(dof_handler.n_dofs());
  }

  // @sect4{The assemble_system function}

  // Here we see the major difference to assembling by hand. Instead of
  // writing loops over cells and faces, the logic is contained in the call to
  // MeshWorker::mesh_loop() and we only need to specify what should happen on
  // each cell, each boundary face, and each interior face. These three tasks
  // are handled by the lambda functions inside the function below.

  template <int dim>
  void AdvectionProblem<dim>::assemble_system()
  {
    using Iterator = typename DoFHandler<dim>::active_cell_iterator;
    const BoundaryValues<dim> boundary_function;

    // This is the function that will be executed for each cell.
    const auto cell_worker = [&](const Iterator   &cell,
                                 ScratchData<dim> &scratch_data,
                                 CopyData         &copy_data) {
      const unsigned int n_dofs =
        scratch_data.fe_values.get_fe().n_dofs_per_cell();
      copy_data.reinit(cell, n_dofs);
      scratch_data.fe_values.reinit(cell);

      const auto &q_points = scratch_data.fe_values.get_quadrature_points();

      const FEValues<dim>       &fe_v = scratch_data.fe_values;
      const std::vector<double> &JxW  = fe_v.get_JxW_values();

      // We solve a homogeneous equation, thus no right hand side shows up in
      // the cell term.  What's left is integrating the matrix entries.
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

    // This is the function called for boundary faces and consists of a normal
    // integration using FEFaceValues. New is the logic to decide if the term
    // goes into the system matrix (outflow) or the right-hand side (inflow).
    const auto boundary_worker = [&](const Iterator     &cell,
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

    // This is the function called on interior faces. The arguments specify
    // cells, face and subface indices (for adaptive refinement). We just pass
    // them along to the reinit() function of FEInterfaceValues.
    const auto face_worker = [&](const Iterator     &cell,
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

    // The following lambda function will handle copying the data from the
    // cell and face assembly into the global matrix and right-hand side.
    //
    // While we would not need an AffineConstraints object, because there are
    // no hanging node constraints in DG discretizations, we use an empty
    // object here as this allows us to use its `copy_local_to_global`
    // functionality.
    const AffineConstraints<double> constraints;

    const auto copier = [&](const CopyData &c) {
      constraints.distribute_local_to_global(c.cell_matrix,
                                             c.cell_rhs,
                                             c.local_dof_indices,
                                             system_matrix,
                                             right_hand_side);

      for (const auto &cdf : c.face_data)
        {
          constraints.distribute_local_to_global(cdf.cell_matrix,
                                                 cdf.joint_dof_indices,
                                                 system_matrix);
        }
    };

    ScratchData<dim> scratch_data(mapping, fe, quadrature, quadrature_face);
    CopyData         copy_data;

    // Here, we finally handle the assembly. We pass in ScratchData and
    // CopyData objects, the lambda functions from above, an specify that we
    // want to assemble interior faces once.
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

  // @sect3{All the rest}
  //
  // For this simple problem we use a standard iterative solver, called GMRES,
  // that creates approximate solutions minimizing the residual in each
  // iterations by adding a new basis vector to the Krylov subspace. This, in
  // combination with a block SSOR preconditioner, that uses the special block
  // matrix structure of system matrices arising from DG discretizations. The
  // size of these blocks are the number of DoFs per cell. Here, we use a SSOR
  // preconditioning as we have not renumbered the DoFs according to the flow
  // field. If the DoFs are renumbered in the downstream direction of the flow,
  // then a block Gauss-Seidel preconditioner (see the PreconditionBlockSOR
  // class with relaxation=1) does a much better job.

  // We create an additional data object for the GMRES solver to increase the
  // maximum number of basis vectors of the Krylov subspace. When this number
  // is reached the GMRES algorithm is restarted using the solution of the
  // previous iteration as the starting approximation. The choice of the
  // number of basis vectors is a trade-off between memory consumption and
  // convergence speed, since a longer basis means minimization over a larger
  // space.
  template <int dim>
  void AdvectionProblem<dim>::solve()
  {
    SolverControl solver_control(1000, 1e-6 * right_hand_side.l2_norm());

    SolverGMRES<Vector<double>>::AdditionalData additional_data;
    additional_data.max_basis_size = 100;
    SolverGMRES<Vector<double>> solver(solver_control, additional_data);

    // Here we create the preconditioner,
    PreconditionBlockSSOR<SparseMatrix<double>> preconditioner;

    // then assign the matrix to it and set the right block size:
    preconditioner.initialize(system_matrix, fe.n_dofs_per_cell());

    // After these preparations we are ready to start the linear solver.
    solver.solve(system_matrix, solution, right_hand_side, preconditioner);

    std::cout << "  Solver converged in " << solver_control.last_step()
              << " iterations." << std::endl;
  }


  // We refine the grid according to a very simple refinement criterion, namely
  // an approximation to the gradient of the solution. As here we consider the
  // DG(1) method (i.e. we use piecewise bilinear shape functions) we could
  // simply compute the gradients on each cell. But we do not want to base our
  // refinement indicator on the gradients on each cell only, but want to base
  // them also on jumps of the discontinuous solution function over faces
  // between neighboring cells. The simplest way of doing that is to compute
  // approximative gradients by difference quotients including the cell under
  // consideration and its neighbors. This is done by the
  // <code>DerivativeApproximation</code> class that computes the approximate
  // gradients in a way similar to the <code>GradientEstimation</code> described
  // in step-9 of this tutorial. In fact, the
  // <code>DerivativeApproximation</code> class was developed following the
  // <code>GradientEstimation</code> class of step-9. Relating to the discussion
  // in step-9, here we consider $h^{1+d/2}|\nabla_h u_h|$. Furthermore we note
  // that we do not consider approximate second derivatives because solutions to
  // the linear advection equation are in general not in $H^2$ but only in $H^1$
  // (or, to be more precise: in $H^1_\beta$, i.e., the space of functions whose
  // derivatives in direction $\beta$ are square integrable).
  template <int dim>
  void AdvectionProblem<dim>::refine_grid()
  {
    // The <code>DerivativeApproximation</code> class computes the gradients to
    // float precision. This is sufficient as they are approximate and serve as
    // refinement indicators only.
    Vector<float> gradient_indicator(triangulation.n_active_cells());

    // Now the approximate gradients are computed
    DerivativeApproximation::approximate_gradient(mapping,
                                                  dof_handler,
                                                  solution,
                                                  gradient_indicator);

    // and they are cell-wise scaled by the factor $h^{1+d/2}$
    unsigned int cell_no = 0;
    for (const auto &cell : dof_handler.active_cell_iterators())
      gradient_indicator(cell_no++) *=
        std::pow(cell->diameter(), 1 + 1.0 * dim / 2);

    // Finally they serve as refinement indicator.
    GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                                    gradient_indicator,
                                                    0.3,
                                                    0.1);

    triangulation.execute_coarsening_and_refinement();
  }


  // The output of this program consists of a vtk file of the adaptively
  // refined grids and the numerical solutions. Finally, we also compute the
  // L-infinity norm of the solution using VectorTools::integrate_difference().
  template <int dim>
  void AdvectionProblem<dim>::output_results(const unsigned int cycle) const
  {
    const std::string filename = "solution-" + std::to_string(cycle) + ".vtk";
    std::cout << "  Writing solution to <" << filename << '>' << std::endl;
    std::ofstream output(filename);

    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "u", DataOut<dim>::type_dof_data);

    data_out.build_patches(mapping);

    data_out.write_vtk(output);

    {
      Vector<float> values(triangulation.n_active_cells());
      VectorTools::integrate_difference(mapping,
                                        dof_handler,
                                        solution,
                                        Functions::ZeroFunction<dim>(),
                                        values,
                                        quadrature,
                                        VectorTools::Linfty_norm);
      const double l_infty =
        VectorTools::compute_global_error(triangulation,
                                          values,
                                          VectorTools::Linfty_norm);
      std::cout << "  L-infinity norm: " << l_infty << std::endl;
    }
  }


  // The following <code>run</code> function is similar to previous examples.
  template <int dim>
  void AdvectionProblem<dim>::run()
  {
    for (unsigned int cycle = 0; cycle < 6; ++cycle)
      {
        std::cout << "Cycle " << cycle << std::endl;

        if (cycle == 0)
          {
            GridGenerator::hyper_cube(triangulation);
            triangulation.refine_global(3);
          }
        else
          refine_grid();

        std::cout << "  Number of active cells:       "
                  << triangulation.n_active_cells() << std::endl;

        setup_system();

        std::cout << "  Number of degrees of freedom: " << dof_handler.n_dofs()
                  << std::endl;

        assemble_system();
        solve();

        output_results(cycle);
      }
  }
} // namespace Step12


// The following <code>main</code> function is similar to previous examples as
// well, and need not be commented on.
int main()
{
  try
    {
      Step12::AdvectionProblem<2> dgmethod;
      dgmethod.run();
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
