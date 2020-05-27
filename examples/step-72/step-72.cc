/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2021 by the deal.II authors
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
 * Authors: Jean-Paul Pelteret,
 *          Wolfgang Bangerth, Colorado State University, 2021.
 * Based on step-15, authored by Sven Wetterauer, University of Heidelberg, 2012
 */


// The majority of this tutorial is an exact replica of step-15. So, in the
// interest of brevity and maintaining a focus on the changes implemented here,
// we will only document what's new and simply indicate which sections of
// code are a repetition of what has come before.


// @sect3{Include files}

// There are a few new header files that have been included in this tutorial.
// The first is the one that provides the declaration of the ParameterAcceptor
// class.
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

// This is the second, which is an all-inclusive header that will allow us
// to incorporate the automatic differentiation (AD) functionality within this
// code.
#include <deal.II/differentiation/ad.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/affine_constraints.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/fe_q.h>

// And the next three provide some multi-threading capability using the generic
// MeshWorker::mesh_loop() framework.
#include <deal.II/meshworker/copy_data.h>
#include <deal.II/meshworker/mesh_loop.h>
#include <deal.II/meshworker/scratch_data.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>


#include <fstream>
#include <iostream>

#include <deal.II/numerics/solution_transfer.h>

// We then open a namespace for this program and import everything from the
// dealii namespace into it, as in previous programs:
namespace Step72
{
  using namespace dealii;

  // @sect3{The <code>MinimalSurfaceProblemParameters</code> class}

  // In this tutorial we will implement three different approaches for
  // assembling the linear system. One mirrors the hand implementation
  // originally provided in step-15, while the other two use the Sacado
  // automatic differentiation library that is provided as a part of the
  // Trilinos framework.
  //
  // To facilitate switching between the three implementations, we have
  // this really basic parameters class that has only two options that are
  // configurable.
  class MinimalSurfaceProblemParameters : public ParameterAcceptor
  {
  public:
    MinimalSurfaceProblemParameters();

    // Selection for the formulation and corresponding AD framework to be used:
    // -  formulation = 0 : Unassisted implementation (full hand linearization).
    // -  formulation = 1 : Automated linearization of the finite element
    //                      residual.
    // -  formulation = 2 : Automated computation of finite element
    //                      residual and linearization using a
    //                      variational formulation.
    unsigned int formulation = 0;

    // The maximum acceptable tolerance for the linear system residual.
    // We will see that the assembly time becomes appreciable once we use
    // the AD framework, so we have increased the tolerance selected in
    // step-15 by one order of magnitude. This way, the computations do
    // not take too long to complete.
    double tolerance = 1e-2;
  };


  MinimalSurfaceProblemParameters::MinimalSurfaceProblemParameters()
    : ParameterAcceptor("Minimal Surface Problem/")
  {
    add_parameter(
      "Formulation", formulation, "", this->prm, Patterns::Integer(0, 2));
    add_parameter("Tolerance", tolerance, "", this->prm, Patterns::Double(0.0));
  }



  // @sect3{The <code>MinimalSurfaceProblem</code> class template}

  // The class template is essentially the same as in step-15.
  // The only functional changes to the class are that:
  // - the run() function now takes in two arguments: one to choose which
  //   assembly approach is to be adopted, and one for the tolerance for
  //   the permissible final residual is, and
  // - there are now three different assembly functions that implement the
  //   three methods of assembling the linear system. We'll provide details
  //   on these later on.

  template <int dim>
  class MinimalSurfaceProblem
  {
  public:
    MinimalSurfaceProblem();

    void run(const int formulation, const double tolerance);

  private:
    void   setup_system(const bool initial_step);
    void   assemble_system_unassisted();
    void   assemble_system_with_residual_linearization();
    void   assemble_system_using_energy_functional();
    void   solve();
    void   refine_mesh();
    void   set_boundary_values();
    double compute_residual(const double alpha) const;
    double determine_step_length() const;
    void   output_results(const unsigned int refinement_cycle) const;

    Triangulation<dim> triangulation;

    DoFHandler<dim> dof_handler;
    FE_Q<dim>       fe;
    QGauss<dim>     quadrature_formula;

    AffineConstraints<double> hanging_node_constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double> current_solution;
    Vector<double> newton_update;
    Vector<double> system_rhs;
  };

  // @sect3{Boundary condition}

  // There are no changes to the boundary conditions applied to the problem.
  template <int dim>
  class BoundaryValues : public Function<dim>
  {
  public:
    virtual double value(const Point<dim> & p,
                         const unsigned int component = 0) const override;
  };


  template <int dim>
  double BoundaryValues<dim>::value(const Point<dim> &p,
                                    const unsigned int /*component*/) const
  {
    return std::sin(2 * numbers::PI * (p[0] + p[1]));
  }


  // @sect3{The <code>MinimalSurfaceProblem</code> class implementation}

  // @sect4{MinimalSurfaceProblem::MinimalSurfaceProblem}

  // There have been no changes made to the class constructor.
  template <int dim>
  MinimalSurfaceProblem<dim>::MinimalSurfaceProblem()
    : dof_handler(triangulation)
    , fe(2)
    , quadrature_formula(fe.degree + 1)
  {}


  // @sect4{MinimalSurfaceProblem::setup_system}

  // There have been no changes made to the function that sets up the class
  // data structures, namely the DoFHandler, the hanging node constraints
  // applied to the problem, and the linear system.
  template <int dim>
  void MinimalSurfaceProblem<dim>::setup_system(const bool initial_step)
  {
    if (initial_step)
      {
        dof_handler.distribute_dofs(fe);
        current_solution.reinit(dof_handler.n_dofs());

        hanging_node_constraints.clear();
        DoFTools::make_hanging_node_constraints(dof_handler,
                                                hanging_node_constraints);
        hanging_node_constraints.close();
      }

    newton_update.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());

    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp);

    hanging_node_constraints.condense(dsp);

    sparsity_pattern.copy_from(dsp);
    system_matrix.reinit(sparsity_pattern);
  }

  // @sect4{Assembling the linear system}

  // @sect5{Manual assembly}

  // The assembly functions are the interesting contributions to this tutorial.
  // The assemble_system_unassisted() method implements exactly the same
  // assembly function as is detailed in step-15, but in this instance we
  // use the MeshWorker::mesh_loop() function to multithread the assembly
  // process. The reason for doing this is quite simple: When using
  // automatic differentiation, we know that there is to be some additional
  // computational overhead incurred. In order to mitigate this performance
  // loss, we'd like to take advantage of as many (easily available)
  // computational resources as possible. The MeshWorker::mesh_loop() concept
  // makes this a relatively straightforward task. At the same time, for the
  // purposes of fair comparison, we need to do the same to the implementation
  // that uses no assistance when computing the residual or its linearization.
  // (The MeshWorker::mesh_loop() function is first discussed in step-12 and
  // step-16, if you'd like to read up on it.)
  //
  // The steps required to implement the multithreading are the same between the
  // three functions, so we'll use the assemble_system_unassisted() function
  // as an opportunity to focus on the multithreading itself.
  template <int dim>
  void MinimalSurfaceProblem<dim>::assemble_system_unassisted()
  {
    system_matrix = 0;
    system_rhs    = 0;

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();

    // The MeshWorker::mesh_loop() expects that we provide two exemplar data
    // structures. The first, `ScratchData`, is to store all large data that
    // is to be reused between threads. The `CopyData` will hold the
    // contributions to the linear system that come from each cell. These
    // independent matrix-vector pairs must be accumulated into the
    // global linear system sequentially. Since we don't need anything
    // on top of what the MeshWorker::ScratchData and MeshWorker::CopyData
    // classes already provide, we use these exact class definitions for
    // our problem. Note that we only require a single instance of a local
    // matrix, local right-hand side vector, and cell degree of freedom index
    // vector -- the MeshWorker::CopyData therefore has `1` for all three
    // of its template arguments.
    using ScratchData = MeshWorker::ScratchData<dim>;
    using CopyData    = MeshWorker::CopyData<1, 1, 1>;

    // We also need to know what type of iterator we'll be working with
    // during assembly. For simplicity, we just ask the compiler to work
    // this out for us using the decltype() specifier, knowing that we'll
    // be iterating over active cells owned by the @p dof_handler.
    using CellIteratorType = decltype(dof_handler.begin_active());

    // Here we initialize the exemplar data structures. Since we know that
    // we need to compute the shape function gradients, weighted Jacobian,
    // and the position of the quadrate points in real space, we pass these
    // flags into the class constructor.
    const ScratchData sample_scratch_data(fe,
                                          quadrature_formula,
                                          update_gradients |
                                            update_quadrature_points |
                                            update_JxW_values);
    const CopyData    sample_copy_data(dofs_per_cell);

    // Now we define a lambda function that will perform the assembly on
    // a single cell. The three arguments are those that will be expected by
    // MeshWorker::mesh_loop(), due to the arguments that we'll pass to that
    // final call. We also capture the @p this pointer, which means that we'll
    // have access to "this" (i.e., the current `MinimalSurfaceProblem<dim>`)
    // class instance, and its private member data (since the lambda function is
    // defined within a MinimalSurfaceProblem<dim> method).
    //
    // At the top of the function, we initialize the data structures
    // that are dependent on the cell for which the work is being
    // performed. Observe that the reinitialization call actually
    // returns an instance to an FEValues object that is initialized
    // and stored within (and, therefore, reused by) the
    // `scratch_data` object.
    //
    // Similarly, we get aliases to the local matrix, local RHS
    // vector, and local cell DoF indices from the `copy_data`
    // instance that MeshWorker::mesh_loop() provides. We then
    // initialize the cell DoF indices, knowing that the local matrix
    // and vector are already correctly sized.
    const auto cell_worker = [this](const CellIteratorType &cell,
                                    ScratchData &           scratch_data,
                                    CopyData &              copy_data) {
      const auto &fe_values = scratch_data.reinit(cell);

      FullMatrix<double> &                  cell_matrix = copy_data.matrices[0];
      Vector<double> &                      cell_rhs    = copy_data.vectors[0];
      std::vector<types::global_dof_index> &local_dof_indices =
        copy_data.local_dof_indices[0];
      cell->get_dof_indices(local_dof_indices);

      // For Newton's method, we require the gradient of the solution at the
      // point about which the problem is being linearized.
      //
      // Once we have that, we can perform assembly for this cell in
      // the usual way.  One minor difference to step-15 is that we've
      // used the (rather convenient) range-based loops to iterate
      // over all quadrature points and degrees-of-freedom.
      std::vector<Tensor<1, dim>> old_solution_gradients(
        fe_values.n_quadrature_points);
      fe_values.get_function_gradients(current_solution,
                                       old_solution_gradients);

      for (const unsigned int q : fe_values.quadrature_point_indices())
        {
          const double coeff =
            1.0 / std::sqrt(1.0 + old_solution_gradients[q] *
                                    old_solution_gradients[q]);

          for (const unsigned int i : fe_values.dof_indices())
            {
              for (const unsigned int j : fe_values.dof_indices())
                cell_matrix(i, j) +=
                  (((fe_values.shape_grad(i, q)      // ((\nabla \phi_i
                     * coeff                         //   * a_n
                     * fe_values.shape_grad(j, q))   //   * \nabla \phi_j)
                    -                                //  -
                    (fe_values.shape_grad(i, q)      //  (\nabla \phi_i
                     * coeff * coeff * coeff         //   * a_n^3
                     * (fe_values.shape_grad(j, q)   //   * (\nabla \phi_j
                        * old_solution_gradients[q]) //      * \nabla u_n)
                     * old_solution_gradients[q]))   //   * \nabla u_n)))
                   * fe_values.JxW(q));              // * dx

              cell_rhs(i) -= (fe_values.shape_grad(i, q)  // \nabla \phi_i
                              * coeff                     // * a_n
                              * old_solution_gradients[q] // * u_n
                              * fe_values.JxW(q));        // * dx
            }
        }
    };

    // The second lambda function that MeshWorker::mesh_loop() requires is
    // one that performs the task of accumulating the local contributions
    // in the global linear system. That is precisely what this one does,
    // and the details of the implementation have been seen before. The
    // primary point to recognize is that the local contributions are stored
    // in the `copy_data` instance that is passed into this function. This
    // `copy_data` has been filled with data during @a some call to the
    // `cell_worker`.
    const auto copier = [dofs_per_cell, this](const CopyData &copy_data) {
      const FullMatrix<double> &cell_matrix = copy_data.matrices[0];
      const Vector<double> &    cell_rhs    = copy_data.vectors[0];
      const std::vector<types::global_dof_index> &local_dof_indices =
        copy_data.local_dof_indices[0];

      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
            system_matrix.add(local_dof_indices[i],
                              local_dof_indices[j],
                              cell_matrix(i, j));

          system_rhs(local_dof_indices[i]) += cell_rhs(i);
        }
    };

    // We have all of the required functions definitions in place, so
    // now we call the MeshWorker::mesh_loop() to perform the actual
    // assembly.  We pass a flag as the last parameter which states
    // that we only want to perform the assembly on the
    // cells. Internally, MeshWorker::mesh_loop() then distributes the
    // available work to different threads, making efficient use of
    // the multiple cores almost all of today's processors have to
    // offer.
    MeshWorker::mesh_loop(dof_handler.active_cell_iterators(),
                          cell_worker,
                          copier,
                          sample_scratch_data,
                          sample_copy_data,
                          MeshWorker::assemble_own_cells);

    // And finally, as is done in step-15, we remove hanging nodes from the
    // system and apply zero boundary values to the linear system that defines
    // the Newton updates $\delta u^n$.
    hanging_node_constraints.condense(system_matrix);
    hanging_node_constraints.condense(system_rhs);

    std::map<types::global_dof_index, double> boundary_values;
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             Functions::ZeroFunction<dim>(),
                                             boundary_values);
    MatrixTools::apply_boundary_values(boundary_values,
                                       system_matrix,
                                       newton_update,
                                       system_rhs);
  }

  // @sect5{Assembly via differentiation of the residual vector}

  // As outlined in the introduction, what we need to do for this
  // second approach is implement the local contributions $F(U)^K$
  // from cell $K$ to the residual vector, and then let the
  // AD machinery deal with how to compute the
  // derivatives $J(U)_{ij}^K=\frac{\partial F(U)^K_i}{\partial U_j}$
  // from it.
  //
  // For the following, recall that
  // @f[
  //   F(U)_i^K \dealcoloneq
  //   \int\limits_K\nabla \varphi_i \cdot \left[ \frac{1}{\sqrt{1+|\nabla
  //   u|^{2}}} \nabla u \right] \, dV ,
  // @f]
  // where $u(\mathbf x)=\sum_j U_j \varphi_j(\mathbf x)$.
  //
  // Let us see how this is implemented in practice:
  template <int dim>
  void MinimalSurfaceProblem<dim>::assemble_system_with_residual_linearization()
  {
    system_matrix = 0;
    system_rhs    = 0;

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();

    using ScratchData      = MeshWorker::ScratchData<dim>;
    using CopyData         = MeshWorker::CopyData<1, 1, 1>;
    using CellIteratorType = decltype(dof_handler.begin_active());

    const ScratchData sample_scratch_data(fe,
                                          quadrature_formula,
                                          update_gradients |
                                            update_quadrature_points |
                                            update_JxW_values);
    const CopyData    sample_copy_data(dofs_per_cell);

    // We'll define up front the AD data structures that we'll be using,
    // utilizing the techniques shown in step-71.
    // In this case, we choose the helper class that will automatically compute
    // the linearization of the finite element residual using Sacado forward
    // automatic differentiation types. These number types can be used to
    // compute first derivatives only. This is exactly what we want, because we
    // know that we'll only be linearizing the residual, which means that we
    // only need to compute first-order derivatives. The return values from the
    // calculations are to be of type `double`.
    //
    // We also need an extractor to retrieve some data related to the field
    // solution to the problem.
    using ADHelper = Differentiation::AD::ResidualLinearization<
      Differentiation::AD::NumberTypes::sacado_dfad,
      double>;
    using ADNumberType = typename ADHelper::ad_type;

    const FEValuesExtractors::Scalar u_fe(0);

    // With this, let us define the lambda function that will be used
    // to compute the cell contributions to the Jacobian matrix and
    // the right hand side:
    const auto cell_worker = [&u_fe, this](const CellIteratorType &cell,
                                           ScratchData &           scratch_data,
                                           CopyData &              copy_data) {
      const auto &       fe_values     = scratch_data.reinit(cell);
      const unsigned int dofs_per_cell = fe_values.get_fe().n_dofs_per_cell();

      FullMatrix<double> &                  cell_matrix = copy_data.matrices[0];
      Vector<double> &                      cell_rhs    = copy_data.vectors[0];
      std::vector<types::global_dof_index> &local_dof_indices =
        copy_data.local_dof_indices[0];
      cell->get_dof_indices(local_dof_indices);

      // We'll now create and initialize an instance of the AD helper class.
      // To do this, we need to specify how many independent variables and
      // dependent variables there are. The independent variables will be the
      // number of local degrees of freedom that our solution vector has,
      // i.e., the number $j$ in the per-element representation of the
      // discretized solution vector
      // $u (\mathbf{x})|_K = \sum\limits_{j} U^K_i \varphi_j(\mathbf{x})$
      // that indicates how many solution coefficients are associated with
      // each finite element. In deal.II, this equals
      // FiniteElement::dofs_per_cell. The number of dependent variables will be
      // the number of entries in the local residual vector that we will be
      // forming. In this particular problem (like many others that employ the
      // [standard Galerkin
      // method](https://en.wikipedia.org/wiki/Galerkin_method)) the number of
      // local solution coefficients matches the number of local residual
      // equations.
      const unsigned int n_independent_variables = local_dof_indices.size();
      const unsigned int n_dependent_variables   = dofs_per_cell;
      ADHelper ad_helper(n_independent_variables, n_dependent_variables);

      // Next we inform the helper of the values of the solution, i.e., the
      // actual values for $U_j$ about which we
      // wish to linearize. As this is done on each element individually, we
      // have to extract the solution coefficients from the global solution
      // vector. In other words, we define all of those coefficients $U_j$
      // where $j$ is a local degree of freedom as the independent variables
      // that enter the computation of the vector $F(U)^{K}$ (the dependent
      // function).
      //
      // Then we get the complete set of degree of freedom values as
      // represented by auto-differentiable numbers. The operations
      // performed with these variables are tracked by the AD library
      // from this point until the object goes out of scope. So it is
      // <em>precisely these variables</em> with respect to which we will
      // compute derivatives of the residual entries.
      ad_helper.register_dof_values(current_solution, local_dof_indices);

      const std::vector<ADNumberType> &dof_values_ad =
        ad_helper.get_sensitive_dof_values();

      // Then we do some problem specific tasks, the first being to
      // compute all values, (spatial) gradients, and the like based on
      // "sensitive" AD degree of freedom values. In this instance we are
      // retrieving the solution gradients at each quadrature point. Observe
      // that the solution gradients are now sensitive
      // to the values of the degrees of freedom as they use the @p ADNumberType
      // as the scalar type and the @p dof_values_ad vector provides the local
      // DoF values.
      std::vector<Tensor<1, dim, ADNumberType>> old_solution_gradients(
        fe_values.n_quadrature_points);
      fe_values[u_fe].get_function_gradients_from_local_dof_values(
        dof_values_ad, old_solution_gradients);

      // The next variable that we declare will store the cell residual vector
      // contributions. This is rather self-explanatory, save for one
      // <b>very important</b> detail:
      // Note that each entry in the vector is hand-initialized with a value
      // of zero. This is a <em>highly recommended</em> practice, as some AD
      // libraries appear not to safely initialize the internal data
      // structures of these number types. Not doing so could lead to some
      // very hard to understand or detect bugs (appreciate that the author
      // of this program mentions this out of, generally bad, experience). So
      // out of an abundance of caution it's worthwhile zeroing the initial
      // value explicitly. After that, apart from a sign change the residual
      // assembly looks much the same as we saw for the cell RHS vector before:
      // We loop over all quadrature points, ensure that the coefficient now
      // encodes its dependence on the (sensitive) finite element DoF values by
      // using the correct `ADNumberType`, and finally we assemble the
      // components of the residual vector. For complete clarity, the finite
      // element shape functions (and their gradients, etc.) as well as the
      // "JxW" values remain scalar
      // valued, but the @p coeff and the  @p old_solution_gradients at each
      // quadrature point are computed in terms of the independent
      // variables.
      std::vector<ADNumberType> residual_ad(n_dependent_variables,
                                            ADNumberType(0.0));
      for (const unsigned int q : fe_values.quadrature_point_indices())
        {
          const ADNumberType coeff =
            1.0 / std::sqrt(1.0 + old_solution_gradients[q] *
                                    old_solution_gradients[q]);

          for (const unsigned int i : fe_values.dof_indices())
            {
              residual_ad[i] += (fe_values.shape_grad(i, q)   // \nabla \phi_i
                                 * coeff                      // * a_n
                                 * old_solution_gradients[q]) // * u_n
                                * fe_values.JxW(q);           // * dx
            }
        }

      // Once we have the full cell residual vector computed, we can register
      // it with the helper class.
      //
      // Thereafter, we compute the residual values (basically,
      // extracting the real values from what we already computed) and
      // their Jacobian (the linearization of each residual component
      // with respect to all cell DoFs) at the evaluation point. For
      // the purposes of assembly into the global linear system, we
      // have to respect the sign difference between the residual and
      // the RHS contribution: For Newton's method, the right hand
      // side vector needs to be equal to the *negative* residual
      // vector.
      ad_helper.register_residual_vector(residual_ad);

      ad_helper.compute_residual(cell_rhs);
      cell_rhs *= -1.0;

      ad_helper.compute_linearization(cell_matrix);
    };

    // The remainder of the function equals what we had previously:
    const auto copier = [dofs_per_cell, this](const CopyData &copy_data) {
      const FullMatrix<double> &cell_matrix = copy_data.matrices[0];
      const Vector<double> &    cell_rhs    = copy_data.vectors[0];
      const std::vector<types::global_dof_index> &local_dof_indices =
        copy_data.local_dof_indices[0];

      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
            system_matrix.add(local_dof_indices[i],
                              local_dof_indices[j],
                              cell_matrix(i, j));

          system_rhs(local_dof_indices[i]) += cell_rhs(i);
        }
    };

    MeshWorker::mesh_loop(dof_handler.active_cell_iterators(),
                          cell_worker,
                          copier,
                          sample_scratch_data,
                          sample_copy_data,
                          MeshWorker::assemble_own_cells);

    hanging_node_constraints.condense(system_matrix);
    hanging_node_constraints.condense(system_rhs);

    std::map<types::global_dof_index, double> boundary_values;
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             Functions::ZeroFunction<dim>(),
                                             boundary_values);
    MatrixTools::apply_boundary_values(boundary_values,
                                       system_matrix,
                                       newton_update,
                                       system_rhs);
  }

  // @sect5{Assembly via differentiation of the energy functional}

  // In this third approach, we compute residual and Jacobian as first
  // and second derivatives of the local energy functional
  // @f[
  //    E\left( U \right)^K
  //     \dealcoloneq \int\limits_{K} \Psi \left( u \right) \, dV
  //     \approx \sum\limits_{q}^{n_{\textrm{q-points}}} \Psi \left( u \left(
  //     \mathbf{X}_{q} \right) \right) \underbrace{\vert J_{q} \vert \times
  //     W_{q}}_{\text{JxW(q)}}
  // @f]
  // with the energy density given by
  // @f[
  //   \Psi \left( u \right) = \sqrt{1+|\nabla u|^{2}} .
  // @f]
  //
  // Let us again see how this is done:
  template <int dim>
  void MinimalSurfaceProblem<dim>::assemble_system_using_energy_functional()
  {
    system_matrix = 0;
    system_rhs    = 0;

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();

    using ScratchData      = MeshWorker::ScratchData<dim>;
    using CopyData         = MeshWorker::CopyData<1, 1, 1>;
    using CellIteratorType = decltype(dof_handler.begin_active());

    const ScratchData sample_scratch_data(fe,
                                          quadrature_formula,
                                          update_gradients |
                                            update_quadrature_points |
                                            update_JxW_values);
    const CopyData    sample_copy_data(dofs_per_cell);

    // In this implementation of the assembly process, we choose the helper
    // class that will automatically compute both the residual and its
    // linearization from the cell contribution to an energy functional using
    // nested Sacado forward automatic differentiation types.
    // The selected number types can be used to compute both first and
    // second derivatives. We require this, as the residual defined as the
    // sensitivity of the potential energy with respect to the DoF values (i.e.
    // its gradient). We'll then need to linearize the residual, implying that
    // second derivatives of the potential energy must be computed. You might
    // want to compare this with the definition of `ADHelper` used int
    // previous function, where we used
    // `Differentiation::AD::ResidualLinearization<Differentiation::AD::NumberTypes::sacado_dfad,double>`.
    using ADHelper = Differentiation::AD::EnergyFunctional<
      Differentiation::AD::NumberTypes::sacado_dfad_dfad,
      double>;
    using ADNumberType = typename ADHelper::ad_type;

    const FEValuesExtractors::Scalar u_fe(0);

    // Let us then again define the lambda function that does the integration on
    // a cell.
    //
    // To initialize an instance of the helper class, we now only require
    // that the number of independent variables (that is, the number
    // of degrees of freedom associated with the element solution vector)
    // are known up front. This is because the second-derivative matrix that
    // results from an energy functional is necessarily square (and also,
    // incidentally, symmetric).
    const auto cell_worker = [&u_fe, this](const CellIteratorType &cell,
                                           ScratchData &           scratch_data,
                                           CopyData &              copy_data) {
      const auto &fe_values = scratch_data.reinit(cell);

      FullMatrix<double> &                  cell_matrix = copy_data.matrices[0];
      Vector<double> &                      cell_rhs    = copy_data.vectors[0];
      std::vector<types::global_dof_index> &local_dof_indices =
        copy_data.local_dof_indices[0];
      cell->get_dof_indices(local_dof_indices);

      const unsigned int n_independent_variables = local_dof_indices.size();
      ADHelper           ad_helper(n_independent_variables);

      // Once more, we register all cell DoFs values with the helper, followed
      // by extracting the "sensitive" variant of these values that are to be
      // used in subsequent operations that must be differentiated -- one of
      // those being the calculation of the solution gradients.
      ad_helper.register_dof_values(current_solution, local_dof_indices);

      const std::vector<ADNumberType> &dof_values_ad =
        ad_helper.get_sensitive_dof_values();

      std::vector<Tensor<1, dim, ADNumberType>> old_solution_gradients(
        fe_values.n_quadrature_points);
      fe_values[u_fe].get_function_gradients_from_local_dof_values(
        dof_values_ad, old_solution_gradients);

      // We next create a variable that stores the cell total energy.
      // Once more we emphasize that we explicitly zero-initialize this value,
      // thereby ensuring the integrity of the data for this starting value.
      //
      // The aim for our approach is then to compute the cell total
      // energy, which is the sum of the internal (due to right hand
      // side functions, typically linear in $U$) and external
      // energies. In this particular case, we have no external
      // energies (e.g., from source terms or Neumann boundary
      // conditions), so we'll focus on the internal energy part.
      //
      // In fact, computing $E(U)^K$ is almost trivial, requiring only
      // the following lines:
      ADNumberType energy_ad = ADNumberType(0.0);
      for (const unsigned int q : fe_values.quadrature_point_indices())
        {
          const ADNumberType psi = std::sqrt(1.0 + old_solution_gradients[q] *
                                                     old_solution_gradients[q]);

          energy_ad += psi * fe_values.JxW(q);
        }

      // After we've computed the total energy on this cell, we'll
      // register it with the helper.  Based on that, we may now
      // compute the desired quantities, namely the residual values
      // and their Jacobian at the evaluation point. As before, the
      // Newton right hand side needs to be the negative of the
      // residual:
      ad_helper.register_energy_functional(energy_ad);

      ad_helper.compute_residual(cell_rhs);
      cell_rhs *= -1.0;

      ad_helper.compute_linearization(cell_matrix);
    };

    // As in the previous two functions, the remainder of the function is as
    // before:
    const auto copier = [dofs_per_cell, this](const CopyData &copy_data) {
      const FullMatrix<double> &cell_matrix = copy_data.matrices[0];
      const Vector<double> &    cell_rhs    = copy_data.vectors[0];
      const std::vector<types::global_dof_index> &local_dof_indices =
        copy_data.local_dof_indices[0];

      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
            system_matrix.add(local_dof_indices[i],
                              local_dof_indices[j],
                              cell_matrix(i, j));

          system_rhs(local_dof_indices[i]) += cell_rhs(i);
        }
    };

    MeshWorker::mesh_loop(dof_handler.active_cell_iterators(),
                          cell_worker,
                          copier,
                          sample_scratch_data,
                          sample_copy_data,
                          MeshWorker::assemble_own_cells);

    hanging_node_constraints.condense(system_matrix);
    hanging_node_constraints.condense(system_rhs);

    std::map<types::global_dof_index, double> boundary_values;
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             Functions::ZeroFunction<dim>(),
                                             boundary_values);
    MatrixTools::apply_boundary_values(boundary_values,
                                       system_matrix,
                                       newton_update,
                                       system_rhs);
  }


  // @sect4{MinimalSurfaceProblem::solve}

  // The solve function is the same as is used in step-15.
  template <int dim>
  void MinimalSurfaceProblem<dim>::solve()
  {
    SolverControl            solver_control(system_rhs.size(),
                                 system_rhs.l2_norm() * 1e-6);
    SolverCG<Vector<double>> solver(solver_control);

    PreconditionSSOR<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);

    solver.solve(system_matrix, newton_update, system_rhs, preconditioner);

    hanging_node_constraints.distribute(newton_update);

    const double alpha = determine_step_length();
    current_solution.add(alpha, newton_update);
  }


  // @sect4{MinimalSurfaceProblem::refine_mesh}

  // Nothing has changed since step-15 with respect to the mesh refinement
  // procedure and transfer of the solution between adapted meshes.
  template <int dim>
  void MinimalSurfaceProblem<dim>::refine_mesh()
  {
    Vector<float> estimated_error_per_cell(triangulation.n_active_cells());

    KellyErrorEstimator<dim>::estimate(
      dof_handler,
      QGauss<dim - 1>(fe.degree + 1),
      std::map<types::boundary_id, const Function<dim> *>(),
      current_solution,
      estimated_error_per_cell);

    GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                                    estimated_error_per_cell,
                                                    0.3,
                                                    0.03);

    triangulation.prepare_coarsening_and_refinement();
    SolutionTransfer<dim> solution_transfer(dof_handler);
    solution_transfer.prepare_for_coarsening_and_refinement(current_solution);
    triangulation.execute_coarsening_and_refinement();

    dof_handler.distribute_dofs(fe);

    Vector<double> tmp(dof_handler.n_dofs());
    solution_transfer.interpolate(current_solution, tmp);
    current_solution = tmp;

    hanging_node_constraints.clear();
    DoFTools::make_hanging_node_constraints(dof_handler,
                                            hanging_node_constraints);
    hanging_node_constraints.close();

    set_boundary_values();

    setup_system(false);
  }



  // @sect4{MinimalSurfaceProblem::set_boundary_values}

  // The choice of boundary conditions remains identical to step-15...
  template <int dim>
  void MinimalSurfaceProblem<dim>::set_boundary_values()
  {
    std::map<types::global_dof_index, double> boundary_values;
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             BoundaryValues<dim>(),
                                             boundary_values);
    for (auto &boundary_value : boundary_values)
      current_solution(boundary_value.first) = boundary_value.second;

    hanging_node_constraints.distribute(current_solution);
  }


  // @sect4{MinimalSurfaceProblem::compute_residual}

  // ... as does the function used to compute the residual during the
  // solution iteration procedure. One could replace this by
  // differentiation of the energy functional if one really wanted,
  // but for simplicity we here simply copy what we already had in
  // step-15.
  template <int dim>
  double MinimalSurfaceProblem<dim>::compute_residual(const double alpha) const
  {
    Vector<double> residual(dof_handler.n_dofs());

    Vector<double> evaluation_point(dof_handler.n_dofs());
    evaluation_point = current_solution;
    evaluation_point.add(alpha, newton_update);

    const QGauss<dim> quadrature_formula(fe.degree + 1);
    FEValues<dim>     fe_values(fe,
                            quadrature_formula,
                            update_gradients | update_quadrature_points |
                              update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    Vector<double>              cell_residual(dofs_per_cell);
    std::vector<Tensor<1, dim>> gradients(n_q_points);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        cell_residual = 0;
        fe_values.reinit(cell);

        fe_values.get_function_gradients(evaluation_point, gradients);

        for (unsigned int q = 0; q < n_q_points; ++q)
          {
            const double coeff =
              1.0 / std::sqrt(1.0 + gradients[q] * gradients[q]);

            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              cell_residual(i) -= (fe_values.shape_grad(i, q) // \nabla \phi_i
                                   * coeff                    // * a_n
                                   * gradients[q]             // * u_n
                                   * fe_values.JxW(q));       // * dx
          }

        cell->get_dof_indices(local_dof_indices);
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          residual(local_dof_indices[i]) += cell_residual(i);
      }

    hanging_node_constraints.condense(residual);

    for (types::global_dof_index i :
         DoFTools::extract_boundary_dofs(dof_handler))
      residual(i) = 0;

    return residual.l2_norm();
  }



  // @sect4{MinimalSurfaceProblem::determine_step_length}

  // The choice of step length (or, under-relaxation factor) for the nonlinear
  // iterations procedure remains fixed at the value chosen and discussed in
  // step-15.
  template <int dim>
  double MinimalSurfaceProblem<dim>::determine_step_length() const
  {
    return 0.1;
  }



  // @sect4{MinimalSurfaceProblem::output_results}

  // This last function to be called from `run()` outputs the current solution
  // (and the Newton update) in graphical form as a VTU file. It is entirely the
  // same as what has been used in previous tutorials.
  template <int dim>
  void MinimalSurfaceProblem<dim>::output_results(
    const unsigned int refinement_cycle) const
  {
    DataOut<dim> data_out;

    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(current_solution, "solution");
    data_out.add_data_vector(newton_update, "update");
    data_out.build_patches();

    const std::string filename =
      "solution-" + Utilities::int_to_string(refinement_cycle, 2) + ".vtu";
    std::ofstream output(filename);
    data_out.write_vtu(output);
  }


  // @sect4{MinimalSurfaceProblem::run}

  // In the run function, most remains the same as was first implemented
  // in step-15. The only observable changes are that we can now choose (via
  // the parameter file) what the final acceptable tolerance for the system
  // residual is, and that we can choose which method of assembly we wish to
  // utilize. To make the second choice clear, we output to the console some
  // message which indicates the selection. Since we're interested in comparing
  // the time taken to assemble for each of the three methods, we've also
  // added a timer that keeps a track of how much time is spent during assembly.
  // We also track the time taken to solve the linear system, so that we can
  // contrast those numbers to the part of the code which would normally take
  // the longest time to execute.
  template <int dim>
  void MinimalSurfaceProblem<dim>::run(const int    formulation,
                                       const double tolerance)
  {
    std::cout << "******** Assembly approach ********" << std::endl;
    const std::array<std::string, 3> method_descriptions = {
      {"Unassisted implementation (full hand linearization).",
       "Automated linearization of the finite element residual.",
       "Automated computation of finite element residual and linearization using a variational formulation."}};
    AssertIndexRange(formulation, method_descriptions.size());
    std::cout << method_descriptions[formulation] << std::endl << std::endl;


    TimerOutput timer(std::cout, TimerOutput::summary, TimerOutput::wall_times);

    GridGenerator::hyper_ball(triangulation);
    triangulation.refine_global(2);

    setup_system(/*first time=*/true);
    set_boundary_values();

    double       last_residual_norm = std::numeric_limits<double>::max();
    unsigned int refinement_cycle   = 0;
    do
      {
        std::cout << "Mesh refinement step " << refinement_cycle << std::endl;

        if (refinement_cycle != 0)
          refine_mesh();

        std::cout << "  Initial residual: " << compute_residual(0) << std::endl;

        for (unsigned int inner_iteration = 0; inner_iteration < 5;
             ++inner_iteration)
          {
            {
              TimerOutput::Scope t(timer, "Assemble");

              if (formulation == 0)
                assemble_system_unassisted();
              else if (formulation == 1)
                assemble_system_with_residual_linearization();
              else if (formulation == 2)
                assemble_system_using_energy_functional();
              else
                AssertThrow(false, ExcNotImplemented());
            }

            last_residual_norm = system_rhs.l2_norm();

            {
              TimerOutput::Scope t(timer, "Solve");
              solve();
            }


            std::cout << "  Residual: " << compute_residual(0) << std::endl;
          }

        output_results(refinement_cycle);

        ++refinement_cycle;
        std::cout << std::endl;
      }
    while (last_residual_norm > tolerance);
  }
} // namespace Step72

// @sect4{The main function}

// Finally the main function. This follows the scheme of most other main
// functions, with two obvious exceptions:
// - We call Utilities::MPI::MPI_InitFinalize in order to set up (via a hidden
//   default parameter) the number of threads using the execution of
//   multithreaded tasks.
// - We also have a few lines dedicates to reading in or initializing the
//   user-defined parameters that will be considered during the execution of the
//   program.
int main(int argc, char *argv[])
{
  try
    {
      using namespace Step72;

      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);

      std::string prm_file;
      if (argc > 1)
        prm_file = argv[1];
      else
        prm_file = "parameters.prm";

      const MinimalSurfaceProblemParameters parameters;
      ParameterAcceptor::initialize(prm_file);

      MinimalSurfaceProblem<2> minimal_surface_problem_2d;
      minimal_surface_problem_2d.run(parameters.formulation,
                                     parameters.tolerance);
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
