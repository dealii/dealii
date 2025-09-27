/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2005 - 2024 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 */


// @sect3{Include files}

// Since this program is only an adaptation of step-4, there is not much new
// stuff in terms of header files. In deal.II, we usually list include files
// in the order base-lac-grid-dofs-fe-numerics, followed by C++ standard
// include files:
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

// The only two new header files that deserve some attention are those for
// the LinearOperator and PackagedOperation classes:
#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/packaged_operation.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <iostream>

// This is the only significant new header, namely the one in which the
// Raviart-Thomas finite element is declared:
#include <deal.II/fe/fe_raviart_thomas.h>

// Finally, as a bonus in this program, we will use a tensorial
// coefficient. Since it may have a spatial dependence, we consider it a
// tensor-valued function. The following include file provides the
// <code>TensorFunction</code> class that offers such functionality:
#include <deal.II/base/tensor_function.h>

// The last step is as in all previous programs: We put all of the code relevant
// to this program into a namespace. (This idea was first introduced in step-7.)
namespace Step20
{
  using namespace dealii;

  // @sect3{The <code>MixedLaplaceProblem</code> class template}

  // Again, since this is an adaptation of step-6, the main class is almost
  // the same as the one in that tutorial program. In terms of member
  // functions, the main differences are that the constructor takes the degree
  // of the Raviart-Thomas element as an argument (and that there is a
  // corresponding member variable to store this value) and the addition of
  // the <code>compute_error</code> function in which, no surprise, we will
  // compute the difference between the exact and the numerical solution to
  // determine convergence of our computations:
  template <int dim>
  class MixedLaplaceProblem
  {
  public:
    MixedLaplaceProblem(const unsigned int degree);
    void run();

  private:
    void make_grid_and_dofs();
    void assemble_system();
    void solve();
    void compute_errors() const;
    void output_results() const;

    const unsigned int degree;

    Triangulation<dim>  triangulation;
    const FESystem<dim> fe;
    DoFHandler<dim>     dof_handler;

    // The second difference is that the sparsity pattern, the system matrix,
    // and solution and right hand side vectors are now blocked. What this
    // means and what one can do with such objects is explained in the
    // introduction to this program as well as further down below when we
    // explain the linear solvers and preconditioners for this problem:
    BlockSparsityPattern      sparsity_pattern;
    BlockSparseMatrix<double> system_matrix;

    BlockVector<double> solution;
    BlockVector<double> system_rhs;
  };


  // @sect3{Right hand side, boundary values, and exact solution}

  // Our next task is to define the right hand side of our problem (i.e., the
  // scalar right hand side for the pressure in the original Laplace
  // equation), boundary values for the pressure, and a function that
  // describes both the pressure and the velocity of the exact solution for
  // later computations of the error. Note that these functions have one, one,
  // and <code>dim+1</code> components, respectively, and that we pass the
  // number of components down to the <code>Function@<dim@></code> base
  // class. For the exact solution, we only declare the function that actually
  // returns the entire solution vector (i.e. all components of it) at
  // once. Here are the respective declarations:
  namespace PrescribedSolution
  {
    constexpr double alpha = 0.3;
    constexpr double beta  = 1;


    template <int dim>
    class RightHandSide : public Function<dim>
    {
    public:
      RightHandSide()
        : Function<dim>(1)
      {}

      virtual double value(const Point<dim>  &p,
                           const unsigned int component = 0) const override;
    };



    template <int dim>
    class PressureBoundaryValues : public Function<dim>
    {
    public:
      PressureBoundaryValues()
        : Function<dim>(1)
      {}

      virtual double value(const Point<dim>  &p,
                           const unsigned int component = 0) const override;
    };


    template <int dim>
    class ExactSolution : public Function<dim>
    {
    public:
      ExactSolution()
        : Function<dim>(dim + 1)
      {}

      virtual void vector_value(const Point<dim> &p,
                                Vector<double>   &value) const override;
    };


    // And then we also have to define these respective functions, of
    // course. Given our discussion in the introduction of how the solution
    // should look, the following computations should be straightforward:
    template <int dim>
    double RightHandSide<dim>::value(const Point<dim> & /*p*/,
                                     const unsigned int /*component*/) const
    {
      return 0;
    }



    template <int dim>
    double
    PressureBoundaryValues<dim>::value(const Point<dim> &p,
                                       const unsigned int /*component*/) const
    {
      return -(alpha * p[0] * p[1] * p[1] / 2 + beta * p[0] -
               alpha * p[0] * p[0] * p[0] / 6);
    }



    template <int dim>
    void ExactSolution<dim>::vector_value(const Point<dim> &p,
                                          Vector<double>   &values) const
    {
      AssertDimension(values.size(), dim + 1);

      values(0) = alpha * p[1] * p[1] / 2 + beta - alpha * p[0] * p[0] / 2;
      values(1) = alpha * p[0] * p[1];
      values(2) = -(alpha * p[0] * p[1] * p[1] / 2 + beta * p[0] -
                    alpha * p[0] * p[0] * p[0] / 6);
    }



    // @sect3{The inverse permeability tensor}

    // In addition to the other equation data, we also want to use a
    // permeability tensor, or better -- because this is all that appears in the
    // weak form -- the inverse of the permeability tensor,
    // <code>KInverse</code>. For the purpose of verifying the exactness of the
    // solution and determining convergence orders, this tensor is more in the
    // way than helpful. We will therefore simply set it to the identity matrix.
    //
    // However, a spatially varying permeability tensor is indispensable in
    // real-life porous media flow simulations, and we would like to use the
    // opportunity to demonstrate the technique to use tensor valued functions.
    //
    // Possibly unsurprisingly, deal.II also has a base class not only for
    // scalar and generally vector-valued functions (the <code>Function</code>
    // base class) but also for functions that return tensors of fixed dimension
    // and rank, the <code>TensorFunction</code> template. Here, the function
    // under consideration returns a dim-by-dim matrix, i.e. a tensor of rank 2
    // and dimension <code>dim</code>. We then choose the template arguments of
    // the base class appropriately.
    //
    // The interface that the <code>TensorFunction</code> class provides is
    // essentially equivalent to the <code>Function</code> class. In particular,
    // there exists a <code>value_list</code> function that takes a list of
    // points at which to evaluate the function, and returns the values of the
    // function in the second argument, a list of tensors:
    template <int dim>
    class KInverse : public TensorFunction<2, dim>
    {
    public:
      KInverse()
        : TensorFunction<2, dim>()
      {}

      virtual void
      value_list(const std::vector<Point<dim>> &points,
                 std::vector<Tensor<2, dim>>   &values) const override;
    };


    // The implementation is less interesting. As in previous examples, we add a
    // check to the beginning of the class to make sure that the sizes of input
    // and output parameters are the same (see step-5 for a discussion of this
    // technique). Then we loop over all evaluation points, and for each one
    // set the output tensor to the identity matrix.
    //
    // There is an oddity at the top of the function (the
    // `(void)points;` statement) that is worth discussing. The values
    // we put into the output `values` array does not actually depend
    // on the `points` arrays of coordinates at which the function is
    // evaluated. In other words, the `points` argument is in fact
    // unused, and we could have just not given it a name if we had
    // wanted. But we want to use the `points` object for checking
    // that the `values` object has the correct size. The problem is
    // that in release mode, `AssertDimension` is defined as a macro
    // that expands to nothing; the compiler will then complain that
    // the `points` object is unused. The idiomatic approach to
    // silencing this warning is to have a statement that evaluates
    // (reads) variable but doesn't actually do anything: That's what
    // `(void)points;` does: It reads from `points`, and then casts
    // the result of the read to `void`, i.e., nothing. This statement
    // is, in other words, completely pointless and implies no actual
    // action except to explain to the compiler that yes, this
    // variable is in fact used even in release mode. (In debug mode,
    // the `AssertDimension` macro expands to something that reads
    // from the variable, and so the funny statement would not be
    // necessary in debug mode.)
    template <int dim>
    void KInverse<dim>::value_list(const std::vector<Point<dim>> &points,
                                   std::vector<Tensor<2, dim>>   &values) const
    {
      (void)points;
      AssertDimension(points.size(), values.size());

      for (auto &value : values)
        value = unit_symmetric_tensor<dim>();
    }
  } // namespace PrescribedSolution



  // @sect3{MixedLaplaceProblem class implementation}

  // @sect4{MixedLaplaceProblem::MixedLaplaceProblem}

  // In the constructor of this class, we first store the value that was
  // passed in concerning the degree of the finite elements we shall use (a
  // degree of zero, for example, means to use RT(0) and DG(0)), and then
  // construct the vector valued element belonging to the space $X_h$ described
  // in the introduction. The rest of the constructor is as in the early
  // tutorial programs.
  //
  // The only thing worth describing here is the constructor call of the
  // <code>fe</code> variable. The <code>FESystem</code> class to which this
  // variable belongs has a number of different constructors that all refer to
  // binding simpler elements together into one larger element. In the present
  // case, we want to couple a single RT(degree) element with a single
  // DQ(degree) element. The constructor to <code>FESystem</code> that does
  // this requires us to specify first the first base element (the
  // <code>FE_RaviartThomas</code> object of given degree) and then the number
  // of copies for this base element, and then similarly the kind and number
  // of <code>FE_DGQ</code> elements. Note that the Raviart-Thomas element
  // already has <code>dim</code> vector components, so that the coupled
  // element will have <code>dim+1</code> vector components, the first
  // <code>dim</code> of which correspond to the velocity variable whereas the
  // last one corresponds to the pressure.
  //
  // It is also worth comparing the way we constructed this element from its
  // base elements, with the way we have done so in step-8: there, we have
  // built it as <code>fe (FE_Q@<dim@>(1), dim)</code>, i.e. we have simply
  // used <code>dim</code> copies of the <code>FE_Q(1)</code> element, one
  // copy for the displacement in each coordinate direction.
  template <int dim>
  MixedLaplaceProblem<dim>::MixedLaplaceProblem(const unsigned int degree)
    : degree(degree)
    , fe(FE_RaviartThomas<dim>(degree), FE_DGQ<dim>(degree))
    , dof_handler(triangulation)
  {}



  // @sect4{MixedLaplaceProblem::make_grid_and_dofs}

  // This next function starts out with well-known functions calls that create
  // and refine a mesh, and then associate degrees of freedom with it:
  template <int dim>
  void MixedLaplaceProblem<dim>::make_grid_and_dofs()
  {
    GridGenerator::hyper_cube(triangulation, -1, 1);
    triangulation.refine_global(5);

    dof_handler.distribute_dofs(fe);

    // However, then things become different. As mentioned in the
    // introduction, we want to subdivide the matrix into blocks corresponding
    // to the two different kinds of variables, velocity and pressure. To this
    // end, we first have to make sure that the indices corresponding to
    // velocities and pressures are not intermingled: First all velocity
    // degrees of freedom, then all pressure DoFs. This way, the global matrix
    // separates nicely into a $2 \times 2$ system. To achieve this, we have to
    // renumber degrees of freedom based on their vector component, an
    // operation that conveniently is already implemented:
    DoFRenumbering::component_wise(dof_handler);

    // The next thing is that we want to figure out the sizes of these blocks
    // so that we can allocate an appropriate amount of space. To this end, we
    // call the DoFTools::count_dofs_per_fe_component() function that
    // counts how many shape functions are non-zero for a particular vector
    // component. We have <code>dim+1</code> vector components, and
    // DoFTools::count_dofs_per_fe_component() will count how many shape
    // functions belong to each of these components.
    //
    // There is one problem here. As described in the documentation of that
    // function, it <i>wants</i> to put the number of $x$-velocity shape
    // functions into <code>dofs_per_component[0]</code>, the number of
    // $y$-velocity shape functions into <code>dofs_per_component[1]</code>
    // (and similar in 3d), and the number of pressure shape functions into
    // <code>dofs_per_component[dim]</code>. But, the Raviart-Thomas element
    // is special in that it is non-@ref GlossPrimitive "primitive", i.e.,
    // for Raviart-Thomas elements all velocity shape functions
    // are nonzero in all components. In other words, the function cannot
    // distinguish between $x$ and $y$ velocity functions because there
    // <i>is</i> no such distinction. It therefore puts the overall number
    // of velocity into each of <code>dofs_per_component[c]</code>,
    // $0\le c\le \text{dim}$. On the other hand, the number
    // of pressure variables equals the number of shape functions that are
    // nonzero in the dim-th component.
    //
    // Using this knowledge, we can get the number of velocity shape
    // functions from any of the first <code>dim</code> elements of
    // <code>dofs_per_component</code>, and then use this below to initialize
    // the vector and matrix block sizes, as well as create output.
    //
    // @note If you find this concept difficult to understand, you may
    // want to consider using the function DoFTools::count_dofs_per_fe_block()
    // instead, as we do in the corresponding piece of code in step-22.
    // You might also want to read up on the difference between
    // @ref GlossBlock "blocks" and @ref GlossComponent "components"
    // in the glossary.
    const std::vector<types::global_dof_index> dofs_per_component =
      DoFTools::count_dofs_per_fe_component(dof_handler);
    const unsigned int n_u = dofs_per_component[0],
                       n_p = dofs_per_component[dim];

    std::cout << "Number of active cells: " << triangulation.n_active_cells()
              << std::endl
              << "Total number of cells: " << triangulation.n_cells()
              << std::endl
              << "Number of degrees of freedom: " << dof_handler.n_dofs()
              << " (" << n_u << '+' << n_p << ')' << std::endl;

    // The next task is to allocate a sparsity pattern for the matrix that we
    // will create. We use a compressed sparsity pattern like in the previous
    // steps, but as <code>system_matrix</code> is a block matrix we use the
    // class <code>BlockDynamicSparsityPattern</code> instead of just
    // <code>DynamicSparsityPattern</code>. This block sparsity pattern has
    // four blocks in a $2 \times 2$ pattern. The blocks' sizes depend on
    // <code>n_u</code> and <code>n_p</code>, which hold the number of velocity
    // and pressure variables.
    const std::vector<types::global_dof_index> block_sizes = {n_u, n_p};
    BlockDynamicSparsityPattern                dsp(block_sizes, block_sizes);
    DoFTools::make_sparsity_pattern(dof_handler, dsp);

    // We use the compressed block sparsity pattern in the same way as the
    // non-block version to create the sparsity pattern and then the system
    // matrix:
    sparsity_pattern.copy_from(dsp);
    system_matrix.reinit(sparsity_pattern);

    // Then we have to resize the solution and right hand side vectors in
    // exactly the same way as the block compressed sparsity pattern:
    solution.reinit(block_sizes);
    system_rhs.reinit(block_sizes);
  }


  // @sect4{MixedLaplaceProblem::assemble_system}

  // Similarly, the function that assembles the linear system has mostly been
  // discussed already in the introduction to this example. At its top, what
  // happens are all the usual steps, with the addition that we do not only
  // allocate quadrature and <code>FEValues</code> objects for the cell terms,
  // but also for face terms. After that, we define the usual abbreviations
  // for variables, and the allocate space for the local matrix and right hand
  // side contributions, and the array that holds the global numbers of the
  // degrees of freedom local to the present cell.
  template <int dim>
  void MixedLaplaceProblem<dim>::assemble_system()
  {
    const QGauss<dim>     quadrature_formula(degree + 2);
    const QGauss<dim - 1> face_quadrature_formula(degree + 2);

    FEValues<dim>     fe_values(fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);
    FEFaceValues<dim> fe_face_values(fe,
                                     face_quadrature_formula,
                                     update_values | update_normal_vectors |
                                       update_quadrature_points |
                                       update_JxW_values);

    const unsigned int dofs_per_cell   = fe.n_dofs_per_cell();
    const unsigned int n_q_points      = quadrature_formula.size();
    const unsigned int n_face_q_points = face_quadrature_formula.size();

    FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     local_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    // The next step is to declare objects that represent the source term,
    // pressure boundary value, and coefficient in the equation. In addition
    // to these objects that represent continuous functions, we also need
    // arrays to hold their values at the quadrature points of individual
    // cells (or faces, for the boundary values). Note that in the case of the
    // coefficient, the array has to be one of matrices.
    const PrescribedSolution::RightHandSide<dim> right_hand_side;
    const PrescribedSolution::PressureBoundaryValues<dim>
                                            pressure_boundary_values;
    const PrescribedSolution::KInverse<dim> k_inverse;

    std::vector<double>         rhs_values(n_q_points);
    std::vector<double>         boundary_values(n_face_q_points);
    std::vector<Tensor<2, dim>> k_inverse_values(n_q_points);

    // Finally, we need a couple of extractors that we will use to get at the
    // velocity and pressure components of vector-valued shape
    // functions. Their function and use is described in detail in the @ref
    // vector_valued report. Essentially, we will use them as subscripts on
    // the FEValues objects below: the FEValues object describes all vector
    // components of shape functions, while after subscription, it will only
    // refer to the velocities (a set of <code>dim</code> components starting
    // at component zero) or the pressure (a scalar component located at
    // position <code>dim</code>):
    const FEValuesExtractors::Vector velocities(0);
    const FEValuesExtractors::Scalar pressure(dim);

    // With all this in place, we can go on with the loop over all cells. The
    // body of this loop has been discussed in the introduction, and will not
    // be commented any further here:
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        fe_values.reinit(cell);

        local_matrix = 0;
        local_rhs    = 0;

        right_hand_side.value_list(fe_values.get_quadrature_points(),
                                   rhs_values);
        k_inverse.value_list(fe_values.get_quadrature_points(),
                             k_inverse_values);

        for (unsigned int q = 0; q < n_q_points; ++q)
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              const Tensor<1, dim> phi_i_u = fe_values[velocities].value(i, q);
              const double div_phi_i_u = fe_values[velocities].divergence(i, q);
              const double phi_i_p     = fe_values[pressure].value(i, q);

              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                {
                  const Tensor<1, dim> phi_j_u =
                    fe_values[velocities].value(j, q);
                  const double div_phi_j_u =
                    fe_values[velocities].divergence(j, q);
                  const double phi_j_p = fe_values[pressure].value(j, q);

                  local_matrix(i, j) +=
                    (phi_i_u * k_inverse_values[q] * phi_j_u //
                     - phi_i_p * div_phi_j_u                 //
                     - div_phi_i_u * phi_j_p)                //
                    * fe_values.JxW(q);
                }

              local_rhs(i) += -phi_i_p * rhs_values[q] * fe_values.JxW(q);
            }

        for (const auto &face : cell->face_iterators())
          if (face->at_boundary())
            {
              fe_face_values.reinit(cell, face);

              pressure_boundary_values.value_list(
                fe_face_values.get_quadrature_points(), boundary_values);

              for (unsigned int q = 0; q < n_face_q_points; ++q)
                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                  local_rhs(i) += -(fe_face_values[velocities].value(i, q) * //
                                    fe_face_values.normal_vector(q) *        //
                                    boundary_values[q] *                     //
                                    fe_face_values.JxW(q));
            }

        // The final step in the loop over all cells is to transfer local
        // contributions into the global matrix and right hand side
        // vector. Note that we use exactly the same interface as in previous
        // examples, although we now use block matrices and vectors instead of
        // the regular ones. In other words, to the outside world, block
        // objects have the same interface as matrices and vectors, but they
        // additionally allow to access individual blocks.
        cell->get_dof_indices(local_dof_indices);
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
            system_matrix.add(local_dof_indices[i],
                              local_dof_indices[j],
                              local_matrix(i, j));
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          system_rhs(local_dof_indices[i]) += local_rhs(i);
      }
  }


  // @sect3{Implementation of linear solvers and preconditioners}

  // The linear solvers and preconditioners we use in this example have
  // been discussed in significant detail already in the introduction. We
  // will therefore not discuss the rationale for our approach here any
  // more, but rather only comment on some remaining implementational
  // aspects.

  // @sect4{MixedLaplace::solve}

  // As already outlined in the introduction, the solve function consists
  // essentially of two steps. First, we have to form the first equation
  // involving the Schur complement and solve for the pressure (component 1
  // of the solution). Then, we can reconstruct the velocities from the
  // second equation (component 0 of the solution).
  template <int dim>
  void MixedLaplaceProblem<dim>::solve()
  {
    // As a first step we declare references to all block components of the
    // matrix, the right hand side and the solution vector that we will
    // need.
    const auto &M = system_matrix.block(0, 0);
    const auto &B = system_matrix.block(0, 1);

    const auto &F = system_rhs.block(0);
    const auto &G = system_rhs.block(1);

    auto &U = solution.block(0);
    auto &P = solution.block(1);

    // Then, we will create corresponding LinearOperator objects and create
    // the <code>op_M_inv</code> operator:
    const auto op_M = linear_operator(M);
    const auto op_B = linear_operator(B);

    ReductionControl         reduction_control_M(2000, 1.0e-18, 1.0e-10);
    SolverCG<Vector<double>> solver_M(reduction_control_M);
    PreconditionJacobi<SparseMatrix<double>> preconditioner_M;

    preconditioner_M.initialize(M);

    const auto op_M_inv = inverse_operator(op_M, solver_M, preconditioner_M);

    // This allows us to declare the Schur complement <code>op_S</code> and
    // the approximate Schur complement <code>op_aS</code>:
    const auto op_S = transpose_operator(op_B) * op_M_inv * op_B;
    const auto op_aS =
      transpose_operator(op_B) * linear_operator(preconditioner_M) * op_B;

    // We now create a preconditioner out of <code>op_aS</code> that
    // applies a fixed number of 30 (inexpensive) CG iterations:
    IterationNumberControl   iteration_number_control_aS(30, 1.e-18);
    SolverCG<Vector<double>> solver_aS(iteration_number_control_aS);

    const auto preconditioner_S =
      inverse_operator(op_aS, solver_aS, PreconditionIdentity());

    // Now on to the first equation. The right hand side of it is
    // $B^TM^{-1}F-G$, which is what we compute in the first few lines. We
    // then solve the first equation with a CG solver and the
    // preconditioner we just declared.
    const auto schur_rhs = transpose_operator(op_B) * op_M_inv * F - G;

    SolverControl            solver_control_S(2000, 1.e-12);
    SolverCG<Vector<double>> solver_S(solver_control_S);

    const auto op_S_inv = inverse_operator(op_S, solver_S, preconditioner_S);

    P = op_S_inv * schur_rhs;

    std::cout << solver_control_S.last_step()
              << " CG Schur complement iterations to obtain convergence."
              << std::endl;

    // After we have the pressure, we can compute the velocity. The equation
    // reads $MU=-BP+F$, and we solve it by first computing the right hand
    // side, and then multiplying it with the object that represents the
    // inverse of the @ref GlossMassMatrix "mass matrix":
    U = op_M_inv * (F - op_B * P);
  }


  // @sect3{MixedLaplaceProblem class implementation (continued)}

  // @sect4{MixedLaplace::compute_errors}

  // After we have dealt with the linear solver and preconditioners, we
  // continue with the implementation of our main class. In particular, the
  // next task is to compute the errors in our numerical solution, in both the
  // pressures as well as velocities.
  //
  // To compute errors in the solution, we have already introduced the
  // <code>VectorTools::integrate_difference</code> function in step-7 and
  // step-11. However, there we only dealt with scalar solutions, whereas here
  // we have a vector-valued solution with components that even denote
  // different quantities and may have different orders of convergence (this
  // isn't the case here, by choice of the used finite elements, but is
  // frequently the case in mixed finite element applications). What we
  // therefore have to do is to `mask' the components that we are interested
  // in. This is easily done: the
  // <code>VectorTools::integrate_difference</code> function takes as one of its
  // arguments a pointer to a weight function (the parameter defaults to the
  // null pointer, meaning unit weights). What we have to do is to pass
  // a function object that equals one in the components we are interested in,
  // and zero in the other ones. For example, to compute the pressure error,
  // we should pass a function that represents the constant vector with a unit
  // value in component <code>dim</code>, whereas for the velocity the
  // constant vector should be one in the first <code>dim</code> components,
  // and zero in the location of the pressure.
  //
  // In deal.II, the <code>ComponentSelectFunction</code> does exactly this:
  // it wants to know how many vector components the function it is to
  // represent should have (in our case this would be <code>dim+1</code>, for
  // the joint velocity-pressure space) and which individual or range of
  // components should be equal to one. We therefore define two such masks at
  // the beginning of the function, following by an object representing the
  // exact solution and a vector in which we will store the cellwise errors as
  // computed by <code>integrate_difference</code>:
  template <int dim>
  void MixedLaplaceProblem<dim>::compute_errors() const
  {
    const ComponentSelectFunction<dim> pressure_mask(dim, dim + 1);
    const ComponentSelectFunction<dim> velocity_mask(std::make_pair(0, dim),
                                                     dim + 1);

    PrescribedSolution::ExactSolution<dim> exact_solution;
    Vector<double> cellwise_errors(triangulation.n_active_cells());

    // As already discussed in step-7, we have to realize that it is
    // impossible to integrate the errors exactly. All we can do is
    // approximate this integral using quadrature. This actually presents a
    // slight twist here: if we naively chose an object of type
    // <code>QGauss@<dim@>(degree+1)</code> as one may be inclined to do (this
    // is what we used for integrating the linear system), one realizes that
    // the error is very small and does not follow the expected convergence
    // curves at all. What is happening is that for the mixed finite elements
    // used here, the Gauss points happen to be superconvergence points in
    // which the pointwise error is much smaller (and converges with higher
    // order) than anywhere else. These are therefore not particularly good
    // points for integration. To avoid this problem, we simply use a
    // trapezoidal rule and iterate it <code>degree+2</code> times in each
    // coordinate direction (again as explained in step-7):
    const QTrapezoid<1>  q_trapez;
    const QIterated<dim> quadrature(q_trapez, degree + 2);

    // With this, we can then let the library compute the errors and output
    // them to the screen:
    VectorTools::integrate_difference(dof_handler,
                                      solution,
                                      exact_solution,
                                      cellwise_errors,
                                      quadrature,
                                      VectorTools::L2_norm,
                                      &pressure_mask);
    const double p_l2_error =
      VectorTools::compute_global_error(triangulation,
                                        cellwise_errors,
                                        VectorTools::L2_norm);

    VectorTools::integrate_difference(dof_handler,
                                      solution,
                                      exact_solution,
                                      cellwise_errors,
                                      quadrature,
                                      VectorTools::L2_norm,
                                      &velocity_mask);
    const double u_l2_error =
      VectorTools::compute_global_error(triangulation,
                                        cellwise_errors,
                                        VectorTools::L2_norm);

    std::cout << "Errors: ||e_p||_L2 = " << p_l2_error
              << ",   ||e_u||_L2 = " << u_l2_error << std::endl;
  }


  // @sect4{MixedLaplace::output_results}

  // The last interesting function is the one in which we generate graphical
  // output. Note that all velocity components get the same solution name
  // "u". Together with using
  // DataComponentInterpretation::component_is_part_of_vector this will
  // cause DataOut<dim>::write_vtu() to generate a vector representation of
  // the individual velocity components, see step-22 or the
  // @ref VVOutput "Generating graphical output"
  // section of the
  // @ref vector_valued
  // topic for more information. Finally, it seems inappropriate for higher
  // order elements to only show a single bilinear quadrilateral per cell in
  // the graphical output. We therefore generate patches of size
  // (degree+1)x(degree+1) to capture the full information content of the
  // solution. See the step-7 tutorial program for more information on this.
  template <int dim>
  void MixedLaplaceProblem<dim>::output_results() const
  {
    std::vector<std::string> solution_names(dim, "u");
    solution_names.emplace_back("p");
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      interpretation(dim,
                     DataComponentInterpretation::component_is_part_of_vector);
    interpretation.push_back(DataComponentInterpretation::component_is_scalar);

    DataOut<dim> data_out;
    data_out.add_data_vector(dof_handler,
                             solution,
                             solution_names,
                             interpretation);

    data_out.build_patches(degree + 1);

    std::ofstream output("solution.vtu");
    data_out.write_vtu(output);
  }



  // @sect4{MixedLaplace::run}

  // This is the final function of our main class. It's only job is to call
  // the other functions in their natural order:
  template <int dim>
  void MixedLaplaceProblem<dim>::run()
  {
    make_grid_and_dofs();
    assemble_system();
    solve();
    compute_errors();
    output_results();
  }
} // namespace Step20


// @sect3{The <code>main</code> function}

// The main function we stole from step-6 instead of step-4. It is almost
// equal to the one in step-6 (apart from the changed class names, of course),
// the only exception is that we pass the degree of the finite element space
// to the constructor of the mixed Laplace problem (here, we use zero-th order
// elements).
int main()
{
  try
    {
      using namespace Step20;

      const unsigned int     fe_degree = 0;
      MixedLaplaceProblem<2> mixed_laplace_problem(fe_degree);
      mixed_laplace_problem.run();
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
