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
 * Authors: Katharina Kormann, Martin Kronbichler, Uppsala University,
 * 2009-2012, updated to MPI version with parallel vectors in 2016
 */


// First include the necessary files from the deal.II library.
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/multigrid/multigrid.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_matrix.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

// This includes the data structures for the efficient implementation of
// matrix-free methods or more generic finite element operators with the class
// MatrixFree.
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>
#include <deal.II/matrix_free/fe_evaluation.h>

#include <iostream>
#include <fstream>


namespace Step37
{
  using namespace dealii;


  // To be efficient, the operations performed in the matrix-free
  // implementation require knowledge of loop lengths at compile time, which
  // are given by the degree of the finite element. Hence, we collect the
  // values of the two template parameters that can be changed at one place in
  // the code. Of course, one could make the degree of the finite element a
  // run-time parameter by compiling the computational kernels for all degrees
  // that are likely (say, between 1 and 6) and selecting the appropriate
  // kernel at run time. Here, we simply choose second order $Q_2$ elements
  // and choose dimension 3 as standard.
  const unsigned int degree_finite_element = 2;
  const unsigned int dimension             = 3;


  // @sect3{Equation data}

  // We define a variable coefficient function for the Poisson problem. It is
  // similar to the function in step-5 but we use the form $a(\mathbf
  // x)=\frac{1}{0.05 + 2\|\bf x\|^2}$ instead of a discontinuous one. It is
  // merely to demonstrate the possibilities of this implementation, rather
  // than making much sense physically. We define the coefficient in the same
  // way as functions in earlier tutorial programs. There is one new function,
  // namely a @p value method with template argument @p number.
  template <int dim>
  class Coefficient : public Function<dim>
  {
  public:
    virtual double value(const Point<dim>  &p,
                         const unsigned int component = 0) const override;

    template <typename number>
    number value(const Point<dim, number> &p,
                 const unsigned int        component = 0) const;
  };



  // This is the new function mentioned above: Evaluate the coefficient for
  // abstract type @p number. It might be just a usual double, but it can also
  // be a somewhat more complicated type that we call VectorizedArray. This
  // data type is essentially a short array of doubles as discussed in the
  // introduction that holds data from several cells. For example, we evaluate
  // the coefficient shown here not on a simple point as usually done, but we
  // hand it a Point<dim,VectorizedArray<double> > point, which is actually a
  // collection of four points in the case of AVX. Do not confuse the entries
  // in VectorizedArray with the different coordinates of the point. Indeed,
  // the data is laid out such that <code>p[0]</code> returns a
  // VectorizedArray, which in turn contains the x-coordinate for the first
  // point and the second point. You may access the coordinates individually
  // using e.g. <code>p[0][j]</code>, j=0,1,2,3, but it is recommended to
  // define operations on a VectorizedArray as much as possible in order to
  // make use of vectorized operations.
  //
  // In the function implementation, we assume that the number type overloads
  // basic arithmetic operations, so we just write the code as usual. The base
  // class function @p value is then computed from the templated function with
  // double type, in order to avoid duplicating code.
  template <int dim>
  template <typename number>
  number Coefficient<dim>::value(const Point<dim, number> &p,
                                 const unsigned int /*component*/) const
  {
    return 1. / (0.05 + 2. * p.square());
  }



  template <int dim>
  double Coefficient<dim>::value(const Point<dim>  &p,
                                 const unsigned int component) const
  {
    return value<double>(p, component);
  }


  // @sect3{Matrix-free implementation}

  // The following class, called <code>LaplaceOperator</code>, implements the
  // differential operator. For all practical purposes, it is a matrix, i.e.,
  // you can ask it for its size (member functions <code>m(), n()</code>) and
  // you can apply it to a vector (the <code>vmult()</code> function). The
  // difference to a real matrix of course lies in the fact that this class
  // does not actually store the <i>elements</i> of the matrix, but only knows
  // how to compute the action of the operator when applied to a vector.
  //
  // The infrastructure describing the matrix size, the initialization from a
  // MatrixFree object, and the various interfaces to matrix-vector products
  // through vmult() and Tvmult() methods, is provided by the class
  // MatrixFreeOperator::Base from which this class derives. The
  // LaplaceOperator class defined here only has to provide a few interfaces,
  // namely the actual action of the operator through the apply_add() method
  // that gets used in the vmult() functions, and a method to compute the
  // diagonal entries of the underlying matrix. We need the diagonal for the
  // definition of the multigrid smoother. Since we consider a problem with
  // variable coefficient, we further implement a method that can fill the
  // coefficient values.
  //
  // Note that the file <code>include/deal.II/matrix_free/operators.h</code>
  // already contains an implementation of the Laplacian through the class
  // MatrixFreeOperators::LaplaceOperator. For educational purposes, the
  // operator is re-implemented in this tutorial program, explaining the
  // ingredients and concepts used there.
  //
  // This program makes use of the data cache for finite element operator
  // application that is integrated in deal.II. This data cache class is
  // called MatrixFree. It contains mapping information (Jacobians) and index
  // relations between local and global degrees of freedom. It also contains
  // constraints like the ones from hanging nodes or Dirichlet boundary
  // conditions. Moreover, it can issue a loop over all cells in %parallel,
  // making sure that only cells are worked on that do not share any degree of
  // freedom (this makes the loop thread-safe when writing into destination
  // vectors). This is a more advanced strategy compared to the WorkStream
  // class described in the @ref threads topic. Of course, to not destroy
  // thread-safety, we have to be careful when writing into class-global
  // structures.
  //
  // The class implementing the Laplace operator has three template arguments,
  // one for the dimension (as many deal.II classes carry), one for the degree
  // of the finite element (which we need to enable efficient computations
  // through the FEEvaluation class), and one for the underlying scalar
  // type. We want to use <code>double</code> numbers (i.e., double precision,
  // 64-bit floating point) for the final matrix, but floats (single
  // precision, 32-bit floating point numbers) for the multigrid level
  // matrices (as that is only a preconditioner, and floats can be processed
  // twice as fast). The class FEEvaluation also takes a template argument for
  // the number of quadrature points in one dimension. In the code below, we
  // hard-code it to <code>fe_degree+1</code>. If we wanted to change it
  // independently of the polynomial degree, we would need to add a template
  // parameter as is done in the MatrixFreeOperators::LaplaceOperator class.
  //
  // As a sidenote, if we implemented several different operations on the same
  // grid and degrees of freedom (like a @ref GlossMassMatrix "mass matrix" and a Laplace matrix), we
  // would define two classes like the current one for each of the operators
  // (derived from the MatrixFreeOperators::Base class), and let both of them
  // refer to the same MatrixFree data cache from the general problem
  // class. The interface through MatrixFreeOperators::Base requires us to
  // only provide a minimal set of functions. This concept allows for writing
  // complex application codes with many matrix-free operations.
  //
  // @note Storing values of type <code>VectorizedArray<number></code>
  // requires care: Here, we use the deal.II table class which is prepared to
  // hold the data with correct alignment. However, storing e.g. an
  // <code>std::vector<VectorizedArray<number> ></code> is not possible with
  // vectorization: A certain alignment of the data with the memory address
  // boundaries is required (essentially, a VectorizedArray that is 32 bytes
  // long in case of AVX needs to start at a memory address that is divisible
  // by 32). The table class (as well as the AlignedVector class it is based
  // on) makes sure that this alignment is respected, whereas std::vector does
  // not in general, which may lead to segmentation faults at strange places
  // for some systems or suboptimal performance for other systems.
  template <int dim, int fe_degree, typename number>
  class LaplaceOperator
    : public MatrixFreeOperators::
        Base<dim, LinearAlgebra::distributed::Vector<number>>
  {
  public:
    using value_type = number;

    LaplaceOperator();

    void clear() override;

    void evaluate_coefficient(const Coefficient<dim> &coefficient_function);

    virtual void compute_diagonal() override;

  private:
    virtual void apply_add(
      LinearAlgebra::distributed::Vector<number>       &dst,
      const LinearAlgebra::distributed::Vector<number> &src) const override;

    void
    local_apply(const MatrixFree<dim, number>                    &data,
                LinearAlgebra::distributed::Vector<number>       &dst,
                const LinearAlgebra::distributed::Vector<number> &src,
                const std::pair<unsigned int, unsigned int> &cell_range) const;

    void local_compute_diagonal(
      FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number> &integrator) const;

    Table<2, VectorizedArray<number>> coefficient;
  };



  // This is the constructor of the @p LaplaceOperator class. All it does is
  // to call the default constructor of the base class
  // MatrixFreeOperators::Base, which in turn is based on the
  // EnableObserverPointer class that asserts that this class is
  // not accessed after going out of scope e.g. in a preconditioner.
  template <int dim, int fe_degree, typename number>
  LaplaceOperator<dim, fe_degree, number>::LaplaceOperator()
    : MatrixFreeOperators::Base<dim,
                                LinearAlgebra::distributed::Vector<number>>()
  {}



  template <int dim, int fe_degree, typename number>
  void LaplaceOperator<dim, fe_degree, number>::clear()
  {
    coefficient.reinit(0, 0);
    MatrixFreeOperators::Base<dim, LinearAlgebra::distributed::Vector<number>>::
      clear();
  }



  // @sect4{Computation of coefficient}

  // To initialize the coefficient, we directly give it the Coefficient class
  // defined above and then select the method
  // <code>coefficient_function.value</code> with vectorized number (which the
  // compiler can deduce from the point data type). The use of the
  // FEEvaluation class (and its template arguments) will be explained below.
  template <int dim, int fe_degree, typename number>
  void LaplaceOperator<dim, fe_degree, number>::evaluate_coefficient(
    const Coefficient<dim> &coefficient_function)
  {
    const unsigned int n_cells = this->data->n_cell_batches();
    FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number> phi(*this->data);

    coefficient.reinit(n_cells, phi.n_q_points);
    for (unsigned int cell = 0; cell < n_cells; ++cell)
      {
        phi.reinit(cell);
        for (const unsigned int q : phi.quadrature_point_indices())
          coefficient(cell, q) =
            coefficient_function.value(phi.quadrature_point(q));
      }
  }



  // @sect4{Local evaluation of Laplace operator}

  // Here comes the main function of this class, the evaluation of the
  // matrix-vector product (or, in general, a finite element operator
  // evaluation). This is done in a function that takes exactly four
  // arguments, the MatrixFree object, the destination and source vectors, and
  // a range of cells that are to be worked on. The method
  // <code>cell_loop</code> in the MatrixFree class will internally call this
  // function with some range of cells that is obtained by checking which
  // cells are possible to work on simultaneously so that write operations do
  // not cause any race condition. Note that the cell range used in the loop
  // is not directly the number of (active) cells in the current mesh, but
  // rather a collection of batches of cells.  In other word, "cell" may be
  // the wrong term to begin with, since FEEvaluation groups data from several
  // cells together. This means that in the loop over quadrature points we are
  // actually seeing a group of quadrature points of several cells as one
  // block. This is done to enable a higher degree of vectorization.  The
  // number of such "cells" or "cell batches" is stored in MatrixFree and can
  // be queried through MatrixFree::n_cell_batches(). Compared to the deal.II
  // cell iterators, in this class all cells are laid out in a plain array
  // with no direct knowledge of level or neighborship relations, which makes
  // it possible to index the cells by unsigned integers.
  //
  // The implementation of the Laplace operator is quite simple: First, we
  // need to create an object FEEvaluation that contains the computational
  // kernels and has data fields to store temporary results (e.g. gradients
  // evaluated on all quadrature points on a collection of a few cells). Note
  // that temporary results do not use a lot of memory, and since we specify
  // template arguments with the element order, the data is stored on the
  // stack (without expensive memory allocation). Usually, one only needs to
  // set two template arguments, the dimension as a first argument and the
  // degree of the finite element as the second argument (this is equal to the
  // number of degrees of freedom per dimension minus one for FE_Q
  // elements). However, here we also want to be able to use float numbers for
  // the multigrid preconditioner, which is the last (fifth) template
  // argument. Therefore, we cannot rely on the default template arguments and
  // must also fill the third and fourth field, consequently. The third
  // argument specifies the number of quadrature points per direction and has
  // a default value equal to the degree of the element plus one. The fourth
  // argument sets the number of components (one can also evaluate
  // vector-valued functions in systems of PDEs, but the default is a scalar
  // element), and finally the last argument sets the number type.
  //
  // Next, we loop over the given cell range and then we continue with the
  // actual implementation: <ol> <li>Tell the FEEvaluation object the (macro)
  // cell we want to work on.  <li>Read in the values of the source vectors
  // (@p read_dof_values), including the resolution of constraints. This
  // stores $u_\mathrm{cell}$ as described in the introduction.  <li>Compute
  // the unit-cell gradient (the evaluation of finite element
  // functions). Since FEEvaluation can combine value computations with
  // gradient computations, it uses a unified interface to all kinds of
  // derivatives of order between zero and two. We only want gradients, no
  // values and no second derivatives, so we set the function arguments to
  // true in the gradient slot (second slot), and to false in the values slot
  // (first slot). There is also a third slot for the Hessian which is
  // false by default, so it needs not be given. Note that the FEEvaluation
  // class internally evaluates shape functions in an efficient way where one
  // dimension is worked on at a time (using the tensor product form of shape
  // functions and quadrature points as mentioned in the introduction). This
  // gives complexity equal to $\mathcal O(d^2 (p+1)^{d+1})$ for polynomial
  // degree $p$ in $d$ dimensions, compared to the naive approach with loops
  // over all local degrees of freedom and quadrature points that is used in
  // FEValues and costs $\mathcal O(d (p+1)^{2d})$.  <li>Next comes the
  // application of the Jacobian transformation, the multiplication by the
  // variable coefficient and the quadrature weight. FEEvaluation has an
  // access function @p get_gradient that applies the Jacobian and returns the
  // gradient in real space. Then, we just need to multiply by the (scalar)
  // coefficient, and let the function @p submit_gradient apply the second
  // Jacobian (for the test function) and the quadrature weight and Jacobian
  // determinant (JxW). Note that the submitted gradient is stored in the same
  // data field as where it is read from in @p get_gradient. Therefore, you
  // need to make sure to not read from the same quadrature point again after
  // having called @p submit_gradient on that particular quadrature point. In
  // general, it is a good idea to copy the result of @p get_gradient when it
  // is used more often than once.  <li>Next follows the summation over
  // quadrature points for all test functions that corresponds to the actual
  // integration step. For the Laplace operator, we just multiply by the
  // gradient, so we call the integrate function with the respective argument
  // set. If you have an equation where you test by both the values of the
  // test functions and the gradients, both template arguments need to be set
  // to true. Calling first the integrate function for values and then
  // gradients in a separate call leads to wrong results, since the second
  // call will internally overwrite the results from the first call. Note that
  // there is no function argument for the second derivative for integrate
  // step.  <li>Eventually, the local contributions in the vector
  // $v_\mathrm{cell}$ as mentioned in the introduction need to be added into
  // the result vector (and constraints are applied). This is done with a call
  // to @p distribute_local_to_global, the same name as the corresponding
  // function in the AffineConstraints (only that we now store the local vector
  // in the FEEvaluation object, as are the indices between local and global
  // degrees of freedom).  </ol>
  template <int dim, int fe_degree, typename number>
  void LaplaceOperator<dim, fe_degree, number>::local_apply(
    const MatrixFree<dim, number>                    &data,
    LinearAlgebra::distributed::Vector<number>       &dst,
    const LinearAlgebra::distributed::Vector<number> &src,
    const std::pair<unsigned int, unsigned int>      &cell_range) const
  {
    FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number> phi(data);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        AssertDimension(coefficient.size(0), data.n_cell_batches());
        AssertDimension(coefficient.size(1), phi.n_q_points);

        phi.reinit(cell);
        phi.read_dof_values(src);
        phi.evaluate(EvaluationFlags::gradients);
        for (const unsigned int q : phi.quadrature_point_indices())
          phi.submit_gradient(coefficient(cell, q) * phi.get_gradient(q), q);
        phi.integrate(EvaluationFlags::gradients);
        phi.distribute_local_to_global(dst);
      }
  }



  // This function implements the loop over all cells for the
  // Base::apply_add() interface. This is done with the @p cell_loop of the
  // MatrixFree class, which takes the operator() of this class with arguments
  // MatrixFree, OutVector, InVector, cell_range. When working with MPI
  // parallelization (but no threading) as is done in this tutorial program,
  // the cell loop corresponds to the following three lines of code:
  //
  // @code
  // src.update_ghost_values();
  // local_apply(*this->data, dst, src, std::make_pair(0U,
  //                                                   data.n_cell_batches()));
  // dst.compress(VectorOperation::add);
  // @endcode
  //
  // Here, the two calls update_ghost_values() and compress() perform the data
  // exchange on processor boundaries for MPI, once for the source vector
  // where we need to read from entries owned by remote processors, and once
  // for the destination vector where we have accumulated parts of the
  // residuals that need to be added to the respective entry of the owner
  // processor. However, MatrixFree::cell_loop does not only abstract away
  // those two calls, but also performs some additional optimizations. On the
  // one hand, it will split the update_ghost_values() and compress() calls in
  // a way to allow for overlapping communication and computation. The
  // local_apply function is then called with three cell ranges representing
  // partitions of the cell range from 0 to MatrixFree::n_cell_batches(). On
  // the other hand, cell_loop also supports thread parallelism in which case
  // the cell ranges are split into smaller chunks and scheduled in an
  // advanced way that avoids access to the same vector entry by several
  // threads. That feature is explained in step-48.
  //
  // Note that after the cell loop, the constrained degrees of freedom need to
  // be touched once more for sensible vmult() operators: Since the assembly
  // loop automatically resolves constraints (just as the
  // AffineConstraints::distribute_local_to_global() call does), it does not
  // compute any contribution for constrained degrees of freedom, leaving the
  // respective entries zero. This would represent a matrix that had empty
  // rows and columns for constrained degrees of freedom. However, iterative
  // solvers like CG only work for non-singular matrices. The easiest way to
  // do that is to set the sub-block of the matrix that corresponds to
  // constrained DoFs to an identity matrix, in which case application of the
  // matrix would simply copy the elements of the right hand side vector into
  // the left hand side. Fortunately, the vmult() implementations
  // MatrixFreeOperators::Base do this automatically for us outside the
  // apply_add() function, so we do not need to take further action here.
  //
  // When using the combination of MatrixFree and FEEvaluation in parallel
  // with MPI, there is one aspect to be careful about &mdash; the indexing
  // used for accessing the vector. For performance reasons, MatrixFree and
  // FEEvaluation are designed to access vectors in MPI-local index space also
  // when working with multiple processors. Working in local index space means
  // that no index translation needs to be performed at the place the vector
  // access happens, apart from the unavoidable indirect addressing. However,
  // local index spaces are ambiguous: While it is standard convention to
  // access the locally owned range of a vector with indices between 0 and the
  // local size, the numbering is not so clear for the ghosted entries and
  // somewhat arbitrary. For the matrix-vector product, only the indices
  // appearing on locally owned cells (plus those referenced via hanging node
  // constraints) are necessary. However, in deal.II we often set all the
  // degrees of freedom on ghosted elements as ghosted vector entries, called
  // the
  // @ref GlossLocallyRelevantDof "locally relevant DoFs described in the glossary".
  // In that case, the MPI-local index of a ghosted vector entry can in
  // general be different in the two possible ghost sets, despite referring
  // to the same global index. To avoid problems, FEEvaluation checks that
  // the partitioning of the vector used for the matrix-vector product does
  // indeed match with the partitioning of the indices in MatrixFree by a
  // check called
  // LinearAlgebra::distributed::Vector::partitioners_are_compatible. To
  // facilitate things, the MatrixFreeOperators::Base class includes a
  // mechanism to fit the ghost set to the correct layout. This happens in the
  // ghost region of the vector, so keep in mind that the ghost region might
  // be modified in both the destination and source vector after a call to a
  // vmult() method. This is legitimate because the ghost region of a
  // distributed deal.II vector is a mutable section and filled on
  // demand. Vectors used in matrix-vector products must not be ghosted upon
  // entry of vmult() functions, so no information gets lost.
  template <int dim, int fe_degree, typename number>
  void LaplaceOperator<dim, fe_degree, number>::apply_add(
    LinearAlgebra::distributed::Vector<number>       &dst,
    const LinearAlgebra::distributed::Vector<number> &src) const
  {
    this->data->cell_loop(&LaplaceOperator::local_apply, this, dst, src);
  }



  // The following function implements the computation of the diagonal of the
  // operator. Computing matrix entries of a matrix-free operator evaluation
  // turns out to be more complicated than evaluating the
  // operator. Fundamentally, we could obtain a matrix representation of the
  // operator by applying the operator on <i>all</i> unit vectors. Of course,
  // that would be very inefficient since we would need to perform <i>n</i>
  // operator evaluations to retrieve the whole matrix. Furthermore, this
  // approach would completely ignore the matrix sparsity. On an individual
  // cell, however, this is the way to go and actually not that inefficient as
  // there usually is a coupling between all degrees of freedom inside the
  // cell.
  //
  // We first initialize the diagonal vector to the correct parallel
  // layout. This vector is encapsulated in a member called
  // inverse_diagonal_entries of type DiagonalMatrix in the base class
  // MatrixFreeOperators::Base. This member is a shared pointer that we first
  // need to initialize and then get the vector representing the diagonal
  // entries in the matrix. As to the actual diagonal computation, we could
  // manually write a cell_loop and invoke a local worker that applies all unit
  // vectors on each cell. Instead, we use MatrixFreeTools::compute_diagonal()
  // to do this for us. Afterwards, we need to set the vector entries
  // subject to Dirichlet boundary conditions to one (either those on the
  // boundary described by the AffineConstraints object inside MatrixFree or
  // the indices at the interface between different grid levels in adaptive
  // multigrid). This is done through the function
  // MatrixFreeOperators::Base::set_constrained_entries_to_one() and matches
  // with the setting in the matrix-vector product provided by the Base
  // operator. Finally, we need to invert the diagonal entries which is the
  // form required by the Chebyshev smoother based on the Jacobi iteration. In
  // the loop, we assert that all entries are non-zero, because they should
  // either have obtained a positive contribution from integrals or be
  // constrained and treated by @p set_constrained_entries_to_one().
  template <int dim, int fe_degree, typename number>
  void LaplaceOperator<dim, fe_degree, number>::compute_diagonal()
  {
    this->inverse_diagonal_entries.reset(
      new DiagonalMatrix<LinearAlgebra::distributed::Vector<number>>());
    LinearAlgebra::distributed::Vector<number> &inverse_diagonal =
      this->inverse_diagonal_entries->get_vector();
    this->data->initialize_dof_vector(inverse_diagonal);

    MatrixFreeTools::compute_diagonal(*this->data,
                                      inverse_diagonal,
                                      &LaplaceOperator::local_compute_diagonal,
                                      this);

    this->set_constrained_entries_to_one(inverse_diagonal);

    for (unsigned int i = 0; i < inverse_diagonal.locally_owned_size(); ++i)
      {
        Assert(inverse_diagonal.local_element(i) > 0.,
               ExcMessage("No diagonal entry in a positive definite operator "
                          "should be zero"));
        inverse_diagonal.local_element(i) =
          1. / inverse_diagonal.local_element(i);
      }
  }



  // In the local compute loop, we compute the diagonal by a loop over all
  // columns in the local matrix and putting the entry 1 in the <i>i</i>th
  // slot and a zero entry in all other slots, i.e., we apply the cell-wise
  // differential operator on one unit vector at a time. The inner part
  // invoking FEEvaluation::evaluate(), the loop over quadrature points, and
  // FEEvaluation::integrate(), is exactly the same as in the local_apply
  // function. Afterwards, we pick out the <i>i</i>th entry of the local
  // result and put it to a temporary storage (as we overwrite all entries in
  // the array behind FEEvaluation::get_dof_value() with the next loop
  // iteration). Finally, the temporary storage is written to the destination
  // vector. Note how we use FEEvaluation::get_dof_value() and
  // FEEvaluation::submit_dof_value() to read and write to the data field that
  // FEEvaluation uses for the integration on the one hand and writes into the
  // global vector on the other hand.
  //
  // Given that we are only interested in the matrix diagonal, we simply throw
  // away all other entries of the local matrix that have been computed along
  // the way. While it might seem wasteful to compute the complete cell matrix
  // and then throw away everything but the diagonal, the integration are so
  // efficient that the computation does not take too much time. Note that the
  // complexity of operator evaluation per element is $\mathcal
  // O((p+1)^{d+1})$ for polynomial degree $k$, so computing the whole matrix
  // costs us $\mathcal O((p+1)^{2d+1})$ operations, not too far away from
  // $\mathcal O((p+1)^{2d})$ complexity for computing the diagonal with
  // FEValues. Since FEEvaluation is also considerably faster due to
  // vectorization and other optimizations, the diagonal computation with this
  // function is actually the fastest (simple) variant. (It would be possible
  // to compute the diagonal with sum factorization techniques in $\mathcal
  // O((p+1)^{d+1})$ operations involving specifically adapted
  // kernels&mdash;but since such kernels are only useful in that particular
  // context and the diagonal computation is typically not on the critical
  // path, they have not been implemented in deal.II.)
  //
  // Note that the code that calls distribute_local_to_global on the vector to
  // accumulate the diagonal entries into the global matrix has some
  // limitations. For operators with hanging node constraints that distribute
  // an integral contribution of a constrained DoF to several other entries
  // inside the distribute_local_to_global call, the vector interface used
  // here does not exactly compute the diagonal entries, but lumps some
  // contributions located on the diagonal of the local matrix that would end
  // up in a off-diagonal position of the global matrix to the diagonal. The
  // result is correct up to discretization accuracy as explained in <a
  // href="http://dx.doi.org/10.4208/cicp.101214.021015a">Kormann (2016),
  // section 5.3</a>, but not mathematically equal. In this tutorial program,
  // no harm can happen because the diagonal is only used for the multigrid
  // level matrices where no hanging node constraints appear.
  template <int dim, int fe_degree, typename number>
  void LaplaceOperator<dim, fe_degree, number>::local_compute_diagonal(
    FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number> &phi) const
  {
    const unsigned int cell = phi.get_current_cell_index();

    phi.evaluate(EvaluationFlags::gradients);

    for (const unsigned int q : phi.quadrature_point_indices())
      {
        phi.submit_gradient(coefficient(cell, q) * phi.get_gradient(q), q);
      }
    phi.integrate(EvaluationFlags::gradients);
  }



  // @sect3{LaplaceProblem class}

  // This class is based on the one in step-16. However, we replaced the
  // SparseMatrix<double> class by our matrix-free implementation, which means
  // that we can also skip the sparsity patterns. Notice that we define the
  // LaplaceOperator class with the degree of finite element as template
  // argument (the value is defined at the top of the file), and that we use
  // float numbers for the multigrid level matrices.
  //
  // The class also has a member variable to keep track of all the detailed
  // timings for setting up the entire chain of data before we actually go
  // about solving the problem. In addition, there is an output stream (that
  // is disabled by default) that can be used to output details for the
  // individual setup operations instead of the summary only that is printed
  // out by default.
  //
  // Since this program is designed to be used with MPI, we also provide the
  // usual @p pcout output stream that only prints the information of the
  // processor with MPI rank 0. The grid used for this programs can either be
  // a distributed triangulation based on p4est (in case deal.II is configured
  // to use p4est), otherwise it is a serial grid that only runs without MPI.
  template <int dim>
  class LaplaceProblem
  {
  public:
    LaplaceProblem();
    void run();

  private:
    void setup_system();
    void assemble_rhs();
    void solve();
    void output_results(const unsigned int cycle) const;

#ifdef DEAL_II_WITH_P4EST
    parallel::distributed::Triangulation<dim> triangulation;
#else
    Triangulation<dim> triangulation;
#endif

    const FE_Q<dim> fe;
    DoFHandler<dim> dof_handler;

    const MappingQ1<dim> mapping;

    AffineConstraints<double> constraints;
    using SystemMatrixType =
      LaplaceOperator<dim, degree_finite_element, double>;
    SystemMatrixType system_matrix;

    MGConstrainedDoFs mg_constrained_dofs;
    using LevelMatrixType = LaplaceOperator<dim, degree_finite_element, float>;
    MGLevelObject<LevelMatrixType> mg_matrices;

    LinearAlgebra::distributed::Vector<double> solution;
    LinearAlgebra::distributed::Vector<double> system_rhs;

    double             setup_time;
    ConditionalOStream pcout;
    ConditionalOStream time_details;
  };



  // When we initialize the finite element, we of course have to use the
  // degree specified at the top of the file as well (otherwise, an exception
  // will be thrown at some point, since the computational kernel defined in
  // the templated LaplaceOperator class and the information from the finite
  // element read out by MatrixFree will not match). The constructor of the
  // triangulation needs to set an additional flag that tells the grid to
  // conform to the 2:1 cell balance over vertices, which is needed for the
  // convergence of the geometric multigrid routines. For the distributed
  // grid, we also need to specifically enable the multigrid hierarchy.
  template <int dim>
  LaplaceProblem<dim>::LaplaceProblem()
#ifdef DEAL_II_WITH_P4EST
    : triangulation(MPI_COMM_WORLD,
                    Triangulation<dim>::limit_level_difference_at_vertices,
                    parallel::distributed::Triangulation<
                      dim>::construct_multigrid_hierarchy)
#else
    : triangulation(Triangulation<dim>::limit_level_difference_at_vertices)
#endif
    , fe(degree_finite_element)
    , dof_handler(triangulation)
    , setup_time(0.)
    , pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    // The LaplaceProblem class holds an additional output stream that
    // collects detailed timings about the setup phase. This stream, called
    // time_details, is disabled by default through the @p false argument
    // specified here. For detailed timings, removing the @p false argument
    // prints all the details.
    , time_details(std::cout,
                   false &&
                     Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
  {}



  // @sect4{LaplaceProblem::setup_system}

  // The setup stage is in analogy to step-16 with relevant changes due to the
  // LaplaceOperator class. The first thing to do is to set up the DoFHandler,
  // including the degrees of freedom for the multigrid levels, and to
  // initialize constraints from hanging nodes and homogeneous Dirichlet
  // conditions. Since we intend to use this programs in %parallel with MPI,
  // we need to make sure that the constraints get to know the locally
  // relevant degrees of freedom, otherwise the storage would explode when
  // using more than a few hundred millions of degrees of freedom, see
  // step-40.

  // Once we have created the multigrid dof_handler and the constraints, we
  // can call the reinit function for the global matrix operator as well as
  // each level of the multigrid scheme. The main action is to set up the
  // <code> MatrixFree </code> instance for the problem. The base class of the
  // <code>LaplaceOperator</code> class, MatrixFreeOperators::Base, is
  // initialized with a shared pointer to MatrixFree object. This way, we can
  // simply create it here and then pass it on to the system matrix and level
  // matrices, respectively. For setting up MatrixFree, we need to activate
  // the update flag in the AdditionalData field of MatrixFree that enables
  // the storage of quadrature point coordinates in real space (by default, it
  // only caches data for gradients (inverse transposed Jacobians) and JxW
  // values). Note that if we call the reinit function without specifying the
  // level (i.e., giving <code>level = numbers::invalid_unsigned_int</code>),
  // MatrixFree constructs a loop over the active cells. In this tutorial, we
  // do not use threads in addition to MPI, which is why we explicitly disable
  // it by setting the MatrixFree::AdditionalData::tasks_parallel_scheme to
  // MatrixFree::AdditionalData::none. Finally, the coefficient is evaluated
  // and vectors are initialized as explained above.
  template <int dim>
  void LaplaceProblem<dim>::setup_system()
  {
    Timer time;
    setup_time = 0;
    {
      system_matrix.clear();
      mg_matrices.clear_elements();

      dof_handler.distribute_dofs(fe);
      dof_handler.distribute_mg_dofs();

      pcout << "Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl;

      constraints.clear();
      constraints.reinit(dof_handler.locally_owned_dofs(),
                         DoFTools::extract_locally_relevant_dofs(dof_handler));
      DoFTools::make_hanging_node_constraints(dof_handler, constraints);
      VectorTools::interpolate_boundary_values(
        mapping, dof_handler, 0, Functions::ZeroFunction<dim>(), constraints);
      constraints.close();
    }
    setup_time += time.wall_time();
    time_details << "Distribute DoFs & B.C.     (CPU/wall) " << time.cpu_time()
                 << "s/" << time.wall_time() << 's' << std::endl;
    time.restart();
    {
      {
        typename MatrixFree<dim, double>::AdditionalData additional_data;
        additional_data.tasks_parallel_scheme =
          MatrixFree<dim, double>::AdditionalData::none;
        additional_data.mapping_update_flags =
          (update_gradients | update_JxW_values | update_quadrature_points);
        std::shared_ptr<MatrixFree<dim, double>> system_mf_storage(
          new MatrixFree<dim, double>());
        system_mf_storage->reinit(mapping,
                                  dof_handler,
                                  constraints,
                                  QGauss<1>(fe.degree + 1),
                                  additional_data);
        system_matrix.initialize(system_mf_storage);
      }

      system_matrix.evaluate_coefficient(Coefficient<dim>());

      system_matrix.initialize_dof_vector(solution);
      system_matrix.initialize_dof_vector(system_rhs);
    }
    setup_time += time.wall_time();
    time_details << "Setup matrix-free system   (CPU/wall) " << time.cpu_time()
                 << "s/" << time.wall_time() << 's' << std::endl;
    time.restart();

    // Next, initialize the matrices for the multigrid method on all the
    // levels. The data structure MGConstrainedDoFs keeps information about
    // the indices subject to boundary conditions as well as the indices on
    // edges between different refinement levels as described in the step-16
    // tutorial program. We then go through the levels of the mesh and
    // construct the constraints and matrices on each level. These follow
    // closely the construction of the system matrix on the original mesh,
    // except the slight difference in naming when accessing information on
    // the levels rather than the active cells.
    {
      const unsigned int nlevels = triangulation.n_global_levels();
      mg_matrices.resize(0, nlevels - 1);

      const std::set<types::boundary_id> dirichlet_boundary_ids = {0};
      mg_constrained_dofs.initialize(dof_handler);
      mg_constrained_dofs.make_zero_boundary_constraints(
        dof_handler, dirichlet_boundary_ids);

      for (unsigned int level = 0; level < nlevels; ++level)
        {
          AffineConstraints<double> level_constraints(
            dof_handler.locally_owned_mg_dofs(level),
            DoFTools::extract_locally_relevant_level_dofs(dof_handler, level));
          for (const types::global_dof_index dof_index :
               mg_constrained_dofs.get_boundary_indices(level))
            level_constraints.constrain_dof_to_zero(dof_index);
          level_constraints.close();

          typename MatrixFree<dim, float>::AdditionalData additional_data;
          additional_data.tasks_parallel_scheme =
            MatrixFree<dim, float>::AdditionalData::none;
          additional_data.mapping_update_flags =
            (update_gradients | update_JxW_values | update_quadrature_points);
          additional_data.mg_level = level;
          std::shared_ptr<MatrixFree<dim, float>> mg_mf_storage_level =
            std::make_shared<MatrixFree<dim, float>>();
          mg_mf_storage_level->reinit(mapping,
                                      dof_handler,
                                      level_constraints,
                                      QGauss<1>(fe.degree + 1),
                                      additional_data);

          mg_matrices[level].initialize(mg_mf_storage_level,
                                        mg_constrained_dofs,
                                        level);
          mg_matrices[level].evaluate_coefficient(Coefficient<dim>());
        }
    }
    setup_time += time.wall_time();
    time_details << "Setup matrix-free levels   (CPU/wall) " << time.cpu_time()
                 << "s/" << time.wall_time() << 's' << std::endl;
  }



  // @sect4{LaplaceProblem::assemble_rhs}

  // The assemble function is very simple since all we have to do is to
  // assemble the right hand side. Thanks to FEEvaluation and all the data
  // cached in the MatrixFree class, which we query from
  // MatrixFreeOperators::Base, this can be done in a few lines. Since this
  // call is not wrapped into a MatrixFree::cell_loop (which would be an
  // alternative), we must not forget to call compress() at the end of the
  // assembly to send all the contributions of the right hand side to the
  // owner of the respective degree of freedom.
  template <int dim>
  void LaplaceProblem<dim>::assemble_rhs()
  {
    Timer time;

    system_rhs = 0;
    FEEvaluation<dim, degree_finite_element> phi(
      *system_matrix.get_matrix_free());
    for (unsigned int cell = 0;
         cell < system_matrix.get_matrix_free()->n_cell_batches();
         ++cell)
      {
        phi.reinit(cell);
        for (const unsigned int q : phi.quadrature_point_indices())
          phi.submit_value(make_vectorized_array<double>(1.0), q);
        phi.integrate(EvaluationFlags::values);
        phi.distribute_local_to_global(system_rhs);
      }
    system_rhs.compress(VectorOperation::add);

    setup_time += time.wall_time();
    time_details << "Assemble right hand side   (CPU/wall) " << time.cpu_time()
                 << "s/" << time.wall_time() << 's' << std::endl;
  }



  // @sect4{LaplaceProblem::solve}

  // The solution process is similar as in step-16. We start with the setup of
  // the transfer. For LinearAlgebra::distributed::Vector, there is a very
  // fast transfer class called MGTransferMatrixFree that does the
  // interpolation between the grid levels with the same fast sum
  // factorization kernels that get also used in FEEvaluation.
  template <int dim>
  void LaplaceProblem<dim>::solve()
  {
    Timer                            time;
    MGTransferMatrixFree<dim, float> mg_transfer(mg_constrained_dofs);
    mg_transfer.build(dof_handler);
    setup_time += time.wall_time();
    time_details << "MG build transfer time     (CPU/wall) " << time.cpu_time()
                 << "s/" << time.wall_time() << "s\n";
    time.restart();

    // As a smoother, this tutorial program uses a Chebyshev iteration instead
    // of SOR in step-16. (SOR would be very difficult to implement because we
    // do not have the matrix elements available explicitly, and it is
    // difficult to make it work efficiently in %parallel.)  The smoother is
    // initialized with our level matrices and the mandatory additional data
    // for the Chebyshev smoother. We use a relatively high degree here (5),
    // since matrix-vector products are comparably cheap. We choose to smooth
    // out a range of $[1.2 \hat{\lambda}_{\max}/15,1.2 \hat{\lambda}_{\max}]$
    // in the smoother where $\hat{\lambda}_{\max}$ is an estimate of the
    // largest eigenvalue (the factor 1.2 is applied inside
    // PreconditionChebyshev). In order to compute that eigenvalue, the
    // Chebyshev initialization performs a few steps of a CG algorithm
    // without preconditioner. Since the highest eigenvalue is usually the
    // easiest one to find and a rough estimate is enough, we choose 10
    // iterations. Finally, we also set the inner preconditioner type in the
    // Chebyshev method which is a Jacobi iteration. This is represented by
    // the DiagonalMatrix class that gets the inverse diagonal entry provided
    // by our LaplaceOperator class.
    //
    // On level zero, we initialize the smoother differently because we want
    // to use the Chebyshev iteration as a solver. PreconditionChebyshev
    // allows the user to switch to solver mode where the number of iterations
    // is internally chosen to the correct value. In the additional data
    // object, this setting is activated by choosing the polynomial degree to
    // @p numbers::invalid_unsigned_int. The algorithm will then attack all
    // eigenvalues between the smallest and largest one in the coarse level
    // matrix. The number of steps in the Chebyshev smoother are chosen such
    // that the Chebyshev convergence estimates guarantee to reduce the
    // residual by the number specified in the variable @p
    // smoothing_range. Note that for solving, @p smoothing_range is a
    // relative tolerance and chosen smaller than one, in this case, we select
    // three orders of magnitude, whereas it is a number larger than 1 when
    // only selected eigenvalues are smoothed.
    //
    // From a computational point of view, the Chebyshev iteration is a very
    // attractive coarse grid solver as long as the coarse size is
    // moderate. This is because the Chebyshev method performs only
    // matrix-vector products and vector updates, which typically parallelize
    // better to the largest cluster size with more than a few tens of
    // thousands of cores than inner product involved in other iterative
    // methods. The former involves only local communication between neighbors
    // in the (coarse) mesh, whereas the latter requires global communication
    // over all processors.
    using SmootherType =
      PreconditionChebyshev<LevelMatrixType,
                            LinearAlgebra::distributed::Vector<float>>;
    mg::SmootherRelaxation<SmootherType,
                           LinearAlgebra::distributed::Vector<float>>
                                                         mg_smoother;
    MGLevelObject<typename SmootherType::AdditionalData> smoother_data;
    smoother_data.resize(0, triangulation.n_global_levels() - 1);
    for (unsigned int level = 0; level < triangulation.n_global_levels();
         ++level)
      {
        if (level > 0)
          {
            smoother_data[level].smoothing_range     = 15.;
            smoother_data[level].degree              = 5;
            smoother_data[level].eig_cg_n_iterations = 10;
          }
        else
          {
            smoother_data[0].smoothing_range = 1e-3;
            smoother_data[0].degree          = numbers::invalid_unsigned_int;
            smoother_data[0].eig_cg_n_iterations = mg_matrices[0].m();
          }
        mg_matrices[level].compute_diagonal();
        smoother_data[level].preconditioner =
          mg_matrices[level].get_matrix_diagonal_inverse();
      }
    mg_smoother.initialize(mg_matrices, smoother_data);

    MGCoarseGridApplySmoother<LinearAlgebra::distributed::Vector<float>>
      mg_coarse;
    mg_coarse.initialize(mg_smoother);

    // The next step is to set up the interface matrices that are needed for the
    // case with hanging nodes. The adaptive multigrid realization in deal.II
    // implements an approach called local smoothing. This means that the
    // smoothing on the finest level only covers the local part of the mesh
    // defined by the fixed (finest) grid level and ignores parts of the
    // computational domain where the terminal cells are coarser than this
    // level. As the method progresses to coarser levels, more and more of the
    // global mesh will be covered. At some coarser level, the whole mesh will
    // be covered. Since all level matrices in the multigrid method cover a
    // single level in the mesh, no hanging nodes appear on the level matrices.
    // At the interface between multigrid levels, homogeneous Dirichlet boundary
    // conditions are set while smoothing. When the residual is transferred to
    // the next coarser level, however, the coupling over the multigrid
    // interface needs to be taken into account. This is done by the so-called
    // interface (or edge) matrices that compute the part of the residual that
    // is missed by the level matrix with
    // homogeneous Dirichlet conditions. We refer to the @ref mg_paper
    // "Multigrid paper by Janssen and Kanschat" for more details.
    //
    // For the implementation of those interface matrices, there is already a
    // pre-defined class MatrixFreeOperators::MGInterfaceOperator that wraps
    // the routines MatrixFreeOperators::Base::vmult_interface_down() and
    // MatrixFreeOperators::Base::vmult_interface_up() in a new class with @p
    // vmult() and @p Tvmult() operations (that were originally written for
    // matrices, hence expecting those names). Note that vmult_interface_down
    // is used during the restriction phase of the multigrid V-cycle, whereas
    // vmult_interface_up is used during the prolongation phase.
    //
    // Once the interface matrix is created, we set up the remaining Multigrid
    // preconditioner infrastructure in complete analogy to step-16 to obtain
    // a @p preconditioner object that can be applied to a matrix.
    mg::Matrix<LinearAlgebra::distributed::Vector<float>> mg_matrix(
      mg_matrices);

    MGLevelObject<MatrixFreeOperators::MGInterfaceOperator<LevelMatrixType>>
      mg_interface_matrices;
    mg_interface_matrices.resize(0, triangulation.n_global_levels() - 1);
    for (unsigned int level = 0; level < triangulation.n_global_levels();
         ++level)
      mg_interface_matrices[level].initialize(mg_matrices[level]);
    mg::Matrix<LinearAlgebra::distributed::Vector<float>> mg_interface(
      mg_interface_matrices);

    Multigrid<LinearAlgebra::distributed::Vector<float>> mg(
      mg_matrix, mg_coarse, mg_transfer, mg_smoother, mg_smoother);
    mg.set_edge_matrices(mg_interface, mg_interface);

    PreconditionMG<dim,
                   LinearAlgebra::distributed::Vector<float>,
                   MGTransferMatrixFree<dim, float>>
      preconditioner(dof_handler, mg, mg_transfer);

    // The setup of the multigrid routines is quite easy and one cannot see
    // any difference in the solve process as compared to step-16. All the
    // magic is hidden behind the implementation of the LaplaceOperator::vmult
    // operation. Note that we print out the solve time and the accumulated
    // setup time through standard out, i.e., in any case, whereas detailed
    // times for the setup operations are only printed in case the flag for
    // detail_times in the constructor is changed.

    SolverControl solver_control(100, 1e-12 * system_rhs.l2_norm());
    SolverCG<LinearAlgebra::distributed::Vector<double>> cg(solver_control);
    setup_time += time.wall_time();
    time_details << "MG build smoother time     (CPU/wall) " << time.cpu_time()
                 << "s/" << time.wall_time() << "s\n";
    pcout << "Total setup time               (wall) " << setup_time << "s\n";

    time.reset();
    time.start();
    constraints.set_zero(solution);
    cg.solve(system_matrix, solution, system_rhs, preconditioner);

    constraints.distribute(solution);

    pcout << "Time solve (" << solver_control.last_step() << " iterations)"
          << (solver_control.last_step() < 10 ? "  " : " ") << "(CPU/wall) "
          << time.cpu_time() << "s/" << time.wall_time() << "s\n";
  }



  // @sect4{LaplaceProblem::output_results}

  // Here is the data output, which is a simplified version of step-5. We use
  // the standard VTU (= compressed VTK) output for each grid produced in the
  // refinement process. In addition, we use a compression algorithm that is
  // optimized for speed rather than disk usage. The default setting (which
  // optimizes for disk usage) makes saving the output take about 4 times as
  // long as running the linear solver, while setting
  // DataOutBase::CompressionLevel to
  // best_speed lowers this to only one fourth the time
  // of the linear solve.
  //
  // We disable the output when the mesh gets too large. A variant of this
  // program has been run on hundreds of thousands MPI ranks with as many as
  // 100 billion grid cells, which is not directly accessible to classical
  // visualization tools.
  template <int dim>
  void LaplaceProblem<dim>::output_results(const unsigned int cycle) const
  {
    Timer time;
    if (triangulation.n_global_active_cells() > 1000000)
      return;

    DataOut<dim> data_out;

    solution.update_ghost_values();
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "solution");
    data_out.build_patches(mapping);

    DataOutBase::VtkFlags flags;
    flags.compression_level = DataOutBase::CompressionLevel::best_speed;
    data_out.set_flags(flags);
    data_out.write_vtu_with_pvtu_record(
      "./", "solution", cycle, MPI_COMM_WORLD, 3);

    time_details << "Time write output          (CPU/wall) " << time.cpu_time()
                 << "s/" << time.wall_time() << "s\n";
  }



  // @sect4{LaplaceProblem::run}

  // The function that runs the program is very similar to the one in
  // step-16. We do few refinement steps in 3d compared to 2d, but that's
  // it.
  //
  // Before we run the program, we output some information about the detected
  // vectorization level as discussed in the introduction.
  template <int dim>
  void LaplaceProblem<dim>::run()
  {
    {
      const unsigned int n_vect_doubles = VectorizedArray<double>::size();
      const unsigned int n_vect_bits    = 8 * sizeof(double) * n_vect_doubles;

      pcout << "Vectorization over " << n_vect_doubles
            << " doubles = " << n_vect_bits << " bits ("
            << Utilities::System::get_current_vectorization_level() << ')'
            << std::endl;
    }

    for (unsigned int cycle = 0; cycle < 9 - dim; ++cycle)
      {
        pcout << "Cycle " << cycle << std::endl;

        if (cycle == 0)
          {
            GridGenerator::hyper_cube(triangulation, 0., 1.);
            triangulation.refine_global(3 - dim);
          }
        triangulation.refine_global(1);
        setup_system();
        assemble_rhs();
        solve();
        output_results(cycle);
        pcout << std::endl;
      };
  }
} // namespace Step37



// @sect3{The <code>main</code> function}

// Apart from the fact that we set up the MPI framework according to step-40,
// there are no surprises in the main function.
int main(int argc, char *argv[])
{
  try
    {
      using namespace Step37;

      Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);

      LaplaceProblem<dimension> laplace_problem;
      laplace_problem.run();
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
