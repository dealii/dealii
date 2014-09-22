/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2009 - 2013 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Authors: Katharina Kormann, Martin Kronbichler, Uppsala University, 2009-2012
 */


// First include the necessary files from the deal.II library.
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/timer.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/multigrid/multigrid.h>
#include <deal.II/multigrid/mg_transfer.h>
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
#include <deal.II/matrix_free/fe_evaluation.h>

#include <fstream>
#include <sstream>


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
  const unsigned int dimension = 3;


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
    Coefficient ()  : Function<dim>() {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;

    template <typename number>
    number value (const Point<dim,number> &p,
                  const unsigned int component = 0) const;

    virtual void value_list (const std::vector<Point<dim> > &points,
                             std::vector<double>            &values,
                             const unsigned int              component = 0) const;
  };



  // This is the new function mentioned above: Evaluate the coefficient for
  // abstract type @p number. It might be just a usual double, but it can also
  // be a somewhat more complicated type that we call VectorizedArray. This
  // data type is essentially a short array of doubles as discussed in the
  // introduction that holds data from several cells. For example, we evaluate
  // the coefficient shown here not on a simple point as usually done, but we
  // hand it a Point<dim,VectorizedArray<double> > point, which is actually a
  // collection of two points in the case of SSE2. Do not confuse the entries
  // in VectorizedArray<double> with the different coordinates of the
  // point. Indeed, the data is laid out such that <code>p[0]</code> returns a
  // VectorizedArray<double>, which in turn contains the x-coordinate for the
  // first point and the second point. You may access the coordinates
  // individually using e.g. <code>p[0][j]</code>, j=0,1, but it is
  // recommended to define operations on a VectorizedArray as much as possible
  // in order to make use of vectorized operations.
  //
  // In the function implementation, we assume that the number type overloads
  // basic arithmetic operations, so we just write the code as usual. The
  // standard functions @p value and value_list that are virtual functions
  // contained in the base class are then computed from the templated function
  // with double type, in order to avoid duplicating code.
  template <int dim>
  template <typename number>
  number Coefficient<dim>::value (const Point<dim,number> &p,
                                  const unsigned int /*component*/) const
  {
    return 1. / (0.05 + 2.*p.square());
  }



  template <int dim>
  double Coefficient<dim>::value (const Point<dim>  &p,
                                  const unsigned int component) const
  {
    return value<double>(p,component);
  }



  template <int dim>
  void Coefficient<dim>::value_list (const std::vector<Point<dim> > &points,
                                     std::vector<double>            &values,
                                     const unsigned int              component) const
  {
    Assert (values.size() == points.size(),
            ExcDimensionMismatch (values.size(), points.size()));
    Assert (component == 0,
            ExcIndexRange (component, 0, 1));

    const unsigned int n_points = points.size();
    for (unsigned int i=0; i<n_points; ++i)
      values[i] = value<double>(points[i],component);
  }



  // @sect3{Matrix-free implementation}

  // The following class, called <code>LaplaceOperator</code>, implements the
  // differential operator. For all practical purposes, it is a matrix, i.e.,
  // you can ask it for its size (member functions <code>m(), n()</code>) and
  // you can apply it to a vector (the various variants of the
  // <code>vmult()</code> function). The difference to a real matrix of course
  // lies in the fact that this class doesn't actually store the
  // <i>elements</i> of the matrix, but only knows how to compute the action
  // of the operator when applied to a vector.

  // In this program, we want to make use of the data cache for finite element
  // operator application that is integrated in deal.II. The main class that
  // collects all data is called MatrixFree. It contains mapping information
  // (Jacobians) and index relations between local and global degrees of
  // freedom. It also contains constraints like the ones from Dirichlet
  // boundary conditions (or hanging nodes, if we had any). Moreover, it can
  // issue a loop over all cells in %parallel, where it makes sure that only
  // cells are worked on that do not share any degree of freedom (this makes
  // the loop thread-safe when writing into destination vectors). This is a
  // more advanced strategy compared to the WorkStream class described in the
  // @ref threads module that serializes operations that might not be
  // thread-safe. Of course, to not destroy thread-safety, we have to be
  // careful when writing into class-global structures.
  //
  // First comes the implementation of the matrix-free class. It provides some
  // standard information we expect for matrices (like returning the
  // dimensions of the matrix), it implements matrix-vector multiplications in
  // several forms (transposed and untransposed), and it provides functions
  // for initializing the structure with data. The class has three template
  // arguments, one for the dimension (as many deal.II classes carry), one for
  // the degree of the finite element (which we need to enable efficient
  // computations through the FEEvaluation class), and one for the underlying
  // scalar type. We want to use <code>double</code> numbers (i.e., double
  // precision, 64-bit floating point) for the final matrix, but floats
  // (single precision, 32-bit floating point numbers) for the multigrid level
  // matrices (as that is only a preconditioner, and floats can be worked with
  // twice as fast).
  //
  // In this class, we store the actual MatrixFree object, the variable
  // coefficient that is evaluated at all quadrature points (so that we don't
  // have to recompute it during matrix-vector products), and a vector that
  // contains the diagonal of the matrix that we need for the multigrid
  // smoother. We choose to let the user provide the diagonal in this program,
  // but we could also integrate a function in this class to evaluate the
  // diagonal. Unfortunately, this forces us to define matrix entries at two
  // places, once when we evaluate the product and once for the diagonal, but
  // the work is still much less than when we compute sparse matrices.
  //
  // As a sidenote, if we implemented several different operations on the same
  // grid and degrees of freedom (like a mass matrix and a Laplace matrix), we
  // would have to have two classes like the current one for each of the
  // operators (maybe with a common base class). However, in that case, we
  // would not store a MatrixFree object in this class to avoid doing the
  // expensive work of precomputing everything MatrixFree stores
  // twice. Rather, we would keep this object in the main class and simply
  // store a reference.
  //
  // @note Note that storing values of type
  // <code>VectorizedArray<number></code> requires care: Here, we use the
  // deal.II table class which is prepared to hold the data with correct
  // alignment. However, storing it in e.g.
  // <code>std::vector<VectorizedArray<number> ></code> is not possible with
  // vectorization: A certain alignment of the data with the memory address
  // boundaries is required (essentially, a VectorizedArray of 16 bytes length
  // as in SSE needs to start at a memory address that is divisible by
  // 16). The table class (as well as the AlignedVector class it is based on)
  // makes sure that this alignment is respected, whereas std::vector can in
  // general not, which may lead to segmentation faults at strange places for
  // some systems or suboptimal performance for other systems.
  template <int dim, int fe_degree, typename number>
  class LaplaceOperator : public Subscriptor
  {
  public:
    LaplaceOperator ();

    void clear();

    void reinit (const DoFHandler<dim>  &dof_handler,
                 const ConstraintMatrix  &constraints,
                 const unsigned int      level = numbers::invalid_unsigned_int);

    unsigned int m () const;
    unsigned int n () const;

    void vmult (Vector<double> &dst,
                const Vector<double> &src) const;
    void Tvmult (Vector<double> &dst,
                 const Vector<double> &src) const;
    void vmult_add (Vector<double> &dst,
                    const Vector<double> &src) const;
    void Tvmult_add (Vector<double> &dst,
                     const Vector<double> &src) const;

    number el (const unsigned int row,
               const unsigned int col) const;
    void set_diagonal (const Vector<number> &diagonal);

    std::size_t memory_consumption () const;

  private:
    void local_apply (const MatrixFree<dim,number>    &data,
                      Vector<double>                      &dst,
                      const Vector<double>                &src,
                      const std::pair<unsigned int,unsigned int> &cell_range) const;

    void evaluate_coefficient(const Coefficient<dim> &function);

    MatrixFree<dim,number>      data;
    Table<2, VectorizedArray<number> > coefficient;

    Vector<number>  diagonal_values;
    bool            diagonal_is_available;
  };



  // This is the constructor of the @p LaplaceOperator class. All it does is
  // to subscribe to the general deal.II @p Subscriptor scheme that makes sure
  // that we do not delete an object of this class as long as it used
  // somewhere else, e.g. in a preconditioner.
  template <int dim, int fe_degree, typename number>
  LaplaceOperator<dim,fe_degree,number>::LaplaceOperator ()
    :
    Subscriptor()
  {}



  // The next functions return the number of rows and columns of the global
  // matrix (i.e. the dimensions of the operator this class represents, the
  // point of this tutorial program was, after all, that we don't actually
  // store the elements of the rows and columns of this operator). Since the
  // matrix is square, the returned numbers are the same. We get the number
  // from the vector partitioner stored in the data field (a partitioner
  // distributes elements of a vector onto a number of different machines if
  // programs are run in %parallel; since this program is written to run on
  // only a single machine, the partitioner will simply say that all elements
  // of the vector -- or, in the current case, all rows and columns of a
  // matrix -- are stored on the current machine).
  template <int dim, int fe_degree, typename number>
  unsigned int
  LaplaceOperator<dim,fe_degree,number>::m () const
  {
    return data.get_vector_partitioner()->size();
  }



  template <int dim, int fe_degree, typename number>
  unsigned int
  LaplaceOperator<dim,fe_degree,number>::n () const
  {
    return data.get_vector_partitioner()->size();
  }



  template <int dim, int fe_degree, typename number>
  void
  LaplaceOperator<dim,fe_degree,number>::clear ()
  {
    data.clear();
    diagonal_is_available = false;
    diagonal_values.reinit(0);
  }


  // @sect4{Initialization}

  // Once we have created the multigrid dof_handler and the constraints, we
  // can call the reinit function for each level of the multigrid routine
  // (and the active cells). The main purpose of the reinit function is to
  // setup the <code> MatrixFree </code> instance for the problem. Also, the
  // coefficient is evaluated. For this, we need to activate the update flag
  // in the AdditionalData field of MatrixFree that enables the storage of
  // quadrature point coordinates in real space (by default, it only caches
  // data for gradients (inverse transposed Jacobians) and JxW values). Note
  // that if we call the reinit function without specifying the level (i.e.,
  // giving <code>level = numbers::invalid_unsigned_int</code>), we have told
  // the class to loop over the active cells.
  //
  // We also set one option regarding task parallelism. We choose to use the
  // @p partition_color strategy, which is based on subdivision of cells into
  // partitions where cells in partition $k$ (or, more precisely, the degrees
  // of freedom on these cells) only interact with cells in partitions $k-1$,
  // $k$, and $k+1$. Within each partition, cells are colored in such a way
  // that cells with the same color do not share degrees of freedom and can,
  // therefore, be worked on at the same time without interference. This
  // determines a task dependency graph that is scheduled by the Intel
  // Threading Building Blocks library. Another option would be the strategy
  // @p partition_partition, which performs better when the grid is more
  // unstructured. We could also manually set the size of chunks that form one
  // task in the scheduling process by setting @p tasks_block_size, but the
  // default strategy to let the function decide works well already.
  //
  // To initialize the coefficient, we directly give it the Coefficient class
  // defined above and then select the method
  // <code>coefficient_function.value</code> with vectorized number (which the
  // compiler can deduce from the point data type). The use of the
  // FEEvaluation class (and its template arguments) will be explained below.
  template <int dim, int fe_degree, typename number>
  void
  LaplaceOperator<dim,fe_degree,number>::reinit (const DoFHandler<dim>  &dof_handler,
                                                 const ConstraintMatrix  &constraints,
                                                 const unsigned int      level)
  {
    typename MatrixFree<dim,number>::AdditionalData additional_data;
    additional_data.tasks_parallel_scheme =
      MatrixFree<dim,number>::AdditionalData::partition_color;
    additional_data.level_mg_handler = level;
    additional_data.mapping_update_flags = (update_gradients | update_JxW_values |
                                            update_quadrature_points);
    data.reinit (dof_handler, constraints, QGauss<1>(fe_degree+1),
                 additional_data);
    evaluate_coefficient(Coefficient<dim>());
  }



  template <int dim, int fe_degree, typename number>
  void
  LaplaceOperator<dim,fe_degree,number>::
  evaluate_coefficient (const Coefficient<dim> &coefficient_function)
  {
    const unsigned int n_cells = data.n_macro_cells();
    FEEvaluation<dim,fe_degree,fe_degree+1,1,number> phi (data);
    coefficient.reinit (n_cells, phi.n_q_points);
    for (unsigned int cell=0; cell<n_cells; ++cell)
      {
        phi.reinit (cell);
        for (unsigned int q=0; q<phi.n_q_points; ++q)
          coefficient(cell,q) =
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
  // not cause any race condition. Note that the total range of cells as
  // visible in this class is usually not equal to the number of (active)
  // cells in the triangulation.  In fact, "cell" may be the wrong term to
  // begin with, since it is rather a collection of quadrature points from
  // several cells, and the MatrixFree class groups the quadrature points of
  // several cells into one block to enable a higher degree of vectorization.
  // The number of such "cells" is stored in MatrixFree and can be queried
  // through MatrixFree::n_macro_cells(). Compared to the
  // deal.II cell iterators, in this class all cells are laid out in a plain
  // array with no direct knowledge of level or neighborship relations, which
  // makes it possible to index the cells by unsigned integers.
  //
  // The implementation of the Laplace operator is quite simple: First, we
  // need to create an object FEEvaluation that contains the computational
  // kernels and has data fields to store temporary results (e.g. gradients
  // evaluated on all quadrature points on a collection of a few cells). Note
  // that temporary results do not use a lot of memory, and since we specify
  // template arguments with the element order, the data is stored on the
  // stack (without expensive memory allocation). Usually, one only needs to
  // set two template arguments, the dimension as first argument and the
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
  // (first slot) and Hessian slot (third slot). Note that the FEEvaluation
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
  // function in the ConstraintMatrix (only that we now store the local vector
  // in the FEEvaluation object, as are the indices between local and global
  // degrees of freedom).  </ol>
  template <int dim, int fe_degree, typename number>
  void
  LaplaceOperator<dim,fe_degree,number>::
  local_apply (const MatrixFree<dim,number>         &data,
               Vector<double>                       &dst,
               const Vector<double>                 &src,
               const std::pair<unsigned int,unsigned int> &cell_range) const
  {
    FEEvaluation<dim,fe_degree,fe_degree+1,1,number> phi (data);

    for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
      {
        phi.reinit (cell);
        phi.read_dof_values(src);
        phi.evaluate (false,true,false);
        for (unsigned int q=0; q<phi.n_q_points; ++q)
          phi.submit_gradient (coefficient(cell,q) *
                               phi.get_gradient(q), q);
        phi.integrate (false,true);
        phi.distribute_local_to_global (dst);
      }
  }



  // @sect4{vmult functions}

  // Now to the @p vmult function that is called externally: In addition to
  // what we do in a @p vmult_add function further down, we set the
  // destination to zero first. The transposed matrix-vector is needed for
  // well-defined multigrid preconditioner operations. Since we solve a
  // Laplace problem, this is the same operation, and we just refer to the
  // vmult operation.
  template <int dim, int fe_degree, typename number>
  void
  LaplaceOperator<dim,fe_degree,number>::vmult (Vector<double>       &dst,
                                                const Vector<double> &src) const
  {
    dst = 0;
    vmult_add (dst, src);
  }



  template <int dim, int fe_degree, typename number>
  void
  LaplaceOperator<dim,fe_degree,number>::Tvmult (Vector<double>       &dst,
                                                 const Vector<double> &src) const
  {
    dst = 0;
    vmult_add (dst,src);
  }



  template <int dim, int fe_degree, typename number>
  void
  LaplaceOperator<dim,fe_degree,number>::Tvmult_add (Vector<double>       &dst,
                                                     const Vector<double> &src) const
  {
    vmult_add (dst,src);
  }



  // This function implements the loop over all cells. This is done with the
  // @p cell_loop of the MatrixFree class, which takes the operator() of this
  // class with arguments MatrixFree, OutVector, InVector, cell_range. Note
  // that we could also use a simple function as local operation in case we
  // had constant coefficients (all we need then is the MatrixFree, the
  // vectors and the cell range), but since the coefficient is stored in a
  // variable of this class, we cannot use that variant here. The cell loop is
  // automatically performed on several threads if multithreading is enabled
  // (this class uses a quite elaborate algorithm to work on cells that do not
  // share any degrees of freedom that could possibly give rise to race
  // conditions, using the dynamic task scheduler of the Intel Threading
  // Building Blocks).
  //
  // After the cell loop, we need to touch the constrained degrees of freedom:
  // Since the assembly loop automatically resolves constraints (just as the
  // ConstraintMatrix::distribute_local_to_global call does), it does not
  // compute any contribution for constrained degrees of freedom. In other
  // words, the entries for constrained DoFs remain zero after the first part
  // of this function, as if the matrix had empty rows and columns for
  // constrained degrees of freedom. On the other hand, iterative solvers like
  // CG only work for non-singular matrices, so we have to modify the
  // operation on constrained DoFs. The easiest way to do that is to pretend
  // that the sub-block of the matrix that corresponds to constrained DoFs is
  // the identity matrix, in which case application of the matrix would simply
  // copy the elements of the right hand side vector into the left hand
  // side. In general, however, one needs to make sure that the diagonal
  // entries of this sub-block are of the same order of magnitude as the
  // diagonal elements of the rest of the matrix.  Here, the domain extent is
  // of unit size, so we can simply choose unit size. If we had domains that
  // are far away from unit size, we would need to choose a number that is
  // close to the size of other diagonal matrix entries, so that these
  // artificial eigenvalues do not change the eigenvalue spectrum (and make
  // convergence with CG more difficult).
  template <int dim, int fe_degree, typename number>
  void
  LaplaceOperator<dim,fe_degree,number>::vmult_add (Vector<double>       &dst,
                                                    const Vector<double> &src) const
  {
    data.cell_loop (&LaplaceOperator::local_apply, this, dst, src);

    const std::vector<unsigned int> &
    constrained_dofs = data.get_constrained_dofs();
    for (unsigned int i=0; i<constrained_dofs.size(); ++i)
      dst(constrained_dofs[i]) += src(constrained_dofs[i]);
  }



  // The next function is used to return entries of the matrix. Since this
  // class is intended not to store the matrix entries, it would make no sense
  // to provide access to all those elements. However, diagonal entries are
  // explicitly needed for the implementation of the Chebyshev smoother that
  // we intend to use in the multigrid preconditioner. This matrix is equipped
  // with a vector that stores the diagonal.
  template <int dim, int fe_degree, typename number>
  number
  LaplaceOperator<dim,fe_degree,number>::el (const unsigned int row,
                                             const unsigned int col) const
  {
    Assert (row == col, ExcNotImplemented());
    Assert (diagonal_is_available == true, ExcNotInitialized());
    return diagonal_values(row);
  }



  // Regarding the calculation of the diagonal, we expect the user to provide
  // a vector with the diagonal entries (and we will compute them in the code
  // below). We only need it for the level matrices of multigrid, not the
  // system matrix (since we only need these diagonals for the multigrid
  // smoother). Since we fill only elements into unconstrained entries, we
  // have to set constrained entries to one in order to avoid the same
  // problems as discussed above.
  template <int dim, int fe_degree, typename number>
  void
  LaplaceOperator<dim,fe_degree,number>::set_diagonal(const Vector<number> &diagonal)
  {
    AssertDimension (m(), diagonal.size());

    diagonal_values = diagonal;

    const std::vector<unsigned int> &
    constrained_dofs = data.get_constrained_dofs();
    for (unsigned int i=0; i<constrained_dofs.size(); ++i)
      diagonal_values(constrained_dofs[i]) = 1.0;

    diagonal_is_available = true;
  }



  // Eventually, we provide a function that calculates how much memory this
  // class uses. We just need to sum up the memory consumption in the
  // MatrixFree object and the memory for storing the other member
  // variables. As a remark: In 3D and for Cartesian meshes, most memory is
  // consumed for storing the vector indices on the local cells (corresponding
  // to local_dof_indices). For general (non-Cartesian) meshes, the cached
  // Jacobian transformation consumes most memory.
  template <int dim, int fe_degree, typename number>
  std::size_t
  LaplaceOperator<dim,fe_degree,number>::memory_consumption () const
  {
    return (data.memory_consumption () +
            coefficient.memory_consumption() +
            diagonal_values.memory_consumption() +
            MemoryConsumption::memory_consumption(diagonal_is_available));
  }



  // @sect3{LaplaceProblem class}

  // This class is based on the one in step-16. However, we replaced the
  // SparseMatrix<double> class by our matrix-free implementation, which means
  // that we can also skip the sparsity patterns. Notice that we define the
  // LaplaceOperator class with the degree of finite element as template
  // argument (the value is defined at the top of the file), and that we use
  // float numbers for the multigrid level matrices.
  //
  // The class also has a member variable to keep track of all the time we
  // spend on setting up the entire chain of data before we actually go about
  // solving the problem. In addition, there is an output stream (that is
  // disabled by default) that can be used to output details for the
  // individual setup operations instead of the summary only that is printed
  // out by default.
  template <int dim>
  class LaplaceProblem
  {
  public:
    LaplaceProblem ();
    void run ();

  private:
    void setup_system ();
    void assemble_system ();
    void assemble_multigrid ();
    void solve ();
    void output_results (const unsigned int cycle) const;

    typedef LaplaceOperator<dim,degree_finite_element,double> SystemMatrixType;
    typedef LaplaceOperator<dim,degree_finite_element,float>  LevelMatrixType;

    Triangulation<dim>               triangulation;
    FE_Q<dim>                        fe;
    DoFHandler<dim>                  dof_handler;
    ConstraintMatrix                 constraints;

    SystemMatrixType                 system_matrix;
    MGLevelObject<LevelMatrixType>   mg_matrices;
    FullMatrix<float>                coarse_matrix;
    MGLevelObject<ConstraintMatrix>  mg_constraints;

    Vector<double>                   solution;
    Vector<double>                   system_rhs;

    double                           setup_time;
    ConditionalOStream               time_details;
  };



  // When we initialize the finite element, we of course have to use the
  // degree specified at the top of the file as well (otherwise, an exception
  // will be thrown at some point, since the computational kernel defined in
  // the templated LaplaceOperator class and the information from the finite
  // element read out by MatrixFree will not match).
  template <int dim>
  LaplaceProblem<dim>::LaplaceProblem ()
    :
    fe (degree_finite_element),
    dof_handler (triangulation),
    time_details (std::cout, false)
  {}



  // @sect4{LaplaceProblem::setup_system}

  // This is the function of step-16 with relevant changes due to the
  // LaplaceOperator class. We do not use adaptive grids, so we do not have to
  // compute edge matrices. Thus, all we do is to implement Dirichlet boundary
  // conditions through the ConstraintMatrix, set up the (one-dimensional)
  // quadrature that should be used by the matrix-free class, and call the
  // initialization functions.
  //
  // In the process, we output data on both the run time of the program as
  // well as on memory consumption, where we output memory data in megabytes
  // (1 million bytes).
  template <int dim>
  void LaplaceProblem<dim>::setup_system ()
  {
    Timer time;
    time.start ();
    setup_time = 0;

    system_matrix.clear();
    mg_matrices.clear();
    mg_constraints.clear();

    dof_handler.distribute_dofs (fe);
    dof_handler.distribute_mg_dofs (fe);

    std::cout << "Number of degrees of freedom: "
              << dof_handler.n_dofs()
              << std::endl;

    constraints.clear();
    VectorTools::interpolate_boundary_values (dof_handler,
                                              0,
                                              ZeroFunction<dim>(),
                                              constraints);
    constraints.close();
    setup_time += time.wall_time();
    time_details << "Distribute DoFs & B.C.     (CPU/wall) "
                 << time() << "s/" << time.wall_time() << "s" << std::endl;
    time.restart();

    system_matrix.reinit (dof_handler, constraints);
    std::cout.precision(4);
    std::cout << "System matrix memory consumption:     "
              << system_matrix.memory_consumption()*1e-6
              << " MB."
              << std::endl;

    solution.reinit (dof_handler.n_dofs());
    system_rhs.reinit (dof_handler.n_dofs());

    setup_time += time.wall_time();
    time_details << "Setup matrix-free system   (CPU/wall) "
                 << time() << "s/" << time.wall_time() << "s" << std::endl;
    time.restart();

    // Next, initialize the matrices for the multigrid method on all the
    // levels. The function MGTools::make_boundary_list returns for each
    // multigrid level which degrees of freedom are located on a Dirichlet
    // boundary; we force these DoFs to have value zero by adding to the
    // ConstraintMatrix object a zero condition by using the command
    // ConstraintMatrix::add_line. Once this is done, we close the
    // ConstraintMatrix on each level so it can be used to read out indices
    // internally in the MatrixFree.
    const unsigned int nlevels = triangulation.n_levels();
    mg_matrices.resize(0, nlevels-1);
    mg_constraints.resize (0, nlevels-1);

    typename FunctionMap<dim>::type dirichlet_boundary;
    ZeroFunction<dim>               homogeneous_dirichlet_bc (1);
    dirichlet_boundary[0] = &homogeneous_dirichlet_bc;
    std::vector<std::set<types::global_dof_index> > boundary_indices(triangulation.n_levels());
    MGTools::make_boundary_list (dof_handler,
                                 dirichlet_boundary,
                                 boundary_indices);
    for (unsigned int level=0; level<nlevels; ++level)
      {
        std::set<types::global_dof_index>::iterator bc_it = boundary_indices[level].begin();
        for ( ; bc_it != boundary_indices[level].end(); ++bc_it)
          mg_constraints[level].add_line(*bc_it);

        mg_constraints[level].close();
        mg_matrices[level].reinit(dof_handler,
                                  mg_constraints[level],
                                  level);
      }
    coarse_matrix.reinit (dof_handler.n_dofs(0),
                          dof_handler.n_dofs(0));
    setup_time += time.wall_time();
    time_details << "Setup matrix-free levels   (CPU/wall) "
                 << time() << "s/" << time.wall_time() << "s" << std::endl;
  }



  // @sect4{LaplaceProblem::assemble_system}

  // The assemble function is significantly reduced compared to step-16. All
  // we need to do is to assemble the right hand side. That is the same as in
  // many other tutorial programs. In the end, we condense the constraints
  // from Dirichlet boundary conditions away from the right hand side.
  template <int dim>
  void LaplaceProblem<dim>::assemble_system ()
  {
    Timer time;
    QGauss<dim>  quadrature_formula(fe.degree+1);
    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values   | update_JxW_values);

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.size();

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                   endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
        cell->get_dof_indices (local_dof_indices);
        fe_values.reinit (cell);
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            double rhs_val = 0;
            for (unsigned int q=0; q<n_q_points; ++q)
              rhs_val += (fe_values.shape_value(i,q) * 1.0 *
                          fe_values.JxW(q));
            system_rhs(local_dof_indices[i]) += rhs_val;
          }
      }
    constraints.condense(system_rhs);
    setup_time += time.wall_time();
    time_details << "Assemble right hand side   (CPU/wall) "
                 << time() << "s/" << time.wall_time() << "s" << std::endl;
  }


  // @sect4{LaplaceProblem::assemble_multigrid}

  // Here is another assemble function. Again, it is simpler than assembling
  // matrices. We need to compute the diagonal of the Laplace matrices on the
  // individual levels, send the final matrices to the LaplaceOperator class,
  // and we need to compute the full matrix on the coarsest level (since that
  // is inverted exactly in the deal.II multigrid implementation).
  template <int dim>
  void LaplaceProblem<dim>::assemble_multigrid ()
  {
    Timer time;
    coarse_matrix = 0;
    QGauss<dim>  quadrature_formula(fe.degree+1);
    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_gradients  | update_inverse_jacobians |
                             update_quadrature_points | update_JxW_values);

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.size();

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    const Coefficient<dim>    coefficient;
    std::vector<double>       coefficient_values (n_q_points);
    FullMatrix<double>        local_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>            local_diagonal (dofs_per_cell);

    const unsigned int n_levels = triangulation.n_levels();
    std::vector<Vector<float> > diagonals (n_levels);
    for (unsigned int level=0; level<n_levels; ++level)
      diagonals[level].reinit (dof_handler.n_dofs(level));

    std::vector<unsigned int> cell_no(triangulation.n_levels());
    typename DoFHandler<dim>::cell_iterator cell = dof_handler.begin(),
                                            endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
        const unsigned int level = cell->level();
        cell->get_mg_dof_indices (local_dof_indices);
        fe_values.reinit (cell);
        coefficient.value_list (fe_values.get_quadrature_points(),
                                coefficient_values);

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            double local_diag = 0;
            for (unsigned int q=0; q<n_q_points; ++q)
              local_diag += ((fe_values.shape_grad(i,q) *
                              fe_values.shape_grad(i,q)) *
                             coefficient_values[q] * fe_values.JxW(q));
            local_diagonal(i) = local_diag;
          }
        mg_constraints[level].distribute_local_to_global(local_diagonal,
                                                         local_dof_indices,
                                                         diagonals[level]);

        if (level == 0)
          {
            local_matrix = 0;

            for (unsigned int i=0; i<dofs_per_cell; ++i)
              for (unsigned int j=0; j<dofs_per_cell; ++j)
                {
                  double add_value = 0;
                  for (unsigned int q=0; q<n_q_points; ++q)
                    add_value += (fe_values.shape_grad(i,q) *
                                  fe_values.shape_grad(j,q) *
                                  coefficient_values[q] *
                                  fe_values.JxW(q));
                  local_matrix(i,j) = add_value;
                }
            mg_constraints[0].distribute_local_to_global (local_matrix,
                                                          local_dof_indices,
                                                          coarse_matrix);
          }
      }

    for (unsigned int level=0; level<n_levels; ++level)
      mg_matrices[level].set_diagonal (diagonals[level]);

    setup_time += time.wall_time();
    time_details << "Assemble MG diagonal       (CPU/wall) "
                 << time() << "s/" << time.wall_time() << "s" << std::endl;
  }



  // @sect4{LaplaceProblem::solve}

  // The solution process again looks like step-16. We now use a Chebyshev
  // smoother instead of SOR (SOR would be very difficult to implement because
  // we do not have the matrix elements available explicitly, and it is
  // difficult to make it work efficiently in %parallel). The multigrid
  // classes provide a simple interface for using the Chebyshev smoother which
  // is defined in a preconditioner class: MGSmootherPrecondition.
  template <int dim>
  void LaplaceProblem<dim>::solve ()
  {
    Timer time;
    MGTransferPrebuilt<Vector<double> > mg_transfer;
    mg_transfer.build_matrices(dof_handler);
    setup_time += time.wall_time();
    time_details << "MG build transfer time     (CPU/wall) " << time()
                 << "s/" << time.wall_time() << "s\n";
    time.restart();

    MGCoarseGridHouseholder<float, Vector<double> > mg_coarse;
    mg_coarse.initialize(coarse_matrix);
    setup_time += time.wall_time();
    time_details << "MG coarse time             (CPU/wall) " << time()
                 << "s/" << time.wall_time() << "s\n";
    time.restart();

    typedef PreconditionChebyshev<LevelMatrixType,Vector<double> > SMOOTHER;
    MGSmootherPrecondition<LevelMatrixType, SMOOTHER, Vector<double> >
    mg_smoother;

    // Then, we initialize the smoother with our level matrices and the
    // mandatory additional data for the Chebyshev smoother. We use quite a
    // high degree here (6), since matrix-vector products are comparably cheap
    // and more parallel than the level-transfer operations. We choose to
    // smooth out a range of $[1.2 \hat{\lambda}_{\max}/10,1.2
    // \hat{\lambda}_{\max}]$ in the smoother where $\hat{\lambda}_{\max}$ is
    // an estimate of the largest eigenvalue. In order to compute that
    // eigenvalue, the Chebyshev initializations performs a few steps of a CG
    // algorithm without preconditioner. Since the highest eigenvalue is
    // usually the easiest one to find and a rough estimate is enough, we
    // choose 10 iterations.
    typename SMOOTHER::AdditionalData smoother_data;
    smoother_data.smoothing_range = 10.;
    smoother_data.degree = 6;
    smoother_data.eig_cg_n_iterations = 10;
    mg_smoother.initialize(mg_matrices, smoother_data);

    MGMatrix<LevelMatrixType, Vector<double> >
    mg_matrix(&mg_matrices);

    Multigrid<Vector<double> > mg(dof_handler,
                                  mg_matrix,
                                  mg_coarse,
                                  mg_transfer,
                                  mg_smoother,
                                  mg_smoother);
    PreconditionMG<dim, Vector<double>,
                   MGTransferPrebuilt<Vector<double> > >
                   preconditioner(dof_handler, mg, mg_transfer);

    // Finally, write out the memory consumption of the Multigrid object (or
    // rather, of its most significant components, since there is no built-in
    // function for the total multigrid object), then create the solver object
    // and solve the system. This is very easy, and we didn't even see any
    // difference in the solve process compared to step-16. The magic is all
    // hidden behind the implementation of the LaplaceOperator::vmult
    // operation. Note that we print out the solve time and the accumulated
    // setup time through standard out, i.e., in any case, whereas detailed
    // times for the setup operations are only printed in case the flag for
    // detail_times in the constructor is changed.
    const std::size_t multigrid_memory
      = (mg_matrices.memory_consumption() +
         mg_transfer.memory_consumption() +
         coarse_matrix.memory_consumption());
    std::cout << "Multigrid objects memory consumption: "
              << multigrid_memory * 1e-6
              << " MB."
              << std::endl;

    SolverControl           solver_control (1000, 1e-12*system_rhs.l2_norm());
    SolverCG<>              cg (solver_control);
    setup_time += time.wall_time();
    time_details << "MG build smoother time     (CPU/wall) " << time()
                 << "s/" << time.wall_time() << "s\n";
    std::cout << "Total setup time               (wall) " << setup_time
              << "s\n";

    time.reset();
    time.start();
    cg.solve (system_matrix, solution, system_rhs,
              preconditioner);


    std::cout << "Time solve ("
              << solver_control.last_step()
              << " iterations)  (CPU/wall) " << time() << "s/"
              << time.wall_time() << "s\n";
  }



  // @sect4{LaplaceProblem::output_results}

  // Here is the data output, which is a simplified version of step-5. We use
  // the standard VTU (= compressed VTK) output for each grid produced in the
  // refinement process.
  template <int dim>
  void LaplaceProblem<dim>::output_results (const unsigned int cycle) const
  {
    DataOut<dim> data_out;

    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (solution, "solution");
    data_out.build_patches ();

    std::ostringstream filename;
    filename << "solution-"
             << cycle
             << ".vtu";

    std::ofstream output (filename.str().c_str());
    data_out.write_vtu (output);
  }



  // @sect4{LaplaceProblem::run}

  // The function that runs the program is very similar to the one in
  // step-16. We make less refinement steps in 3D compared to 2D, but that's
  // it.
  template <int dim>
  void LaplaceProblem<dim>::run ()
  {
    for (unsigned int cycle=0; cycle<9-dim; ++cycle)
      {
        std::cout << "Cycle " << cycle << std::endl;

        if (cycle == 0)
          {
            GridGenerator::hyper_cube (triangulation, 0., 1.);
            triangulation.refine_global (3-dim);
          }
        triangulation.refine_global (1);
        setup_system ();
        assemble_system ();
        assemble_multigrid ();
        solve ();
        output_results (cycle);
        std::cout << std::endl;
      };
  }
}



// @sect3{The <code>main</code> function}

// This is as in most other programs.
int main ()
{
  try
    {
      using namespace Step37;

      deallog.depth_console(0);
      LaplaceProblem<dimension> laplace_problem;
      laplace_problem.run ();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
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
      std::cerr << std::endl << std::endl
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
