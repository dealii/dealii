/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2005 - 2013 by the deal.II authors
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
 * Author: Wolfgang Bangerth, Texas A&M University, 2005, 2006
 */


// @sect3{Include files}

// Since this program is only an adaptation of step-4, there is not much new
// stuff in terms of header files. In deal.II, we usually list include files
// in the order base-lac-grid-dofs-fe-numerics, followed by C++ standard
// include files:
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
// For our Schur complement solver, we need two new objects. One is a matrix
// object which acts as the inverse of a matrix by calling an iterative
// solver.
#include <deal.II/lac/iterative_inverse.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
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

// The last step is as in all previous programs:
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
    MixedLaplaceProblem (const unsigned int degree);
    void run ();

  private:
    void make_grid_and_dofs ();
    void assemble_system ();
    void solve ();
    void compute_errors () const;
    void output_results () const;

    const unsigned int   degree;

    Triangulation<dim>   triangulation;
    FESystem<dim>        fe;
    DoFHandler<dim>      dof_handler;

    // The second difference is that the sparsity pattern, the system matrix,
    // and solution and right hand side vectors are now blocked. What this
    // means and what one can do with such objects is explained in the
    // introduction to this program as well as further down below when we
    // explain the linear solvers and preconditioners for this problem:
    BlockSparsityPattern      sparsity_pattern;
    BlockSparseMatrix<double> system_matrix;

    BlockVector<double>       solution;
    BlockVector<double>       system_rhs;
  };


  // @sect3{Right hand side, boundary values, and exact solution}

  // Our next task is to define the right hand side of our problem (i.e., the
  // scalar right hand side for the pressure in the original Laplace
  // equation), boundary values for the pressure, as well as a function that
  // describes both the pressure and the velocity of the exact solution for
  // later computations of the error. Note that these functions have one, one,
  // and <code>dim+1</code> components, respectively, and that we pass the
  // number of components down to the <code>Function@<dim@></code> base
  // class. For the exact solution, we only declare the function that actually
  // returns the entire solution vector (i.e. all components of it) at
  // once. Here are the respective declarations:
  template <int dim>
  class RightHandSide : public Function<dim>
  {
  public:
    RightHandSide () : Function<dim>(1) {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
  };



  template <int dim>
  class PressureBoundaryValues : public Function<dim>
  {
  public:
    PressureBoundaryValues () : Function<dim>(1) {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
  };


  template <int dim>
  class ExactSolution : public Function<dim>
  {
  public:
    ExactSolution () : Function<dim>(dim+1) {}

    virtual void vector_value (const Point<dim> &p,
                               Vector<double>   &value) const;
  };


  // And then we also have to define these respective functions, of
  // course. Given our discussion in the introduction of how the solution
  // should look like, the following computations should be straightforward:
  template <int dim>
  double RightHandSide<dim>::value (const Point<dim>  &/*p*/,
                                    const unsigned int /*component*/) const
  {
    return 0;
  }



  template <int dim>
  double PressureBoundaryValues<dim>::value (const Point<dim>  &p,
                                             const unsigned int /*component*/) const
  {
    const double alpha = 0.3;
    const double beta = 1;
    return -(alpha*p[0]*p[1]*p[1]/2 + beta*p[0] - alpha*p[0]*p[0]*p[0]/6);
  }



  template <int dim>
  void
  ExactSolution<dim>::vector_value (const Point<dim> &p,
                                    Vector<double>   &values) const
  {
    Assert (values.size() == dim+1,
            ExcDimensionMismatch (values.size(), dim+1));

    const double alpha = 0.3;
    const double beta = 1;

    values(0) = alpha*p[1]*p[1]/2 + beta - alpha*p[0]*p[0]/2;
    values(1) = alpha*p[0]*p[1];
    values(2) = -(alpha*p[0]*p[1]*p[1]/2 + beta*p[0] - alpha*p[0]*p[0]*p[0]/6);
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
  // Possibly unsurprising, deal.II also has a base class not only for scalar
  // and generally vector-valued functions (the <code>Function</code> base
  // class) but also for functions that return tensors of fixed dimension and
  // rank, the <code>TensorFunction</code> template. Here, the function under
  // consideration returns a dim-by-dim matrix, i.e. a tensor of rank 2 and
  // dimension <code>dim</code>. We then choose the template arguments of the
  // base class appropriately.
  //
  // The interface that the <code>TensorFunction</code> class provides is
  // essentially equivalent to the <code>Function</code> class. In particular,
  // there exists a <code>value_list</code> function that takes a list of
  // points at which to evaluate the function, and returns the values of the
  // function in the second argument, a list of tensors:
  template <int dim>
  class KInverse : public TensorFunction<2,dim>
  {
  public:
    KInverse () : TensorFunction<2,dim>() {}

    virtual void value_list (const std::vector<Point<dim> > &points,
                             std::vector<Tensor<2,dim> >    &values) const;
  };


  // The implementation is less interesting. As in previous examples, we add a
  // check to the beginning of the class to make sure that the sizes of input
  // and output parameters are the same (see step-5 for a discussion of this
  // technique). Then we loop over all evaluation points, and for each one
  // first clear the output tensor and then set all its diagonal elements to
  // one (i.e. fill the tensor with the identity matrix):
  template <int dim>
  void
  KInverse<dim>::value_list (const std::vector<Point<dim> > &points,
                             std::vector<Tensor<2,dim> >    &values) const
  {
    Assert (points.size() == values.size(),
            ExcDimensionMismatch (points.size(), values.size()));

    for (unsigned int p=0; p<points.size(); ++p)
      {
        values[p].clear ();

        for (unsigned int d=0; d<dim; ++d)
          values[p][d][d] = 1.;
      }
  }



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
  MixedLaplaceProblem<dim>::MixedLaplaceProblem (const unsigned int degree)
    :
    degree (degree),
    fe (FE_RaviartThomas<dim>(degree), 1,
        FE_DGQ<dim>(degree), 1),
    dof_handler (triangulation)
  {}



  // @sect4{MixedLaplaceProblem::make_grid_and_dofs}

  // This next function starts out with well-known functions calls that create
  // and refine a mesh, and then associate degrees of freedom with it:
  template <int dim>
  void MixedLaplaceProblem<dim>::make_grid_and_dofs ()
  {
    GridGenerator::hyper_cube (triangulation, -1, 1);
    triangulation.refine_global (3);

    dof_handler.distribute_dofs (fe);

    // However, then things become different. As mentioned in the
    // introduction, we want to subdivide the matrix into blocks corresponding
    // to the two different kinds of variables, velocity and pressure. To this
    // end, we first have to make sure that the indices corresponding to
    // velocities and pressures are not intermingled: First all velocity
    // degrees of freedom, then all pressure DoFs. This way, the global matrix
    // separates nicely into a $2 \times 2$ system. To achieve this, we have to
    // renumber degrees of freedom base on their vector component, an
    // operation that conveniently is already implemented:
    DoFRenumbering::component_wise (dof_handler);

    // The next thing is that we want to figure out the sizes of these blocks,
    // so that we can allocate an appropriate amount of space. To this end, we
    // call the <code>DoFTools::count_dofs_per_component</code> function that
    // counts how many shape functions are non-zero for a particular vector
    // component. We have <code>dim+1</code> vector components, and we have to
    // use the knowledge that for Raviart-Thomas elements all shape functions
    // are nonzero in all components. In other words, the number of velocity
    // shape functions equals the number of overall shape functions that are
    // nonzero in the zeroth vector component. On the other hand, the number
    // of pressure variables equals the number of shape functions that are
    // nonzero in the dim-th component. Let us compute these numbers and then
    // create some nice output with that:
    std::vector<types::global_dof_index> dofs_per_component (dim+1);
    DoFTools::count_dofs_per_component (dof_handler, dofs_per_component);
    const unsigned int n_u = dofs_per_component[0],
                       n_p = dofs_per_component[dim];

    std::cout << "Number of active cells: "
              << triangulation.n_active_cells()
              << std::endl
              << "Total number of cells: "
              << triangulation.n_cells()
              << std::endl
              << "Number of degrees of freedom: "
              << dof_handler.n_dofs()
              << " (" << n_u << '+' << n_p << ')'
              << std::endl;

    // The next task is to allocate a sparsity pattern for the matrix that we
    // will create. The way this works is that we first obtain a guess for the
    // maximal number of nonzero entries per row (this could be done more
    // efficiently in this case, but we only want to solve relatively small
    // problems for which this is not so important). In the second step, we
    // allocate a $2 \times 2$ block pattern and then reinitialize each of the blocks
    // to its correct size using the <code>n_u</code> and <code>n_p</code>
    // variables defined above that hold the number of velocity and pressure
    // variables. In this second step, we only operate on the individual
    // blocks of the system. In the third step, we therefore have to instruct
    // the overlying block system to update its knowledge about the sizes of
    // the blocks it manages; this happens with the
    // <code>sparsity_pattern.collect_sizes()</code> call:
    const unsigned int
    n_couplings = dof_handler.max_couplings_between_dofs();

    sparsity_pattern.reinit (2,2);
    sparsity_pattern.block(0,0).reinit (n_u, n_u, n_couplings);
    sparsity_pattern.block(1,0).reinit (n_p, n_u, n_couplings);
    sparsity_pattern.block(0,1).reinit (n_u, n_p, n_couplings);
    sparsity_pattern.block(1,1).reinit (n_p, n_p, n_couplings);
    sparsity_pattern.collect_sizes();

    // Now that the sparsity pattern and its blocks have the correct sizes, we
    // actually need to construct the content of this pattern, and as usual
    // compress it, before we also initialize a block matrix with this block
    // sparsity pattern:
    DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
    sparsity_pattern.compress();

    system_matrix.reinit (sparsity_pattern);

    // Then we have to resize the solution and right hand side vectors in
    // exactly the same way:
    solution.reinit (2);
    solution.block(0).reinit (n_u);
    solution.block(1).reinit (n_p);
    solution.collect_sizes ();

    system_rhs.reinit (2);
    system_rhs.block(0).reinit (n_u);
    system_rhs.block(1).reinit (n_p);
    system_rhs.collect_sizes ();
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
  void MixedLaplaceProblem<dim>::assemble_system ()
  {
    QGauss<dim>   quadrature_formula(degree+2);
    QGauss<dim-1> face_quadrature_formula(degree+2);

    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values    | update_gradients |
                             update_quadrature_points  | update_JxW_values);
    FEFaceValues<dim> fe_face_values (fe, face_quadrature_formula,
                                      update_values    | update_normal_vectors |
                                      update_quadrature_points  | update_JxW_values);

    const unsigned int   dofs_per_cell   = fe.dofs_per_cell;
    const unsigned int   n_q_points      = quadrature_formula.size();
    const unsigned int   n_face_q_points = face_quadrature_formula.size();

    FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>       local_rhs (dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    // The next step is to declare objects that represent the source term,
    // pressure boundary value, and coefficient in the equation. In addition
    // to these objects that represent continuous functions, we also need
    // arrays to hold their values at the quadrature points of individual
    // cells (or faces, for the boundary values). Note that in the case of the
    // coefficient, the array has to be one of matrices.
    const RightHandSide<dim>          right_hand_side;
    const PressureBoundaryValues<dim> pressure_boundary_values;
    const KInverse<dim>               k_inverse;

    std::vector<double> rhs_values (n_q_points);
    std::vector<double> boundary_values (n_face_q_points);
    std::vector<Tensor<2,dim> > k_inverse_values (n_q_points);

    // Finally, we need a couple of extractors that we will use to get at the
    // velocity and pressure components of vector-valued shape
    // functions. Their function and use is described in detail in the @ref
    // vector_valued report. Essentially, we will use them as subscripts on
    // the FEValues objects below: the FEValues object describes all vector
    // components of shape functions, while after subscription, it will only
    // refer to the velocities (a set of <code>dim</code> components starting
    // at component zero) or the pressure (a scalar component located at
    // position <code>dim</code>):
    const FEValuesExtractors::Vector velocities (0);
    const FEValuesExtractors::Scalar pressure (dim);

    // With all this in place, we can go on with the loop over all cells. The
    // body of this loop has been discussed in the introduction, and will not
    // be commented any further here:
    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
        fe_values.reinit (cell);
        local_matrix = 0;
        local_rhs = 0;

        right_hand_side.value_list (fe_values.get_quadrature_points(),
                                    rhs_values);
        k_inverse.value_list (fe_values.get_quadrature_points(),
                              k_inverse_values);

        for (unsigned int q=0; q<n_q_points; ++q)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
              const Tensor<1,dim> phi_i_u     = fe_values[velocities].value (i, q);
              const double        div_phi_i_u = fe_values[velocities].divergence (i, q);
              const double        phi_i_p     = fe_values[pressure].value (i, q);

              for (unsigned int j=0; j<dofs_per_cell; ++j)
                {
                  const Tensor<1,dim> phi_j_u     = fe_values[velocities].value (j, q);
                  const double        div_phi_j_u = fe_values[velocities].divergence (j, q);
                  const double        phi_j_p     = fe_values[pressure].value (j, q);

                  local_matrix(i,j) += (phi_i_u * k_inverse_values[q] * phi_j_u
                                        - div_phi_i_u * phi_j_p
                                        - phi_i_p * div_phi_j_u)
                                       * fe_values.JxW(q);
                }

              local_rhs(i) += -phi_i_p *
                              rhs_values[q] *
                              fe_values.JxW(q);
            }

        for (unsigned int face_no=0;
             face_no<GeometryInfo<dim>::faces_per_cell;
             ++face_no)
          if (cell->at_boundary(face_no))
            {
              fe_face_values.reinit (cell, face_no);

              pressure_boundary_values
              .value_list (fe_face_values.get_quadrature_points(),
                           boundary_values);

              for (unsigned int q=0; q<n_face_q_points; ++q)
                for (unsigned int i=0; i<dofs_per_cell; ++i)
                  local_rhs(i) += -(fe_face_values[velocities].value (i, q) *
                                    fe_face_values.normal_vector(q) *
                                    boundary_values[q] *
                                    fe_face_values.JxW(q));
            }

        // The final step in the loop over all cells is to transfer local
        // contributions into the global matrix and right hand side
        // vector. Note that we use exactly the same interface as in previous
        // examples, although we now use block matrices and vectors instead of
        // the regular ones. In other words, to the outside world, block
        // objects have the same interface as matrices and vectors, but they
        // additionally allow to access individual blocks.
        cell->get_dof_indices (local_dof_indices);
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            system_matrix.add (local_dof_indices[i],
                               local_dof_indices[j],
                               local_matrix(i,j));
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          system_rhs(local_dof_indices[i]) += local_rhs(i);
      }
  }


  // @sect3{Linear solvers and preconditioners}

  // The linear solvers and preconditioners we use in this example have been
  // discussed in significant detail already in the introduction. We will
  // therefore not discuss the rationale for these classes here any more, but
  // rather only comment on implementational aspects.


  // @sect4{The <code>SchurComplement</code> class template}

  // The next class is the Schur complement class. Its rationale has also been
  // discussed in length in the introduction. The only things we would like to
  // note is that the class, too, is derived from the <code>Subscriptor</code>
  // class and that as mentioned above it stores pointers to the entire block
  // matrix and the inverse of the mass matrix block using
  // <code>SmartPointer</code> objects.
  //
  // The <code>vmult</code> function requires two temporary vectors that we do
  // not want to re-allocate and free every time we call this function. Since
  // here, we have full control over the use of these vectors (unlike above,
  // where a class called by the <code>vmult</code> function required these
  // vectors, not the <code>vmult</code> function itself), we allocate them
  // directly, rather than going through the <code>VectorMemory</code>
  // mechanism. However, again, these member variables do not carry any state
  // between successive calls to the member functions of this class (i.e., we
  // never care what values they were set to the last time a member function
  // was called), we mark these vectors as <code>mutable</code>.
  //
  // The rest of the (short) implementation of this class is straightforward
  // if you know the order of matrix-vector multiplications performed by the
  // <code>vmult</code> function:
  class SchurComplement : public Subscriptor
  {
  public:
    SchurComplement (const BlockSparseMatrix<double> &A,
                     const IterativeInverse<Vector<double> > &Minv);

    void vmult (Vector<double>       &dst,
                const Vector<double> &src) const;

  private:
    const SmartPointer<const BlockSparseMatrix<double> > system_matrix;
    const SmartPointer<const IterativeInverse<Vector<double> > > m_inverse;

    mutable Vector<double> tmp1, tmp2;
  };


  SchurComplement::SchurComplement (const BlockSparseMatrix<double> &A,
                                    const IterativeInverse<Vector<double> > &Minv)
    :
    system_matrix (&A),
    m_inverse (&Minv),
    tmp1 (A.block(0,0).m()),
    tmp2 (A.block(0,0).m())
  {}


  void SchurComplement::vmult (Vector<double>       &dst,
                               const Vector<double> &src) const
  {
    system_matrix->block(0,1).vmult (tmp1, src);
    m_inverse->vmult (tmp2, tmp1);
    system_matrix->block(1,0).vmult (dst, tmp2);
  }


  // @sect4{The <code>ApproximateSchurComplement</code> class template}

  // The third component of our solver and preconditioner system is the class
  // that approximates the Schur complement so we can form a an InverseIterate
  // object that approximates the inverse of the Schur complement. It follows
  // the same pattern as the Schur complement class, with the only exception
  // that we do not multiply with the inverse mass matrix in
  // <code>vmult</code>, but rather just do a single Jacobi
  // step. Consequently, the class also does not have to store a pointer to an
  // inverse mass matrix object.
  //
  // We will later use this class as a template argument to the
  // IterativeInverse class which will in turn want to use it as a
  // template argument for the PointerMatrix class. The latter class
  // has a function that requires us to also write a function that
  // provides the product with the transpose of the matrix this object
  // represents. As a consequence, in the code below, we also
  // implement a <tt>Tvmult</tt> function here that represents the
  // product of the transpose matrix with a vector. It is easy to see
  // how this needs to be implemented here: since the matrix is
  // symmetric, we can as well call <code>vmult</code> wherever the
  // product with the transpose matrix is required. (Note, however,
  // that even though we implement this function here, there will in
  // fact not be any need for it as long as we use SolverCG as the
  // solver since that solver does not ever call the function that
  // provides this operation.)
  class ApproximateSchurComplement : public Subscriptor
  {
  public:
    ApproximateSchurComplement (const BlockSparseMatrix<double> &A);

    void vmult (Vector<double>       &dst,
                const Vector<double> &src) const;
    void Tvmult (Vector<double>       &dst,
                 const Vector<double> &src) const;

  private:
    const SmartPointer<const BlockSparseMatrix<double> > system_matrix;

    mutable Vector<double> tmp1, tmp2;
  };


  ApproximateSchurComplement::ApproximateSchurComplement (const BlockSparseMatrix<double> &A)
    :
    system_matrix (&A),
    tmp1 (A.block(0,0).m()),
    tmp2 (A.block(0,0).m())
  {}


  void ApproximateSchurComplement::vmult (Vector<double>       &dst,
                                          const Vector<double> &src) const
  {
    system_matrix->block(0,1).vmult (tmp1, src);
    system_matrix->block(0,0).precondition_Jacobi (tmp2, tmp1);
    system_matrix->block(1,0).vmult (dst, tmp2);
  }


  void ApproximateSchurComplement::Tvmult (Vector<double>       &dst,
                                           const Vector<double> &src) const
  {
    vmult (dst, src);
  }



  // @sect4{MixedLaplace::solve}

  // After all these preparations, we can finally write the function that
  // actually solves the linear problem. We will go through the two parts it
  // has that each solve one of the two equations, the first one for the
  // pressure (component 1 of the solution), then the velocities (component 0
  // of the solution). Both parts need an object representing the inverse mass
  // matrix and an auxiliary vector, and we therefore declare these objects at
  // the beginning of this function.
  template <int dim>
  void MixedLaplaceProblem<dim>::solve ()
  {
    PreconditionIdentity identity;
    IterativeInverse<Vector<double> > m_inverse;
    m_inverse.initialize(system_matrix.block(0,0), identity);
    m_inverse.solver.select("cg");
    static ReductionControl inner_control(1000, 0., 1.e-13);
    m_inverse.solver.set_control(inner_control);

    Vector<double> tmp (solution.block(0).size());

    // Now on to the first equation. The right hand side of it is $B^TM^{-1}F-G$,
    // which is what we compute in the first few lines. We then declare the
    // objects representing the Schur complement, its approximation, and the
    // inverse of the approximation. Finally, we declare a solver object and
    // hand off all these matrices and vectors to it to compute block 1 (the
    // pressure) of the solution:
    {
      Vector<double> schur_rhs (solution.block(1).size());

      m_inverse.vmult (tmp, system_rhs.block(0));
      system_matrix.block(1,0).vmult (schur_rhs, tmp);
      schur_rhs -= system_rhs.block(1);


      SchurComplement
      schur_complement (system_matrix, m_inverse);

      ApproximateSchurComplement
      approximate_schur_complement (system_matrix);

      IterativeInverse<Vector<double> >
      preconditioner;
      preconditioner.initialize(approximate_schur_complement, identity);
      preconditioner.solver.select("cg");
      preconditioner.solver.set_control(inner_control);


      SolverControl solver_control (solution.block(1).size(),
                                    1e-12*schur_rhs.l2_norm());
      SolverCG<>    cg (solver_control);

      cg.solve (schur_complement, solution.block(1), schur_rhs,
                preconditioner);

      std::cout << solver_control.last_step()
                << " CG Schur complement iterations to obtain convergence."
                << std::endl;
    }

    // After we have the pressure, we can compute the velocity. The equation
    // reads $MU=-BP+F$, and we solve it by first computing the right hand
    // side, and then multiplying it with the object that represents the
    // inverse of the mass matrix:
    {
      system_matrix.block(0,1).vmult (tmp, solution.block(1));
      tmp *= -1;
      tmp += system_rhs.block(0);

      m_inverse.vmult (solution.block(0), tmp);
    }
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
  // <code>VectorTools::integrate_difference</code> function takes as its last
  // argument a pointer to a weight function (the parameter defaults to the
  // null pointer, meaning unit weights). What we simply have to do is to pass
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
  void MixedLaplaceProblem<dim>::compute_errors () const
  {
    const ComponentSelectFunction<dim>
    pressure_mask (dim, dim+1);
    const ComponentSelectFunction<dim>
    velocity_mask(std::make_pair(0, dim), dim+1);

    ExactSolution<dim> exact_solution;
    Vector<double> cellwise_errors (triangulation.n_active_cells());

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
    QTrapez<1>     q_trapez;
    QIterated<dim> quadrature (q_trapez, degree+2);

    // With this, we can then let the library compute the errors and output
    // them to the screen:
    VectorTools::integrate_difference (dof_handler, solution, exact_solution,
                                       cellwise_errors, quadrature,
                                       VectorTools::L2_norm,
                                       &pressure_mask);
    const double p_l2_error = cellwise_errors.l2_norm();

    VectorTools::integrate_difference (dof_handler, solution, exact_solution,
                                       cellwise_errors, quadrature,
                                       VectorTools::L2_norm,
                                       &velocity_mask);
    const double u_l2_error = cellwise_errors.l2_norm();

    std::cout << "Errors: ||e_p||_L2 = " << p_l2_error
              << ",   ||e_u||_L2 = " << u_l2_error
              << std::endl;
  }


  // @sect4{MixedLaplace::output_results}

  // The last interesting function is the one in which we generate graphical
  // output. Everything here looks obvious and familiar. Note how we construct
  // unique names for all the solution variables at the beginning, like we did
  // in step-8 and other programs later on. The only thing worth mentioning is
  // that for higher order elements, in seems inappropriate to only show a
  // single bilinear quadrilateral per cell in the graphical output. We
  // therefore generate patches of size (degree+1)x(degree+1) to capture the
  // full information content of the solution. See the step-7 tutorial program
  // for more information on this.
  //
  // Note that we output the <code>dim+1</code> components of the solution
  // vector as a collection of individual scalars here. Most visualization
  // programs will then only offer to visualize them individually, rather than
  // allowing us to plot the flow field as a vector field. However, as
  // explained in the corresponding function of step-22 or the @ref VVOutput
  // "Generating graphical output" section of the @ref vector_valued module,
  // instructing the DataOut class to identify components of the FESystem
  // object as elements of a <code>dim</code>-dimensional vector is not
  // actually very difficult and will then allow us to show results as vector
  // plots. We skip this here for simplicity and refer to the links above for
  // more information.
  template <int dim>
  void MixedLaplaceProblem<dim>::output_results () const
  {
    std::vector<std::string> solution_names;
    switch (dim)
      {
      case 2:
        solution_names.push_back ("u");
        solution_names.push_back ("v");
        solution_names.push_back ("p");
        break;

      case 3:
        solution_names.push_back ("u");
        solution_names.push_back ("v");
        solution_names.push_back ("w");
        solution_names.push_back ("p");
        break;

      default:
        Assert (false, ExcNotImplemented());
      }


    DataOut<dim> data_out;

    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (solution, solution_names);

    data_out.build_patches (degree+1);

    std::ofstream output ("solution.gmv");
    data_out.write_gmv (output);
  }



  // @sect4{MixedLaplace::run}

  // This is the final function of our main class. It's only job is to call
  // the other functions in their natural order:
  template <int dim>
  void MixedLaplaceProblem<dim>::run ()
  {
    make_grid_and_dofs();
    assemble_system ();
    solve ();
    compute_errors ();
    output_results ();
  }
}


// @sect3{The <code>main</code> function}

// The main function we stole from step-6 instead of step-4. It is almost
// equal to the one in step-6 (apart from the changed class names, of course),
// the only exception is that we pass the degree of the finite element space
// to the constructor of the mixed Laplace problem (here, we use zero-th order
// elements).
int main ()
{
  try
    {
      using namespace dealii;
      using namespace Step20;

      deallog.depth_console (0);

      MixedLaplaceProblem<2> mixed_laplace_problem(0);
      mixed_laplace_problem.run ();
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
