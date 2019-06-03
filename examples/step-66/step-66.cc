/* ---------------------------------------------------------------------
 * 
 * Copyright (C) 2003 - 2018 by the deal.II authors
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
 * Authors: Fabian Castelli, Karlsruhe Institute of Technology (KIT)
 *          Timo Heister, University of Utah
 */



// First we include the typical headers of the deal.II library needed for this
// tutroial:
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_c1.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>

// In particular, we need to include the headers for the matrix-free framework:
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>
#include <deal.II/matrix_free/fe_evaluation.h>

// And since we want to use a geometric multigrid preconditioner, we need also
// the multilevel headers:
#include <deal.II/multigrid/multigrid.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_constrained_dofs.h>


// Finally some common C++ headers for in and output:
#include <iostream>
#include <fstream>



namespace Step66
{
  using namespace dealii;
  
  
  
  // As in the previous tutorials on the matrix-free framework we define the
  // degree of the finite element space and the space dimension as constants
  // variables. For this example we solve a tow dimensional problem and use the
  // fourth order Lagrangian finite element space.
  const unsigned int degree_finite_element = 4;
  const unsigned int dimension             = 2;
  
  
  
  // sect3{Matrix-free operators}
  
  // In the beginning we define the two matrix-free operator for the Jacobian
  // and the residual. As guideline we follow the tutorials step-37 and step-48,
  // where the precise interface of the MatrixFreeOperators::Base class was
  // extensively documented.
  
  
  // @setc4{JacobianOperator}
  
  // Since we want to use the Jacobian as system matrix and pass it to the
  // linear solver as well as to the multilevel preconditioner classes we derive
  // the JacobianOperator class from the MatrixFreeOperators::Base class, such
  // that we have already the right interface. The two functions we need to
  // reimplement from the base class are the apply_add and the compute_diagonal
  // function. To allow preconditioning with float precision we define the
  // number type as template argument.
  //
  // As mentioned already in the introduction, for the Jacobian $F'$ of the
  // $n+1$-th Newton step we need to evaluate it at the last Newton step
  // $u_h^n$. To get the information of the last Newton step $u_h^n$ we do
  // pretty much the same as we in step-37, where we stored the values of a
  // coefficient function in a table once before we use the matrix-free
  // operator. Instead of a function evaluate_coefficient, we here implement a
  // function evaluate_newton_step.
  //
  // As private member functions of the JacobianOperator we implement the
  // local_apply and the local_compute_diagonal function, which we call in the
  // cell_loop for computing the matrix-vector product or the diagonal.
  template <int dim, int fe_degree, typename number>
  class JacobianOperator : public MatrixFreeOperators::Base<dim,LinearAlgebra::distributed::Vector<number> >
  {
  public:
    using value_type = number;
    
    JacobianOperator();
    
    virtual
    void
    clear() override;
    
    void
    evaluate_newton_step(const LinearAlgebra::distributed::Vector<number> &src);
    
    virtual
    void
    compute_diagonal() override;
    
  private:
    virtual
    void
    apply_add(LinearAlgebra::distributed::Vector<number> &dst, const LinearAlgebra::distributed::Vector<number> &src) const override;
    
    void
    local_apply(const MatrixFree<dim,number> &data, LinearAlgebra::distributed::Vector<number> &dst, const LinearAlgebra::distributed::Vector<number> &src, const std::pair<unsigned int,unsigned int> &cell_range) const;
    
    void
    local_compute_diagonal(const MatrixFree<dim,number> &data, LinearAlgebra::distributed::Vector<number> &dst, const unsigned int &dummy, const std::pair<unsigned int,unsigned int> &cell_range) const;
    
    Table<2, VectorizedArray<number> > nonlinear_values;
  };
  
  
  
  // The constructor of the JacobainOperator just calls the constructor of the
  // base class MatrixFreeOperators::Base, which is itself derived from the
  // Subscriptor class.
  template <int dim, int fe_degree, typename number>
  JacobianOperator<dim,fe_degree,number>::JacobianOperator()
  :
  MatrixFreeOperators::Base<dim,LinearAlgebra::distributed::Vector<number> >()
  {}
  
  
  
  // The clear function resets out the table holding the values for the
  // nonlinearity and call the clear function of the base class.
  template <int dim, int fe_degree, typename number>
  void
  JacobianOperator<dim,fe_degree,number>::clear()
  {
    nonlinear_values.reinit(0, 0);
    MatrixFreeOperators::Base<dim,LinearAlgebra::distributed::Vector<number> >::clear();
  }
  
  
  
  // The following evlauate_newton_step function is based on the
  // evlauate_coefficient function from step-37, however, it does not evaluate a
  // function object, it evaluates a vector representing a finite element
  // function, namely the last Newton step needed for the Jacobian. Therefore we
  // set up a FEEvaluateion object and evaluate the finite element function in
  // the quadrature points with the gather_evaluate function. Remember if we
  // only loop over cells this function is just a wrapper around functions
  // read_dof_values and evaluate. Since we store the evaluated values of the
  // finite element function in a table we do not have to call integrate in
  // combination with distribute_local_to_global or integrate scatter.
  //
  // This will work well and in the local_apply function we can use the values
  // stored in the table to apply the matrix-vector product, however, we can
  // also optimize the implementation of the Jacobian at this stage. We can
  // directly evaluate the nonlinear function std::exp(u_h^n[q]) and store these
  // values in the table. This saves in each call of the vmult function all
  // evaluations of the nonlinearity.
  template <int dim, int fe_degree, typename number>
  void
  JacobianOperator<dim,fe_degree,number>::evaluate_newton_step(const LinearAlgebra::distributed::Vector<number> &src)
  {
    const unsigned int n_cells = this->data->n_macro_cells();
    FEEvaluation<dim,fe_degree,fe_degree+1,1,number> phi(*this->data);
    
    nonlinear_values.reinit(n_cells, phi.n_q_points);
    for(unsigned int cell=0; cell<n_cells; ++cell)
    {
      phi.reinit(cell);
      phi.gather_evaluate(src, true, false);
      
      for(unsigned int q=0; q<phi.n_q_points; ++q)
        nonlinear_values(cell, q) = std::exp(phi.get_value(q));
    }
  }
  
  
  
  // Now in the local_apply function, which actually implements the local action
  // of the system matrix on a cell, we can use the information of the last
  // Newton step stored in the table nonlinear_values. The rest of this function
  // as basically the same in step-37. We setup the FEEvaluation object, gather
  // and evaluate the values and gradients of the input vector, submit the
  // values and gradients according to the form of the Jacobian and finally call
  // integrate_scatter to distribute the local contributions into the global
  // vector dst.
  template <int dim, int fe_degree, typename number>
  void
  JacobianOperator<dim,fe_degree,number>::local_apply(const MatrixFree<dim,number> &data, LinearAlgebra::distributed::Vector<number> &dst, const LinearAlgebra::distributed::Vector<number> &src, const std::pair<unsigned int,unsigned int> &cell_range) const
  {
    FEEvaluation<dim,fe_degree,fe_degree+1,1,number> phi(data);
    
    for(unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
    {
      AssertDimension(nonlinear_values.size(0), data.n_macro_cells());
      AssertDimension(nonlinear_values.size(1), phi.n_q_points);
      
      phi.reinit(cell);
      phi.gather_evaluate(src, true, true);
      
      for(unsigned int q=0; q<phi.n_q_points; ++q)
      {
        phi.submit_value(-nonlinear_values(cell,q)*phi.get_value(q), q);
        phi.submit_gradient(phi.get_gradient(q), q);
      }
      
      phi.integrate_scatter(true, true, dst);
    }
  }
  
  
  
  // Next we use cell_loop of the MatrixFree class to perform the actual loop
  // over all cells computing the cell contribution to the matrix-vector
  // product.
  template <int dim, int fe_degree, typename number>
  void
  JacobianOperator<dim,fe_degree,number>::apply_add(LinearAlgebra::distributed::Vector<number> &dst, const LinearAlgebra::distributed::Vector<number> &src) const
  {
    this->data->cell_loop(&JacobianOperator::local_apply, this, dst, src, true);
  }
  
  
  
  // The missing two functions of the JacobianOperator calculate the diagonal
  // entries of the Jacobian. The only difference compared to step-37 is the
  // calculation of the cell contribution in the local_compute_diagonal
  // function. However, therefore we only have to extend and change the
  // arguments for the submit functions in the loop over all quadrature points
  // and this can be done according to the local_apply function. So no further
  // comments to these two functions schould be necessary.
  template <int dim, int fe_degree, typename number>
  void
  JacobianOperator<dim,fe_degree,number>::compute_diagonal()
  {
    this->inverse_diagonal_entries.reset(
      new DiagonalMatrix<LinearAlgebra::distributed::Vector<number> >());
    LinearAlgebra::distributed::Vector<number> &inverse_diagonal =
    this->inverse_diagonal_entries->get_vector();
    this->data->initialize_dof_vector(inverse_diagonal);
    
    unsigned int dummy = 0;
    
    this->data->cell_loop(&JacobianOperator::local_compute_diagonal,
                          this,
                          inverse_diagonal,
                          dummy);
    
    this->set_constrained_entries_to_one(inverse_diagonal);
    
    for(unsigned int i=0; i<inverse_diagonal.local_size(); ++i)
    {
      Assert(inverse_diagonal.local_element(i) > 0.,
             ExcMessage("No diagonal entry in a positive definite operator "
             "should be zero"));
      inverse_diagonal.local_element(i) =
      1./inverse_diagonal.local_element(i);
    }
  }
  
  
  
  template <int dim, int fe_degree, typename number>
  void
  JacobianOperator<dim,fe_degree,number>::local_compute_diagonal(const MatrixFree<dim,number> &data, LinearAlgebra::distributed::Vector<number> &dst, const unsigned int &, const std::pair<unsigned int,unsigned int> &cell_range) const
  {
    FEEvaluation<dim,fe_degree,fe_degree+1,1,number> phi(data);
    
    AlignedVector<VectorizedArray<number> > diagonal(phi.dofs_per_cell);
    
    for(unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
    {
      AssertDimension(nonlinear_values.size(0), data.n_macro_cells());
      AssertDimension(nonlinear_values.size(1), phi.n_q_points);
      
      phi.reinit(cell);
      for(unsigned int i=0; i<phi.dofs_per_cell; ++i)
      {
        for(unsigned int j=0; j<phi.dofs_per_cell; ++j)
          phi.submit_dof_value(VectorizedArray<number>(), j);
        phi.submit_dof_value(make_vectorized_array<number>(1.), i);
        
        phi.evaluate(true, true);
        for(unsigned int q=0; q<phi.n_q_points; ++q)
        {
          phi.submit_value(-nonlinear_values(cell, q)*phi.get_value(q), q);
          phi.submit_gradient(phi.get_gradient(q), q);
        }
        phi.integrate(true, true);
        diagonal[i] = phi.get_dof_value(i);
      }
      for(unsigned int i=0; i<phi.dofs_per_cell; ++i)
        phi.submit_dof_value(diagonal[i], i);
      phi.distribute_local_to_global(dst);
    }
  }

