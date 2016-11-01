// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2015 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


#ifndef dealii__matrix_free_operators_h
#define dealii__matrix_free_operators_h


#include <deal.II/base/exceptions.h>
#include <deal.II/base/vectorization.h>
#include <deal.II/base/subscriptor.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/fe_evaluation.h>


DEAL_II_NAMESPACE_OPEN


namespace MatrixFreeOperators
{
  /**
   * Abstract base class for matrix-free operators which can be used both at
   * the finest mesh or at a certain level in geometric multigrid.
   *
   * A derived class has to implement apply_add() method as well as
   * compute_diagonal() to fill the protected member inverse_diagonal_entries.
   * In case of a non-symmetric operator, Tapply_add() should be additionally
   * implemented.
   *
   * @author Denis Davydov, 2016
   */
  template <int dim, typename Number = double>
  class Base : public Subscriptor
  {
  public:
    /**
     * Number typedef.
     */
    typedef Number number_type;

    /**
     * size_type needed for preconditioner classes.
     */
    typedef typename LinearAlgebra::distributed::Vector<Number>::size_type size_type;

    /**
     * Default constructor.
     */
    Base ();

    /**
     * Virtual destructor.
     */
    virtual ~Base();

    /**
     * Release all memory and return to a state just like after having called
     * the default constructor.
     */
    virtual void clear();

    /**
     * Initialize operator on fine scale.
     */
    void initialize (const MatrixFree<dim,Number> &data);

    /**
     * Initialize operator on a level @p level.
     */
    void initialize (const MatrixFree<dim,Number> &data,
                     const MGConstrainedDoFs &mg_constrained_dofs,
                     const unsigned int level);

    /**
     * Return the dimension of the codomain (or range) space.
     */
    size_type m () const;

    /**
     * Return the dimension of the domain space.
     */
    size_type n () const;

    /**
     * vmult operator for interface.
     */
    void vmult_interface_down(LinearAlgebra::distributed::Vector<Number> &dst,
                              const LinearAlgebra::distributed::Vector<Number> &src) const;

    /**
     * vmult operator for interface.
     */
    void vmult_interface_up(LinearAlgebra::distributed::Vector<Number> &dst,
                            const LinearAlgebra::distributed::Vector<Number> &src) const;

    /**
     * Matrix-vector multiplication.
     */
    void vmult (LinearAlgebra::distributed::Vector<Number> &dst,
                const LinearAlgebra::distributed::Vector<Number> &src) const;

    /**
     * Transpose matrix-vector multiplication.
     */
    void Tvmult (LinearAlgebra::distributed::Vector<Number> &dst,
                 const LinearAlgebra::distributed::Vector<Number> &src) const;

    /**
     * Adding Matrix-vector multiplication.
     */
    void vmult_add (LinearAlgebra::distributed::Vector<Number> &dst,
                    const LinearAlgebra::distributed::Vector<Number> &src) const;

    /**
     * Adding transpose matrix-vector multiplication.
     */
    void Tvmult_add (LinearAlgebra::distributed::Vector<Number> &dst,
                     const LinearAlgebra::distributed::Vector<Number> &src) const;

    /**
     * Returns the value of the matrix entry (row,col). In matrix-free context
     * this function is valid only for row==col when diagonal is initialized.
     */
    Number el (const unsigned int row,
               const unsigned int col) const;

    /**
     * Determine an estimate for the memory consumption (in bytes) of this object.
     */
    virtual std::size_t memory_consumption () const;

    /**
     * A wrapper for initialize_dof_vector() of MatrixFree object.
     */
    void initialize_dof_vector (LinearAlgebra::distributed::Vector<Number> &vec) const;

    /**
     * Compute diagonal of this operator.
     *
     * A derived class needs to implement this function and resize and fill
     * the protected member inverse_diagonal_entries accordingly.
     */
    virtual void compute_diagonal () = 0;

    /**
     * Get read access to the inverse diagonal of this operator.
     */
    const LinearAlgebra::distributed::Vector<Number> &get_matrix_diagonal_inverse() const;

    /**
     * Apply the Jacobi preconditioner, which multiplies every element of the
     * <tt>src</tt> vector by the inverse of the respective diagonal element and
     * multiplies the result with the relaxation factor <tt>omega</tt>.
     */
    void precondition_Jacobi(LinearAlgebra::distributed::Vector<Number> &dst,
                             const LinearAlgebra::distributed::Vector<Number> &src,
                             const Number omega) const;

  protected:

    /**
     * Set constrained entries (both from hanging nodes and edge constraints)
     * of @p dst to one.
     */
    void set_constrained_entries_to_one (LinearAlgebra::distributed::Vector<Number> &dst) const;

    /**
     * Apply operator to @p src and add result in @p dst.
     */
    virtual void apply_add(LinearAlgebra::distributed::Vector<Number> &dst,
                           const LinearAlgebra::distributed::Vector<Number> &src) const = 0;

    /**
     * Apply transpose operator to @p src and add result in @p dst.
     *
     * Default implementation is to call apply_add().
     */
    virtual void Tapply_add(LinearAlgebra::distributed::Vector<Number> &dst,
                            const LinearAlgebra::distributed::Vector<Number> &src) const;

    /**
     * MatrixFree object to be used with this operator.
     */
    SmartPointer<const MatrixFree<dim,Number>, Base<dim,Number> > data;

    /**
     * A vector to store inverse of diagonal elements.
     */
    LinearAlgebra::distributed::Vector<Number> inverse_diagonal_entries;

  private:

    /**
     * Indices of DoFs on edge in case the operator is used in GMG context.
     */
    std::vector<unsigned int> edge_constrained_indices;

    /**
     * Auxiliary vector.
     */
    mutable std::vector<std::pair<Number,Number> > edge_constrained_values;

    /**
     * A flag which determines whether or not this operator has interface
     * matrices in GMG context.
     */
    bool have_interface_matrices;

    /**
     * Function which implements vmult_add (@p transpose = false) and
     * Tvmult_add (@p transpose = true).
     */
    void mult_add (LinearAlgebra::distributed::Vector<Number> &dst,
                   const LinearAlgebra::distributed::Vector<Number> &src,
                   const bool transpose) const;
  };



  /**
   * This class implements the operation of the action of the inverse of a
   * mass matrix on an element for the special case of an evaluation object
   * with as many quadrature points as there are cell degrees of freedom. It
   * uses algorithms from FEEvaluation and produces the exact mass matrix for
   * DGQ elements. This algorithm uses tensor products of inverse 1D shape
   * matrices over quadrature points, so the inverse operation is exactly as
   * expensive as applying the forward operator on each cell. Of course, for
   * continuous finite elements this operation does not produce the inverse of
   * a mass operation as the coupling between the elements cannot be
   * considered by this operation.
   *
   * The equation may contain variable coefficients, so the user is required
   * to provide an array for the inverse of the local coefficient (this class
   * provide a helper method 'fill_inverse_JxW_values' to get the inverse of a
   * constant-coefficient operator).
   *
   * @author Martin Kronbichler, 2014
   */
  template <int dim, int fe_degree, int n_components = 1, typename Number = double>
  class CellwiseInverseMassMatrix
  {
  public:
    /**
     * Constructor. Initializes the shape information from the ShapeInfo field
     * in the class FEEval.
     */
    CellwiseInverseMassMatrix (const FEEvaluationBase<dim,n_components,Number> &fe_eval);

    /**
     * Applies the inverse mass matrix operation on an input array. It is
     * assumed that the passed input and output arrays are of correct size,
     * namely FEEval::dofs_per_cell * n_components long. The inverse of the
     * local coefficient (also containing the inverse JxW values) must be
     * passed as first argument. Passing more than one component in the
     * coefficient is allowed.
     */
    void apply(const AlignedVector<VectorizedArray<Number> > &inverse_coefficient,
               const unsigned int             n_actual_components,
               const VectorizedArray<Number> *in_array,
               VectorizedArray<Number>       *out_array) const;

    /**
     * Fills the given array with the inverse of the JxW values, i.e., a mass
     * matrix with coefficient 1. Non-unit coefficients must be multiplied (in
     * inverse form) to this array.
     */
    void fill_inverse_JxW_values(AlignedVector<VectorizedArray<Number> > &inverse_jxw) const;

  private:
    /**
     * A reference to the FEEvaluation object for getting the JxW_values.
     */
    const FEEvaluationBase<dim,n_components,Number> &fe_eval;

    /**
     * A structure to hold inverse shape functions
     */
    AlignedVector<VectorizedArray<Number> > inverse_shape;
  };



  /**
   * This class implements the operation of the action of a mass matrix.
   *
   * @author Daniel Arndt, 2016
   */
  template <int dim, int fe_degree, int n_components = 1, typename Number = double>
  class MassOperator : public Base<dim, Number>
  {
  public:

    /**
     * Constructor.
     */
    MassOperator ();

    /**
     * For preconditioning, we store a lumped mass matrix at the diagonal entries.
     */
    virtual void compute_diagonal ();

  private:
    /**
     * Applies the mass matrix operation on an input vector. It is
     * assumed that the passed input and output vector are correctly initialized
     * using initialize_dof_vector().
     */
    virtual void apply_add (LinearAlgebra::distributed::Vector<Number>       &dst,
                            const LinearAlgebra::distributed::Vector<Number> &src) const;

    /**
     * For this operator, there is just a cell contribution.
     */
    void local_apply_cell (const MatrixFree<dim,Number>                     &data,
                           LinearAlgebra::distributed::Vector<Number>       &dst,
                           const LinearAlgebra::distributed::Vector<Number> &src,
                           const std::pair<unsigned int,unsigned int>  &cell_range) const;
  };



  /**
   * This class implements the operation of the action of a Laplace matrix,
   * namely $ L_{ij} = \int_\Omega c(\mathbf x) \mathbf \nabla N_i(\mathbf x) \cdot \mathbf \nabla N_j(\mathbf x)\,d \mathbf x$,
   * where $c(\mathbf x)$ is the scalar heterogeneity coefficient.
   *
   * @author Denis Davydov, 2016
   */
  template <int dim, int fe_degree, int n_q_points_1d = fe_degree+1, int n_components = 1, typename Number = double>
  class LaplaceOperator : public Base<dim, Number>
  {
  public:

    /**
     * Constructor.
     */
    LaplaceOperator ();

    /**
     * The diagonal is approximated by computing a local diagonal matrix per element
     * and distributing it to the global diagonal. This will lead to wrong results
     * on element with hanging nodes but is still an acceptable approximation
     * to be used in preconditioners.
     */
    virtual void compute_diagonal ();

    /**
     * Set the heterogeneous scalar coefficient @p scalar_coefficient to be used at
     * the quadrature points. The Table should be of correct size, consistent
     * with the total number of quadrature points in <code>dim</code>-dimensions,
     * controlled by the @p n_q_points_1d template parameter. Here,
     * <code>(*scalar_coefficient)(cell,q)</code> corresponds to the value of the
     * coefficient, where <code>cell</code> is an index into a set of cell
     * batches as administered by the MatrixFree framework (which does not work
     * on individual cells, but instead of batches of cells at once), and
     * <code>q</code> is the number of the quadrature point within this batch.
     *
     * Such tables can be initialized by
     * @code
     * std_cxx11::shared_ptr<Table<2, VectorizedArray<double> > > coefficient;
     * coefficient = std_cxx11::make_shared<Table<2, VectorizedArray<double> > >();
     * {
     *   FEEvaluation<dim,fe_degree,n_q_points_1d,1,double> fe_eval(mf_data);
     *   const unsigned int n_cells = mf_data.n_macro_cells();
     *   const unsigned int n_q_points = fe_eval.n_q_points;
     *   coefficient->reinit(n_cells, n_q_points);
     *   for (unsigned int cell=0; cell<n_cells; ++cell)
     *     {
     *       fe_eval.reinit(cell);
     *       for (unsigned int q=0; q<n_q_points; ++q)
     *         (*coefficient)(cell,q) = function.value(fe_eval.quadrature_point(q));
     *     }
     * }
     * @endcode
     * where <code>mf_data</code> is a MatrixFree object and <code>function</code>
     * is a function which provides the following method
     * <code>VectorizedArray<double> value(const Point<dim, VectorizedArray<double> > &p_vec)</code>.
     *
     * If this function is not called, the coefficient is assumed to be unity.
     *
     * The argument to this function is a shared pointer to such a table. The
     * class stores the shared pointer to this table, not a deep copy
     * and uses it to form the Laplace matrix. Consequently, you can update the
     * table and re-use the current object to obtain the action of a Laplace
     * matrix with this updated coefficient. Alternatively, if the table values
     * are only to be filled once, the original shared pointer can also go out
     * of scope in user code and the clear() command or destructor of this class
     * will delete the table.
     */
    void set_coefficient(const std_cxx11::shared_ptr<Table<2, VectorizedArray<Number> > > &scalar_coefficient );

    virtual void clear();

    /**
     * Read/Write access to coefficients to be used in Laplace operator.
     *
     * The function will throw an error if coefficients are not previously set
     * by set_coefficient() function.
     */
    std_cxx11::shared_ptr< Table<2, VectorizedArray<Number> > > get_coefficient();

  private:
    /**
     * Applies the laplace matrix operation on an input vector. It is
     * assumed that the passed input and output vector are correctly initialized
     * using initialize_dof_vector().
     */
    virtual void apply_add (LinearAlgebra::distributed::Vector<Number>       &dst,
                            const LinearAlgebra::distributed::Vector<Number> &src) const;

    /**
     * Applies the Laplace operator on a cell.
     */
    void local_apply_cell (const MatrixFree<dim,Number>                     &data,
                           LinearAlgebra::distributed::Vector<Number>       &dst,
                           const LinearAlgebra::distributed::Vector<Number> &src,
                           const std::pair<unsigned int,unsigned int>  &cell_range) const;

    /**
     * Apply diagonal part of the Laplace operator on a cell.
     */
    void local_diagonal_cell (const MatrixFree<dim,Number>                &data,
                              LinearAlgebra::distributed::Vector<Number>  &dst,
                              const unsigned int &,
                              const std::pair<unsigned int,unsigned int>  &cell_range) const;

    /**
     * Apply Laplace operator on a cell @p cell.
     */
    void do_operation_on_cell(FEEvaluation<dim,fe_degree,n_q_points_1d,n_components,Number> &phi,
                              const unsigned int cell) const;

    /**
     * User-provided heterogeneity coefficient.
     */
    std_cxx11::shared_ptr< Table<2, VectorizedArray<Number> > > scalar_coefficient;
  };



  // ------------------------------------ inline functions ---------------------

  template <int dim, int fe_degree, int n_components, typename Number>
  inline
  CellwiseInverseMassMatrix<dim,fe_degree,n_components,Number>
  ::CellwiseInverseMassMatrix (const FEEvaluationBase<dim,n_components,Number> &fe_eval)
    :
    fe_eval (fe_eval)
  {
    FullMatrix<double> shapes_1d(fe_degree+1, fe_degree+1);
    for (unsigned int i=0, c=0; i<shapes_1d.m(); ++i)
      for (unsigned int j=0; j<shapes_1d.n(); ++j, ++c)
        shapes_1d(i,j) = fe_eval.get_shape_info().shape_values_number[c];
    shapes_1d.gauss_jordan();
    const unsigned int stride = (fe_degree+2)/2;
    inverse_shape.resize(stride*(fe_degree+1));
    for (unsigned int i=0; i<stride; ++i)
      for (unsigned int q=0; q<(fe_degree+2)/2; ++q)
        {
          inverse_shape[i*stride+q] =
            0.5 * (shapes_1d(i,q) + shapes_1d(i,fe_degree-q));
          inverse_shape[(fe_degree-i)*stride+q] =
            0.5 * (shapes_1d(i,q) - shapes_1d(i,fe_degree-q));
        }
    if (fe_degree % 2 == 0)
      for (unsigned int q=0; q<(fe_degree+2)/2; ++q)
        inverse_shape[fe_degree/2*stride+q] = shapes_1d(fe_degree/2,q);
  }



  template <int dim, int fe_degree, int n_components, typename Number>
  inline
  void
  CellwiseInverseMassMatrix<dim,fe_degree,n_components,Number>
  ::fill_inverse_JxW_values(AlignedVector<VectorizedArray<Number> > &inverse_jxw) const
  {
    const unsigned int dofs_per_cell = Utilities::fixed_int_power<fe_degree+1,dim>::value;
    Assert(inverse_jxw.size() > 0 &&
           inverse_jxw.size() % dofs_per_cell == 0,
           ExcMessage("Expected diagonal to be a multiple of scalar dof per cells"));

    // temporarily reduce size of inverse_jxw to dofs_per_cell to get JxW values
    // from fe_eval (will not reallocate any memory)
    const unsigned int previous_size = inverse_jxw.size();
    inverse_jxw.resize(dofs_per_cell);
    fe_eval.fill_JxW_values(inverse_jxw);

    // invert
    inverse_jxw.resize_fast(previous_size);
    for (unsigned int q=0; q<dofs_per_cell; ++q)
      inverse_jxw[q] = 1./inverse_jxw[q];
    // copy values to rest of vector
    for (unsigned int q=dofs_per_cell; q<previous_size; )
      for (unsigned int i=0; i<dofs_per_cell; ++i, ++q)
        inverse_jxw[q] = inverse_jxw[i];
  }



  template <int dim, int fe_degree, int n_components, typename Number>
  inline
  void
  CellwiseInverseMassMatrix<dim,fe_degree,n_components,Number>
  ::apply(const AlignedVector<VectorizedArray<Number> > &inverse_coefficients,
          const unsigned int             n_actual_components,
          const VectorizedArray<Number> *in_array,
          VectorizedArray<Number>       *out_array) const
  {
    const unsigned int dofs_per_cell = Utilities::fixed_int_power<fe_degree+1,dim>::value;
    Assert(inverse_coefficients.size() > 0 &&
           inverse_coefficients.size() % dofs_per_cell == 0,
           ExcMessage("Expected diagonal to be a multiple of scalar dof per cells"));
    if (inverse_coefficients.size() != dofs_per_cell)
      AssertDimension(n_actual_components * dofs_per_cell, inverse_coefficients.size());

    Assert(dim == 2 || dim == 3, ExcNotImplemented());

    internal::EvaluatorTensorProduct<internal::evaluate_evenodd,dim,fe_degree,
             fe_degree+1, VectorizedArray<Number> >
             evaluator(inverse_shape, inverse_shape, inverse_shape);

    const unsigned int shift_coefficient =
      inverse_coefficients.size() > dofs_per_cell ? dofs_per_cell : 0;
    const VectorizedArray<Number> *inv_coefficient = &inverse_coefficients[0];
    VectorizedArray<Number> temp_data_field[dofs_per_cell];
    for (unsigned int d=0; d<n_actual_components; ++d)
      {
        const VectorizedArray<Number> *in = in_array+d*dofs_per_cell;
        VectorizedArray<Number> *out = out_array+d*dofs_per_cell;
        // Need to select 'apply' method with hessian slot because values
        // assume symmetries that do not exist in the inverse shapes
        evaluator.template hessians<0,false,false> (in, temp_data_field);
        evaluator.template hessians<1,false,false> (temp_data_field, out);

        if (dim == 3)
          {
            evaluator.template hessians<2,false,false> (out, temp_data_field);
            for (unsigned int q=0; q<dofs_per_cell; ++q)
              temp_data_field[q] *= inv_coefficient[q];
            evaluator.template hessians<2,true,false> (temp_data_field, out);
          }
        else if (dim == 2)
          for (unsigned int q=0; q<dofs_per_cell; ++q)
            out[q] *= inv_coefficient[q];

        evaluator.template hessians<1,true,false>(out, temp_data_field);
        evaluator.template hessians<0,true,false>(temp_data_field, out);

        inv_coefficient += shift_coefficient;
      }
  }

  //----------------- Base operator -----------------------------
  template <int dim, typename Number>
  Base<dim,Number>::~Base ()
  {
  }



  template <int dim, typename Number>
  Base<dim,Number>::Base ()
    :
    Subscriptor(),
    data(NULL),
    have_interface_matrices(false)
  {
  }



  template <int dim, typename Number>
  typename Base<dim,Number>::size_type
  Base<dim,Number>::m () const
  {
    Assert(data != NULL,
           ExcNotInitialized());
    return data->get_vector_partitioner()->size();
  }



  template <int dim, typename Number>
  typename Base<dim,Number>::size_type
  Base<dim,Number>::n () const
  {
    return m();
  }



  template <int dim, typename Number>
  void
  Base<dim,Number>::clear ()
  {
    data = NULL;
    inverse_diagonal_entries.reinit(0);
  }



  template <int dim, typename Number>
  Number
  Base<dim,Number>::el (const unsigned int row,
                        const unsigned int col) const
  {
    (void) col;
    Assert (row == col, ExcNotImplemented());
    Assert (inverse_diagonal_entries.size() > 0, ExcNotInitialized());
    return 1.0/inverse_diagonal_entries(row);
  }



  template <int dim, typename Number>
  void
  Base<dim,Number>::initialize_dof_vector (LinearAlgebra::distributed::Vector<Number> &vec) const
  {
    Assert(data != NULL,
           ExcNotInitialized());
    if (!vec.partitioners_are_compatible(*data->get_dof_info(0).vector_partitioner))
      data->initialize_dof_vector(vec);
    Assert(vec.partitioners_are_globally_compatible(*data->get_dof_info(0).vector_partitioner),
           ExcInternalError());
  }



  template <int dim, typename Number>
  void
  Base<dim,Number>::
  initialize (const MatrixFree<dim,Number> &data_)
  {
    data =  SmartPointer<const MatrixFree<dim,Number>, Base<dim,Number> >(&data_,typeid(*this).name());
    edge_constrained_indices.clear();
    have_interface_matrices = false;
  }



  template <int dim, typename Number>
  void
  Base<dim,Number>::
  initialize (const MatrixFree<dim,Number> &data_,
              const MGConstrainedDoFs &mg_constrained_dofs,
              const unsigned int level)
  {
    AssertThrow (level != numbers::invalid_unsigned_int,
                 ExcMessage("level is not set"));

    data =  SmartPointer<const MatrixFree<dim,Number>, Base<dim,Number> >(&data_,typeid(*this).name());

    // setup edge_constrained indices
    std::vector<types::global_dof_index> interface_indices;
    mg_constrained_dofs.get_refinement_edge_indices(level).fill_index_vector(interface_indices);
    edge_constrained_indices.clear();
    edge_constrained_indices.reserve(interface_indices.size());
    edge_constrained_values.resize(interface_indices.size());
    const IndexSet &locally_owned = data->get_dof_handler().locally_owned_mg_dofs(level);
    for (unsigned int i=0; i<interface_indices.size(); ++i)
      if (locally_owned.is_element(interface_indices[i]))
        edge_constrained_indices.push_back(locally_owned.index_within_set(interface_indices[i]));
    have_interface_matrices = Utilities::MPI::max((unsigned int)edge_constrained_indices.size(),
                                                  data_.get_vector_partitioner()->get_mpi_communicator()) > 0;
  }



  template <int dim, typename Number>
  void
  Base<dim,Number>::set_constrained_entries_to_one (LinearAlgebra::distributed::Vector<Number> &dst) const
  {
    const std::vector<unsigned int> &
    constrained_dofs = data->get_constrained_dofs();
    for (unsigned int i=0; i<constrained_dofs.size(); ++i)
      dst.local_element(constrained_dofs[i]) = 1.;
    for (unsigned int i=0; i<edge_constrained_indices.size(); ++i)
      dst.local_element(edge_constrained_indices[i]) = 1.;
  }



  template <int dim, typename Number>
  void
  Base<dim,Number>::vmult (LinearAlgebra::distributed::Vector<Number>       &dst,
                           const LinearAlgebra::distributed::Vector<Number> &src) const
  {
    dst = 0;
    vmult_add (dst, src);
  }



  template <int dim, typename Number>
  void
  Base<dim,Number>::vmult_add (LinearAlgebra::distributed::Vector<Number> &dst,
                               const LinearAlgebra::distributed::Vector<Number> &src) const
  {
    mult_add (dst, src, false);
  }



  template <int dim, typename Number>
  void
  Base<dim,Number>::Tvmult_add (LinearAlgebra::distributed::Vector<Number> &dst,
                                const LinearAlgebra::distributed::Vector<Number> &src) const
  {
    mult_add (dst, src, true);
  }



  template <int dim, typename Number>
  void
  Base<dim,Number>::mult_add (LinearAlgebra::distributed::Vector<Number> &dst,
                              const LinearAlgebra::distributed::Vector<Number> &src,
                              const bool transpose) const
  {
    Assert(src.partitioners_are_globally_compatible(*data->get_dof_info(0).vector_partitioner), ExcInternalError());
    Assert(dst.partitioners_are_globally_compatible(*data->get_dof_info(0).vector_partitioner), ExcInternalError());

    // set zero Dirichlet values on the input vector (and remember the src and
    // dst values because we need to reset them at the end)
    for (unsigned int i=0; i<edge_constrained_indices.size(); ++i)
      {
        edge_constrained_values[i] =
          std::pair<Number,Number>(src.local_element(edge_constrained_indices[i]),
                                   dst.local_element(edge_constrained_indices[i]));
        const_cast<LinearAlgebra::distributed::Vector<Number>&>(src).local_element(edge_constrained_indices[i]) = 0.;
      }

    if (transpose)
      Tapply_add(dst,src);
    else
      apply_add(dst,src);

    const std::vector<unsigned int> &
    constrained_dofs = data->get_constrained_dofs();
    for (unsigned int i=0; i<constrained_dofs.size(); ++i)
      dst.local_element(constrained_dofs[i]) += src.local_element(constrained_dofs[i]);

    // reset edge constrained values, multiply by unit matrix and add into
    // destination
    for (unsigned int i=0; i<edge_constrained_indices.size(); ++i)
      {
        const_cast<LinearAlgebra::distributed::Vector<Number>&>(src).local_element(edge_constrained_indices[i]) = edge_constrained_values[i].first;
        dst.local_element(edge_constrained_indices[i]) = edge_constrained_values[i].second + edge_constrained_values[i].first;
      }
  }



  template <int dim, typename Number>
  void
  Base<dim,Number>::
  vmult_interface_down(LinearAlgebra::distributed::Vector<Number> &dst,
                       const LinearAlgebra::distributed::Vector<Number> &src) const
  {
    Assert(src.partitioners_are_globally_compatible(*data->get_dof_info(0).vector_partitioner), ExcInternalError());
    Assert(dst.partitioners_are_globally_compatible(*data->get_dof_info(0).vector_partitioner), ExcInternalError());

    dst = 0;

    if (!have_interface_matrices)
      return;

    // set zero Dirichlet values on the input vector (and remember the src and
    // dst values because we need to reset them at the end)
    for (unsigned int i=0; i<edge_constrained_indices.size(); ++i)
      {
        edge_constrained_values[i] =
          std::pair<Number,Number>(src.local_element(edge_constrained_indices[i]),
                                   dst.local_element(edge_constrained_indices[i]));
        const_cast<LinearAlgebra::distributed::Vector<Number>&>(src).local_element(edge_constrained_indices[i]) = 0.;
      }

    apply_add(dst,src);

    unsigned int c=0;
    for (unsigned int i=0; i<edge_constrained_indices.size(); ++i)
      {
        for ( ; c<edge_constrained_indices[i]; ++c)
          dst.local_element(c) = 0.;
        ++c;

        // reset the src values
        const_cast<LinearAlgebra::distributed::Vector<Number>&>(src).local_element(edge_constrained_indices[i]) = edge_constrained_values[i].first;
      }
    for ( ; c<dst.local_size(); ++c)
      dst.local_element(c) = 0.;
  }



  template <int dim, typename Number>
  void
  Base<dim,Number>::
  vmult_interface_up(LinearAlgebra::distributed::Vector<Number> &dst,
                     const LinearAlgebra::distributed::Vector<Number> &src) const
  {
    Assert(src.partitioners_are_globally_compatible(*data->get_dof_info(0).vector_partitioner), ExcInternalError());
    Assert(dst.partitioners_are_globally_compatible(*data->get_dof_info(0).vector_partitioner), ExcInternalError());

    dst = 0;

    if (!have_interface_matrices)
      return;

    LinearAlgebra::distributed::Vector<Number> src_cpy (src);
    unsigned int c=0;
    for (unsigned int i=0; i<edge_constrained_indices.size(); ++i)
      {
        for ( ; c<edge_constrained_indices[i]; ++c)
          src_cpy.local_element(c) = 0.;
        ++c;
      }
    for ( ; c<src_cpy.local_size(); ++c)
      src_cpy.local_element(c) = 0.;

    apply_add(dst,src_cpy);

    for (unsigned int i=0; i<edge_constrained_indices.size(); ++i)
      {
        dst.local_element(edge_constrained_indices[i]) = 0.;
      }
  }



  template <int dim, typename Number>
  void
  Base<dim,Number>::Tvmult (LinearAlgebra::distributed::Vector<Number>       &dst,
                            const LinearAlgebra::distributed::Vector<Number> &src) const
  {
    dst = 0;
    Tvmult_add (dst,src);
  }



  template <int dim, typename Number>
  std::size_t
  Base<dim,Number>::memory_consumption () const
  {
    return inverse_diagonal_entries.memory_consumption();
  }



  template <int dim, typename Number>
  const LinearAlgebra::distributed::Vector<Number> &
  Base<dim,Number>::get_matrix_diagonal_inverse() const
  {
    Assert(inverse_diagonal_entries.size() > 0, ExcNotInitialized());
    return inverse_diagonal_entries;
  }



  template <int dim, typename Number>
  void
  Base<dim,Number>::Tapply_add(LinearAlgebra::distributed::Vector<Number> &dst,
                               const LinearAlgebra::distributed::Vector<Number> &src) const
  {
    apply_add(dst,src);
  }



  template <int dim, typename Number>
  void
  Base<dim,Number>::precondition_Jacobi(LinearAlgebra::distributed::Vector<Number> &dst,
                                        const LinearAlgebra::distributed::Vector<Number> &src,
                                        const Number omega) const
  {
    Assert(inverse_diagonal_entries.size() > 0, ExcNotInitialized());

    dst = src;
    dst.scale(inverse_diagonal_entries);
    dst*= omega;
  }



  //-----------------------------MassOperator----------------------------------

  template <int dim, int fe_degree, int n_components, typename Number>
  MassOperator<dim, fe_degree, n_components, Number>::
  MassOperator ()
    :
    Base<dim, Number>()
  {}



  template <int dim, int fe_degree, int n_components, typename Number>
  void
  MassOperator<dim, fe_degree, n_components, Number>::
  compute_diagonal()
  {
    Assert((Base<dim, Number>::data != NULL), ExcNotInitialized());

    LinearAlgebra::distributed::Vector<Number> ones;
    Base<dim, Number>::initialize_dof_vector(Base<dim, Number>::inverse_diagonal_entries);
    Base<dim, Number>::initialize_dof_vector(ones);
    ones = 1.;
    ones.update_ghost_values();
    apply_add(Base<dim, Number>::inverse_diagonal_entries, ones);

    const std::vector<unsigned int> &constrained_dofs
      = Base<dim, Number>::data->get_constrained_dofs();
    for (unsigned int i=0; i< constrained_dofs.size(); ++i)
      Base<dim, Number>::inverse_diagonal_entries.local_element(constrained_dofs[i]) = 1.;

    const unsigned int local_size = Base<dim, Number>::inverse_diagonal_entries.local_size();
    for (unsigned int i=0; i<local_size; ++i)
      Base<dim, Number>::inverse_diagonal_entries.local_element(i)
        =1./Base<dim, Number>::inverse_diagonal_entries.local_element(i);

    Base<dim, Number>::inverse_diagonal_entries.compress(VectorOperation::insert);
    Base<dim, Number>::inverse_diagonal_entries.update_ghost_values();
  }



  template <int dim, int fe_degree, int n_components, typename Number>
  void
  MassOperator<dim, fe_degree, n_components, Number>::
  apply_add (LinearAlgebra::distributed::Vector<Number>       &dst,
             const LinearAlgebra::distributed::Vector<Number> &src) const
  {
    Base<dim, Number>::data->cell_loop (&MassOperator::local_apply_cell,
                                        this, dst, src);
  }



  template <int dim, int fe_degree, int n_components, typename Number>
  void
  MassOperator<dim, fe_degree, n_components, Number>::
  local_apply_cell (const MatrixFree<dim,Number>                     &data,
                    LinearAlgebra::distributed::Vector<Number>       &dst,
                    const LinearAlgebra::distributed::Vector<Number> &src,
                    const std::pair<unsigned int,unsigned int>  &cell_range) const
  {
    FEEvaluation<dim, fe_degree, fe_degree+1, n_components, Number> phi(data);
    for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
      {
        phi.reinit (cell);
        phi.read_dof_values(src);
        phi.evaluate (true,false,false);
        for (unsigned int q=0; q<phi.n_q_points; ++q)
          phi.submit_value (phi.get_value(q), q);
        phi.integrate (true,false);
        phi.distribute_local_to_global (dst);
      }
  }


  //-----------------------------LaplaceOperator----------------------------------

  template <int dim, int fe_degree, int n_q_points_1d, int n_components, typename Number>
  LaplaceOperator<dim, fe_degree, n_q_points_1d, n_components, Number>::
  LaplaceOperator ()
    :
    Base<dim, Number>()
  {
  }



  template <int dim, int fe_degree, int n_q_points_1d, int n_components, typename Number>
  void
  LaplaceOperator<dim, fe_degree, n_q_points_1d, n_components, Number>::
  clear ()
  {
    Base<dim, Number>::clear();
    scalar_coefficient = NULL;
  }



  template <int dim, int fe_degree, int n_q_points_1d, int n_components, typename Number>
  void
  LaplaceOperator<dim, fe_degree, n_q_points_1d, n_components, Number>::
  set_coefficient(const std_cxx11::shared_ptr<Table<2, VectorizedArray<Number> > > &scalar_coefficient_ )
  {
    scalar_coefficient = scalar_coefficient_;
  }



  template <int dim, int fe_degree, int n_q_points_1d, int n_components, typename Number>
  std_cxx11::shared_ptr< Table<2, VectorizedArray<Number> > >
  LaplaceOperator<dim, fe_degree, n_q_points_1d, n_components, Number>::
  get_coefficient()
  {
    Assert (scalar_coefficient.get(),
            ExcNotInitialized());
    return scalar_coefficient;
  }



  template <int dim, int fe_degree, int n_q_points_1d, int n_components, typename Number>
  void
  LaplaceOperator<dim, fe_degree, n_q_points_1d, n_components, Number>::
  compute_diagonal()
  {
    Assert((Base<dim, Number>::data != NULL), ExcNotInitialized());

    unsigned int dummy = 0;
    LinearAlgebra::distributed::Vector<Number> &inverse_diagonal_entries = Base<dim,Number>::inverse_diagonal_entries;
    this->initialize_dof_vector(inverse_diagonal_entries);
    Base<dim,Number>::
    data->cell_loop (&LaplaceOperator::local_diagonal_cell,
                     this, inverse_diagonal_entries, dummy);

    this->set_constrained_entries_to_one(inverse_diagonal_entries);

    for (unsigned int i=0; i<inverse_diagonal_entries.local_size(); ++i)
      if (std::abs(inverse_diagonal_entries.local_element(i)) > std::sqrt(std::numeric_limits<Number>::epsilon()))
        inverse_diagonal_entries.local_element(i) = 1./inverse_diagonal_entries.local_element(i);
      else
        inverse_diagonal_entries.local_element(i) = 1.;

    Base<dim, Number>::inverse_diagonal_entries.compress(VectorOperation::insert);
    Base<dim, Number>::inverse_diagonal_entries.update_ghost_values();
  }



  template <int dim, int fe_degree, int n_q_points_1d, int n_components, typename Number>
  void
  LaplaceOperator<dim, fe_degree, n_q_points_1d, n_components, Number>::
  apply_add (LinearAlgebra::distributed::Vector<Number>       &dst,
             const LinearAlgebra::distributed::Vector<Number> &src) const
  {
    Base<dim, Number>::data->cell_loop (&LaplaceOperator::local_apply_cell,
                                        this, dst, src);
  }

  namespace
  {
    template<typename Number>
    bool
    non_negative(const VectorizedArray<Number> &n)
    {
      for (unsigned int v=0; v<VectorizedArray<Number>::n_array_elements; ++v)
        if (n[v] < 0.)
          return false;

      return true;
    }
  }



  template <int dim, int fe_degree, int n_q_points_1d, int n_components, typename Number>
  void
  LaplaceOperator<dim, fe_degree, n_q_points_1d, n_components, Number>::
  do_operation_on_cell(FEEvaluation<dim,fe_degree,n_q_points_1d,n_components,Number> &phi,
                       const unsigned int cell) const
  {
    phi.evaluate (false,true,false);
    if (scalar_coefficient.get())
      {
        for (unsigned int q=0; q<phi.n_q_points; ++q)
          {
            Assert (non_negative((*scalar_coefficient)(cell,q)),
                    ExcMessage("Coefficient must be non-negative"));
            phi.submit_gradient ((*scalar_coefficient)(cell,q)*phi.get_gradient(q), q);
          }
      }
    else
      {
        for (unsigned int q=0; q<phi.n_q_points; ++q)
          {
            phi.submit_gradient (phi.get_gradient(q), q);
          }
      }
    phi.integrate (false,true);
  }




  template <int dim, int fe_degree, int n_q_points_1d, int n_components, typename Number>
  void
  LaplaceOperator<dim, fe_degree, n_q_points_1d, n_components, Number>::
  local_apply_cell (const MatrixFree<dim,Number>                     &data,
                    LinearAlgebra::distributed::Vector<Number>       &dst,
                    const LinearAlgebra::distributed::Vector<Number> &src,
                    const std::pair<unsigned int,unsigned int>  &cell_range) const
  {
    FEEvaluation<dim,fe_degree,n_q_points_1d,n_components,Number> phi (data);
    for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
      {
        phi.reinit (cell);
        phi.read_dof_values(src);
        do_operation_on_cell(phi,cell);
        phi.distribute_local_to_global (dst);
      }
  }


  template <int dim, int fe_degree, int n_q_points_1d, int n_components, typename Number>
  void
  LaplaceOperator<dim, fe_degree, n_q_points_1d, n_components, Number>::
  local_diagonal_cell (const MatrixFree<dim,Number>                     &data,
                       LinearAlgebra::distributed::Vector<Number>       &dst,
                       const unsigned int &,
                       const std::pair<unsigned int,unsigned int>       &cell_range) const
  {
    FEEvaluation<dim,fe_degree,n_q_points_1d,n_components,Number> phi (data);
    for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
      {
        phi.reinit (cell);
        VectorizedArray<Number> local_diagonal_vector[phi.tensor_dofs_per_cell];
        for (unsigned int i=0; i<phi.dofs_per_cell; ++i)
          {
            for (unsigned int j=0; j<phi.dofs_per_cell; ++j)
              phi.begin_dof_values()[j] = VectorizedArray<Number>();
            phi.begin_dof_values()[i] = 1.;
            do_operation_on_cell(phi,cell);
            local_diagonal_vector[i] = phi.begin_dof_values()[i];
          }
        for (unsigned int i=0; i<phi.tensor_dofs_per_cell; ++i)
          phi.begin_dof_values()[i] = local_diagonal_vector[i];
        phi.distribute_local_to_global (dst);
      }
  }


} // end of namespace MatrixFreeOperators


DEAL_II_NAMESPACE_CLOSE

#endif
