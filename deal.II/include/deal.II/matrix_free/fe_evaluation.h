//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2011, 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#ifndef __deal2__matrix_free_fe_evaluation_h
#define __deal2__matrix_free_fe_evaluation_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/vectorization.h>
#include <deal.II/matrix_free/matrix_free.h>


DEAL_II_NAMESPACE_OPEN

namespace parallel
{
  namespace distributed
  {
    template <typename> class Vector;
  }
}



namespace internal
{
  DeclException0 (ExcAccessToUninitializedField);
}



/**
 * This is the base class for the FEEvaluation classes. This class is a base
 * class and needs usually not be called in user code. Use one of the derived
 * classes instead. It implements access functions to vectors for the @p
 * read_dof_values, @p set_dof_values, and @p distributed_local_to_global
 * functions, as well as the @p reinit method.
 *
 * This class has five template arguments:
 *
 * @param dim Dimension in which this class is to be used
 *
 * @param dofs_per_cell Number of degrees of freedom of the FE per cell,
 *                  usually (fe_degree+1)^dim for elements based on a tensor
 *                  product
 *
 * @param n_q_points Number of points in the quadrature formula, usually
 *                  (fe_degree+1)^dim for tensor-product quadrature formulas
 *
 * @param n_components Number of vector components when solving a system of
 *                  PDEs. If the same operation is applied to several
 *                  components of a PDE (e.g. a vector Laplace equation), they
 *                  can be applied simultaneously with one call (and often
 *                  more efficiently)
 *
 * @param Number Number format, usually @p double or @p float
 *
 * @author Katharina Kormann and Martin Kronbichler, 2010, 2011
 */
template <int dim, int dofs_per_cell_, int n_q_points_,
          int n_components, typename Number>
class FEEvaluationBase
{
public:
  typedef VectorizedArray<Number> vector_t;
  static const std::size_t  n_vectors =
    VectorizedArray<Number>::n_array_elements;
  static const unsigned int dofs_per_cell = dofs_per_cell_;
  static const unsigned int n_q_points    = n_q_points_;

                                /**
                                 * Constructor. Takes all data stored in
                                 * MatrixFree. If applied to problems with
                                 * more than one finite element or more than
                                 * one quadrature formula selected during
                                 * construction of @p matrix_free, @p
                                 * fe_no and @p quad_no allow to select the
                                 * appropriate components.
                                 */
  FEEvaluationBase (const MatrixFree<dim,Number> &matrix_free,
                    const unsigned int                fe_no   = 0,
                    const unsigned int                quad_no = 0);

                                /**
                                 * Initializes the operation pointer to the
                                 * current cell. Unlike the FEValues::reinit
                                 * function, where the information related to
                                 * a particular cell is generated in the
                                 * reinit call, this function is very cheap
                                 * since all data is pre-computed in @p
                                 * matrix_free, and only a few indices
                                 * have to be set appropriately.
                                 */
  void reinit (const unsigned int cell);

                                /**
                                 * For the transformation information stored
                                 * in MappingInfo, this function returns the
                                 * index which belongs to the current cell as
                                 * specified in @p reinit. Note that
                                 * MappingInfo has different fields for
                                 * Cartesian cells, cells with linear mapping
                                 * and with general mappings, so in order to
                                 * access the correct data, this interface
                                 * must be used together with get_cell_type.
                                 */
  unsigned int get_cell_data_number() const;

                                /**
                                 * Returns the type of the cell the @p reinit
                                 * function has been called for. 0 means
                                 * Cartesian cells (which allows for
                                 * considerable data compression), 1 means
                                 * cells with linear mappings, and 2 means
                                 * general cells without any compressed
                                 * storage applied.
                                 */
  unsigned int get_cell_type() const;

                                /**
                                 * Returns a read-only pointer to the first
                                 * field of function values on quadrature
                                 * points. First come the function values on
                                 * all quadrature points for the first
                                 * component, then all values for the second
                                 * component, and so on. This is related to
                                 * the internal data structures used in this
                                 * class. The raw data after a call to @p
                                 * evaluate only contains unit cell
                                 * operations, so possible transformations,
                                 * quadrature weights etc. must be applied
                                 * manually. In general, it is safer to use
                                 * the get_value() function instead, which
                                 * does all the transformation internally.
                                 */
  const vector_t * begin_values () const;

                                /**
                                 * Returns a read and write pointer to the
                                 * first field of function values on
                                 * quadrature points. First come the function
                                 * values on all quadrature points for the
                                 * first component, then all values for the
                                 * second component, and so on. This is
                                 * related to the internal data structures
                                 * used in this class. The raw data after a
                                 * call to @p evaluate only contains unit
                                 * cell operations, so possible
                                 * transformations, quadrature weights
                                 * etc. must be applied manually. In general,
                                 * it is safer to use the get_value() function
                                 * instead, which does all the transformation
                                 * internally.
                                 */
  vector_t * begin_values ();

                                /**
                                 * Returns a read-only pointer to the first
                                 * field of function gradients on quadrature
                                 * points. First comes the x-component of the
                                 * gradient for the first component on all
                                 * quadrature points, then the y-component,
                                 * and so on. Next comes the x-component of
                                 * the second component, and so on. This is
                                 * related to the internal data structures
                                 * used in this class. The raw data after a
                                 * call to @p evaluate only contains unit
                                 * cell operations, so possible
                                 * transformations, quadrature weights
                                 * etc. must be applied manually. In general,
                                 * it is safer to use the get_gradient() function
                                 * instead, which does all the transformation
                                 * internally.
                                 */
  const vector_t * begin_gradients () const;

                                /**
                                 * Returns a read and write pointer to the
                                 * first field of function gradients on
                                 * quadrature points. First comes the
                                 * x-component of the gradient for the first
                                 * component on all quadrature points, then
                                 * the y-component, and so on. Next comes the
                                 * x-component of the second component, and so
                                 * on. This is related to the internal data
                                 * structures used in this class. The raw data
                                 * after a call to @p evaluate only
                                 * contains unit cell operations, so possible
                                 * transformations, quadrature weights
                                 * etc. must be applied manually. In general,
                                 * it is safer to use the get_gradient()
                                 * function instead, which does all the
                                 * transformation internally.
                                 */
  vector_t * begin_gradients ();

                                /**
                                 * Returns a read-only pointer to the first
                                 * field of function hessians on quadrature
                                 * points. First comes the xx-component of the
                                 * hessian for the first component on all
                                 * quadrature points, then the yy-component,
                                 * zz-component in (3D), then the
                                 * xy-component, and so on. Next comes the
                                 * xx-component of the second component, and
                                 * so on. This is related to the internal data
                                 * structures used in this class. The raw data
                                 * after a call to @p evaluate only
                                 * contains unit cell operations, so possible
                                 * transformations, quadrature weights
                                 * etc. must be applied manually. In general,
                                 * it is safer to use the get_laplacian() or
                                 * get_hessian() functions instead, which does
                                 * all the transformation internally.
                                 */
  const vector_t * begin_hessians () const;

                                /**
                                 * Returns a read and write pointer to the
                                 * first field of function hessians on
                                 * quadrature points. First comes the
                                 * xx-component of the hessian for the first
                                 * component on all quadrature points, then
                                 * the yy-component, zz-component in (3D),
                                 * then the xy-component, and so on. Next
                                 * comes the xx-component of the second
                                 * component, and so on. This is related to
                                 * the internal data structures used in this
                                 * class. The raw data after a call to @p
                                 * evaluate only contains unit cell
                                 * operations, so possible transformations,
                                 * quadrature weights etc. must be applied
                                 * manually. In general, it is safer to use
                                 * the get_laplacian() or get_hessian()
                                 * functions instead, which does all the
                                 * transformation internally.
                                 */
  vector_t * begin_hessians ();

                                /**
                                 * For the vector @p src, read out the values
                                 * on the degrees of freedom of the current
                                 * cell, and store them internally. Similar
                                 * functionality as the function
                                 * DoFAccessor::get_interpolated_dof_values
                                 * when no constraints are present, but it
                                 * also includes constraints from hanging
                                 * nodes, so one can see it as a similar
                                 * function to
                                 * ConstraintMatrix::read_dof_values as
                                 * well. Note that if vectorization is
                                 * enabled, the DoF values for several cells
                                 * are set.
                                 *
                                 * If some constraints on the vector are
                                 * inhomogeneous, use the function
                                 * read_dof_values_plain instead and provide
                                 * the vector with useful data also in
                                 * constrained positions by calling
                                 * ConstraintMatrix::distribute. When
                                 * accessing vector entries during the
                                 * solution of linear systems, the temporary
                                 * solution should always have homogeneous
                                 * constraints and this method is the correct
                                 * one.
                                 */
  template <typename VectorType>
  void read_dof_values (const VectorType &src);

                                /**
                                 * For a collection of several vector @p src,
                                 * read out the values on the degrees of
                                 * freedom of the current cell for @p
                                 * n_components (template argument), starting
                                 * at @p first_index, and store them
                                 * internally. Similar functionality as the
                                 * function ConstraintMatrix::read_dof_values.
                                 * Note that if vectorization is enabled, the
                                 * DoF values for several cells are set.
                                 */
  template <typename VectorType>
  void read_dof_values (const std::vector<VectorType> &src,
                       const unsigned int             first_index=0);

                                /**
                                 * Reads data from several vectors. Same as
                                 * other function with std::vector, but
                                 * accepts a vector of pointers to vectors.
                                 */
  template <typename VectorType>
  void read_dof_values (const std::vector<VectorType*> &src,
                       const unsigned int              first_index=0);

                                /**
                                 * For a collection of several vector @p src,
                                 * read out the values on the degrees of
                                 * freedom of the current cell for @p
                                 * n_components (template argument), and store
                                 * them internally. Similar functionality as
                                 * the function
                                 * ConstraintMatrix::read_dof_values. Note
                                 * that if vectorization is enabled, the DoF
                                 * values for several cells are set.
                                 */
  template<typename VectorType>
  void read_dof_values (const VectorType * src_data[]);

                                /**
                                 * For the vector @p src, read out the values
                                 * on the degrees of freedom of the current
                                 * cell, and store them internally. Similar
                                 * functionality as the function
                                 * DoFAccessor::get_interpolated_dof_values. As
                                 * opposed to the read_dof_values function,
                                 * this function reads out the plain entries
                                 * from vectors, without taking stored
                                 * constraints into account. This way of
                                 * access is appropriate when the constraints
                                 * have been distributed on the vector by a
                                 * call to ConstraintMatrix::distribute
                                 * previously. This function is also necessary
                                 * when inhomogeneous constraints are to be
                                 * used, as MatrixFree can only handle
                                 * homogeneous constraints. Note that if
                                 * vectorization is enabled, the DoF values
                                 * for several cells are set.
                                 */
  template <typename VectorType>
  void read_dof_values_plain (const VectorType &src);

                                /**
                                 * For a collection of several vector @p src,
                                 * read out the values on the degrees of
                                 * freedom of the current cell for @p
                                 * n_components (template argument), starting
                                 * at @p first_index, and store them
                                 * internally. Similar functionality as the
                                 * function DoFAccessor::read_dof_values.
                                 * Note that if vectorization is enabled, the
                                 * DoF values for several cells are set.
                                 */
  template <typename VectorType>
  void read_dof_values_plain (const std::vector<VectorType> &src,
                             const unsigned int             first_index=0);

                                /**
                                 * Reads data from several vectors without
                                 * resolving constraints. Same as other
                                 * function with std::vector, but accepts a
                                 * vector of pointers to vectors.
                                 */
  template <typename VectorType>
  void read_dof_values_plain (const std::vector<VectorType*> &src,
                             const unsigned int              first_index=0);

                                /**
                                 * For a collection of several vector @p src,
                                 * read out the values on the degrees of
                                 * freedom of the current cell for @p
                                 * n_components (template argument), and store
                                 * them internally. Similar functionality as
                                 * the function
                                 * DoFAccessor::read_dof_values. Note
                                 * that if vectorization is enabled, the DoF
                                 * values for several cells are set.
                                 */
  template<typename VectorType>
  void read_dof_values_plain (const VectorType * src_data[]);

                                /**
                                 * Takes the values stored internally on dof
                                 * values of the current cell and sums them
                                 * into the vector @p dst. The function also
                                 * applies constraints during the write
                                 * operation. The functionality is hence
                                 * similar to the function
                                 * ConstraintMatrix::distribute_local_to_global.
                                 * Note that if vectorization is enabled, the
                                 * DoF values for several cells are used.
                                 */
  template<typename VectorType>
  void distribute_local_to_global (VectorType &dst) const;

                                /**
                                 * Takes the values stored internally on dof
                                 * values of the current cell for a
                                 * vector-valued problem consisting of @p
                                 * n_components (template argument) and sums
                                 * them into the collection of vectors vector
                                 * @p dst, starting at index @p
                                 * first_index. The function also applies
                                 * constraints during the write operation. The
                                 * functionality is hence similar to the
                                 * function
                                 * ConstraintMatrix::distribute_local_to_global.
                                 * Note that if vectorization is enabled, the
                                 * DoF values for several cells are used.
                                 */
  template<typename VectorType>
  void distribute_local_to_global (std::vector<VectorType> &dst,
                                   const unsigned int       first_index=0) const;

                                /**
                                 * Writes data to several vectors. Same as
                                 * other function with std::vector, but
                                 * accepts a vector of pointers to vectors.
                                 */
  template<typename VectorType>
  void distribute_local_to_global (std::vector<VectorType*> &dst,
                                   const unsigned int       first_index=0) const;

                                /**
                                 * Takes the values stored internally on dof
                                 * values of the current cell for a
                                 * vector-valued problem consisting of @p
                                 * n_components (template argument) and sums
                                 * them into the collection of vectors vector
                                 * @p dst. The function also applies
                                 * constraints during the write operation. The
                                 * functionality is hence similar to the
                                 * function
                                 * ConstraintMatrix::distribute_local_to_global.
                                 * Note that if vectorization is enabled, the
                                 * DoF values for several cells are used.
                                 */
  template<typename VectorType>
  void distribute_local_to_global (VectorType * dst_data[]) const;

                                /**
                                 * Takes the values stored internally on dof
                                 * values of the current cell and sums them
                                 * into the vector @p dst. The function also
                                 * applies constraints during the write
                                 * operation. The functionality is hence
                                 * similar to the function
                                 * ConstraintMatrix::distribute_local_to_global.
                                 * Note that if vectorization is enabled, the
                                 * DoF values for several cells are used.
                                 */
  template<typename VectorType>
  void set_dof_values (VectorType &dst) const;

                                /**
                                 * Takes the values stored internally on dof
                                 * values of the current cell for a
                                 * vector-valued problem consisting of @p
                                 * n_components (template argument) and sums
                                 * them into the collection of vectors vector
                                 * @p dst, starting at index @p
                                 * first_index. The function also applies
                                 * constraints during the write operation. The
                                 * functionality is hence similar to the
                                 * function
                                 * ConstraintMatrix::distribute_local_to_global.
                                 * Note that if vectorization is enabled, the
                                 * DoF values for several cells are used.
                                 */
  template<typename VectorType>
  void set_dof_values (std::vector<VectorType> &dst,
                       const unsigned int       first_index=0) const;

                                /**
                                 * Writes data to several vectors. Same as
                                 * other function with std::vector, but
                                 * accepts a vector of pointers to vectors.
                                 */
  template<typename VectorType>
  void set_dof_values (std::vector<VectorType*> &dst,
                       const unsigned int        first_index=0) const;

                                /**
                                 * Takes the values stored internally on dof
                                 * values of the current cell for a
                                 * vector-valued problem consisting of @p
                                 * n_components (template argument) and sums
                                 * them into the collection of vectors vector
                                 * @p dst. The function also applies
                                 * constraints during the write operation. The
                                 * functionality is hence similar to the
                                 * function
                                 * ConstraintMatrix::distribute_local_to_global.
                                 * Note that if vectorization is enabled, the
                                 * DoF values for several cells are used.
                                 */
  template<typename VectorType>
  void set_dof_values (VectorType * dst_data[]) const;

                                /**
                                 * Returns the value stored for the local
                                 * degree of freedom with index @p dof. If the
                                 * object is vector-valued, a vector-valued
                                 * return argument is given. Note that when
                                 * vectorization is enabled, values from
                                 * several cells are grouped together. If @p
                                 * set_dof_values was called last, the value
                                 * corresponds to the one set there. If @p
                                 * integrate was called last, it instead
                                 * corresponds to the value of the integrated
                                 * function with the test function of the
                                 * given index.
                                 */
  Tensor<1,n_components,vector_t>
  get_dof_value (unsigned int dof) const;

                                /**
                                 * Write a value to the field containing the
                                 * degrees of freedom with component @p
                                 * dof. Access to the same field as through @p
                                 * get_dof_value.
                                 */
  void submit_dof_value (Tensor<1,n_components,vector_t> val_in,
                         unsigned int dof);

                                /**
                                 * Returns the value of a finite
                                 * element function at quadrature
                                 * point number @p q_point after a
                                 * call to @p evaluate(true,...), or
                                 * the value that has been stored
                                 * there with a call to @p
                                 * submit_value. If the object is
                                 * vector-valued, a vector-valued
                                 * return argument is given. Note that
                                 * when vectorization is enabled,
                                 * values from several cells are
                                 * grouped together.
                                 */
  Tensor<1,n_components,vector_t>
  get_value (unsigned int q_point) const;

                                /**
                                 * Write a value to the field containing the
                                 * values on quadrature points with component
                                 * @p q_point. Access to the same field as
                                 * through @p get_value. If applied before the
                                 * function @p integrate(true,...) is
                                 * called, this specifies the value which is
                                 * tested by all basis function on the current
                                 * cell and integrated over.
                                 */
  void submit_value (Tensor<1,n_components,vector_t> val_in,
                     unsigned int q_point);

                                /**
                                 * Returns the gradient of a finite element
                                 * function at quadrature point number @p
                                 * q_point after a call to @p
                                 * evaluate(...,true,...), or the value
                                 * that has been stored there with a call to
                                 * @p submit_gradient.
                                 */
  Tensor<1,n_components,Tensor<1,dim,vector_t> >
  get_gradient (unsigned int q_point) const;

                                /**
                                 * Write a gradient to the field containing
                                 * the values on quadrature points with
                                 * component @p q_point. Access to the same
                                 * field as through @p get_gradient. If
                                 * applied before the function @p
                                 * integrate(...,true) is called,
                                 * this specifies the gradient which is tested
                                 * by all basis function gradients on the
                                 * current cell and integrated over.
                                 */
  void submit_gradient(Tensor<1,n_components,Tensor<1,dim,vector_t> >grad_in,
                       unsigned int q_point);

                                /**
                                 * Returns the Hessian of a finite element
                                 * function at quadrature point number @p
                                 * q_point after a call to @p
                                 * evaluate(...,true). If only the
                                 * diagonal or even the trace of the Hessian,
                                 * the Laplacian, is needed, use the other
                                 * functions below.
                                 */
  Tensor<1,n_components,Tensor<2,dim,vector_t> >
  get_hessian (unsigned int q_point) const;

                                /**
                                 * Returns the diagonal of the Hessian of a
                                 * finite element function at quadrature point
                                 * number @p q_point after a call to @p
                                 * evaluate(...,true).
                                 */
  Tensor<1,n_components,Tensor<1,dim,vector_t> >
  get_hessian_diagonal (unsigned int q_point) const;

                                /**
                                 * Returns the Laplacian of a finite element
                                 * function at quadrature point number @p
                                 * q_point after a call to @p
                                 * evaluate(...,true).
                                 */
  Tensor<1,n_components,vector_t>
  get_laplacian (unsigned int q_point) const;

                                /**
                                 * Takes values on quadrature points,
                                 * multiplies by the Jacobian determinant and
                                 * quadrature weights (JxW) and sums the
                                 * values for all quadrature points on the
                                 * cell. The result is a scalar, representing
                                 * the integral over the function over the
                                 * cell. If a vector-element is used, the
                                 * resulting components are still
                                 * separated. Moreover, if vectorization is
                                 * enabled, the integral values of several
                                 * cells are represented together.
                                 */
  Tensor<1,n_components,vector_t>
  integrate_value ();

                                /**
                                 * Stores a reference to the underlying data.
                                 */
  const MatrixFree<dim,Number>   &matrix_info;

                                /**
                                 * Stores a reference to the underlying DoF
                                 * indices and constraint description for the
                                 * component specified at construction. Also
                                 * contained in matrix_info, but it simplifies
                                 * code if we store a reference to it.
                                 */
  const internal::MatrixFreeFunctions::DoFInfo      &dof_info;

                                /**
                                 * Stores the constraints weights that
                                 * supplement DoFInfo. Also contained in
                                 * matrix_info, but it simplifies code if we
                                 * store a reference to it.
                                 */
  const internal::MatrixFreeFunctions::CompressedMatrix<Number> &constraint_pool;

                                /**
                                 * Stores a reference to the underlying
                                 * transformation data from unit to real cells
                                 * for the given quadrature formula specified
                                 * at construction.  Also contained in
                                 * matrix_info, but it simplifies code if we
                                 * store a reference to it.
                                 */
  const internal::MatrixFreeFunctions::MappingInfo<dim,Number> &mapping_info;

                                /**
                                 * Stores the active fe index for this class
                                 * for efficient indexing in the hp case.
                                 */
  const unsigned int active_fe_index;

                                /**
                                 * Stores the active quadrature index for this
                                 * class for efficient indexing in the hp
                                 * case.
                                 */
  const unsigned int active_quad_index;

                                /**
                                 * Stores a reference to the unit cell data,
                                 * i.e., values, gradients and Hessians in 1D
                                 * at the quadrature points that constitute
                                 * the tensor product. Also contained in
                                 * matrix_info, but it simplifies code if we
                                 * store a reference to it.
                                 */
  const internal::MatrixFreeFunctions::FEEvaluationData<Number> &data;

protected:
                                /**
                                 * Internal data fields that store the
                                 * values. Since all array lengths are known
                                 * at compile time and since they are rarely
                                 * more than a few kilobytes, allocate them on
                                 * the stack. This makes it possible to
                                 * cheaply set up a FEEvaluation object and
                                 * write thread-safe programs by letting each
                                 * thread own a private object of this type.
                                 */
  vector_t values_dofs[n_components][dofs_per_cell>0?dofs_per_cell:1];
  vector_t values_quad[n_components][n_q_points>0?n_q_points:1];
  vector_t gradients_quad[n_components][dim][n_q_points>0?n_q_points:1];
  vector_t hessians_quad[n_components][(dim*(dim+1))/2][n_q_points>0?n_q_points:1];

                                /**
                                 * Stores the indices of the current cell.
                                 */
  const unsigned int quad_no;
  const unsigned int n_fe_components;
  unsigned int cell;
  unsigned int cell_type;
  unsigned int cell_data_number;
  bool         at_irregular_cell;
  unsigned int n_irreg_components_filled;

                                /**
                                 * A pointer to the Cartesian Jacobian
                                 * information of the present cell. Only set
                                 * to a useful value if on a Cartesian cell,
                                 * otherwise zero.
                                 */
  const Tensor<1,dim,vector_t> * cartesian;

                                /**
                                 * A pointer to the Jacobian information of
                                 * the present cell. Only set to a useful
                                 * value if on a non-Cartesian cell.
                                 */
  const Tensor<2,dim,vector_t> * jacobian;

                                /**
                                 * A pointer to the Jacobian determinant of
                                 * the present cell. If on a Cartesian cell or
                                 * on a cell with constant Jacobian, this is
                                 * just the Jacobian determinant, otherwise
                                 * the Jacobian determinant times the
                                 * quadrature weight.
                                 */
  const vector_t * J_value;

                                /**
                                 * A pointer to the quadrature weights of the
                                 * underlying quadrature formula.
                                 */
  const vector_t * quadrature_weights;

                                /**
                                 * A pointer to the quadrature points on the
                                 * present cell.
                                 */
  const Point<dim,vector_t> * quadrature_points;

                                /**
                                 * A pointer to the diagonal part of the
                                 * Jacobian gradient on the present
                                 * cell. Only set to a useful value if on a
                                 * general cell with non-constant Jacobian.
                                 */
  const Tensor<2,dim,vector_t> * jacobian_grad;

                                /**
                                 * A pointer to the upper diagonal part of the
                                 * Jacobian gradient on the present cell. Only
                                 * set to a useful value if on a general cell
                                 * with non-constant Jacobian.
                                 */
  const Tensor<1,(dim>1?dim*(dim-1)/2:1),Tensor<1,dim,vector_t> > * jacobian_grad_upper;

                                /**
                                 * Debug information to track whether we
                                 * uninitialized fields are accessed.
                                 */
  bool     dof_values_initialized;
  bool     values_quad_initialized;
  bool     gradients_quad_initialized;
  bool     hessians_quad_initialized;
  bool     values_quad_submitted;
  bool     gradients_quad_submitted;
};



/**
 * This class provides access to the data fields of the FEEvaluation
 * classes. Generic access is achieved through the base class, and
 * specializations for scalar and vector-valued elements are defined
 * separately.
 *
 * @author Katharina Kormann and Martin Kronbichler, 2010, 2011
 */
template <int dim, int dofs_per_cell_, int n_q_points_,
          int n_components, typename Number>
class FEEvaluationAccess :
  public FEEvaluationBase<dim,dofs_per_cell_,n_q_points_,n_components,Number>
{
 public:
  typedef VectorizedArray<Number> vector_t;
  typedef Tensor<1,n_components,vector_t> value_type;
  typedef Tensor<1,n_components,Tensor<1,dim,vector_t> > gradient_type;
  static const std::size_t  n_vectors = VectorizedArray<Number>::n_array_elements;
  static const unsigned int dofs_per_cell = dofs_per_cell_;
  static const unsigned int n_q_points    = n_q_points_;
  typedef FEEvaluationBase<dim,dofs_per_cell_,n_q_points_,n_components,
                           Number> BaseClass;

                                /**
                                 * Constructor. Takes all data stored in
                                 * MatrixFree. If applied to problems with
                                 * more than one finite element or more than
                                 * one quadrature formula selected during
                                 * construction of @p matrix_free, @p
                                 * fe_no and @p quad_no allow to select the
                                 * appropriate components.
                                 */
  FEEvaluationAccess (const MatrixFree<dim,Number> &matrix_free,
                      const unsigned int                fe_no   = 0,
                      const unsigned int                quad_no = 0);
};




/**
 * This class provides access to the data fields of the FEEvaluation
 * classes. Partial specialization for scalar fields that defines access with
 * simple data fields, i.e., scalars for the values and Tensor<1,dim> for the
 * gradients.
 *
 * @author Katharina Kormann and Martin Kronbichler, 2010, 2011
 */
template <int dim, int dofs_per_cell_, int n_q_points_, typename Number>
class FEEvaluationAccess<dim,dofs_per_cell_,n_q_points_,1,Number> :
  public FEEvaluationBase<dim,dofs_per_cell_,n_q_points_,1,Number>
{
 public:
  typedef Number                            number_type;
  typedef VectorizedArray<Number> vector_t;
  typedef VectorizedArray<Number> value_type;
  typedef Tensor<1,dim,vector_t>            gradient_type;
  static const unsigned int dimension = dim;
  static const std::size_t  n_vectors = VectorizedArray<Number>::n_array_elements;
  static const unsigned int dofs_per_cell = dofs_per_cell_;
  static const unsigned int n_q_points    = n_q_points_;
  typedef FEEvaluationBase<dim,dofs_per_cell_,n_q_points_,1,Number> BaseClass;

                                /**
                                 * Constructor. Takes all data stored in
                                 * MatrixFree. If applied to problems with
                                 * more than one finite element or more than
                                 * one quadrature formula selected during
                                 * construction of @p matrix_free, @p
                                 * fe_no and @p quad_no allow to select the
                                 * appropriate components.
                                 */
  FEEvaluationAccess (const MatrixFree<dim,Number> &matrix_free,
                      const unsigned int                fe_no   = 0,
                      const unsigned int                quad_no = 0);

                                /**
                                 * Returns the value stored for the local
                                 * degree of freedom with index @p dof. If the
                                 * object is vector-valued, a vector-valued
                                 * return argument is given. Note that when
                                 * vectorization is enabled, values from
                                 * several cells are grouped together. If @p
                                 * set_dof_values was called last, the value
                                 * corresponds to the one set there. If @p
                                 * integrate was called last, it instead
                                 * corresponds to the value of the integrated
                                 * function with the test function of the
                                 * given index.
                                 */
  vector_t
  get_dof_value (unsigned int dof) const;

                                /**
                                 * Write a value to the field containing the
                                 * degrees of freedom with component @p
                                 * dof. Access to the same field as through @p
                                 * get_dof_value.
                                 */
  void submit_dof_value (vector_t     val_in,
                         unsigned int dof);

                                /**
                                 * Returns the value of a finite element
                                 * function at quadrature point number @p
                                 * q_point after a call to @p
                                 * evaluate(true,...), or the value that
                                 * has been stored there with a call to @p
                                 * submit_value. If the object is
                                 * vector-valued, a vector-valued return
                                 * argument is given. Note that when
                                 * vectorization is enabled, values from
                                 * several cells are grouped together.
                                 */
  vector_t
  get_value (unsigned int q_point) const;

                                /**
                                 * Write a value to the field
                                 * containing the values on quadrature
                                 * points with component @p
                                 * q_point. Access to the same field
                                 * as through @p get_value. If applied
                                 * before the function @p
                                 * integrate(true,...) is called, this
                                 * specifies the value which is tested
                                 * by all basis function on the
                                 * current cell and integrated over.
                                 */
  void submit_value (vector_t     val_in,
                     unsigned int q_point);

                                /**
                                 * Returns the gradient of a finite
                                 * element function at quadrature
                                 * point number @p q_point after a
                                 * call to @p evaluate(...,true,...),
                                 * or the value that has been stored
                                 * there with a call to @p
                                 * submit_gradient.
                                 */
  gradient_type
  get_gradient (unsigned int q_point) const;

                                /**
                                 * Write a gradient to the field
                                 * containing the values on quadrature
                                 * points with component @p
                                 * q_point. Access to the same field
                                 * as through @p get_gradient. If
                                 * applied before the function @p
                                 * integrate(...,true) is called, this
                                 * specifies the gradient which is
                                 * tested by all basis function
                                 * gradients on the current cell and
                                 * integrated over.
                                 */
  void submit_gradient(gradient_type grad_in,
                       unsigned int  q_point);

                                /**
                                 * Returns the Hessian of a finite
                                 * element function at quadrature
                                 * point number @p q_point after a
                                 * call to @p evaluate(...,true). If
                                 * only the diagonal part of the
                                 * Hessian or its trace, the
                                 * Laplacian, are needed, use the
                                 * respective functions below.
                                 */
  Tensor<2,dim,vector_t>
  get_hessian (unsigned int q_point) const;

                                /**
                                 * Returns the diagonal of the Hessian
                                 * of a finite element function at
                                 * quadrature point number @p q_point
                                 * after a call to @p
                                 * evaluate(...,true).
                                 */
  gradient_type
  get_hessian_diagonal (unsigned int q_point) const;

                                /**
                                 * Returns the Laplacian of a finite
                                 * element function at quadrature
                                 * point number @p q_point after a
                                 * call to @p evaluate(...,true).
                                 */
  value_type
  get_laplacian (unsigned int q_point) const;

                                /**
                                 * Takes values on quadrature points,
                                 * multiplies by the Jacobian determinant and
                                 * quadrature weights (JxW) and sums the
                                 * values for all quadrature points on the
                                 * cell. The result is a scalar, representing
                                 * the integral over the function over the
                                 * cell. If a vector-element is used, the
                                 * resulting components are still
                                 * separated. Moreover, if vectorization is
                                 * enabled, the integral values of several
                                 * cells are represented together.
                                 */
  value_type
  integrate_value ();
};



/**
 * This class provides access to the data fields of the FEEvaluation
 * classes. Partial specialization for fields with as many components as the
 * underlying space dimension, i.e., values are of type Tensor<1,dim> and
 * gradients of type Tensor<2,dim>. Provides some additional functions for
 * access, like the symmetric gradient and divergence.
 *
 * @author Katharina Kormann and Martin Kronbichler, 2010, 2011
 */
template <int dim, int dofs_per_cell_, int n_q_points_, typename Number>
class FEEvaluationAccess<dim,dofs_per_cell_,n_q_points_,dim,Number> :
  public FEEvaluationBase<dim,dofs_per_cell_,n_q_points_,dim,Number>
{
 public:
  typedef VectorizedArray<Number> vector_t;
  typedef Tensor<1,dim,vector_t>            value_type;
  typedef Tensor<2,dim,vector_t>            gradient_type;
  typedef SymmetricTensor<2,dim,vector_t>   sym_gradient_type;
  typedef Tensor<1,dim==2?1:dim,vector_t>   curl_type;

  static const std::size_t  n_vectors = VectorizedArray<Number>::n_array_elements;
  static const unsigned int dofs_per_cell = dofs_per_cell_;
  static const unsigned int n_q_points    = n_q_points_;
  typedef FEEvaluationBase<dim,dofs_per_cell_,n_q_points_,dim,Number> BaseClass;

                                /**
                                 * Constructor. Takes all data stored in
                                 * MatrixFree. If applied to problems with
                                 * more than one finite element or more than
                                 * one quadrature formula selected during
                                 * construction of @p matrix_free, @p
                                 * fe_no and @p quad_no allow to select the
                                 * appropriate components.
                                 */
  FEEvaluationAccess (const MatrixFree<dim,Number> &matrix_free,
                      const unsigned int                fe_no   = 0,
                      const unsigned int                quad_no = 0);

                                /**
                                 * Returns the gradient of a finite element
                                 * function at quadrature point number @p
                                 * q_point after a call to @p
                                 * evaluate(...,true,...).
                                 */
  gradient_type
  get_gradient (unsigned int q_point) const;

                                /**
                                 * Returns the divergence of a vector-valued
                                 * finite element at quadrature point number
                                 * @p q_point after a call to @p
                                 * evaluate(...,true,...).
                                 */
  vector_t
  get_divergence (unsigned int q_point) const;

                                /**
                                 * Returns the symmetric gradient of a
                                 * vector-valued finite element at
                                 * quadrature point number @p q_point
                                 * after a call to @p
                                 * evaluate(...,true,...). It
                                 * corresponds to <tt>0.5
                                 * (grad+grad<sup>T</sup>)</tt>.
                                 */
  sym_gradient_type
  get_symmetric_gradient (unsigned int q_point) const;

                                /**
                                 * Returns the curl of the vector field,
                                 * $nabla \times v$ after a call to @p
                                 * evaluate(...,true,...).
                                 */
  curl_type
  get_curl (unsigned int q_point) const;

                                /**
                                 * Returns the Hessian of a finite
                                 * element function at quadrature
                                 * point number @p q_point after a
                                 * call to @p evaluate(...,true). If
                                 * only the diagonal of the Hessian or
                                 * its trace, the Laplacian, is
                                 * needed, use the respective
                                 * functions.
                                 */
  Tensor<3,dim,vector_t>
  get_hessian (unsigned int q_point) const;

                                /**
                                 * Returns the diagonal of the Hessian
                                 * of a finite element function at
                                 * quadrature point number @p q_point
                                 * after a call to @p
                                 * evaluate(...,true).
                                 */
  gradient_type
  get_hessian_diagonal (unsigned int q_point) const;

                                /**
                                 * Write a gradient to the field containing
                                 * the values on quadrature points with
                                 * component @p q_point. Access to the same
                                 * field as through @p get_gradient. If
                                 * applied before the function @p
                                 * integrate(...,true) is called,
                                 * this specifies the gradient which is tested
                                 * by all basis function gradients on the
                                 * current cell and integrated over.
                                 */
  void submit_gradient(gradient_type grad_in,
                       unsigned int  q_point);

                                /**
                                 * Write a gradient to the field containing
                                 * the values on quadrature points with
                                 * component @p q_point. This function is an
                                 * alternative to the other submit_gradient
                                 * function when using a system of fixed
                                 * number of equations which happens to
                                 * coincide with the dimension for some
                                 * dimensions, but not all. To allow for
                                 * dimension-independent programming, this
                                 * function can be used instead.
                                 */
  void submit_gradient(Tensor<1,dim,Tensor<1,dim,vector_t> > grad_in,
                       unsigned int                          q_point);

                                /**
                                 * Write a gradient to the field containing
                                 * the values on quadrature points with
                                 * component @p q_point. Access to the same
                                 * field as through @p get_gradient. If
                                 * applied before the function @p
                                 * integrate(...,true) is called,
                                 * this specifies the gradient which is tested
                                 * by all basis function gradients on the
                                 * current cell and integrated over.
                                 */
  void submit_symmetric_gradient(sym_gradient_type grad_in,
                                 unsigned int      q_point);

                                /**
                                 * Write the components of a curl containing
                                 * the values on quadrature point @p
                                 * q_point. Access to the same data field as
                                 * through @p get_gradient.
                                 */
  void submit_curl (curl_type    curl_in,
                    unsigned int q_point);
};



/**
 * The class that provides all functions necessary to evaluate functions at
 * quadrature points and cell integrations. In functionality, this class is
 * similar to FEValues<dim>, however, it includes a lot of specialized
 * functions that make it much faster (between 5 and 500 times as fast,
 * depending on the polynomial order). Access to the data fields is provided
 * through functionality in the class FEEvaluationAccess.
 *
 * This class is designed for general local finite element operations based on
 * tensor products of 1D polynomials and quadrature points. Often, there are
 * some symmetries or zeros in the unit cell data that allow for a more
 * efficient operator application. FEEvaluation is specialized to standard
 * FE_Q/FE_DGQ elements and quadrature points symmetric around 0.5 (like Gauss
 * quadrature), and hence the most common situation. FEEvaluationGL is a
 * specialization for elements where quadrature formula and support points are
 * chosen so that a orthogonal relation between the values holds. This is the
 * case for FE_Q based on Gauss-Lobatto support points with Gauss-Lobatto
 * quadrature formula of the same order.
 *
 * This class has five template arguments:
 *
 * @param dim Dimension in which this class is to be used
 *
 * @param n_dofs_1d Number of degrees of freedom of the FE in 1D, usually
 *                   fe_degree+1, for elements based on a tensor product
 *
 * @param n_q_points_1d Number of points in the quadrature formula in 1D,
 *                   usually chosen as fe_degree+1
 *
 * @param n_components Number of vector components when solving a system of
 *                  PDEs. If the same operation is applied to several
 *                  components of a PDE (e.g. a vector Laplace equation), they
 *                  can be applied simultaneously with one call (and often
 *                  more efficiently)
 *
 * @param Number Number format, usually @p double or @p float
 *
 * @author Katharina Kormann and Martin Kronbichler, 2010, 2011
 */
template <int dim, int n_dofs_1d, int n_q_points_1d=n_dofs_1d,
          int n_components=1, typename Number=double >
class FEEvaluationGeneral :
  public FEEvaluationAccess<dim,
                            (n_dofs_1d*(dim>1?n_dofs_1d:1)*(dim>2?n_dofs_1d:1)),
                            (n_q_points_1d*(dim>1?n_q_points_1d:1)*(dim>2?n_q_points_1d:1)),
                            n_components,Number>
{
 public:
  typedef VectorizedArray<Number> vector_t;
  static const std::size_t  n_vectors = VectorizedArray<Number>::n_array_elements;
  typedef FEEvaluationAccess<dim,(n_dofs_1d*(dim>1?n_dofs_1d:1)*
                                  (dim>2?n_dofs_1d:1)),
    (n_q_points_1d*(dim>1?n_q_points_1d:1)*(dim>2?n_q_points_1d:1)),
    n_components, Number> BaseClass;
  static const unsigned int dofs_per_cell = BaseClass::dofs_per_cell;
  static const unsigned int n_q_points    = BaseClass::n_q_points;

                                /**
                                 * Constructor. Takes all data stored in
                                 * MatrixFree. If applied to problems with
                                 * more than one finite element or more than
                                 * one quadrature formula selected during
                                 * construction of @p matrix_free, @p
                                 * fe_no and @p quad_no allow to select the
                                 * appropriate components.
                                 */
  FEEvaluationGeneral (const MatrixFree<dim,Number> &matrix_free,
                       const unsigned int                fe_no   = 0,
                       const unsigned int                quad_no = 0);

                                /**
                                 * Evaluates the function values, the
                                 * gradients, and the Laplacians of the FE
                                 * function given at the DoF values in the
                                 * input vector at the quadrature points.  The
                                 * function arguments specify which parts
                                 * shall actually be computed. Needs to be
                                 * called before the functions @p get_value(),
                                 * @p get_gradient() or @p get_laplacian
                                 * return useful information.
                                 */
  void evaluate (bool evaluate_val, bool evaluate_grad, 
                 bool evaluate_hess=false);

                                /**
                                 * This function takes the values and/or
                                 * gradients that are stored on quadrature
                                 * points, tests them by all the basis
                                 * functions/gradients on the cell and
                                 * performs the cell integration. The two
                                 * function arguments @p integrate_val and @p
                                 * integrate_grad are used to enable/disable
                                 * some of values or gradients.
                                 */
  void integrate (bool integrate_val, bool integrate_grad);

                                /**
                                 * Returns the q-th quadrature point stored in
                                 * MappingInfo.
                                 */
  Point<dim,vector_t> quadrature_point (const unsigned int q_point) const;

protected:

                                /**
                                 * Internal function that applies the shape
                                 * function data of the tensor product in a
                                 * given coordinate direction (first template
                                 * argument), from polynomials to values on
                                 * quadrature points (second flag set to true)
                                 * or in an integration loop from values on
                                 * quadrature points to values tested by
                                 * different test function (second flag set to
                                 * false), and if the result is to be added to
                                 * some previous results or not.
                                 */
  template <int direction, bool dof_to_quad, bool add>
  void apply_tensor_prod (const vector_t * shape_data,
                          const vector_t in [],
                          vector_t out []);
};



/**
 * The class that provides all functions necessary to evaluate functions at
 * quadrature points and cell integrations. In functionality, this class is
 * similar to FEValues<dim>, however, it includes a lot of specialized
 * functions that make it much faster (between 5 and 500, depending on the
 * polynomial order).
 *
 * This class is a specialization of FEEvaluationGeneral designed for standard
 * FE_Q or FE_DGQ elements and quadrature points symmetric around 0.5 (like
 * Gauss quadrature), and hence the most common situation.
 *
 * This class has five template arguments:
 *
 * @param dim Dimension in which this class is to be used
 *
 * @param n_dofs_1d Number of degrees of freedom of the FE in 1D, usually
 *                   fe_degree+1, for elements based on a tensor product
 *
 * @param n_q_points_1d Number of points in the quadrature formula in 1D,
 *                   usually chosen as fe_degree+1
 *
 * @param n_components Number of vector components when solving a system of
 *                  PDEs. If the same operation is applied to several
 *                  components of a PDE (e.g. a vector Laplace equation), they
 *                  can be applied simultaneously with one call (and often
 *                  more efficiently)
 *
 * @param Number Number format, usually @p double or @p float
 *
 * @author Katharina Kormann and Martin Kronbichler, 2010, 2011
 */
template <int dim, int n_dofs_1d, int n_q_points_1d=n_dofs_1d,
          int n_components=1, typename Number=double >
class FEEvaluation :
  public FEEvaluationGeneral<dim,n_dofs_1d,n_q_points_1d,n_components,Number>
{
 public:
  typedef VectorizedArray<Number> vector_t;
  static const std::size_t  n_vectors = VectorizedArray<Number>::n_array_elements;
  typedef FEEvaluationGeneral<dim,n_dofs_1d,n_q_points_1d,n_components,Number> BaseClass;
  static const unsigned int dofs_per_cell = BaseClass::dofs_per_cell;
  static const unsigned int n_q_points    = BaseClass::n_q_points;

                                /**
                                 * Constructor. Takes all data stored in
                                 * MatrixFree. If applied to problems with
                                 * more than one finite element or more than
                                 * one quadrature formula selected during
                                 * construction of @p matrix_free, @p
                                 * fe_no and @p quad_no allow to select the
                                 * appropriate components.
                                 */
  FEEvaluation (const MatrixFree<dim,Number> &matrix_free,
                const unsigned int                fe_no   = 0,
                const unsigned int                quad_no = 0);

                                /**
                                 * Evaluates the function values, the
                                 * gradients, and the Laplacians of the FE
                                 * function given at the DoF values in the
                                 * input vector at the quadrature points on
                                 * the unit cell.  The function arguments
                                 * specify which parts shall actually be
                                 * computed. Needs to be called before the
                                 * functions @p get_value(), @p get_gradient()
                                 * or @p get_laplacian give useful information
                                 * (unless these values have been set
                                 * manually).
                                 */
  void evaluate (bool evaluate_val, bool evaluate_grad, 
                 bool evaluate_hess=false);

                                /**
                                 * This function takes the values and/or
                                 * gradients that are stored on quadrature
                                 * points, tests them by all the basis
                                 * functions/gradients on the cell and
                                 * performs the cell integration. The two
                                 * function arguments @p integrate_val and @p
                                 * integrate_grad are used to enable/disable
                                 * some of values or gradients.
                                 */
  void integrate (bool integrate_val, bool integrate_grad);

protected:
                                /**
                                 * Internal function that applies the function
                                 * values of the tensor product in a given
                                 * coordinate direction (first template
                                 * argument), from polynomials to values on
                                 * quadrature points (second flag set to true)
                                 * or in an integration loop from values on
                                 * quadrature points to values tested by
                                 * different test function (second flag set to
                                 * false), and if the result is to be added to
                                 * previous content in the data fields or
                                 * not.
                                 */
  template <int direction, bool dof_to_quad, bool add>
  void apply_values (const vector_t in [], vector_t out []);

                                /**
                                 * Internal function that applies the gradient
                                 * operation of the tensor product in a given
                                 * coordinate direction (first template
                                 * argument), from polynomials to values on
                                 * quadrature points (second flag set to true)
                                 * or in an integration loop from values on
                                 * quadrature points to values tested by
                                 * different test function (second flag set to
                                 * false), and if the result is to be added to
                                 * previous content in the data fields or
                                 * not.
                                 */
  template <int direction, bool dof_to_quad, bool add>
  void apply_gradients (const vector_t in [], vector_t out []);

                                /**
                                 * Internal function that applies the second
                                 * derivative operation (Hessian) of the
                                 * tensor product in a given coordinate
                                 * direction (first template argument), from
                                 * polynomials to values on quadrature points
                                 * (second flag set to true) or in an
                                 * integration loop from values on quadrature
                                 * points to values tested by different test
                                 * function (second flag set to false), and if
                                 * the result is to be added to previous
                                 * content in the data fields or not.
                                 */
  template <int direction, bool dof_to_quad, bool add>
  void apply_hessians (const vector_t in [], vector_t out []);
};



/**
 * The class that provides all functions necessary to evaluate functions at
 * quadrature points and cell integrations. In functionality, this class is
 * similar to FEValues<dim>, however, it includes a lot of specialized
 * functions that make it much faster (between 5 and 500, depending on the
 * polynomial order).
 *
 * This class is a specialization of FEEvaluation for elements where
 * quadrature formula and support points are chosen so that a orthonormal
 * relation between the values holds. This is the case for FE_Q based on
 * Gauss-Lobatto support points with Gauss-Lobatto quadrature formula of the
 * same order (QGaussLobatto). In that case, application of values is trivial
 * (as the transformation is the identity matrix), and application of
 * gradients is considerably simpler (since all value applications in
 * directions other than the gradient direction are again identity
 * operations).
 *
 * This class has five template arguments:
 *
 * @param dim Dimension in which this class is to be used
 *
 * @param n_dofs_1d Number of degrees of freedom of the FE in 1D, usually
 *                  fe_degree+1, for elements based on a tensor product
 *
 * @param n_q_points_1d Number of points in the quadrature formula in 1D,
 *                   usually chosen as fe_degree+1
 *
 * @param n_components Number of vector components when solving a system of
 *                  PDEs. If the same operation is applied to several
 *                  components of a PDE (e.g. a vector Laplace equation), they
 *                  can be applied simultaneously with one call (and often
 *                  more efficiently)
 *
 * @param Number Number format, usually @p double or @p float
 *
 * @author Katharina Kormann and Martin Kronbichler, 2010, 2011
 */
template <int dim, int n_points_1d, int n_components=1, typename Number=double >
class FEEvaluationGL :
  public FEEvaluation<dim,n_points_1d,n_points_1d,n_components,Number>
{
 public:
  typedef VectorizedArray<Number> vector_t;
  static const std::size_t  n_vectors = VectorizedArray<Number>::n_array_elements;
  typedef FEEvaluation<dim,n_points_1d,n_points_1d,n_components,Number> BaseClass;
  static const unsigned int dofs_per_cell = BaseClass::dofs_per_cell;
  static const unsigned int n_q_points    = BaseClass::n_q_points;

                                /**
                                 * Constructor. Takes all data stored in
                                 * MatrixFree. If applied to problems with
                                 * more than one finite element or more than
                                 * one quadrature formula selected during
                                 * construction of @p matrix_free, @p
                                 * fe_no and @p quad_no allow to select the
                                 * appropriate components.
                                 */
  FEEvaluationGL (const MatrixFree<dim,Number> &matrix_free,
                    const unsigned int                fe_no   = 0,
                    const unsigned int                quad_no = 0);

                                /**
                                 * Evaluates the function values, the
                                 * gradients, and the Hessians of the FE
                                 * function given at the DoF values in the
                                 * input vector at the quadrature points of
                                 * the unit cell. The function arguments
                                 * specify which parts shall actually be
                                 * computed. Needs to be called before the
                                 * functions @p get_value(), @p get_gradient()
                                 * or @p get_laplacian give useful information
                                 * (unless these values have been set
                                 * manually).
                                 */
  void evaluate (bool evaluate_val, bool evaluate_grad, 
                 bool evaluate_lapl=false);

                                /**
                                 * This function takes the values and/or
                                 * gradients that are stored on quadrature
                                 * points, tests them by all the basis
                                 * functions/gradients on the cell and
                                 * performs the cell integration. The two
                                 * function arguments @p integrate_val and @p
                                 * integrate_grad are used to enable/disable
                                 * some of values or gradients.
                                 */
  void integrate (bool integrate_val, bool integrate_grad);

protected:
                                /**
                                 * Internal function that applies the gradient
                                 * operation of the tensor product in a given
                                 * coordinate direction (first template
                                 * argument), from polynomials to values on
                                 * quadrature points (second flag set to true)
                                 * or in an integration loop from values on
                                 * quadrature points to values tested by
                                 * different test function (second flag set to
                                 * false), and if the result is to be added to
                                 * some previous results or not.
                                 */
  template <int direction, bool dof_to_quad, bool add>
  void apply_gradients (const vector_t in [], vector_t out []);
};




/*----------------------- Inline functions ----------------------------------*/

#ifndef DOXYGEN


/*----------------------- FEEvaluationBase -------------------------------*/

template <int dim, int dofs_per_cell_, int n_q_points_,
          int n_components, typename Number>
inline
FEEvaluationBase<dim,dofs_per_cell_,n_q_points_,n_components,Number>::
FEEvaluationBase (const MatrixFree<dim,Number> &data_in,
                  const unsigned int fe_no_in,
                  const unsigned int quad_no_in)
  :
  matrix_info        (data_in),
  dof_info           (data_in.get_dof_info(fe_no_in)),
  constraint_pool    (data_in.get_constraint_pool()),
  mapping_info       (data_in.get_mapping_info()),
  active_fe_index    (dof_info.fe_index_from_dofs_per_cell
                      (dofs_per_cell_ * dof_info.n_components)),
  active_quad_index  (mapping_info.
                      mapping_data_gen[quad_no_in].
                      quad_index_from_n_q_points(n_q_points_)),
  data               (data_in.get_fe_evaluation
                      (fe_no_in, quad_no_in, active_fe_index,
                       active_quad_index)),
  quad_no            (quad_no_in),
  n_fe_components    (dof_info.n_components),
  cell               (numbers::invalid_unsigned_int),
  cell_type          (numbers::invalid_unsigned_int),
  cartesian          (0),
  jacobian           (0),
  J_value            (0),
  quadrature_weights (mapping_info.mapping_data_gen[quad_no].
                      quadrature_weights[active_quad_index].begin()),
  quadrature_points  (0),
  jacobian_grad      (0),
  jacobian_grad_upper(0)
{
  Assert (matrix_info.indices_initialized() == true,
          ExcNotInitialized());
  Assert (matrix_info.mapping_initialized() == true,
          ExcNotInitialized());
  AssertDimension (matrix_info.get_size_info().n_vectors, n_vectors);
  Assert (n_fe_components == 1 ||
          n_components == n_fe_components,
          ExcMessage ("The underlying FE is vector-valued. In this case, the "
                      "template argument n_components must be a the same "
                      "as the number of underlying vector components."));


                                // do not check for correct dimensions of data
                                // fields here, should be done in derived
                                // classes
}



template <int dim, int dofs_per_cell_, int n_q_points_,
          int n_components, typename Number>
inline
void
FEEvaluationBase<dim,dofs_per_cell_,n_q_points_,n_components,Number>::
reinit (const unsigned int cell_in)
{
  AssertIndexRange (cell_in, dof_info.row_starts.size()-1);
  AssertDimension (((dof_info.cell_active_fe_index.size() > 0) ?
                    dof_info.cell_active_fe_index[cell_in] : 0),
                   active_fe_index);
  cell = cell_in;
  cell_type = mapping_info.get_cell_type(cell);
  cell_data_number = mapping_info.get_cell_data_index(cell);

  if (mapping_info.quadrature_points_initialized == true)
    {
      AssertIndexRange (cell_data_number, mapping_info.
                        mapping_data_gen[quad_no].rowstart_q_points.size());
      const unsigned int index = mapping_info.mapping_data_gen[quad_no].
        rowstart_q_points[cell];
      AssertIndexRange (index, mapping_info.mapping_data_gen[quad_no].
                        quadrature_points.size());
      quadrature_points =
        &mapping_info.mapping_data_gen[quad_no].quadrature_points[index];
    }

  if (cell_type == 0)
    {
      cartesian = &mapping_info.cartesian[cell_data_number].first;
      J_value   = &mapping_info.cartesian[cell_data_number].second;
    }
  else if (cell_type == 1)
    {
      jacobian  = &mapping_info.linear[cell_data_number].first;
      J_value   = &mapping_info.linear[cell_data_number].second;
    }
  else
    {
      const unsigned int rowstart = mapping_info.
        mapping_data_gen[quad_no].rowstart_jacobians[cell_data_number];
      AssertIndexRange (rowstart, mapping_info.
                        mapping_data_gen[quad_no].jacobians.size());
      jacobian =
        &mapping_info.mapping_data_gen[quad_no].jacobians[rowstart];
      if (mapping_info.JxW_values_initialized == true)
        {
          AssertIndexRange (rowstart, mapping_info.
                            mapping_data_gen[quad_no].JxW_values.size());
          J_value = &(mapping_info.mapping_data_gen[quad_no].
                      JxW_values[rowstart]);
        }
      if (mapping_info.second_derivatives_initialized == true)
        {
          AssertIndexRange(rowstart, mapping_info.
                           mapping_data_gen[quad_no].jacobians_grad_diag.size());
          jacobian_grad = &mapping_info.mapping_data_gen[quad_no].
            jacobians_grad_diag[rowstart];
          AssertIndexRange(rowstart, mapping_info.
                           mapping_data_gen[quad_no].jacobians_grad_upper.size());
          jacobian_grad_upper = &mapping_info.mapping_data_gen[quad_no].
            jacobians_grad_upper[rowstart];
        }
    }

  n_irreg_components_filled =
    std_cxx1x::get<2>(dof_info.row_starts[cell_in]);
  at_irregular_cell = n_irreg_components_filled > 0;
#ifdef DEBUG
  dof_values_initialized      = false;
  values_quad_initialized     = false;
  gradients_quad_initialized  = false;
  hessians_quad_initialized   = false;
#endif
}



template <int dim, int dofs_per_cell_, int n_q_points_,
          int n_components, typename Number>
inline
unsigned int
FEEvaluationBase<dim,dofs_per_cell_,n_q_points_,n_components,Number>::
get_cell_data_number () const
{
  Assert (cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  return cell_data_number;
}



template <int dim, int dofs_per_cell_, int n_q_points_,
          int n_components, typename Number>
inline
unsigned int
FEEvaluationBase<dim,dofs_per_cell_,n_q_points_,n_components,Number>::
get_cell_type () const
{
  Assert (cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  return cell_type;
}



template <int dim, int dofs_per_cell_, int n_q_points_,
          int n_components, typename Number>
inline
const VectorizedArray<Number> *
FEEvaluationBase<dim,dofs_per_cell_,n_q_points_,n_components,Number>::
begin_values () const
{
  Assert (values_quad_initialized || values_quad_submitted,
          ExcNotInitialized());
  return &values_quad[0][0];
}



template <int dim, int dofs_per_cell_, int n_q_points_,
          int n_components, typename Number>
inline
VectorizedArray<Number> *
FEEvaluationBase<dim,dofs_per_cell_,n_q_points_,n_components,Number>::
begin_values ()
{
#ifdef DEBUG
  values_quad_submitted = true;
#endif
  return &values_quad[0][0];
}



template <int dim, int dofs_per_cell_, int n_q_points_,
          int n_components, typename Number>
inline
const VectorizedArray<Number> *
FEEvaluationBase<dim,dofs_per_cell_,n_q_points_,n_components,Number>::
begin_gradients () const
{
  Assert (gradients_quad_initialized || gradients_quad_submitted,
          ExcNotInitialized());
  return &gradients_quad[0][0][0];
}



template <int dim, int dofs_per_cell_, int n_q_points_,
          int n_components, typename Number>
inline
VectorizedArray<Number> *
FEEvaluationBase<dim,dofs_per_cell_,n_q_points_,n_components,Number>::
begin_gradients ()
{
#ifdef DEBUG
  gradients_quad_submitted = true;
#endif
  return &gradients_quad[0][0][0];
}



template <int dim, int dofs_per_cell_, int n_q_points_,
          int n_components, typename Number>
inline
const VectorizedArray<Number> *
FEEvaluationBase<dim,dofs_per_cell_,n_q_points_,n_components,Number>::
begin_hessians () const
{
  Assert (hessians_quad_initialized, ExcNotInitialized());
  return &hessians_quad[0][0][0];
}



template <int dim, int dofs_per_cell_, int n_q_points_,
          int n_components, typename Number>
inline
VectorizedArray<Number> *
FEEvaluationBase<dim,dofs_per_cell_,n_q_points_,n_components,Number>::
begin_hessians ()
{
  return &hessians_quad[0][0][0];
}



namespace internal
{
                                // write access to generic vectors that have
                                // operator ().
  template <typename VectorType>
  inline
  typename VectorType::value_type &
  vector_access (VectorType         &vec,
                 const unsigned int  entry)
  {
    return vec(entry);
  }



                                // read access to generic vectors that have
                                // operator ().
  template <typename VectorType>
  inline
  typename VectorType::value_type
  vector_access (const VectorType   &vec,
                 const unsigned int  entry)
  {
    return vec(entry);
  }



                                // write access to distributed MPI vectors
                                // that have operator [] to access data in
                                // local index space, which is what we use in
                                // DoFInfo and hence in read_dof_values etc.
  template <typename Number>
  inline
  Number &
  vector_access (parallel::distributed::Vector<Number> &vec,
                 const unsigned int                     entry)
  {
    return vec.local_element(entry);
  }



                                // read access to distributed MPI vectors that
                                // have operator [] to access data in local
                                // index space, which is what we use in
                                // DoFInfo and hence in read_dof_values etc.
  template <typename Number>
  inline
  Number
  vector_access (const parallel::distributed::Vector<Number> &vec,
                 const unsigned int                           entry)
  {
    return vec.local_element(entry);
  }



                                // this is to make sure that the parallel
                                // partitioning in the
                                // parallel::distributed::Vector is really the
                                // same as stored in MatrixFree
  template <typename VectorType>
  inline
  void check_vector_compatibility (const VectorType &vec,
                                   const internal::MatrixFreeFunctions::DoFInfo      &dof_info)
  {
    AssertDimension (vec.size(),
                     dof_info.vector_partitioner->size());
  }

  template <typename Number>
  inline
  void check_vector_compatibility (const parallel::distributed::Vector<Number> &vec,
                                   const internal::MatrixFreeFunctions::DoFInfo      &dof_info)
  {
    Assert (vec.partitioners_are_compatible(*dof_info.vector_partitioner),
            ExcMessage("The parallel layout of the given vector is not "
                       "compatible with the parallel partitioning in MatrixFree. "
                       "Use MatrixFree::initialize_dof_vector to get a "
                       "compatible vector."));
  }
}



template <int dim, int dofs_per_cell_, int n_q_points_,
          int n_components, typename Number>
template<typename VectorType>
inline
void
FEEvaluationBase<dim,dofs_per_cell_,n_q_points_,n_components,Number>::
read_dof_values (const VectorType &src)
{
  AssertDimension (n_components, n_fe_components);
                                // only need one component, but to avoid
                                // compiler warnings, use n_components copies
                                // here (but these will not be used)
  const VectorType * src_data[n_components];
  for (unsigned int d=0; d<n_components; ++d)
    src_data[d] = &src;
  read_dof_values (src_data);
}



template <int dim, int dofs_per_cell_, int n_q_points_,
          int n_components, typename Number>
template<typename VectorType>
inline
void
FEEvaluationBase<dim,dofs_per_cell_,n_q_points_,n_components,Number>::
read_dof_values (const std::vector<VectorType> &src,
                const unsigned int             first_index)
{
  AssertIndexRange (first_index, src.size());
  Assert (n_fe_components == 1, ExcNotImplemented());
  Assert ((n_fe_components == 1 ?
           (first_index+n_components <= src.size()) : true),
          ExcIndexRange (first_index + n_components, 0, src.size()));
  const VectorType * src_data [n_components];
  for (unsigned int comp=0; comp<n_components; ++comp)
    src_data[comp] = &src[comp+first_index];
  read_dof_values (src_data);
}



template <int dim, int dofs_per_cell_, int n_q_points_,
          int n_components, typename Number>
template<typename VectorType>
inline
void
FEEvaluationBase<dim,dofs_per_cell_,n_q_points_,n_components,Number>::
read_dof_values (const std::vector<VectorType*> &src,
                const unsigned int              first_index)
{
  AssertIndexRange (first_index, src.size());
  Assert (n_fe_components == 1, ExcNotImplemented());
  Assert ((n_fe_components == 1 ?
           (first_index+n_components <= src.size()) : true),
          ExcIndexRange (first_index + n_components, 0, src.size()));
  const VectorType * src_data [n_components];
  for (unsigned int comp=0; comp<n_components; ++comp)
    src_data[comp] = src[comp+first_index];
  read_dof_values (src_data);
}



template <int dim, int dofs_per_cell_, int n_q_points_,
          int n_components, typename Number>
template<typename VectorType>
inline
void
FEEvaluationBase<dim,dofs_per_cell_,n_q_points_,n_components,Number>::
read_dof_values (const VectorType * src[])
{
  Assert (cell != numbers::invalid_unsigned_int, ExcNotInitialized());

                                // loop over all local dofs. ind_local holds
                                // local number on cell, index iterates over
                                // the elements of index_local_to_global and
                                // dof_indices points to the global indices
                                // stored in index_local_to_global
  const unsigned int * dof_indices = dof_info.begin_indices(cell);
  const std::pair<unsigned short,unsigned short> * indicators =
    dof_info.begin_indicators(cell);
  const std::pair<unsigned short,unsigned short> * indicators_end =
    dof_info.end_indicators(cell);
  unsigned int ind_local = 0;

                                // scalar case (or case when all components
                                // have the same degrees of freedom and sit on
                                // a different vector each)
  if (n_fe_components == 1)
    {
      for (unsigned int comp=0; comp<n_components; ++comp)
        internal::check_vector_compatibility (*src[comp], dof_info);
      Number * local_src_number [n_components];
      for (unsigned int comp=0; comp<n_components; ++comp)
        local_src_number[comp] = reinterpret_cast<Number*>(values_dofs[comp]);

                                // standard case where there are sufficiently
                                // many cells to fill all vectors
      if (at_irregular_cell == false)
        {
                                // check whether there is any constraint on
                                // the current cell
          if (indicators != indicators_end)
            {
              for ( ; indicators != indicators_end; ++indicators)
                {
                                // run through values up to next constraint
                  for (unsigned int j=0; j<indicators->first; ++j)
                    for (unsigned int comp=0; comp<n_components; ++comp)
                      local_src_number[comp][ind_local+j] =
                        internal::vector_access (*src[comp], dof_indices[j]);
                  ind_local += indicators->first;
                  dof_indices   += indicators->first;

                                // constrained case: build the local value as
                                // a linear combination of the global value
                                // according to constraints
                  Number value [n_components];
                  for (unsigned int comp=0; comp<n_components; ++comp)
                    value[comp] = 0;
                  const Number * data_val =
                    constraint_pool.begin(indicators->second);
                  const unsigned int row_length =
                    constraint_pool.row_length(indicators->second);
                  for (unsigned int k=0; k<row_length; ++k)
                    for (unsigned int comp=0; comp<n_components; ++comp)
                      value[comp] +=
                        (internal::vector_access (*src[comp], dof_indices[k]) *
                         data_val[k]);
                  dof_indices += row_length;

                  for (unsigned int comp=0; comp<n_components; ++comp)
                    local_src_number[comp][ind_local] = value[comp];
                  ind_local++;
                }

                                // get the dof values past the last
                                // constraint
              for(; ind_local<n_vectors*dofs_per_cell; ++dof_indices, ++ind_local)
                {
                  for (unsigned int comp=0; comp<n_components; ++comp)
                    local_src_number[comp][ind_local] =
                      internal::vector_access (*src[comp], *dof_indices);
                }
            }
          else
            {
                                // no constraint at all: loop bounds are
                                // known, compiler can unroll without checks
              AssertDimension (dof_info.end_indices(cell)-dof_indices,
                               n_vectors*dofs_per_cell);
              for (unsigned int j=0; j<dofs_per_cell*n_vectors; ++j)
                for (unsigned int comp=0; comp<n_components; ++comp)
                  local_src_number[comp][j] =
                    internal::vector_access (*src[comp], dof_indices[j]);
            }
        }

                                // non-standard case: need to fill in zeros
                                // for those components that are not present
                                // (a bit more expensive), but there is not
                                // more than one such cell
      else
        {
          Assert (n_irreg_components_filled > 0, ExcInternalError());
          for ( ; indicators != indicators_end; ++indicators)
            {
              for(unsigned int j=0; j<indicators->first; ++j)
                {
                                // non-constrained case: copy the data from
                                // the global vector, src, to the local one,
                                // local_src.
                  for (unsigned int comp=0; comp<n_components; ++comp)
                    local_src_number[comp][ind_local] =
                      internal::vector_access (*src[comp], dof_indices[j]);

                                // here we jump over all the components that
                                // are artificial
                  ++ind_local;
                  while (ind_local % n_vectors >= n_irreg_components_filled)
                    {
                      for (unsigned int comp=0; comp<n_components; ++comp)
                        local_src_number[comp][ind_local] = 0.;
                      ++ind_local;
                    }
                }
              dof_indices += indicators->first;

                                // constrained case: build the local value as
                                // a linear combination of the global value
                                // according to constraint
              Number value [n_components];
              for (unsigned int comp=0; comp<n_components; ++comp)
                value[comp] = 0.;
              const Number * data_val =
                constraint_pool.begin(indicators->second);
              const unsigned int row_length =
                constraint_pool.row_length(indicators->second);
              for (unsigned int k=0; k<row_length; ++k)
                for (unsigned int comp=0; comp<n_components; ++comp)
                  value[comp] +=
                    internal::vector_access (*src[comp], dof_indices[k]) * data_val[k];
              dof_indices += row_length;
              for (unsigned int comp=0; comp<n_components; ++comp)
                local_src_number[comp][ind_local] = value[comp];
              ind_local++;
              while (ind_local % n_vectors >= n_irreg_components_filled)
                {
                  for (unsigned int comp=0; comp<n_components; ++comp)
                    local_src_number[comp][ind_local] = 0.;
                  ++ind_local;
                }
            }
          for(; ind_local<n_vectors*dofs_per_cell; ++dof_indices)
            {
              Assert (dof_indices != dof_info.end_indices(cell),
                      ExcInternalError());

                                // non-constrained case: copy the data from
                                // the global vector, src, to the local one,
                                // local_dst.
              for (unsigned int comp=0; comp<n_components; ++comp)
                local_src_number[comp][ind_local] =
                  internal::vector_access (*src[comp], *dof_indices);
              ++ind_local;
              while (ind_local % n_vectors >= n_irreg_components_filled)
                {
                  for (unsigned int comp=0; comp<n_components; ++comp)
                    local_src_number[comp][ind_local] = 0.;
                  ++ind_local;
                }
            }
        }
    }
  else
                                // case with vector-valued finite elements
                                // where all components are included in one
                                // single vector. Assumption: first come all
                                // entries to the first component, then all
                                // entries to the second one, and so on. This
                                // is ensured by the way MatrixFree reads
                                // out the indices.
    {
      internal::check_vector_compatibility (*src[0], dof_info);
      Assert (n_fe_components == n_components, ExcNotImplemented());
      const unsigned int total_dofs_per_cell =
        dofs_per_cell * n_vectors * n_components;
      Number * local_src_number = reinterpret_cast<Number*>(values_dofs[0]);
      if (at_irregular_cell == false)
        {
                                // check whether there is any constraint on
                                // the current cell
          if (indicators != indicators_end)
            {
              for ( ; indicators != indicators_end; ++indicators)
                {
                                // run through values up to next constraint
                  for (unsigned int j=0; j<indicators->first; ++j)
                      local_src_number[ind_local+j] =
                        internal::vector_access (*src[0], dof_indices[j]);
                  ind_local += indicators->first;
                  dof_indices   += indicators->first;

                                // constrained case: build the local value as
                                // a linear combination of the global value
                                // according to constraints
                  Number value = 0;
                  const Number * data_val =
                    constraint_pool.begin(indicators->second);
                  const unsigned int row_length =
                    constraint_pool.row_length(indicators->second);
                  for (unsigned int k=0; k<row_length; ++k)
                    value +=
                      (internal::vector_access (*src[0], dof_indices[k]) *
                       data_val[k]);
                  dof_indices += row_length;

                  local_src_number[ind_local] = value;
                  ind_local++;
                }

                                // get the dof values past the last
                                // constraint
              for(; ind_local<total_dofs_per_cell; ++dof_indices, ++ind_local)
                local_src_number[ind_local] =
                  internal::vector_access (*src[0], *dof_indices);
              Assert (dof_indices == dof_info.end_indices(cell),
                      ExcInternalError());
            }
          else
            {
                                // no constraint at all: loop bounds are
                                // known, compiler can unroll without checks
              AssertDimension (dof_info.end_indices(cell)-dof_indices,
                               static_cast<int>(total_dofs_per_cell));
              for (unsigned int j=0; j<total_dofs_per_cell; ++j)
                local_src_number[j] =
                  internal::vector_access (*src[0], dof_indices[j]);
            }
        }

                                // non-standard case: need to fill in zeros
                                // for those components that are not present
                                // (a bit more expensive), but there is not
                                // more than one such cell
      else
        {
          Assert (n_irreg_components_filled > 0, ExcInternalError());
          for ( ; indicators != indicators_end; ++indicators)
            {
              for(unsigned int j=0; j<indicators->first; ++j)
                {
                                // non-constrained case: copy the data from
                                // the global vector, src, to the local one,
                                // local_src.
                  local_src_number[ind_local] =
                    internal::vector_access (*src[0], dof_indices[j]);

                                // here we jump over all the components that
                                // are artificial
                  ++ind_local;
                  while (ind_local % n_vectors >= n_irreg_components_filled)
                    {
                      local_src_number[ind_local] = 0.;
                      ++ind_local;
                    }
                }
              dof_indices += indicators->first;

                                // constrained case: build the local value as
                                // a linear combination of the global value
                                // according to constraint
              Number value = 0;
              const Number * data_val =
                constraint_pool.begin(indicators->second);
              const unsigned int row_length =
                constraint_pool.row_length(indicators->second);
              for (unsigned int k=0; k<row_length; ++k)
                value +=
                  internal::vector_access (*src[0], dof_indices[k]) * data_val[k];
              dof_indices += row_length;
              local_src_number[ind_local] = value;
              ind_local++;
              while (ind_local % n_vectors >= n_irreg_components_filled)
                {
                  local_src_number[ind_local] = 0.;
                  ++ind_local;
                }
            }
          for(; ind_local<total_dofs_per_cell; ++dof_indices)
            {
              Assert (dof_indices != dof_info.end_indices(cell),
                      ExcInternalError());

                                // non-constrained case: copy the data from
                                // the global vector, src, to the local one,
                                // local_dst.
              local_src_number[ind_local] =
                internal::vector_access (*src[0], *dof_indices);
              ++ind_local;
              while (ind_local % n_vectors >= n_irreg_components_filled)
                {
                  local_src_number[ind_local] = 0.;
                  ++ind_local;
                }
            }
        }
    }

#ifdef DEBUG
  dof_values_initialized = true;
#endif
}






template <int dim, int dofs_per_cell_, int n_q_points_,
          int n_components, typename Number>
template<typename VectorType>
inline
void
FEEvaluationBase<dim,dofs_per_cell_,n_q_points_,n_components,Number>::
read_dof_values_plain (const VectorType &src)
{
  AssertDimension (n_components, n_fe_components);
                                // only need one component, but to avoid
                                // compiler warnings, use n_components copies
                                // here (but these will not be used)
  const VectorType * src_data[n_components];
  for (unsigned int d=0; d<n_components; ++d)
    src_data[d] = &src;
  read_dof_values_plain (src_data);
}



template <int dim, int dofs_per_cell_, int n_q_points_,
          int n_components, typename Number>
template<typename VectorType>
inline
void
FEEvaluationBase<dim,dofs_per_cell_,n_q_points_,n_components,Number>::
read_dof_values_plain (const std::vector<VectorType> &src,
                      const unsigned int             first_index)
{
  AssertIndexRange (first_index, src.size());
  Assert (n_fe_components == 1, ExcNotImplemented());
  Assert ((n_fe_components == 1 ?
           (first_index+n_components <= src.size()) : true),
          ExcIndexRange (first_index + n_components, 0, src.size()));
  const VectorType * src_data [n_components];
  for (unsigned int comp=0; comp<n_components; ++comp)
    src_data[comp] = &src[comp+first_index];
  read_dof_values_plain (src_data);
}



template <int dim, int dofs_per_cell_, int n_q_points_,
          int n_components, typename Number>
template<typename VectorType>
inline
void
FEEvaluationBase<dim,dofs_per_cell_,n_q_points_,n_components,Number>::
read_dof_values_plain (const std::vector<VectorType*> &src,
                      const unsigned int              first_index)
{
  AssertIndexRange (first_index, src.size());
  Assert (n_fe_components == 1, ExcNotImplemented());
  Assert ((n_fe_components == 1 ?
           (first_index+n_components <= src.size()) : true),
          ExcIndexRange (first_index + n_components, 0, src.size()));
  const VectorType * src_data [n_components];
  for (unsigned int comp=0; comp<n_components; ++comp)
    src_data[comp] = src[comp+first_index];
  read_dof_values_plain (src_data);
}



template <int dim, int dofs_per_cell_, int n_q_points_,
          int n_components, typename Number>
template<typename VectorType>
inline
void
FEEvaluationBase<dim,dofs_per_cell_,n_q_points_,n_components,Number>::
read_dof_values_plain (const VectorType * src[])
{
  Assert (cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  Assert (dof_info.store_plain_indices == true, ExcNotInitialized());

                                // loop over all local dofs. ind_local holds
                                // local number on cell, index iterates over
                                // the elements of index_local_to_global and
                                // dof_indices points to the global indices
                                // stored in index_local_to_global
  const unsigned int * dof_indices = dof_info.begin_indices_plain(cell);

                                // scalar case (or case when all components
                                // have the same degrees of freedom and sit on
                                // a different vector each)
  if (n_fe_components == 1)
    {
      for (unsigned int comp=0; comp<n_components; ++comp)
        internal::check_vector_compatibility (*src[comp], dof_info);
      Number * local_src_number [n_components];
      for (unsigned int comp=0; comp<n_components; ++comp)
        local_src_number[comp] = reinterpret_cast<Number*>(values_dofs[comp]);

                                // standard case where there are sufficiently
                                // many cells to fill all vectors
      if (at_irregular_cell == false)
        {
          for (unsigned int j=0; j<dofs_per_cell*n_vectors; ++j)
            for (unsigned int comp=0; comp<n_components; ++comp)
              local_src_number[comp][j] =
                internal::vector_access (*src[comp], dof_indices[j]);
        }

                                // non-standard case: need to fill in zeros
                                // for those components that are not present
                                // (a bit more expensive), but there is not
                                // more than one such cell
      else
        {
          Assert (n_irreg_components_filled > 0, ExcInternalError());
          for(unsigned int ind_local=0; ind_local<n_vectors*dofs_per_cell;
              ++dof_indices)
            {
                                // non-constrained case: copy the data from
                                // the global vector, src, to the local one,
                                // local_dst.
              for (unsigned int comp=0; comp<n_components; ++comp)
                local_src_number[comp][ind_local] =
                  internal::vector_access (*src[comp], *dof_indices);
              ++ind_local;
              while (ind_local % n_vectors >= n_irreg_components_filled)
                {
                  for (unsigned int comp=0; comp<n_components; ++comp)
                    local_src_number[comp][ind_local] = 0.;
                  ++ind_local;
                }
            }
        }
    }
  else
                                // case with vector-valued finite elements
                                // where all components are included in one
                                // single vector. Assumption: first come all
                                // entries to the first component, then all
                                // entries to the second one, and so on. This
                                // is ensured by the way MatrixFree reads
                                // out the indices.
    {
      internal::check_vector_compatibility (*src[0], dof_info);
      Assert (n_fe_components == n_components, ExcNotImplemented());
      const unsigned int total_dofs_per_cell =
        dofs_per_cell * n_vectors * n_components;
      Number * local_src_number = reinterpret_cast<Number*>(values_dofs[0]);
      if (at_irregular_cell == false)
        {
          for (unsigned int j=0; j<total_dofs_per_cell; ++j)
            local_src_number[j] =
              internal::vector_access (*src[0], dof_indices[j]);
        }

                                // non-standard case: need to fill in zeros
                                // for those components that are not present
                                // (a bit more expensive), but there is not
                                // more than one such cell
      else
        {
          Assert (n_irreg_components_filled > 0, ExcInternalError());
          for(unsigned int ind_local=0; ind_local<total_dofs_per_cell; ++dof_indices)
            {
                                // non-constrained case: copy the data from
                                // the global vector, src, to the local one,
                                // local_dst.
              local_src_number[ind_local] =
                internal::vector_access (*src[0], *dof_indices);
              ++ind_local;
              while (ind_local % n_vectors >= n_irreg_components_filled)
                {
                  local_src_number[ind_local] = 0.;
                  ++ind_local;
                }
            }
        }
    }

#ifdef DEBUG
  dof_values_initialized = true;
#endif
}



template <int dim, int dofs_per_cell_, int n_q_points_,
          int n_components, typename Number>
template<typename VectorType>
inline
void
FEEvaluationBase<dim,dofs_per_cell_,n_q_points_,n_components,Number>::
distribute_local_to_global (VectorType &dst) const
{
  AssertDimension (n_components, n_fe_components);
                                // only need one component, but to avoid
                                // compiler warnings, use n_components copies
                                // here (but these will not be used)
  VectorType * dst_data [n_components];
  for (unsigned int d=0; d<n_components; ++d)
    dst_data[d] = &dst;
  distribute_local_to_global (dst_data);
}



template <int dim, int dofs_per_cell_, int n_q_points_,
          int n_components, typename Number>
template<typename VectorType>
inline
void
FEEvaluationBase<dim,dofs_per_cell_,n_q_points_,n_components,Number>::
distribute_local_to_global (std::vector<VectorType>  &dst,
                            const unsigned int        first_index) const
{
  AssertIndexRange (first_index, dst.size());
  Assert (n_fe_components == 1, ExcNotImplemented());
  Assert ((n_fe_components == 1 ?
           (first_index+n_components <= dst.size()) : true),
          ExcIndexRange (first_index + n_components, 0, dst.size()));

  VectorType * dst_data [n_components];
  for (unsigned int comp=0; comp<n_components; ++comp)
    dst_data[comp] = &dst[comp+first_index];
  distribute_local_to_global (dst_data);
}



template <int dim, int dofs_per_cell_, int n_q_points_,
          int n_components, typename Number>
template<typename VectorType>
inline
void
FEEvaluationBase<dim,dofs_per_cell_,n_q_points_,n_components,Number>::
distribute_local_to_global (std::vector<VectorType*>  &dst,
                            const unsigned int         first_index) const
{
  AssertIndexRange (first_index, dst.size());
  Assert (n_fe_components == 1, ExcNotImplemented());
  Assert ((n_fe_components == 1 ?
           (first_index+n_components <= dst.size()) : true),
          ExcIndexRange (first_index + n_components, 0, dst.size()));

  VectorType * dst_data [n_components];
  for (unsigned int comp=0; comp<n_components; ++comp)
    dst_data[comp] = dst[comp+first_index];
  distribute_local_to_global (dst_data);
}



template <int dim, int dofs_per_cell_, int n_q_points_,
          int n_components, typename Number>
template<typename VectorType>
inline
void
FEEvaluationBase<dim,dofs_per_cell_,n_q_points_,n_components,Number>::
distribute_local_to_global (VectorType * dst[]) const
{
  Assert (cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  Assert (dof_values_initialized==true,
          internal::ExcAccessToUninitializedField());

                                // loop over all local dofs. ind_local holds
                                // local number on cell, index iterates over
                                // the elements of index_local_to_global and
                                // dof_indices points to the global indices
                                // stored in index_local_to_global
  const unsigned int * dof_indices = dof_info.begin_indices(cell);
  const std::pair<unsigned short,unsigned short> * indicators =
    dof_info.begin_indicators(cell);
  const std::pair<unsigned short,unsigned short> * indicators_end =
    dof_info.end_indicators(cell);
  unsigned int ind_local = 0;

                                // scalar case (or case when all components
                                // have the same degrees of freedom and sit on
                                // a different vector each)
  if (n_fe_components == 1)
    {
      for (unsigned int comp=0; comp<n_components; ++comp)
        internal::check_vector_compatibility (*dst[comp], dof_info);

      const Number * local_dst_number [n_components];
      for (unsigned int comp=0; comp<n_components; ++comp)
        local_dst_number[comp] =
          reinterpret_cast<const Number*>(values_dofs[comp]);
      if (at_irregular_cell == false)
        {
                                // check whether there is no constraint at all
          if (indicators != indicators_end)
            {
                                // run from one constraint to the next
              for ( ; indicators != indicators_end; ++indicators)
                {
                                // distribute values up to the constraint
                                // (values not constrained)
                  for (unsigned int j=0; j<indicators->first; ++j)
                    for (unsigned int comp=0; comp<n_components; ++comp)
                      internal::vector_access (*dst[comp], dof_indices[j])
                        += local_dst_number[comp][ind_local+j];
                  dof_indices += indicators->first;
                  ind_local   += indicators->first;

                                // constrained case: build the local value as
                                // a linear combination of the global value
                                // according to constraint
                  const Number * data_val =
                    constraint_pool.begin(indicators->second);
                  const unsigned int row_length =
                    constraint_pool.row_length(indicators->second);
                  for (unsigned int k=0; k<row_length; ++k)
                    for (unsigned int comp=0; comp<n_components; ++comp)
                      internal::vector_access (*dst[comp], dof_indices[k])
                        += local_dst_number[comp][ind_local] * data_val[k];
                  dof_indices += row_length;
                  ++ind_local;
                }
                                // distribute values after the last constraint
                                // (values not constrained)
              for(; ind_local<dofs_per_cell*n_vectors; ++dof_indices, ++ind_local)
                for (unsigned int comp=0; comp<n_components; ++comp)
                  internal::vector_access (*dst[comp], *dof_indices)
                    += local_dst_number[comp][ind_local];
            }
                                // no constraint at all: loop bounds are
                                // known, compiler can unroll without checks
          else
            {
              AssertDimension (dof_info.end_indices(cell)-dof_indices,
                               n_vectors * dofs_per_cell);
              for (unsigned int j=0; j<dofs_per_cell*n_vectors; ++j)
                for (unsigned int comp=0; comp<n_components; ++comp)
                  internal::vector_access (*dst[comp], dof_indices[j])
                    += local_dst_number[comp][j];
            }
          return;
        }

                                // irregular case
      Assert (n_irreg_components_filled > 0, ExcInternalError());
      for ( ; indicators != indicators_end; ++indicators)
        {
          for(unsigned int j=0; j<indicators->first; ++j)
            {
                                // non-constrained case: copy the data from
                                // the local vector to the global one.
              for (unsigned int comp=0; comp<n_components; ++comp)
                internal::vector_access (*dst[comp], dof_indices[j])
                  += local_dst_number[comp][ind_local];
              ++ind_local;
              if (ind_local % n_vectors == n_irreg_components_filled)
                ind_local += n_vectors-n_irreg_components_filled;
            }
          dof_indices += indicators->first;

                                // constrained case: distribute according to
                                // the constraint
          const Number * data_val =
            constraint_pool.begin(indicators->second);
          const unsigned int row_length =
            constraint_pool.row_length(indicators->second);
          for (unsigned int k=0; k<row_length; ++k)
            {
              for (unsigned int comp=0; comp<n_components; ++comp)
                internal::vector_access (*dst[comp], dof_indices[k])
                  += local_dst_number[comp][ind_local] * data_val[k];
            }
          dof_indices += row_length;
          ++ind_local;
          if (ind_local % n_vectors == n_irreg_components_filled)
            ind_local += n_vectors-n_irreg_components_filled;
        }
      for(; ind_local<dofs_per_cell*n_vectors; ++dof_indices)
        {
          Assert (dof_indices != dof_info.end_indices(cell),
                  ExcInternalError());

                                // non-constrained case
          for (unsigned int comp=0; comp<n_components; ++comp)
            internal::vector_access (*dst[comp], *dof_indices)
              += local_dst_number[comp][ind_local];
          ++ind_local;
          if (ind_local % n_vectors == n_irreg_components_filled)
            ind_local += n_vectors-n_irreg_components_filled;
        }
    }
  else
                                // case with vector-valued finite elements
                                // where all components are included in one
                                // single vector. Assumption: first come all
                                // entries to the first component, then all
                                // entries to the second one, and so on. This
                                // is ensured by the way MatrixFree reads
                                // out the indices.
    {
      internal::check_vector_compatibility (*dst[0], dof_info);
      Assert (n_fe_components == n_components, ExcNotImplemented());
      const unsigned int total_dofs_per_cell =
        dofs_per_cell * n_vectors * n_components;
      const Number * local_dst_number =
        reinterpret_cast<const Number*>(values_dofs[0]);
      if (at_irregular_cell == false)
        {
                                // check whether there is no constraint at all
          if (indicators != indicators_end)
            {
                                // run from one constraint to the next
              for ( ; indicators != indicators_end; ++indicators)
                {
                                // distribute values up to the constraint
                                // (values not constrained)
                  for (unsigned int j=0; j<indicators->first; ++j)
                    internal::vector_access (*dst[0], dof_indices[j])
                      += local_dst_number[ind_local+j];
                  dof_indices += indicators->first;
                  ind_local   += indicators->first;

                                // constrained case: build the local value as
                                // a linear combination of the global value
                                // according to constraint
                  const Number * data_val =
                    constraint_pool.begin(indicators->second);
                  const unsigned int row_length =
                    constraint_pool.row_length(indicators->second);
                  for (unsigned int k=0; k<row_length; ++k)
                    internal::vector_access (*dst[0], dof_indices[k])
                      += local_dst_number[ind_local] * data_val[k];
                  dof_indices += row_length;
                  ++ind_local;
                }
                                // distribute values after the last constraint
                                // (values not constrained)
              for(; ind_local<total_dofs_per_cell; ++dof_indices, ++ind_local)
                internal::vector_access (*dst[0], *dof_indices)
                  += local_dst_number[ind_local];
            }
                                // no constraint at all: loop bounds are
                                // known, compiler can unroll without checks
          else
            {
              AssertDimension (dof_info.end_indices(cell)-dof_indices,
                               static_cast<int>(total_dofs_per_cell));
              for (unsigned int j=0; j<total_dofs_per_cell; ++j)
                internal::vector_access (*dst[0], dof_indices[j])
                  += local_dst_number[j];
            }
          return;
        }

                                // irregular case
      Assert (n_irreg_components_filled > 0, ExcInternalError());
      for ( ; indicators != indicators_end; ++indicators)
        {
          for(unsigned int j=0; j<indicators->first; ++j)
            {
                                // non-constrained case: copy the data from
                                // the local vector to the global one.
              internal::vector_access (*dst[0], dof_indices[j])
                += local_dst_number[ind_local];
              ++ind_local;
              if (ind_local % n_vectors == n_irreg_components_filled)
                ind_local += n_vectors-n_irreg_components_filled;
            }
          dof_indices += indicators->first;

                                // constrained case: distribute according to
                                // the constraint
          const Number * data_val =
            constraint_pool.begin(indicators->second);
          const unsigned int row_length =
            constraint_pool.row_length(indicators->second);
          for (unsigned int k=0; k<row_length; ++k)
            {
              internal::vector_access (*dst[0], dof_indices[k])
                += local_dst_number[ind_local] * data_val[k];
            }
          dof_indices += row_length;
          ++ind_local;
          if (ind_local % n_vectors == n_irreg_components_filled)
            ind_local += n_vectors-n_irreg_components_filled;
        }
      for(; ind_local<total_dofs_per_cell; ++dof_indices)
        {
          Assert (dof_indices != dof_info.end_indices(cell),
                  ExcInternalError());

                                // non-constrained case
          internal::vector_access (*dst[0], *dof_indices)
              += local_dst_number[ind_local];
          ++ind_local;
          if (ind_local % n_vectors == n_irreg_components_filled)
            ind_local += n_vectors-n_irreg_components_filled;
        }
      Assert (dof_indices == dof_info.end_indices(cell),
              ExcInternalError());
    }
}



template <int dim, int dofs_per_cell_, int n_q_points_,
          int n_components, typename Number>
template<typename VectorType>
inline
void
FEEvaluationBase<dim,dofs_per_cell_,n_q_points_,n_components,Number>::
set_dof_values (VectorType &dst) const
{
  AssertDimension (n_components, n_fe_components);
                                // only need one component, but to avoid
                                // compiler warnings, use n_components copies
                                // here (but these will not be used)
  VectorType * dst_data [n_components];
  for (unsigned int d=0; d<n_components; ++d)
    dst_data[d] = &dst;
  set_dof_values (dst_data);
}



template <int dim, int dofs_per_cell_, int n_q_points_,
          int n_components, typename Number>
template<typename VectorType>
inline
void
FEEvaluationBase<dim,dofs_per_cell_,n_q_points_,n_components,Number>::
set_dof_values (std::vector<VectorType>  &dst,
                const unsigned int        first_index) const
{
  AssertIndexRange (first_index, dst.size());
  Assert (n_fe_components == 1, ExcNotImplemented());
  Assert ((n_fe_components == 1 ?
           (first_index+n_components <= dst.size()) : true),
          ExcIndexRange (first_index + n_components, 0, dst.size()));

  VectorType * dst_data [n_components];
  for (unsigned int comp=0; comp<n_components; ++comp)
    dst_data[comp] = &dst[comp+first_index];
  set_dof_values (dst_data);
}



template <int dim, int dofs_per_cell_, int n_q_points_,
          int n_components, typename Number>
template<typename VectorType>
inline
void
FEEvaluationBase<dim,dofs_per_cell_,n_q_points_,n_components,Number>::
set_dof_values (std::vector<VectorType*>  &dst,
                const unsigned int         first_index) const
{
  AssertIndexRange (first_index, dst.size());
  Assert (n_fe_components == 1, ExcNotImplemented());
  Assert ((n_fe_components == 1 ?
           (first_index+n_components <= dst.size()) : true),
          ExcIndexRange (first_index + n_components, 0, dst.size()));

  VectorType * dst_data [n_components];
  for (unsigned int comp=0; comp<n_components; ++comp)
    dst_data[comp] = dst[comp+first_index];
  set_dof_values (dst_data);
}



template <int dim, int dofs_per_cell_, int n_q_points_,
          int n_components, typename Number>
template<typename VectorType>
inline
void
FEEvaluationBase<dim,dofs_per_cell_,n_q_points_,n_components,Number>::
set_dof_values (VectorType * dst[]) const
{
  Assert (cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  Assert (dof_values_initialized==true,
          internal::ExcAccessToUninitializedField());

                                // loop over all local dofs. ind_local holds
                                // local number on cell, index iterates over
                                // the elements of index_local_to_global and
                                // glob_indices points to the global indices
                                // stored in index_local_to_global
  const unsigned int * dof_indices = dof_info.begin_indices(cell);
  const std::pair<unsigned short,unsigned short> * indicators =
    dof_info.begin_indicators(cell);
  const std::pair<unsigned short,unsigned short> * indicators_end =
    dof_info.end_indicators(cell);
  unsigned int ind_local = 0;

  if (n_fe_components == 1)
    {
      for (unsigned int comp=0; comp<n_components; ++comp)
        AssertDimension (dst[comp]->size(),
                         dof_info.vector_partitioner->size());

      const Number * local_dst_number [n_components];
      for (unsigned int comp=0; comp<n_components; ++comp)
        local_dst_number[comp] =
          reinterpret_cast<const Number*>(values_dofs[comp]);
      if (at_irregular_cell == false)
        {
                                // check whether there is no constraint at all
          if (indicators != indicators_end)
            {
                                // run from one constraint to the next
              for ( ; indicators != indicators_end; ++indicators)
                {
                                // distribute values up to the constraint
                                // (values not constrained)
                  for (unsigned int j=0; j<indicators->first; ++j)
                    for (unsigned int comp=0; comp<n_components; ++comp)
                      internal::vector_access (*dst[comp], dof_indices[j])
                        = local_dst_number[comp][ind_local+j];
                  dof_indices += indicators->first;
                  ind_local   += indicators->first;

                                // jump over constraints
                  const unsigned int row_length =
                    constraint_pool.row_length(indicators->second);
                  dof_indices += row_length;
                  ++ind_local;
                }
                                // distribute values after the last constraint
                                // (values not constrained)
              for(; ind_local<dofs_per_cell*n_vectors; ++dof_indices, ++ind_local)
                for (unsigned int comp=0; comp<n_components; ++comp)
                  internal::vector_access (*dst[comp], *dof_indices)
                    = local_dst_number[comp][ind_local];
            }
                                // no constraint at all: loop bounds are
                                // known, compiler can unroll without checks
          else
            {
              AssertDimension (dof_info.end_indices(cell)-dof_indices,
                               n_vectors * dofs_per_cell);
              for (unsigned int j=0; j<dofs_per_cell*n_vectors; ++j)
                for (unsigned int comp=0; comp<n_components; ++comp)
                  internal::vector_access (*dst[comp], dof_indices[j])
                    = local_dst_number[comp][j];
            }
          return;
        }

                                // irregular case
      Assert (n_irreg_components_filled > 0, ExcInternalError());
      for ( ; indicators != indicators_end; ++indicators)
        {
          for(unsigned int j=0; j<indicators->first; ++j)
            {
                                // non-constrained case
              for (unsigned int comp=0; comp<n_components; ++comp)
                internal::vector_access (*dst[comp], dof_indices[j])
                  = local_dst_number[comp][ind_local];
              ++ind_local;
              if (ind_local % n_vectors == n_irreg_components_filled)
                ind_local += n_vectors-n_irreg_components_filled;
            }
          dof_indices += indicators->first;

                                // jump over constraint
          const unsigned int row_length =
            constraint_pool.row_length(indicators->second);
          dof_indices += row_length;
          ++ind_local;
          if (ind_local % n_vectors == n_irreg_components_filled)
            ind_local += n_vectors-n_irreg_components_filled;
        }
      for(; ind_local<dofs_per_cell*n_vectors; ++dof_indices)
        {
          Assert (dof_indices != dof_info.end_indices(cell),
                  ExcInternalError());
                                // non-constrained case
          for (unsigned int comp=0; comp<n_components; ++comp)
            internal::vector_access (*dst[comp], *dof_indices)
              = local_dst_number[comp][ind_local];
          ++ind_local;
          if (ind_local % n_vectors == n_irreg_components_filled)
            ind_local += n_vectors-n_irreg_components_filled;
        }
    }
  else
                                // case with vector-valued finite elements
                                // where all components are included in one
                                // single vector. Assumption: first come all
                                // entries to the first component, then all
                                // entries to the second one, and so on. This
                                // is ensured by the way MatrixFree reads
                                // out the indices.
    {
      AssertDimension (dst[0]->size(),
                       dof_info.vector_partitioner->size());
      Assert (n_fe_components == n_components, ExcNotImplemented());
      const unsigned int total_dofs_per_cell =
        dofs_per_cell * n_vectors * n_components;
      const Number * local_dst_number =
        reinterpret_cast<const Number*>(values_dofs[0]);

      if (at_irregular_cell == false)
        {
                                // check whether there is no constraint at all
          if (indicators != indicators_end)
            {
                                // run from one constraint to the next
              for ( ; indicators != indicators_end; ++indicators)
                {
                                // distribute values up to the constraint
                                // (values not constrained)
                  for (unsigned int j=0; j<indicators->first; ++j)
                    internal::vector_access (*dst[0], dof_indices[j])
                      = local_dst_number[ind_local+j];
                  dof_indices += indicators->first;
                  ind_local   += indicators->first;

                                // jump over constraints
                  const unsigned int row_length =
                    constraint_pool.row_length(indicators->second);
                  dof_indices += row_length;
                  ++ind_local;
                }
                                // distribute values after the last constraint
                                // (values not constrained)
              for(; ind_local<total_dofs_per_cell; ++dof_indices, ++ind_local)
                internal::vector_access (*dst[0], *dof_indices)
                  = local_dst_number[ind_local];
            }
                                // no constraint at all: loop bounds are
                                // known, compiler can unroll without checks
          else
            {
              AssertDimension (dof_info.end_indices(cell)-dof_indices,
                               total_dofs_per_cell);
              for (unsigned int j=0; j<total_dofs_per_cell; ++j)
                internal::vector_access (*dst[0], dof_indices[j])
                  = local_dst_number[j];
            }
          return;
        }

                                // irregular case
      Assert (n_irreg_components_filled > 0, ExcInternalError());
      for ( ; indicators != indicators_end; ++indicators)
        {
          for(unsigned int j=0; j<indicators->first; ++j)
            {
                                // non-constrained case
              internal::vector_access (*dst[0], dof_indices[j])
                = local_dst_number[ind_local];
              ++ind_local;
              if (ind_local % n_vectors == n_irreg_components_filled)
                ind_local += n_vectors-n_irreg_components_filled;
            }
          dof_indices += indicators->first;

                                // jump over constraint
          const unsigned int row_length =
            constraint_pool.row_length(indicators->second);
          dof_indices += row_length;
          ++ind_local;
          if (ind_local % n_vectors == n_irreg_components_filled)
            ind_local += n_vectors-n_irreg_components_filled;
        }
      for(; ind_local<total_dofs_per_cell; ++dof_indices)
        {
          Assert (dof_indices != dof_info.end_indices(cell),
                  ExcInternalError());
                                // non-constrained case
          internal::vector_access (*dst[0], *dof_indices)
            = local_dst_number[ind_local];
          ++ind_local;
          if (ind_local % n_vectors == n_irreg_components_filled)
            ind_local += n_vectors-n_irreg_components_filled;
        }
      Assert (dof_indices == dof_info.end_indices (cell),
              ExcInternalError());
    }
}



// ------------------------------ access to data fields ---------------------

template <int dim, int dofs_per_cell_, int n_q_points_,
          int n_components, typename Number>
inline
Tensor<1,n_components,VectorizedArray<Number> >
FEEvaluationBase<dim,dofs_per_cell_,n_q_points_,n_components,Number>::
get_dof_value (unsigned int dof) const
{
  AssertIndexRange (dof, dofs_per_cell);
  Tensor<1,n_components,vector_t> return_value (false);
  for(unsigned int comp=0;comp<n_components;comp++)
    return_value[comp] = this->values_dofs[comp][dof];
  return return_value;
}



template <int dim, int dofs_per_cell_, int n_q_points_,
          int n_components, typename Number>
inline
Tensor<1,n_components,VectorizedArray<Number> >
FEEvaluationBase<dim,dofs_per_cell_,n_q_points_,n_components,Number>::
get_value (unsigned int q_point) const
{
  Assert (this->values_quad_initialized==true,
          internal::ExcAccessToUninitializedField());
  AssertIndexRange (q_point, n_q_points);
  Tensor<1,n_components,vector_t> return_value (false);
  for(unsigned int comp=0;comp<n_components;comp++)
    return_value[comp] = this->values_quad[comp][q_point];
  return return_value;
}



template <int dim, int dofs_per_cell_, int n_q_points_,
          int n_components, typename Number>
inline
Tensor<1,n_components,Tensor<1,dim,VectorizedArray<Number> > >
FEEvaluationBase<dim,dofs_per_cell_,n_q_points_,n_components,Number>::
get_gradient (unsigned int q_point) const
{
  Assert (this->gradients_quad_initialized==true,
          internal::ExcAccessToUninitializedField());
  AssertIndexRange (q_point, n_q_points);

  Tensor<1,n_components,Tensor<1,dim,vector_t> > grad_out (false);

                                // Cartesian cell
  if (this->cell_type == 0)
    {
      for (unsigned int comp=0;comp<n_components;comp++)
        for (unsigned int d=0; d<dim; ++d)
          grad_out[comp][d] = (this->gradients_quad[comp][d][q_point] *
                               cartesian[0][d]);
    }
                                // cell with general Jacobian
  else if (this->cell_type == 2)
    {
      for(unsigned int comp=0;comp<n_components;comp++)
        {
          for (unsigned int d=0; d<dim; ++d)
            {
              grad_out[comp][d] = (jacobian[q_point][d][0] *
                                   this->gradients_quad[comp][0][q_point]);
              for (unsigned e=1; e<dim; ++e)
                grad_out[comp][d] += (jacobian[q_point][d][e] *
                                      this->gradients_quad[comp][e][q_point]);
            }
        }
    }
                                // cell with general Jacobian, but constant
                                // within the cell
  else // if (this->cell_type == 1)
    {
      for(unsigned int comp=0;comp<n_components;comp++)
        {
          for (unsigned int d=0; d<dim; ++d)
            {
              grad_out[comp][d] = (jacobian[0][d][0] *
                                   this->gradients_quad[comp][0][q_point]);
              for (unsigned e=1; e<dim; ++e)
                grad_out[comp][d] += (jacobian[0][d][e] *
                                      this->gradients_quad[comp][e][q_point]);
            }
        }
    }
  return grad_out;
}



template <int dim, int dofs_per_cell_, int n_q_points_,
          int n_components, typename Number>
inline
Tensor<1,n_components,Tensor<2,dim,VectorizedArray<Number> > >
FEEvaluationBase<dim,dofs_per_cell_,n_q_points_,n_components,Number>::
get_hessian (unsigned int q_point) const
{
  Assert (this->hessians_quad_initialized==true,
          internal::ExcAccessToUninitializedField());
  AssertIndexRange (q_point, n_q_points);

  Tensor<2,dim,vector_t> hessian_out [n_components];

                                // Cartesian cell
  if (this->cell_type == 0)
    {
      const Tensor<1,dim,vector_t> &jac = cartesian[0];
      for (unsigned int comp=0;comp<n_components;comp++)
        for (unsigned int d=0; d<dim; ++d)
          {
            hessian_out[comp][d][d] = (this->hessians_quad[comp][d][q_point] *
                                       jac[d] * jac[d]);
            switch (dim)
              {
              case 1: break;
              case 2:
                hessian_out[comp][0][1] = (this->hessians_quad[comp][2][q_point] *
                                           jac[0] * jac[1]);
                break;
              case 3:
                hessian_out[comp][0][1] = (this->hessians_quad[comp][3][q_point] *
                                           jac[0] * jac[1]);
                hessian_out[comp][0][2] = (this->hessians_quad[comp][4][q_point] *
                                           jac[0] * jac[2]);
                hessian_out[comp][1][2] = (this->hessians_quad[comp][5][q_point] *
                                           jac[1] * jac[2]);
                break;
              default: Assert (false, ExcNotImplemented());
              }
            for (unsigned int e=d+1; e<dim; ++e)
              hessian_out[comp][e][d] = hessian_out[comp][d][e];
          }
    }
                                // cell with general Jacobian
  else if (this->cell_type == 2)
    {
      Assert (this->mapping_info.second_derivatives_initialized == true,
              ExcNotInitialized());
      const Tensor<2,dim,vector_t> & jac = jacobian[q_point];
      const Tensor<2,dim,vector_t> & jac_grad = jacobian_grad[q_point];
      const typename internal::MatrixFreeFunctions::MappingInfo<dim,Number>::tensorUT
        & jac_grad_UT = jacobian_grad_upper[q_point];
      for(unsigned int comp=0;comp<n_components;comp++)
        {
                                // compute laplacian before the gradient
                                // because it needs to access unscaled
                                // gradient data
          vector_t tmp[dim][dim];

                                // compute tmp = hess_unit(u) * J^T. do this
                                // manually because we do not store the lower
                                // diagonal because of symmetry
          for (unsigned int d=0; d<dim; ++d)
            {
              switch (dim)
                {
                case 1:
                  tmp[0][0] = jac[0][0] * this->hessians_quad[comp][0][q_point];
                  break;
                case 2:
                  tmp[0][d] = (jac[d][0] * this->hessians_quad[comp][0][q_point] +
                               jac[d][1] * this->hessians_quad[comp][2][q_point]);
                  tmp[1][d] = (jac[d][0] * this->hessians_quad[comp][2][q_point] +
                               jac[d][1] * this->hessians_quad[comp][1][q_point]);
                  break;
                case 3:
                  tmp[0][d] = (jac[d][0] * this->hessians_quad[comp][0][q_point] +
                               jac[d][1] * this->hessians_quad[comp][3][q_point] +
                               jac[d][2] * this->hessians_quad[comp][4][q_point]);
                  tmp[1][d] = (jac[d][0] * this->hessians_quad[comp][3][q_point] +
                               jac[d][1] * this->hessians_quad[comp][1][q_point] +
                               jac[d][2] * this->hessians_quad[comp][5][q_point]);
                  tmp[2][d] = (jac[d][0] * this->hessians_quad[comp][4][q_point] +
                               jac[d][1] * this->hessians_quad[comp][5][q_point] +
                               jac[d][2] * this->hessians_quad[comp][2][q_point]);
                  break;
                default: Assert (false, ExcNotImplemented());
                }
            }
                                // compute first part of hessian,
                                // J * tmp = J * hess_unit(u) * J^T
          for (unsigned int d=0; d<dim; ++d)
            for (unsigned int e=d; e<dim; ++e)
              {
                hessian_out[comp][d][e] = jac[d][0] * tmp[0][e];
                for (unsigned int f=1; f<dim; ++f)
                  hessian_out[comp][d][e] += jac[d][f] * tmp[f][e];
              }

                                // add diagonal part of J' * grad(u)
          for (unsigned int d=0; d<dim; ++d)
            for (unsigned int e=0; e<dim; ++e)
              hessian_out[comp][d][d] += (jac_grad[d][e] *
                                          this->gradients_quad[comp][e][q_point]);

                                // add off-diagonal part of J' * grad(u)
          for (unsigned int d=0, count=0; d<dim; ++d)
            for (unsigned int e=d+1; e<dim; ++e, ++count)
              for (unsigned int f=0; f<dim; ++f)
                hessian_out[comp][d][e] += (jac_grad_UT[count][f] *
                                            this->gradients_quad[comp][f][q_point]);

                                // take symmetric part
          for (unsigned int d=0; d<dim; ++d)
            for (unsigned int e=d+1; e<dim; ++e)
              hessian_out[comp][e][d] = hessian_out[comp][d][e];
        }
    }
                                // cell with general Jacobian, but constant
                                // within the cell
  else // if (this->cell_type == 1)
    {
      const Tensor<2,dim,vector_t> &jac = jacobian[0];
      for(unsigned int comp=0;comp<n_components;comp++)
        {
                                // compute laplacian before the gradient
                                // because it needs to access unscaled
                                // gradient data
          vector_t tmp[dim][dim];

                                // compute tmp = hess_unit(u) * J^T. do this
                                // manually because we do not store the lower
                                // diagonal because of symmetry
          for (unsigned int d=0; d<dim; ++d)
            {
              switch (dim)
                {
                case 1:
                  tmp[0][0] = jac[0][0] * this->hessians_quad[comp][0][q_point];
                  break;
                case 2:
                  tmp[0][d] = (jac[d][0] * this->hessians_quad[comp][0][q_point] +
                               jac[d][1] * this->hessians_quad[comp][2][q_point]);
                  tmp[1][d] = (jac[d][0] * this->hessians_quad[comp][2][q_point] +
                               jac[d][1] * this->hessians_quad[comp][1][q_point]);
                  break;
                case 3:
                  tmp[0][d] = (jac[d][0] * this->hessians_quad[comp][0][q_point] +
                               jac[d][1] * this->hessians_quad[comp][3][q_point] +
                               jac[d][2] * this->hessians_quad[comp][4][q_point]);
                  tmp[1][d] = (jac[d][0] * this->hessians_quad[comp][3][q_point] +
                               jac[d][1] * this->hessians_quad[comp][1][q_point] +
                               jac[d][2] * this->hessians_quad[comp][5][q_point]);
                  tmp[2][d] = (jac[d][0] * this->hessians_quad[comp][4][q_point] +
                               jac[d][1] * this->hessians_quad[comp][5][q_point] +
                               jac[d][2] * this->hessians_quad[comp][2][q_point]);
                  break;
                default: Assert (false, ExcNotImplemented());
                }
            }
                                // compute first part of hessian,
                                // J * tmp = J * hess_unit(u) * J^T
          for (unsigned int d=0; d<dim; ++d)
            for (unsigned int e=d; e<dim; ++e)
              {
                hessian_out[comp][d][e] = jac[d][0] * tmp[0][e];
                for (unsigned int f=1; f<dim; ++f)
                  hessian_out[comp][d][e] += jac[d][f] * tmp[f][e];
              }

                                // no J' * grad(u) part here because the
                                // Jacobian is constant throughout the cell
                                // and hence, its derivative is zero

                                // take symmetric part
          for (unsigned int d=0; d<dim; ++d)
            for (unsigned int e=d+1; e<dim; ++e)
              hessian_out[comp][e][d] = hessian_out[comp][d][e];
        }
    }
  return Tensor<1,n_components,Tensor<2,dim,vector_t> >(hessian_out);
}



template <int dim, int dofs_per_cell_, int n_q_points_,
          int n_components, typename Number>
inline
Tensor<1,n_components,Tensor<1,dim,VectorizedArray<Number> > >
FEEvaluationBase<dim,dofs_per_cell_,n_q_points_,n_components,Number>::
get_hessian_diagonal (unsigned int q_point) const
{
  Assert (this->hessians_quad_initialized==true,
          internal::ExcAccessToUninitializedField());
  AssertIndexRange (q_point, n_q_points);

  Tensor<1,n_components,Tensor<1,dim,vector_t> > hessian_out (false);

                                // Cartesian cell
  if (this->cell_type == 0)
    {
      const Tensor<1,dim,vector_t> &jac = cartesian[0];
      for (unsigned int comp=0;comp<n_components;comp++)
        for (unsigned int d=0; d<dim; ++d)
          hessian_out[comp][d] = (this->hessians_quad[comp][d][q_point] *
                                  jac[d] * jac[d]);
    }
                                // cell with general Jacobian
  else if (this->cell_type == 2)
    {
      Assert (this->mapping_info.second_derivatives_initialized == true,
              ExcNotInitialized());
      const Tensor<2,dim,vector_t> &jac = jacobian[q_point];
      const Tensor<2,dim,vector_t> &jac_grad = jacobian_grad[q_point];
      for(unsigned int comp=0;comp<n_components;comp++)
        {
                                // compute laplacian before the gradient
                                // because it needs to access unscaled
                                // gradient data
          vector_t tmp[dim][dim];

                                // compute tmp = hess_unit(u) * J^T. do this
                                // manually because we do not store the lower
                                // diagonal because of symmetry
          for (unsigned int d=0; d<dim; ++d)
            {
              switch (dim)
                {
                case 1:
                  tmp[0][0] = jac[0][0] * this->hessians_quad[comp][0][q_point];
                  break;
                case 2:
                  tmp[0][d] = (jac[d][0] * this->hessians_quad[comp][0][q_point] +
                               jac[d][1] * this->hessians_quad[comp][2][q_point]);
                  tmp[1][d] = (jac[d][0] * this->hessians_quad[comp][2][q_point] +
                               jac[d][1] * this->hessians_quad[comp][1][q_point]);
                  break;
                case 3:
                  tmp[0][d] = (jac[d][0] * this->hessians_quad[comp][0][q_point] +
                               jac[d][1] * this->hessians_quad[comp][3][q_point] +
                               jac[d][2] * this->hessians_quad[comp][4][q_point]);
                  tmp[1][d] = (jac[d][0] * this->hessians_quad[comp][3][q_point] +
                               jac[d][1] * this->hessians_quad[comp][1][q_point] +
                               jac[d][2] * this->hessians_quad[comp][5][q_point]);
                  tmp[2][d] = (jac[d][0] * this->hessians_quad[comp][4][q_point] +
                               jac[d][1] * this->hessians_quad[comp][5][q_point] +
                               jac[d][2] * this->hessians_quad[comp][2][q_point]);
                  break;
                default: Assert (false, ExcNotImplemented());
                }
            }
                                // compute only the trace part of hessian,
                                // J * tmp = J * hess_unit(u) * J^T
          for (unsigned int d=0; d<dim; ++d)
            {
              hessian_out[comp][d] = jac[d][0] * tmp[0][d];
              for (unsigned int f=1; f<dim; ++f)
                hessian_out[comp][d] += jac[d][f] * tmp[f][d];
            }

          for (unsigned int d=0; d<dim; ++d)
            for (unsigned int e=0; e<dim; ++e)
              hessian_out[comp][d] += (jac_grad[d][e] *
                                      this->gradients_quad[comp][e][q_point]);
        }
    }
                                // cell with general Jacobian, but constant
                                // within the cell
  else // if (this->cell_type == 1)
    {
      const Tensor<2,dim,vector_t> & jac = jacobian[0];
      for(unsigned int comp=0;comp<n_components;comp++)
        {
                                // compute laplacian before the gradient
                                // because it needs to access unscaled
                                // gradient data
          vector_t tmp[dim][dim];

                                // compute tmp = hess_unit(u) * J^T. do this
                                // manually because we do not store the lower
                                // diagonal because of symmetry
          for (unsigned int d=0; d<dim; ++d)
            {
              switch (dim)
                {
                case 1:
                  tmp[0][0] = jac[0][0] * this->hessians_quad[comp][0][q_point];
                  break;
                case 2:
                  tmp[0][d] = (jac[d][0] * this->hessians_quad[comp][0][q_point] +
                               jac[d][1] * this->hessians_quad[comp][2][q_point]);
                  tmp[1][d] = (jac[d][0] * this->hessians_quad[comp][2][q_point] +
                               jac[d][1] * this->hessians_quad[comp][1][q_point]);
                  break;
                case 3:
                  tmp[0][d] = (jac[d][0] * this->hessians_quad[comp][0][q_point] +
                               jac[d][1] * this->hessians_quad[comp][3][q_point] +
                               jac[d][2] * this->hessians_quad[comp][4][q_point]);
                  tmp[1][d] = (jac[d][0] * this->hessians_quad[comp][3][q_point] +
                               jac[d][1] * this->hessians_quad[comp][1][q_point] +
                               jac[d][2] * this->hessians_quad[comp][5][q_point]);
                  tmp[2][d] = (jac[d][0] * this->hessians_quad[comp][4][q_point] +
                               jac[d][1] * this->hessians_quad[comp][5][q_point] +
                               jac[d][2] * this->hessians_quad[comp][2][q_point]);
                  break;
                default: Assert (false, ExcNotImplemented());
                }
            }
                                // compute only the trace part of hessian,
                                // J * tmp = J * hess_unit(u) * J^T
          for (unsigned int d=0; d<dim; ++d)
            {
              hessian_out[comp][d] = jac[d][0] * tmp[0][d];
              for (unsigned int f=1; f<dim; ++f)
                hessian_out[comp][d] += jac[d][f] * tmp[f][d];
            }
        }
    }
  return hessian_out;
}



template <int dim, int dofs_per_cell_, int n_q_points_,
          int n_components, typename Number>
inline
Tensor<1,n_components,VectorizedArray<Number> >
FEEvaluationBase<dim,dofs_per_cell_,n_q_points_,n_components,Number>::
get_laplacian (unsigned int q_point) const
{
  Assert (this->hessians_quad_initialized==true,
          internal::ExcAccessToUninitializedField());
  AssertIndexRange (q_point, n_q_points);
  Tensor<1,n_components,vector_t> laplacian_out (false);
  const Tensor<1,n_components,Tensor<1,dim,vector_t> > hess_diag
    = get_hessian_diagonal(q_point);
  for (unsigned int comp=0; comp<n_components; ++comp)
    {
      laplacian_out[comp] = hess_diag[comp][0];
      for (unsigned int d=1; d<dim; ++d)
        laplacian_out[comp] += hess_diag[comp][d];
    }
  return laplacian_out;
}



template <int dim, int dofs_per_cell_, int n_q_points_,
          int n_components, typename Number>
inline
void
FEEvaluationBase<dim,dofs_per_cell_,n_q_points_,n_components,Number>::
submit_dof_value (Tensor<1,n_components,VectorizedArray<Number> > val_in,
                  unsigned int dof)
{
#ifdef DEBUG
  this->dof_values_initialized = true;
#endif
  AssertIndexRange (dof, dofs_per_cell);
  for (unsigned int comp=0;comp<n_components;comp++)
    this->values_dofs[comp][dof] = val_in[comp];
}



template <int dim, int dofs_per_cell_, int n_q_points_,
          int n_components, typename Number>
inline
void
FEEvaluationBase<dim,dofs_per_cell_,n_q_points_,n_components,Number>::
submit_value (Tensor<1,n_components,VectorizedArray<Number> > val_in,
              unsigned int q_point)
{
#ifdef DEBUG
  Assert (this->cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  AssertIndexRange (q_point, n_q_points);
  this->values_quad_submitted = true;
#endif
  if (this->cell_type == 2)
    {
      const vector_t JxW = J_value[q_point];
      for (unsigned int comp=0; comp<n_components; ++comp)
        this->values_quad[comp][q_point] = val_in[comp] * JxW;
    }
  else //if (this->cell_type < 2)
    {
      const vector_t JxW = J_value[0] * quadrature_weights[q_point];
      for (unsigned int comp=0; comp<n_components; ++comp)
        this->values_quad[comp][q_point] = val_in[comp] * JxW;
    }
}



template <int dim, int dofs_per_cell_, int n_q_points_,
          int n_components, typename Number>
inline
void
FEEvaluationBase<dim,dofs_per_cell_,n_q_points_,n_components,Number>::
submit_gradient (Tensor<1,n_components,
                        Tensor<1,dim,VectorizedArray<Number> > > grad_in,
                 unsigned int q_point)
{
#ifdef DEBUG
  Assert (this->cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  AssertIndexRange (q_point, n_q_points);
  this->gradients_quad_submitted = true;
#endif
  if (this->cell_type == 0)
    {
      const vector_t JxW = J_value[0] * quadrature_weights[q_point];
      for (unsigned int comp=0;comp<n_components;comp++)
        for (unsigned int d=0; d<dim; ++d)
          this->gradients_quad[comp][d][q_point] = (grad_in[comp][d] *
                                                    cartesian[0][d] * JxW);
    }
  else if (this->cell_type == 2)
    {
      for (unsigned int comp=0; comp<n_components; ++comp)
        for (unsigned int d=0; d<dim; ++d)
          {
            vector_t new_val = jacobian[q_point][0][d] * grad_in[comp][0];
            for (unsigned e=1; e<dim; ++e)
              new_val += jacobian[q_point][e][d] * grad_in[comp][e];
            this->gradients_quad[comp][d][q_point] = new_val * J_value[q_point];
          }
    }
  else //if (this->cell_type == 1)
    {
      const vector_t JxW = J_value[0] * quadrature_weights[q_point];
      for (unsigned int comp=0; comp<n_components; ++comp)
        for (unsigned int d=0; d<dim; ++d)
          {
            vector_t new_val = jacobian[0][0][d] * grad_in[comp][0];
            for (unsigned e=1; e<dim; ++e)
              new_val += jacobian[0][e][d] * grad_in[comp][e];
            this->gradients_quad[comp][d][q_point] = new_val * JxW;
          }
    }
}



template <int dim, int dofs_per_cell_, int n_q_points_,
          int n_components, typename Number>
inline
Tensor<1,n_components,VectorizedArray<Number> >
FEEvaluationBase<dim,dofs_per_cell_,n_q_points_,n_components,Number>::
integrate_value ()
{
#ifdef DEBUG
  Assert (this->cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  Assert (this->values_quad_submitted == true,
          internal::ExcAccessToUninitializedField());
#endif
  Tensor<1,n_components,vector_t> return_value (false);
  for (unsigned int comp=0; comp<n_components; ++comp)
    return_value[comp] = this->values_quad[comp][0];
  for (unsigned int q=0; q<n_q_points; ++q)
    for (unsigned int comp=0; comp<n_components; ++comp)
      return_value[comp] += this->values_quad[comp][q];
  return (return_value);
}



/*----------------------- FEEvaluationAccess -------------------------------*/


template <int dim, int dofs_per_cell_, int n_q_points_,
          int n_components, typename Number>
inline
FEEvaluationAccess<dim,dofs_per_cell_,n_q_points_,n_components,Number>::
FEEvaluationAccess (const MatrixFree<dim,Number> &data_in,
                    const unsigned int fe_no,
                    const unsigned int quad_no_in)
  :
  FEEvaluationBase <dim,dofs_per_cell_,n_q_points_,n_components,Number>
  (data_in, fe_no, quad_no_in)
{}




/*-------------------- FEEvaluationAccess scalar --------------------------*/


template <int dim, int dofs_per_cell_,  int n_q_points_, typename Number>
inline
FEEvaluationAccess<dim,dofs_per_cell_,n_q_points_,1,Number>::
FEEvaluationAccess (const MatrixFree<dim,Number> &data_in,
                    const unsigned int fe_no,
                    const unsigned int quad_no_in)
  :
  FEEvaluationBase <dim,dofs_per_cell_,n_q_points_,1,Number>
  (data_in, fe_no, quad_no_in)
{}



template <int dim, int dofs_per_cell_,  int n_q_points_, typename Number>
inline
VectorizedArray<Number>
FEEvaluationAccess<dim,dofs_per_cell_,n_q_points_,1,Number>::
get_dof_value (unsigned int dof) const
{
  AssertIndexRange (dof, dofs_per_cell);
  return this->values_dofs[0][dof];
}



template <int dim, int dofs_per_cell_,  int n_q_points_, typename Number>
inline
VectorizedArray<Number>
FEEvaluationAccess<dim,dofs_per_cell_,n_q_points_,1,Number>::
get_value (unsigned int q_point) const
{
  Assert (this->values_quad_initialized==true,
          internal::ExcAccessToUninitializedField());
  AssertIndexRange (q_point, n_q_points);
  return this->values_quad[0][q_point];
}



template <int dim, int dofs_per_cell_,  int n_q_points_, typename Number>
inline
Tensor<1,dim,VectorizedArray<Number> >
FEEvaluationAccess<dim,dofs_per_cell_,n_q_points_,1,Number>::
get_gradient (unsigned int q_point) const
{
                                // could use the base class gradient, but that
                                // involves too many inefficient
                                // initialization

  Assert (this->gradients_quad_initialized==true,
          internal::ExcAccessToUninitializedField());
  AssertIndexRange (q_point, n_q_points);

  Tensor<1,dim,vector_t> grad_out (false);

                                // Cartesian cell
  if (this->cell_type == 0)
    {
      for (unsigned int d=0; d<dim; ++d)
        grad_out[d] = (this->gradients_quad[0][d][q_point] *
                       this->cartesian[0][d]);
    }
                                // cell with general Jacobian
  else if (this->cell_type == 2)
    {
      for (unsigned int d=0; d<dim; ++d)
        {
          grad_out[d] = (this->jacobian[q_point][d][0] *
                         this->gradients_quad[0][0][q_point]);
          for (unsigned e=1; e<dim; ++e)
            grad_out[d] += (this->jacobian[q_point][d][e] *
                            this->gradients_quad[0][e][q_point]);
        }
    }
                                // cell with general Jacobian, but constant
                                // within the cell
  else // if (this->cell_type == 1)
    {
      for (unsigned int d=0; d<dim; ++d)
        {
          grad_out[d] = (this->jacobian[0][d][0] *
                         this->gradients_quad[0][0][q_point]);
          for (unsigned e=1; e<dim; ++e)
            grad_out[d] += (this->jacobian[0][d][e] *
                            this->gradients_quad[0][e][q_point]);
        }
    }
  return grad_out;
}



template <int dim, int dofs_per_cell_,  int n_q_points_, typename Number>
inline
Tensor<2,dim,VectorizedArray<Number> >
FEEvaluationAccess<dim,dofs_per_cell_,n_q_points_,1,Number>::
get_hessian (unsigned int q_point) const
{
  return BaseClass::get_hessian(q_point)[0];
}



template <int dim, int dofs_per_cell_,  int n_q_points_, typename Number>
inline
Tensor<1,dim,VectorizedArray<Number> >
FEEvaluationAccess<dim,dofs_per_cell_,n_q_points_,1,Number>::
get_hessian_diagonal (unsigned int q_point) const
{
  return BaseClass::get_hessian_diagonal(q_point)[0];
}



template <int dim, int dofs_per_cell_,  int n_q_points_, typename Number>
inline
VectorizedArray<Number>
FEEvaluationAccess<dim,dofs_per_cell_,n_q_points_,1,Number>::
get_laplacian (unsigned int q_point) const
{
  return BaseClass::get_laplacian(q_point)[0];
}



template <int dim, int dofs_per_cell_,  int n_q_points_, typename Number>
inline
void
FEEvaluationAccess<dim,dofs_per_cell_,n_q_points_,1,Number>::
submit_dof_value (VectorizedArray<Number> val_in,
                  unsigned int dof)
{
#ifdef DEBUG
  this->dof_values_initialized = true;
  AssertIndexRange (dof, dofs_per_cell);
#endif
  this->values_dofs[0][dof] = val_in;
}



template <int dim, int dofs_per_cell_,  int n_q_points_, typename Number>
inline
void
FEEvaluationAccess<dim,dofs_per_cell_,n_q_points_,1,Number>::
submit_value (VectorizedArray<Number> val_in,
              unsigned int q_point)
{
#ifdef DEBUG
  Assert (this->cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  AssertIndexRange (q_point, n_q_points);
  this->values_quad_submitted = true;
#endif
  if (this->cell_type == 2)
    {
      const vector_t JxW = this->J_value[q_point];
      this->values_quad[0][q_point] = val_in * JxW;
    }
  else //if (this->cell_type < 2)
    {
      const vector_t JxW = this->J_value[0] * this->quadrature_weights[q_point];
      this->values_quad[0][q_point] = val_in * JxW;
    }
}



template <int dim, int dofs_per_cell_,  int n_q_points_, typename Number>
inline
void
FEEvaluationAccess<dim,dofs_per_cell_,n_q_points_,1,Number>::
submit_gradient (Tensor<1,dim,VectorizedArray<Number> > grad_in,
                 unsigned int q_point)
{
#ifdef DEBUG
  Assert (this->cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  AssertIndexRange (q_point, n_q_points);
  this->gradients_quad_submitted = true;
#endif
  if (this->cell_type == 0)
    {
      const vector_t JxW = this->J_value[0] * this->quadrature_weights[q_point];
      for (unsigned int d=0; d<dim; ++d)
        this->gradients_quad[0][d][q_point] = (grad_in[d] *
                                               this->cartesian[0][d] * JxW);
    }
  else if (this->cell_type == 2)
    {
      for (unsigned int d=0; d<dim; ++d)
        {
          vector_t new_val = this->jacobian[q_point][0][d] * grad_in[0];
          for (unsigned e=1; e<dim; ++e)
            new_val += this->jacobian[q_point][e][d] * grad_in[e];
          this->gradients_quad[0][d][q_point] = new_val * this->J_value[q_point];
        }
    }
  else //if (this->cell_type == 1)
    {
      const vector_t JxW = this->J_value[0] * this->quadrature_weights[q_point];
      for (unsigned int d=0; d<dim; ++d)
        {
          vector_t new_val = this->jacobian[0][0][d] * grad_in[0];
          for (unsigned e=1; e<dim; ++e)
            new_val += this->jacobian[0][e][d] * grad_in[e];
          this->gradients_quad[0][d][q_point] = new_val * JxW;
        }
    }
}



template <int dim, int dofs_per_cell_,  int n_q_points_, typename Number>
inline
VectorizedArray<Number>
FEEvaluationAccess<dim,dofs_per_cell_,n_q_points_,1,Number>::
integrate_value ()
{
  return BaseClass::integrate_value()[0];
}




/*----------------- FEEvaluationAccess vector-valued ----------------------*/


template <int dim, int dofs_per_cell_,  int n_q_points_, typename Number>
inline
FEEvaluationAccess<dim,dofs_per_cell_,n_q_points_,dim,Number>::
FEEvaluationAccess (const MatrixFree<dim,Number> &data_in,
                    const unsigned int fe_no,
                    const unsigned int quad_no_in)
  :
  FEEvaluationBase <dim,dofs_per_cell_,n_q_points_,dim,Number>
  (data_in, fe_no, quad_no_in)
{}



template <int dim, int dofs_per_cell_,  int n_q_points_, typename Number>
inline
Tensor<2,dim,VectorizedArray<Number> >
FEEvaluationAccess<dim,dofs_per_cell_,n_q_points_,dim,Number>::
get_gradient (unsigned int q_point) const
{
  return BaseClass::get_gradient (q_point);
}



template <int dim, int dofs_per_cell_,  int n_q_points_, typename Number>
inline
VectorizedArray<Number>
FEEvaluationAccess<dim,dofs_per_cell_,n_q_points_,dim,Number>::
get_divergence (unsigned int q_point) const
{
  Assert (this->gradients_quad_initialized==true,
          internal::ExcAccessToUninitializedField());
  AssertIndexRange (q_point, n_q_points);

  vector_t divergence;

                                // Cartesian cell
  if (this->cell_type == 0)
    {
      divergence = (this->gradients_quad[0][0][q_point] *
                    this->cartesian[0][0]);
      for (unsigned int d=1; d<dim; ++d)
        divergence += (this->gradients_quad[d][d][q_point] *
                       this->cartesian[0][d]);
    }
                                // cell with general Jacobian
  else if (this->cell_type == 2)
    {
      divergence = (this->jacobian[q_point][0][0] *
                    this->gradients_quad[0][0][q_point]);
      for (unsigned e=1; e<dim; ++e)
        divergence += (this->jacobian[q_point][0][e] *
                       this->gradients_quad[0][e][q_point]);
      for (unsigned int d=1; d<dim; ++d)
        for (unsigned e=0; e<dim; ++e)
          divergence += (this->jacobian[q_point][d][e] *
                         this->gradients_quad[d][e][q_point]);
    }
                                // cell with general Jacobian, but constant
                                // within the cell
  else // if (this->cell_type == 1)
    {
      divergence = (this->jacobian[0][0][0] *
                    this->gradients_quad[0][0][q_point]);
      for (unsigned e=1; e<dim; ++e)
        divergence += (this->jacobian[0][0][e] *
                       this->gradients_quad[0][e][q_point]);
      for (unsigned int d=1; d<dim; ++d)
        for (unsigned e=0; e<dim; ++e)
          divergence += (this->jacobian[0][d][e] *
                         this->gradients_quad[d][e][q_point]);
    }
  return divergence;
}



template <int dim, int dofs_per_cell_,  int n_q_points_, typename Number>
inline
SymmetricTensor<2,dim,VectorizedArray<Number> >
FEEvaluationAccess<dim,dofs_per_cell_,n_q_points_,dim,Number>::
get_symmetric_gradient (unsigned int q_point) const
{
                                // copy from generic function into
                                // dim-specialization function
  const Tensor<2,dim,vector_t> grad = get_gradient(q_point);
  vector_t symmetrized [(dim*dim+dim)/2];
  vector_t half = make_vectorized_array (0.5);
  for (unsigned int d=0; d<dim; ++d)
    symmetrized[d] = grad[d][d];
  switch (dim)
    {
    case 1:
      break;
    case 2:
      symmetrized[2] = grad[0][1] + grad[1][0];
      symmetrized[2] *= half;
      break;
    case 3:
      symmetrized[3] = grad[0][1] + grad[1][0];
      symmetrized[3] *= half;
      symmetrized[4] = grad[0][2] + grad[2][0];
      symmetrized[4] *= half;
      symmetrized[5] = grad[1][2] + grad[2][1];
      symmetrized[5] *= half;
      break;
    default:
      Assert (false, ExcNotImplemented());
    }
  return SymmetricTensor<2,dim,vector_t> (symmetrized);
}



template <int dim, int dofs_per_cell_,  int n_q_points_, typename Number>
inline
typename FEEvaluationAccess<dim,dofs_per_cell_,n_q_points_,dim,Number>::curl_type
FEEvaluationAccess<dim,dofs_per_cell_,n_q_points_,dim,Number>::
get_curl (unsigned int q_point) const
{
                                // copy from generic function into
                                // dim-specialization function
  const Tensor<2,dim,vector_t> grad = get_gradient(q_point);
  curl_type curl;
  switch (dim)
    {
    case 1:
      Assert (false,
              ExcMessage("Computing the curl in 1d is not a useful operation"));
      break;
    case 2:
      curl[0] = grad[1][0] - grad[0][1];
      break;
    case 3:
      curl[0] = grad[2][1] - grad[1][2];
      curl[1] = grad[0][2] - grad[2][0];
      curl[2] = grad[1][0] - grad[0][1];
      break;
    default:
      Assert (false, ExcNotImplemented());
    }
  return curl;
}



template <int dim, int dofs_per_cell_,  int n_q_points_, typename Number>
inline
Tensor<2,dim,VectorizedArray<Number> >
FEEvaluationAccess<dim,dofs_per_cell_,n_q_points_,dim,Number>::
get_hessian_diagonal (unsigned int q_point) const
{
  Assert (this->hessians_quad_initialized==true,
          internal::ExcAccessToUninitializedField());
  AssertIndexRange (q_point, n_q_points);

  return BaseClass::get_hessian_diagonal (q_point);
}



template <int dim, int dofs_per_cell_,  int n_q_points_, typename Number>
inline
Tensor<3,dim,VectorizedArray<Number> >
FEEvaluationAccess<dim,dofs_per_cell_,n_q_points_,dim,Number>::
get_hessian (unsigned int q_point) const
{
  Assert (this->hessians_quad_initialized==true,
          internal::ExcAccessToUninitializedField());
  AssertIndexRange (q_point, n_q_points);
  return BaseClass::get_hessian(q_point);
}



template <int dim, int dofs_per_cell_,  int n_q_points_, typename Number>
inline
void
FEEvaluationAccess<dim,dofs_per_cell_,n_q_points_,dim,Number>::
submit_gradient (Tensor<2,dim,VectorizedArray<Number> > grad_in,
                 unsigned int q_point)
{
  BaseClass::submit_gradient (grad_in, q_point);
}



template <int dim, int dofs_per_cell_,  int n_q_points_, typename Number>
inline
void
FEEvaluationAccess<dim,dofs_per_cell_,n_q_points_,dim,Number>::
submit_gradient (Tensor<1,dim,Tensor<1,dim,VectorizedArray<Number> > > grad_in,
                 unsigned int q_point)
{
  BaseClass::submit_gradient(grad_in, q_point);
}



template <int dim, int dofs_per_cell_,  int n_q_points_, typename Number>
inline
void
FEEvaluationAccess<dim,dofs_per_cell_,n_q_points_,dim,Number>::
submit_symmetric_gradient (SymmetricTensor<2,dim,VectorizedArray<Number> >
                           sym_grad,
                           unsigned int q_point)
{
                                // could have used base class operator, but
                                // that involves some overhead which is
                                // inefficient. it is nice to have the
                                // symmetric tensor because that saves some
                                // operations
#ifdef DEBUG
  Assert (this->cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  AssertIndexRange (q_point, n_q_points);
  this->gradients_quad_submitted = true;
#endif
  if (this->cell_type == 0)
    {
      const vector_t JxW = this->J_value[0] * this->quadrature_weights[q_point];
      for (unsigned int d=0; d<dim; ++d)
        this->gradients_quad[d][d][q_point] = (sym_grad.access_raw_entry(d) *
                                               JxW *
                                               this->cartesian[0][d]);
      for (unsigned int e=0, counter=dim; e<dim; ++e)
        for (unsigned int d=e+1; d<dim; ++d, ++counter)
          {
            const vector_t value = sym_grad.access_raw_entry(counter) * JxW;
            this->gradients_quad[e][d][q_point] = (value *
                                                   this->cartesian[0][d]);
            this->gradients_quad[d][e][q_point] = (value *
                                                   this->cartesian[0][e]);
          }
    }
  else if (this->cell_type == 2)
    {
      vector_t weighted [dim][dim];
      {
        const vector_t JxW = this->J_value[q_point];
        for (unsigned int i=0; i<dim; ++i)
          weighted[i][i] = sym_grad.access_raw_entry(i) * JxW;
        for (unsigned int i=0, counter=dim; i<dim; ++i)
          for (unsigned int j=i+1; j<dim; ++j, ++counter)
            {
              const vector_t value = sym_grad.access_raw_entry(counter) * JxW;
              weighted[i][j] = value;
              weighted[j][i] = value;
            }
      }
      for (unsigned int comp=0; comp<dim; ++comp)
        for (unsigned int d=0; d<dim; ++d)
          {
            vector_t new_val = this->jacobian[q_point][0][d] * weighted[comp][0];
            for (unsigned e=1; e<dim; ++e)
              new_val += this->jacobian[q_point][e][d] * weighted[comp][e];
            this->gradients_quad[comp][d][q_point] = new_val;
          }
    }
  else //if (this->cell_type == 1)
    {
      vector_t weighted [dim][dim];
      {
        const vector_t JxW = (this->J_value[0] *
                              this->quadrature_weights[q_point]);
        for (unsigned int i=0; i<dim; ++i)
          weighted[i][i] = sym_grad.access_raw_entry(i) * JxW;
        for (unsigned int i=0, counter=dim; i<dim; ++i)
          for (unsigned int j=i+1; j<dim; ++j, ++counter)
            {
              const vector_t value = sym_grad.access_raw_entry(counter) * JxW;
              weighted[i][j] = value;
              weighted[j][i] = value;
            }
      }
      for (unsigned int comp=0; comp<dim; ++comp)
        for (unsigned int d=0; d<dim; ++d)
          {
            vector_t new_val = this->jacobian[q_point][0][d] * weighted[comp][0];
            for (unsigned e=1; e<dim; ++e)
              new_val += this->jacobian[q_point][e][d] * weighted[comp][e];
            this->gradients_quad[comp][d][q_point] = new_val;
          }
    }
}



template <int dim, int dofs_per_cell_,  int n_q_points_, typename Number>
inline
void
FEEvaluationAccess<dim,dofs_per_cell_,n_q_points_,dim,Number>::
submit_curl (curl_type    curl,
             unsigned int q_point)
{
  Tensor<2,dim,vector_t> grad;
  switch (dim)
    {
    case 1:
      Assert (false,
              ExcMessage("Testing by the curl in 1d is not a useful operation"));
      break;
    case 2:
      grad[1][0] = curl[0];
      grad[0][1] = -curl[0];
      break;
    case 3:
      grad[2][1] = curl[0];
      grad[1][2] = -curl[0];
      grad[0][2] = curl[1];
      grad[2][0] = -curl[1];
      grad[1][0] = curl[2];
      grad[0][1] = -curl[2];
      break;
    default:
      Assert (false, ExcNotImplemented());
    }
  submit_gradient (grad, q_point);
}



/*----------------------- FEEvaluationGeneral -------------------------------*/

template <int dim, int n_dofs_1d,  int n_q_points_1d, int n_components,
          typename Number>
inline
FEEvaluationGeneral<dim,n_dofs_1d,n_q_points_1d,n_components,Number>::
FEEvaluationGeneral (const MatrixFree<dim,Number> &data_in,
                     const unsigned int fe_no,
                     const unsigned int quad_no_in)
  :
  BaseClass (data_in, fe_no, quad_no_in)
{
#ifdef DEBUG
                                // print error message when the dimensions do
                                // not match. Propose a possible fix
  if (dofs_per_cell != this->data.dofs_per_cell ||
      n_q_points != this->data.n_q_points)
    {
      std::string message =
        "-------------------------------------------------------\n";
      message += "Illegal arguments in constructor/wrong template arguments!\n";
      message += "    Called -->   FEEvaluation<dim,";
      message += Utilities::int_to_string(n_dofs_1d) + ",";
      message += Utilities::int_to_string(n_q_points_1d) + ",Number>(data, ";
      message += Utilities::int_to_string(fe_no) + ", ";
      message += Utilities::int_to_string(quad_no_in) + ")\n";

                                // check whether some other vector component
                                // has the correct number of points
      unsigned int proposed_dof_comp = numbers::invalid_unsigned_int,
        proposed_quad_comp = numbers::invalid_unsigned_int;
      for (unsigned int no=0; no<this->matrix_info.n_components(); ++no)
        if (this->matrix_info.get_dof_info(no).dofs_per_cell[0] ==
            dofs_per_cell)
          {
            proposed_dof_comp = no;
            break;
          }
      for (unsigned int no=0; no<this->mapping_info.mapping_data_gen.size(); ++no)
        if (this->mapping_info.mapping_data_gen[no].n_q_points[this->active_quad_index] == n_q_points)
          {
            proposed_quad_comp = no;
            break;
          }
      if (proposed_dof_comp  != numbers::invalid_unsigned_int &&
          proposed_quad_comp != numbers::invalid_unsigned_int)
        {
          message += "Wrong vector component selection:\n";
          message += "    Did you mean FEEvaluation<dim,Number,";
          message += Utilities::int_to_string(n_dofs_1d) + ",";
          message += Utilities::int_to_string(n_q_points_1d) + ">(data, ";
          message += Utilities::int_to_string(proposed_dof_comp) + ", ";
          message += Utilities::int_to_string(proposed_quad_comp) + ")?\n";
          std::string correct_pos;
          if (proposed_dof_comp != fe_no)
            correct_pos = " ^ ";
          else
            correct_pos = "   ";
          if (proposed_quad_comp != quad_no_in)
            correct_pos += " ^\n";
          else
            correct_pos += "  \n";
          message += "                                                   " + correct_pos;
        }
                                // ok, did not find the numbers specified by
                                // the template arguments in the given
                                // list. Suggest correct template arguments
      const unsigned int proposed_n_dofs_1d = static_cast<unsigned int>(std::pow(1.001*this->data.dofs_per_cell,1./dim));
      const unsigned int proposed_n_q_points_1d = static_cast<unsigned int>(std::pow(1.001*this->data.n_q_points,1./dim));
      message += "Wrong template arguments:\n";
      message += "    Did you mean FEEvaluation<dim,";
      message += Utilities::int_to_string(proposed_n_dofs_1d) + ",";
      message += Utilities::int_to_string(proposed_n_q_points_1d);
      message += ",Number>(data, ";
      message += Utilities::int_to_string(fe_no) + ", ";
      message += Utilities::int_to_string(quad_no_in) + ")?\n";
      std::string correct_pos;
      if (proposed_n_dofs_1d != n_dofs_1d)
        correct_pos = " ^";
      else
        correct_pos = "  ";
      if (proposed_n_q_points_1d != n_q_points_1d)
        correct_pos += " ^\n";
      else
        correct_pos += "  \n";
      message += "                                 " + correct_pos;

      Assert (dofs_per_cell == this->data.dofs_per_cell &&
              n_q_points == this->data.n_q_points,
              ExcMessage(message));
    }
  AssertDimension (n_q_points,
                   this->mapping_info.mapping_data_gen[this->quad_no].
                   n_q_points[this->active_quad_index]);
  AssertDimension (dofs_per_cell * this->n_fe_components,
                   this->dof_info.dofs_per_cell[this->active_fe_index]);
#endif
}



template <int dim, int n_dofs_1d,  int n_q_points_1d, int n_components,
          typename Number>
inline
void
FEEvaluationGeneral<dim,n_dofs_1d,n_q_points_1d,n_components,Number>::
evaluate (bool evaluate_val, bool evaluate_grad, bool evaluate_lapl)
{
  Assert (this->cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  Assert (this->dof_values_initialized == true,
          internal::ExcAccessToUninitializedField());

  const vector_t * val  = this->data.shape_values.begin();
  const vector_t * grad = this->data.shape_gradients.begin();
  const vector_t * hess = this->data.shape_hessians.begin();

  for(unsigned int comp=0;comp<n_components;comp++)
    {
      vector_t temp1[n_dofs_1d > n_q_points_1d ? dofs_per_cell : n_q_points];
      vector_t temp2[n_dofs_1d > n_q_points_1d ? dofs_per_cell : n_q_points];

      if (dim == 3)
        {
          if (evaluate_grad == true)
            {
              // grad x
              apply_tensor_prod<0,true,false> (grad, this->values_dofs[comp], temp1);
              apply_tensor_prod<1,true,false> (val, temp1, temp2);
              apply_tensor_prod<2,true,false> (val, temp2, this->gradients_quad[comp][0]);
            }

          if (evaluate_lapl == true)
            {
              // grad xz
              if (evaluate_grad == false)
                {
                  apply_tensor_prod<0,true,false> (grad, this->values_dofs[comp], temp1);
                  apply_tensor_prod<1,true,false> (val, temp1, temp2);
                }
              apply_tensor_prod<2,true,false>(grad, temp2, this->hessians_quad[comp][4]);

                  // grad xy
              apply_tensor_prod<1,true,false>(grad, temp1, temp2);
              apply_tensor_prod<2,true,false> (val, temp2, this->hessians_quad[comp][3]);

              // grad xx
              apply_tensor_prod<0,true,false> (hess, this->values_dofs[comp], temp1);
              apply_tensor_prod<1,true,false> (val, temp1, temp2);
              apply_tensor_prod<2,true,false> (val, temp2, this->hessians_quad[comp][0]);
            }

          // grad y
          apply_tensor_prod<0,true,false> (val, this->values_dofs[comp], temp1);
          if (evaluate_grad == true)
            {
              apply_tensor_prod<1,true,false> (grad, temp1, temp2);
              apply_tensor_prod<2,true,false> (val, temp2, this->gradients_quad[comp][1]);
            }

          if (evaluate_lapl == true)
            {
              // grad yz
              if (evaluate_grad == false)
                apply_tensor_prod<1,true,false> (grad, temp1, temp2);
              apply_tensor_prod<2,true,false> (grad, temp2, this->hessians_quad[comp][5]);

              // grad yy
              apply_tensor_prod<1,true,false> (hess, temp1, temp2);
              apply_tensor_prod<2,true,false> (val, temp2, this->hessians_quad[comp][1]);
            }

          // grad z: can use the values applied in x direction stored in temp1
          apply_tensor_prod<1,true,false> (val, temp1, temp2);
          if (evaluate_grad == true)
            apply_tensor_prod<2,true,false> (grad, temp2, this->gradients_quad[comp][2]);

          // grad zz: can use the values applied in x and y direction stored in temp2
          if (evaluate_lapl == true)
            apply_tensor_prod<2,true,false> (hess, temp2, this->hessians_quad[comp][2]);

          // val: can use the values applied in x & y direction stored in temp2
          if (evaluate_val == true)
            apply_tensor_prod<2,true,false> (val, temp2, this->values_quad[comp]);
        }
      else if (dim == 2)
        {
          // grad x
          if (evaluate_grad == true)
            {
              apply_tensor_prod<0,true,false> (grad, this->values_dofs[comp], temp1);
              apply_tensor_prod<1,true,false> (val, temp1, this->gradients_quad[comp][0]);
            }
          if (evaluate_lapl == true)
            {
              // grad xy
              if (evaluate_grad == false)
                apply_tensor_prod<0,true,false> (grad, this->values_dofs[comp], temp1);
              apply_tensor_prod<1,true,false> (grad, temp1, this->hessians_quad[comp][2]);

              // grad xx
              apply_tensor_prod<0,true,false> (hess, this->values_dofs[comp], temp1);
              apply_tensor_prod<1,true,false> (val, temp1, this->hessians_quad[comp][0]);
            }

          // grad y
          apply_tensor_prod<0,true,false> (val, this->values_dofs[comp], temp1);
          if (evaluate_grad == true)
            apply_tensor_prod<1,true,false> (grad, temp1, this->gradients_quad[comp][1]);

          // grad yy
          if (evaluate_lapl == true)
            apply_tensor_prod<1,true,false> (hess, temp1, this->hessians_quad[comp][1]);

          // val: can use values applied in x
          if (evaluate_val == true)
            apply_tensor_prod<1,true,false> (val, temp1, this->values_quad[comp]);
        }
      else if (dim == 1)
        {
          if (evaluate_val == true)
            apply_tensor_prod<0,true,false> (val, this->values_dofs[comp],
                                             this->values_quad[comp]);
          if (evaluate_grad == true)
            apply_tensor_prod<0,true,false> (grad, this->values_dofs[comp],
                                             this->gradients_quad[comp][0]);
          if (evaluate_lapl == true)
            apply_tensor_prod<0,true,false> (hess, this->values_dofs[comp],
                                             this->hessians_quad[comp][0]);
        }
    }

#ifdef DEBUG
  if (evaluate_val == true)
    this->values_quad_initialized = true;
  if (evaluate_grad == true)
    this->gradients_quad_initialized = true;
  if (evaluate_lapl == true)
    this->hessians_quad_initialized  = true;
#endif
}



template <int dim, int n_dofs_1d,  int n_q_points_1d, int n_components,
          typename Number>
inline
void
FEEvaluationGeneral<dim,n_dofs_1d,n_q_points_1d,n_components,Number>::
integrate (bool integrate_val,bool integrate_grad)
{
#ifdef DEBUG
  Assert (this->cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  if (integrate_val == true)
    Assert (this->values_quad_submitted == true,
            internal::ExcAccessToUninitializedField());
  if (integrate_grad == true)
    Assert (this->gradients_quad_submitted == true,
            internal::ExcAccessToUninitializedField());
#endif

  const vector_t * val  = this->data.shape_values.begin();
  const vector_t * grad = this->data.shape_gradients.begin();

  for(unsigned int comp=0;comp<n_components;comp++)
    {
      vector_t temp1[n_dofs_1d > n_q_points_1d ? dofs_per_cell : n_q_points];
      vector_t temp2[n_dofs_1d > n_q_points_1d ? dofs_per_cell : n_q_points];

      if (dim == 3)
        {
          if (integrate_val == true)
            {
              // val
              apply_tensor_prod<0,false,false> (val, this->values_quad[comp], temp1);
            }
          if (integrate_grad == true)
            {
              // grad x: can sum to temporary value in temp1
              if (integrate_val == true)
                apply_tensor_prod<0,false,true> 
                  (grad, this->gradients_quad[comp][0],temp1);
              else
                apply_tensor_prod<0,false,false> 
                  (grad, this->gradients_quad[comp][0],temp1);
            }
          apply_tensor_prod<1,false,false> (val, temp1, temp2);
          if (integrate_grad == true)
            {
              // grad y: can sum to temporary x value in temp2
              apply_tensor_prod<0,false,false> (val, this->gradients_quad[comp][1], temp1);
              apply_tensor_prod<1,false,true> (grad, temp1, temp2);
            }
          apply_tensor_prod<2,false,false> (val, temp2, this->values_dofs[comp]);
          if (integrate_grad == true)
            {
              // grad z: can sum to temporary x and y value in output
              apply_tensor_prod<0,false,false> (val, this->gradients_quad[comp][2], temp1);
              apply_tensor_prod<1,false,false> (val, temp1, temp2);
              apply_tensor_prod<2,false,true> (grad, temp2, this->values_dofs[comp]);
            }
        }
      else if (dim == 2)
        {
          // val
          if (integrate_val == true)
            apply_tensor_prod<0,false,false> (val, this->values_quad[comp], temp1);
          if (integrate_grad == true)
            {
              //grad x
              if (integrate_val == true)
                apply_tensor_prod<0,false,true> 
                  (grad, this->gradients_quad[comp][0],temp1);
              else
                apply_tensor_prod<0,false,false> 
                  (grad, this->gradients_quad[comp][0],temp1);
            }
          apply_tensor_prod<1,false,false> (val, temp1, this->values_dofs[comp]);
          if (integrate_grad == true)
            {
              // grad y
              apply_tensor_prod<0,false,false> (grad, this->gradients_quad[comp][1], temp1);
              apply_tensor_prod<1,false,true> (val, temp1, this->values_dofs[comp]);
            }
        }
      else if (dim == 1)
        {
          if (integrate_grad == true)
            apply_tensor_prod<0,false,false> (grad, this->gradients_quad[comp][0],
                                              this->values_dofs[comp]);
          if (integrate_val == true)
            {
              if (integrate_grad == true)
                apply_tensor_prod<0,false,true> (val, this->values_quad[comp],
                                                 this->values_dofs[comp]);
              else
                apply_tensor_prod<0,false,false> (val, this->values_quad[comp],
                                                  this->values_dofs[comp]);
            }
        }
    }

#ifdef DEBUG
  this->dof_values_initialized = true;
#endif
}



template <int dim, int n_dofs_1d,  int n_q_points_1d, int n_components,
          typename Number>
inline
Point<dim,VectorizedArray<Number> >
FEEvaluationGeneral<dim,n_dofs_1d,n_q_points_1d,n_components,Number>::
quadrature_point (const unsigned int q) const
{
  Assert (this->mapping_info.quadrature_points_initialized == true,
          ExcNotInitialized());
  AssertIndexRange (q, n_q_points);

                                // Cartesian mesh: not all quadrature points
                                // are stored, only the diagonal. Hence, need
                                // to find the tensor product index and
                                // retrieve the value from that
  if (this->cell_type == 0)
    {
      Point<dim,vector_t> point (false);
      switch (dim)
        {
        case 1:
          return this->quadrature_points[q];
        case 2:
          point[0] = this->quadrature_points[q%n_q_points_1d][0];
          point[1] = this->quadrature_points[q/n_q_points_1d][1];
          return point;
        case 3:
          point[0] = this->quadrature_points[q%n_q_points_1d][0];
          point[1] = this->quadrature_points[(q/n_q_points_1d)%n_q_points_1d][1];
          point[2] = this->quadrature_points[q/(n_q_points_1d*n_q_points_1d)][2];
          return point;
        default:
          Assert (false, ExcNotImplemented());
          return point;
        }
    }
                                // all other cases: just return the respective
                                // data as it is fully stored
  else
    return this->quadrature_points[q];
}


                                // General tensor product application for up
                                // to three spatial dimensions. Does not
                                // assume any symmetry in the shape values
                                // field
template <int dim, int n_dofs_1d,  int n_q_points_1d, int n_components,
          typename Number>
template <int direction, bool dof_to_quad, bool add>
inline
void
FEEvaluationGeneral<dim,n_dofs_1d,n_q_points_1d,n_components,Number>::
apply_tensor_prod (const VectorizedArray<Number>*shape_data,
                   const VectorizedArray<Number> input [],
                   VectorizedArray<Number>       output [])
{
  AssertIndexRange (direction, dim);
  const int mm     = dof_to_quad ? n_dofs_1d : n_q_points_1d,
            nn     = dof_to_quad ? n_q_points_1d : n_dofs_1d;

  const int n_blocks1 = (dim > 1 ? (direction > 0 ? nn : mm) : 1);
  const int n_blocks2 = (dim > 2 ? (direction > 1 ? nn : mm) : 1);
  const int stride    = ((direction > 0 ? nn : 1 ) *
                         (direction > 1 ? nn : 1));

  const vector_t * in = &input[0];
  vector_t * out = &output[0];
  for (int i2=0; i2<n_blocks2; ++i2)
  {
    for (int i1=0; i1<n_blocks1; ++i1)
    {
      for (int col=0; col<nn; ++col)
        {
          vector_t val0;
          if (dof_to_quad == true)
            val0 = shape_data[col];
          else
            val0 = shape_data[col*n_q_points_1d];
          vector_t res0 = val0 * in[0];
          for (int ind=1; ind<mm; ++ind)
            {
              if (dof_to_quad == true)
                val0 = shape_data[ind*n_q_points_1d+col];
              else
                val0 = shape_data[col*n_q_points_1d+ind];
              res0 += val0 * in[stride*ind];
            }
          if (add == false)
            out[stride*col]         = res0;
          else
            out[stride*col]        += res0;
        }

                                // increment: in regular case, just go to the
                                // next point in x-direction. If we are at the
                                // end of one chunk in x-dir, need to jump
                                // over to the next layer in z-direction
      switch (direction)
        {
        case 0:
          in += mm;
          out += nn;
          break;
        case 1:
        case 2:
          ++in;
          ++out;
          break;
        default:
          Assert (false, ExcNotImplemented());
        }
    }
    if (direction == 1)
      {
        in += nn*(mm-1);
        out += nn*(nn-1);
      }
  }
}


/*----------------------- FEEvaluation -------------------------------*/


template <int dim, int n_dofs_1d,  int n_q_points_1d, int n_components,
          typename Number>
inline
FEEvaluation<dim,n_dofs_1d,n_q_points_1d,n_components,Number>::
FEEvaluation (const MatrixFree<dim,Number> &data_in,
              const unsigned int fe_no,
              const unsigned int quad_no)
  :
  BaseClass (data_in, fe_no, quad_no)
{
                                // check whether element is appropriate
#ifdef DEBUG
  const double zero_tol =
    types_are_equal<Number,double>::value==true?1e-8:1e-7;
  std::string error_message = "FEEvaluation not appropriate.\n";
  error_message += "  It assumes symmetry of quadrature points w.r.t. 0.5 \n";
  error_message += " and the basis functions starting from left and right.\n";
  error_message += "Try FEEvaluationGeneral<...> instead!";

                                // symmetry for values
  for (unsigned int i=0; i<(n_dofs_1d+1)/2; ++i)
    for (unsigned int j=0; j<n_q_points_1d; ++j)
      Assert (std::fabs(this->data.shape_values[i*n_q_points_1d+j][0] -
                        this->data.shape_values[(n_dofs_1d-i)*n_q_points_1d
                                                -j-1][0]) < zero_tol,
              ExcMessage(error_message));

                                // shape values should be zero at for all
                                // basis functions except for one where they
                                // are one in the middle
  if (n_q_points_1d%2 == 1 && n_dofs_1d%2 == 1)
    {
      for (int i=0; i<static_cast<int>(n_dofs_1d/2); ++i)
        Assert (std::fabs(this->data.shape_values[i*n_q_points_1d+
                                                  n_q_points_1d/2][0]) < zero_tol,
                ExcMessage(error_message));
      Assert (std::fabs(this->data.shape_values[(n_dofs_1d/2)*n_q_points_1d+
                                                n_q_points_1d/2][0]-1.)< zero_tol,
              ExcMessage(error_message));
    }

                                // skew-symmetry for gradient, zero of middle
                                // basis function in middle quadrature point
  for (unsigned int i=0; i<(n_dofs_1d+1)/2; ++i)
    for (unsigned int j=0; j<n_q_points_1d; ++j)
      Assert (std::fabs(this->data.shape_gradients[i*n_q_points_1d+j][0] +
                        this->data.shape_gradients[(n_dofs_1d-i)*n_q_points_1d-
                                                  j-1][0]) < zero_tol,
              ExcMessage(error_message));
  if (n_dofs_1d%2 == 1 && n_q_points_1d%2 == 1)
    Assert (std::fabs(this->data.shape_gradients[(n_dofs_1d/2)*n_q_points_1d+
                                                (n_q_points_1d/2)][0]) < zero_tol,
            ExcMessage(error_message));


                                // symmetry for Laplacian
  for (unsigned int i=0; i<(n_dofs_1d+1)/2; ++i)
    for (unsigned int j=0; j<n_q_points_1d; ++j)
      Assert (std::fabs(this->data.shape_hessians[i*n_q_points_1d+j][0] -
                        this->data.shape_hessians[(n_dofs_1d-i)*n_q_points_1d-
                                                   j-1][0]) < zero_tol,
              ExcMessage(error_message));
#endif
}



template <int dim, int n_dofs_1d,  int n_q_points_1d, int n_components,
          typename Number>
inline
void
FEEvaluation<dim,n_dofs_1d,n_q_points_1d,n_components,Number>::
evaluate (bool evaluate_val, bool evaluate_grad, bool evaluate_lapl)
{
  Assert (this->cell != numbers::invalid_unsigned_int,
          ExcNotInitialized());
  Assert (this->dof_values_initialized == true,
          internal::ExcAccessToUninitializedField());

  for(unsigned int comp=0;comp<n_components;comp++)
    {
      vector_t temp1[n_dofs_1d > n_q_points_1d ? dofs_per_cell : n_q_points];
      vector_t temp2[n_dofs_1d > n_q_points_1d ? dofs_per_cell : n_q_points];

      if (dim == 3)
        {
          if (evaluate_grad == true)
            {
              // grad x
              apply_gradients<0,true,false> (this->values_dofs[comp], temp1);
              apply_values<1,true,false> (temp1, temp2);
              apply_values<2,true,false> (temp2, this->gradients_quad[comp][0]);
            }

          if (evaluate_lapl == true)
            {
              // grad xz
              if (evaluate_grad == false)
                {
                  apply_gradients<0,true,false> (this->values_dofs[comp], temp1);
                  apply_values<1,true,false> (temp1, temp2);
                }
              apply_gradients<2,true,false>(temp2, this->hessians_quad[comp][4]);

                  // grad xy
              apply_gradients<1,true,false>(temp1, temp2);
              apply_values<2,true,false> (temp2, this->hessians_quad[comp][3]);

              // grad xx
              apply_hessians<0,true,false> (this->values_dofs[comp], temp1);
              apply_values<1,true,false> (temp1, temp2);
              apply_values<2,true,false> (temp2, this->hessians_quad[comp][0]);
            }

          // grad y
          apply_values<0,true,false> (this->values_dofs[comp], temp1);
          if (evaluate_grad == true)
            {
              apply_gradients<1,true,false> (temp1, temp2);
              apply_values<2,true,false> (temp2, this->gradients_quad[comp][1]);
            }

          if (evaluate_lapl == true)
            {
              // grad yz
              if (evaluate_grad == false)
                apply_gradients<1,true,false> (temp1, temp2);
              apply_gradients<2,true,false> (temp2, this->hessians_quad[comp][5]);

              // grad yy
              apply_hessians<1,true,false> (temp1, temp2);
              apply_values<2,true,false> (temp2, this->hessians_quad[comp][1]);
            }

          // grad z: can use the values applied in x direction stored in temp1
          apply_values<1,true,false> (temp1, temp2);
          if (evaluate_grad == true)
            apply_gradients<2,true,false> (temp2, this->gradients_quad[comp][2]);

          // grad zz: can use the values applied in x and y direction stored in temp2
          if (evaluate_lapl == true)
            apply_hessians<2,true,false> (temp2, this->hessians_quad[comp][2]);

          // val: can use the values applied in x & y direction stored in temp2
          if (evaluate_val == true)
            apply_values<2,true,false> (temp2, this->values_quad[comp]);
        }
      else if (dim == 2)
        {
          // grad x
          if (evaluate_grad == true)
            {
              apply_gradients<0,true,false> (this->values_dofs[comp], temp1);
              apply_values<1,true,false> (temp1, this->gradients_quad[comp][0]);
            }
          if (evaluate_lapl == true)
            {
              // grad xy
              if (evaluate_grad == false)
                apply_gradients<0,true,false> (this->values_dofs[comp], temp1);
              apply_gradients<1,true,false> (temp1, this->hessians_quad[comp][2]);

              // grad xx
              apply_hessians<0,true,false> (this->values_dofs[comp], temp1);
              apply_values<1,true,false> (temp1, this->hessians_quad[comp][0]);
            }

          // grad y
          apply_values<0,true,false> (this->values_dofs[comp], temp1);
          if (evaluate_grad == true)
            apply_gradients<1,true,false> (temp1, this->gradients_quad[comp][1]);

          // grad yy
          if (evaluate_lapl == true)
            apply_hessians<1,true,false> (temp1, this->hessians_quad[comp][1]);

          // val: can use values applied in x
          if (evaluate_val == true)
            apply_values<1,true,false> (temp1, this->values_quad[comp]);
        }
      else if (dim == 1)
        {
          if (evaluate_val == true)
            apply_values<0,true,false> (this->values_dofs[comp],
                                        this->values_quad[comp]);
          if (evaluate_grad == true)
            apply_gradients<0,true,false> (this->values_dofs[comp],
                                           this->gradients_quad[comp][0]);
          if (evaluate_lapl == true)
            apply_hessians<0,true,false> (this->values_dofs[comp],
                                            this->hessians_quad[comp][0]);
        }
    }

#ifdef DEBUG
  if (evaluate_val == true)
    this->values_quad_initialized = true;
  if (evaluate_grad == true)
    this->gradients_quad_initialized = true;
  if (evaluate_lapl == true)
    this->hessians_quad_initialized  = true;
#endif
}



template <int dim, int n_dofs_1d,  int n_q_points_1d, int n_components,
          typename Number>
inline
void
FEEvaluation<dim,n_dofs_1d,n_q_points_1d,n_components,Number>::
integrate (bool integrate_val,bool integrate_grad)
{
#ifdef DEBUG
  Assert (this->cell != numbers::invalid_unsigned_int,
          ExcNotInitialized());
  if (integrate_val == true)
    Assert (this->values_quad_submitted == true,
            internal::ExcAccessToUninitializedField());
  if (integrate_grad == true)
    Assert (this->gradients_quad_submitted == true,
            internal::ExcAccessToUninitializedField());
#endif

  for(unsigned int comp=0;comp<n_components;comp++)
    {
      vector_t temp1[n_dofs_1d > n_q_points_1d ? dofs_per_cell : n_q_points];
      vector_t temp2[n_dofs_1d > n_q_points_1d ? dofs_per_cell : n_q_points];

      if (dim == 3)
        {
          if (integrate_val == true)
            {
              // val
              apply_values<0,false,false> (this->values_quad[comp], temp1);
            }
          if (integrate_grad == true)
            {
              // grad x: can sum to temporary value in temp1
              if (integrate_val == true)
                apply_gradients<0,false,true> (this->gradients_quad[comp][0],
                                               temp1);
              else
                apply_gradients<0,false,false> (this->gradients_quad[comp][0],
                                                temp1);
            }
          apply_values<1,false,false> (temp1, temp2);
          if (integrate_grad == true)
            {
              // grad y: can sum to temporary x value in temp2
              apply_values<0,false,false> (this->gradients_quad[comp][1], temp1);
              apply_gradients<1,false,true> (temp1, temp2);
            }
          apply_values<2,false,false> (temp2, this->values_dofs[comp]);
          if (integrate_grad == true)
            {
              // grad z: can sum to temporary x and y value in output
              apply_values<0,false,false> (this->gradients_quad[comp][2], temp1);
              apply_values<1,false,false> (temp1, temp2);
              apply_gradients<2,false,true> (temp2, this->values_dofs[comp]);
            }
        }
      else if (dim == 2)
        {
          // val
          if (integrate_val == true)
            apply_values<0,false,false> (this->values_quad[comp], temp1);
          if (integrate_grad == true)
            {
              //grad x
              if (integrate_val == true)
                apply_gradients<0,false,true> (this->gradients_quad[comp][0],
                                               temp1);
              else
                apply_gradients<0,false,false> (this->gradients_quad[comp][0],
                                               temp1);
            }
          apply_values<1,false,false> (temp1, this->values_dofs[comp]);
          if (integrate_grad == true)
            {
              // grad y
              apply_values<0,false,false> (this->gradients_quad[comp][1], temp1);
              apply_gradients<1,false,true> (temp1, this->values_dofs[comp]);
            }
        }
      else if (dim == 1)
        {
          if (integrate_grad == true)
            apply_gradients<0,false,false> (this->gradients_quad[comp][0],
                                            this->values_dofs[comp]);
          if (integrate_val == true)
            {
              if (integrate_grad == true)
                apply_values<0,false,true> (this->values_quad[comp],
                                            this->values_dofs[comp]);
              else
                apply_values<0,false,false> (this->values_quad[comp],
                                             this->values_dofs[comp]);
            }
        }
    }
#ifdef DEBUG
  this->dof_values_initialized = true;
#endif
}



// ----------------- optimized implementation tensor product symmetric case

template <int dim, int n_dofs_1d,  int n_q_points_1d, int n_components,
          typename Number>
template <int direction, bool dof_to_quad, bool add>
inline
void
FEEvaluation<dim,n_dofs_1d,n_q_points_1d,n_components,Number>::
apply_values (const VectorizedArray<Number> input [],
              VectorizedArray<Number>       output [])
{
  AssertIndexRange (direction, dim);
  const int mm     = dof_to_quad ? n_dofs_1d : n_q_points_1d,
            nn     = dof_to_quad ? n_q_points_1d : n_dofs_1d;
  const int n_cols = nn / 2;
  const int mid    = mm / 2;

  const int n_blocks1 = (dim > 1 ? (direction > 0 ? nn : mm) : 1);
  const int n_blocks2 = (dim > 2 ? (direction > 1 ? nn : mm) : 1);
  const int stride    = ((direction > 0 ? nn : 1 ) *
                         (direction > 1 ? nn : 1));

  const vector_t * in = &input[0];
  vector_t * out = &output[0];
  for (int i2=0; i2<n_blocks2; ++i2)
  {
    for (int i1=0; i1<n_blocks1; ++i1)
    {
      for (int col=0; col<n_cols; ++col)
        {
          vector_t val0, val1, res0, res1;
          if (dof_to_quad == true)
            {
              val0 = this->data.shape_values[col];
              val1 = this->data.shape_values[nn-1-col];
            }
          else
            {
              val0 = this->data.shape_values[col*n_q_points_1d];
              val1 = this->data.shape_values[(col+1)*n_q_points_1d-1];
            }
          if (mid > 0)
            {
              res0 = val0 * in[0];
              res1 = val1 * in[0];
              res0 += val1 * in[stride*(mm-1)];
              res1 += val0 * in[stride*(mm-1)];
              for (int ind=1; ind<mid; ++ind)
                {
                  if (dof_to_quad == true)
                    {
                      val0 = this->data.shape_values[ind*n_q_points_1d+col];
                      val1 = this->data.shape_values[ind*n_q_points_1d+nn-1-col];
                    }
                  else
                    {
                      val0 = this->data.shape_values[col*n_q_points_1d+ind];
                      val1 = this->data.shape_values[(col+1)*n_q_points_1d-1-ind];
                    }
                  res0 += val0 * in[stride*ind];
                  res1 += val1 * in[stride*ind];
                  res0 += val1 * in[stride*(mm-1-ind)];
                  res1 += val0 * in[stride*(mm-1-ind)];
                }
            }
          else
            res0 = res1 = vector_t();
          if (dof_to_quad == true)
            {
              if (mm % 2 == 1)
                {
                  val0 = this->data.shape_values[mid*n_q_points_1d+col];
                  val1 = val0 * in[stride*mid];
                  res0 += val1;
                  res1 += val1;
                }
            }
          else
            {
              if (mm % 2 == 1 && nn % 2 == 0)
                {
                  val0 = this->data.shape_values[col*n_q_points_1d+mid];
                  val1 = val0 * in[stride*mid];
                  res0 += val1;
                  res1 += val1;
                }
            }
          if (add == false)
            {
              out[stride*col]         = res0;
              out[stride*(nn-1-col)]  = res1;
            }
          else
            {
              out[stride*col]        += res0;
              out[stride*(nn-1-col)] += res1;
            }
        }
      if ( dof_to_quad == true && nn%2==1 && mm%2==1 )
        {
          if (add==false)
            out[stride*n_cols]  = in[stride*mid];
          else
            out[stride*n_cols] += in[stride*mid];
        }
      else if (dof_to_quad == true && nn%2==1)
        {
          vector_t res0;
          vector_t val0  = this->data.shape_values[n_cols];
          if (mid > 0)
            {
              res0  = in[0] + in[stride*(mm-1)];
              res0 *= val0;
              for (int ind=1; ind<mid; ++ind)
                {
                  val0  = this->data.shape_values[ind*n_q_points_1d+n_cols];
                  vector_t val1  = in[stride*ind] + in[stride*(mm-1-ind)];
                  val1 *= val0;
                  res0 += val1;
                }
            }
          else
            res0 = vector_t();
          if (mm % 2 == 1)
            {
              val0  = this->data.shape_values[mid*n_q_points_1d+n_cols];
              res0 += val0 * in[stride*mid];
            }
          if (add == false)
            out[stride*n_cols]  = res0;
          else
            out[stride*n_cols] += res0;
        }
      else if (dof_to_quad == false && nn%2 == 1)
        {
          vector_t res0;
          if (mid > 0)
            {
              vector_t val0 = this->data.shape_values[n_cols*n_q_points_1d];
              res0 = in[0] + in[stride*(mm-1)];
              res0 *= val0;
              for (int ind=1; ind<mid; ++ind)
                {
                  val0  = this->data.shape_values[n_cols*n_q_points_1d+ind];
                  vector_t val1 = in[stride*ind] + in[stride*(mm-1-ind)];
                  val1 *= val0;
                  res0 += val1;
                }
              if (mm % 2)
                res0 += in[stride*mid];
            }
          else
            res0 = in[0];
          if (add == false)
            out[stride*n_cols]  = res0;
          else
            out[stride*n_cols] += res0;
        }

                                // increment: in regular case, just go to the
                                // next point in x-direction. If we are at the
                                // end of one chunk in x-dir, need to jump
                                // over to the next layer in z-direction
      switch (direction)
        {
        case 0:
          in += mm;
          out += nn;
          break;
        case 1:
        case 2:
          ++in;
          ++out;
          break;
        default:
          Assert (false, ExcNotImplemented());
        }
    }
    if (direction == 1)
      {
        in += nn*(mm-1);
        out += nn*(nn-1);
      }
  }
}



template <int dim, int n_dofs_1d,  int n_q_points_1d, int n_components,
          typename Number>
template <int direction, bool dof_to_quad, bool add>
inline
void
FEEvaluation<dim,n_dofs_1d,n_q_points_1d,n_components,Number>::
apply_gradients (const VectorizedArray<Number> input [],
                 VectorizedArray<Number>       output [])
{
  AssertIndexRange (direction, dim);
  const int mm     = dof_to_quad ? n_dofs_1d : n_q_points_1d,
            nn     = dof_to_quad ? n_q_points_1d : n_dofs_1d;
  const int n_cols = nn / 2;
  const int mid    = mm / 2;

  const int n_blocks1 = (dim > 1 ? (direction > 0 ? nn : mm) : 1);
  const int n_blocks2 = (dim > 2 ? (direction > 1 ? nn : mm) : 1);
  const int stride    = ((direction > 0 ? nn : 1 ) *
                         (direction > 1 ? nn : 1));

  const vector_t * in = &input[0];
  vector_t * out = &output[0];
  for (int i2=0; i2<n_blocks2; ++i2)
  {
    for (int i1=0; i1<n_blocks1; ++i1)
    {
      for (int col=0; col<n_cols; ++col)
        {
          vector_t val0, val1, res0, res1;
          if (dof_to_quad == true)
            {
              val0 = this->data.shape_gradients[col];
              val1 = this->data.shape_gradients[nn-1-col];
            }
          else
            {
              val0 = this->data.shape_gradients[col*n_q_points_1d];
              val1 = this->data.shape_gradients[(nn-col-1)*n_q_points_1d];
            }
          if (mid > 0)
            {
              res0 = val0 * in[0];
              res1 = val1 * in[0];
              res0 -= val1 * in[stride*(mm-1)];
              res1 -= val0 * in[stride*(mm-1)];
              for (int ind=1; ind<mid; ++ind)
                {
                  if (dof_to_quad == true)
                    {
                      val0 = this->data.shape_gradients[ind*n_q_points_1d+col];
                      val1 = this->data.shape_gradients[ind*n_q_points_1d+nn-1-col];
                    }
                  else
                    {
                      val0 = this->data.shape_gradients[col*n_q_points_1d+ind];
                      val1 = this->data.shape_gradients[(nn-col-1)*n_q_points_1d+ind];
                    }
                  res0 += val0 * in[stride*ind];
                  res1 += val1 * in[stride*ind];
                  res0 -= val1 * in[stride*(mm-1-ind)];
                  res1 -= val0 * in[stride*(mm-1-ind)];
                }
            }
          else
            res0 = res1 = vector_t();
          if (mm % 2 == 1)
            {
              if (dof_to_quad == true)
                val0 = this->data.shape_gradients[mid*n_q_points_1d+col];
              else
                val0 = this->data.shape_gradients[col*n_q_points_1d+mid];
              val1 = val0 * in[stride*mid];
              res0 += val1;
              res1 -= val1;
            }
          if (add == false)
            {
              out[stride*col]         = res0;
              out[stride*(nn-1-col)]  = res1;
            }
          else
            {
              out[stride*col]        += res0;
              out[stride*(nn-1-col)] += res1;
            }
        }
      if ( nn%2 == 1 )
        {
          vector_t val0, res0;
          if (dof_to_quad == true)
            val0 = this->data.shape_gradients[n_cols];
          else
            val0 = this->data.shape_gradients[n_cols*n_q_points_1d];
          res0  = in[0] - in[stride*(mm-1)];
          res0 *= val0;
          for (int ind=1; ind<mid; ++ind)
            {
              if (dof_to_quad == true)
                val0 = this->data.shape_gradients[ind*n_q_points_1d+n_cols];
              else
                val0 = this->data.shape_gradients[n_cols*n_q_points_1d+ind];
              vector_t val1  = in[stride*ind] - in[stride*(mm-1-ind)];
              val1 *= val0;
              res0 += val1;
            }
          if (add == false)
            out[stride*n_cols]  = res0;
          else
            out[stride*n_cols] += res0;
        }

                                // increment: in regular case, just go to the
                                // next point in x-direction. for y-part in 3D
                                // and if we are at the end of one chunk in
                                // x-dir, need to jump over to the next layer
                                // in z-direction
      switch (direction)
        {
        case 0:
          in += mm;
          out += nn;
          break;
        case 1:
        case 2:
          ++in;
          ++out;
          break;
        default:
          Assert (false, ExcNotImplemented());
        }
    }

    if (direction == 1)
      {
        in  += nn * (mm-1);
        out += nn * (nn-1);
      }
  }
}



                                // Laplacian operator application. Very
                                // similar to value application because the
                                // same symmetry relations hold. However, it
                                // is not possible to omit some values that
                                // are zero for the values
template <int dim, int n_dofs_1d,  int n_q_points_1d, int n_components,
          typename Number>
template <int direction, bool dof_to_quad, bool add>
inline
void
FEEvaluation<dim,n_dofs_1d,n_q_points_1d,n_components,Number>::
apply_hessians (const VectorizedArray<Number> input [],
                  VectorizedArray<Number>       output [])
{
  AssertIndexRange (direction, dim);
  const int mm     = dof_to_quad ? n_dofs_1d : n_q_points_1d,
            nn     = dof_to_quad ? n_q_points_1d : n_dofs_1d;
  const int n_cols = nn / 2;
  const int mid    = mm / 2;

  const int n_blocks1 = (dim > 1 ? (direction > 0 ? nn : mm) : 1);
  const int n_blocks2 = (dim > 2 ? (direction > 1 ? nn : mm) : 1);
  const int stride    = ((direction > 0 ? nn : 1 ) *
                         (direction > 1 ? nn : 1));

  const vector_t * in = &input[0];
  vector_t * out = &output[0];
  for (int i2=0; i2<n_blocks2; ++i2)
  {
    for (int i1=0; i1<n_blocks1; ++i1)
    {
      for (int col=0; col<n_cols; ++col)
        {
          vector_t val0, val1, res0, res1;
          if (dof_to_quad == true)
            {
              val0 = this->data.shape_hessians[col];
              val1 = this->data.shape_hessians[nn-1-col];
            }
          else
            {
              val0 = this->data.shape_hessians[col*n_q_points_1d];
              val1 = this->data.shape_hessians[(col+1)*n_q_points_1d-1];
            }
          if (mid > 0)
            {
              res0 = val0 * in[0];
              res1 = val1 * in[0];
              res0 += val1 * in[stride*(mm-1)];
              res1 += val0 * in[stride*(mm-1)];
              for (int ind=1; ind<mid; ++ind)
                {
                  if (dof_to_quad == true)
                    {
                      val0 = this->data.shape_hessians[ind*n_q_points_1d+col];
                      val1 = this->data.shape_hessians[ind*n_q_points_1d+nn-1-col];
                    }
                  else
                    {
                      val0 = this->data.shape_hessians[col*n_q_points_1d+ind];
                      val1 = this->data.shape_hessians[(col+1)*n_q_points_1d-1-ind];
                    }
                  res0 += val0 * in[stride*ind];
                  res1 += val1 * in[stride*ind];
                  res0 += val1 * in[stride*(mm-1-ind)];
                  res1 += val0 * in[stride*(mm-1-ind)];
                }
            }
          else
            res0 = res1 = vector_t();
          if (mm % 2 == 1)
            {
              if (dof_to_quad == true)
                val0 = this->data.shape_hessians[mid*n_q_points_1d+col];
              else
                val0 = this->data.shape_hessians[col*n_q_points_1d+mid];
              val1 = val0 * in[stride*mid];
              res0 += val1;
              res1 += val1;
            }
          if (add == false)
            {
              out[stride*col]         = res0;
              out[stride*(nn-1-col)]  = res1;
            }
          else
            {
              out[stride*col]        += res0;
              out[stride*(nn-1-col)] += res1;
            }
        }
      if ( nn%2 == 1 )
        {
          vector_t val0, res0;
          if (dof_to_quad == true)
            val0 = this->data.shape_hessians[n_cols];
          else
            val0 = this->data.shape_hessians[n_cols*n_q_points_1d];
          if (mid > 0)
            {
              res0  = in[0] + in[stride*(mm-1)];
              res0 *= val0;
              for (int ind=1; ind<mid; ++ind)
                {
                  if (dof_to_quad == true)
                    val0 = this->data.shape_hessians[ind*n_q_points_1d+n_cols];
                  else
                    val0 = this->data.shape_hessians[n_cols*n_q_points_1d+ind];
                  vector_t val1  = in[stride*ind] + in[stride*(mm-1-ind)];
                  val1 *= val0;
                  res0 += val1;
                }
            }
          else
            res0 = vector_t();
          if (mm % 2 == 1)
            {
              if (dof_to_quad == true)
                val0 = this->data.shape_hessians[mid*n_q_points_1d+n_cols];
              else
                val0 = this->data.shape_hessians[n_cols*n_q_points_1d+mid];
              res0 += val0 * in[stride*mid];
            }
          if (add == false)
            out[stride*n_cols]  = res0;
          else
            out[stride*n_cols] += res0;
        }

                                // increment: in regular case, just go to the
                                // next point in x-direction. If we are at the
                                // end of one chunk in x-dir, need to jump
                                // over to the next layer in z-direction
      switch (direction)
        {
        case 0:
          in += mm;
          out += nn;
          break;
        case 1:
        case 2:
          ++in;
          ++out;
          break;
        default:
          Assert (false, ExcNotImplemented());
        }
    }
    if (direction == 1)
      {
        in += nn*(mm-1);
        out += nn*(nn-1);
      }
  }
}


/*----------------------- FEEvaluationGL -------------------------------*/


template <int dim, int n_points_1d, int n_components, typename Number>
inline
FEEvaluationGL<dim,n_points_1d,n_components,Number>::
FEEvaluationGL (const MatrixFree<dim,Number> &data_in,
                  const unsigned int fe_no,
                  const unsigned int quad_no)
  :
  BaseClass (data_in, fe_no, quad_no)
{
#ifdef DEBUG
  std::string error_mess = "FEEvaluationGL not appropriate. It assumes:\n";
  error_mess += "   - identity operation for shape values\n";
  error_mess += "   - zero diagonal at interior points for gradients\n";
  error_mess += "Try FEEvaluation<...> instead!";

  const double zero_tol =
    types_are_equal<Number,double>::value==true?1e-9:1e-7;

  for (unsigned int i=0; i<n_points_1d; ++i)
    for (unsigned int j=0; j<n_points_1d; ++j)
      if (i!=j)
        {
          Assert (std::fabs(this->data.shape_values[i*n_points_1d+j][0])<zero_tol,
                  ExcMessage (error_mess.c_str()));
        }
      else
        {
          Assert (std::fabs(this->data.shape_values[i*n_points_1d+
                                                    j][0]-1.)<zero_tol,
                  ExcMessage (error_mess.c_str()));
        }
  for (unsigned int i=1; i<n_points_1d-1; ++i)
    Assert (std::fabs(this->data.shape_gradients[i*n_points_1d+i][0])<zero_tol,
            ExcMessage (error_mess.c_str()));
#endif
}



template <int dim, int n_points_1d, int n_components, typename Number>
inline
void
FEEvaluationGL<dim,n_points_1d,n_components,Number>::
evaluate (bool evaluate_val,bool evaluate_grad,bool evaluate_lapl)
{
  Assert (this->cell != numbers::invalid_unsigned_int,
          ExcNotInitialized());
  Assert (this->dof_values_initialized == true,
          internal::ExcAccessToUninitializedField());

  if (evaluate_val == true)
    {
      std::memcpy (&this->values_quad[0][0], &this->values_dofs[0][0],
                   dofs_per_cell * n_components *
                   sizeof (this->values_dofs[0][0]));
#ifdef DEBUG
      this->values_quad_initialized = true;
#endif
    }
  if (evaluate_grad == true)
    {
      for(unsigned int comp=0;comp<n_components;comp++)
        {
          if (dim == 3)
            {
              // grad x
              apply_gradients<0,true,false> (this->values_dofs[comp],
                                             this->gradients_quad[comp][0]);
              // grad y
              apply_gradients<1,true,false> (this->values_dofs[comp],
                                             this->gradients_quad[comp][1]);
              // grad y
              apply_gradients<2,true,false> (this->values_dofs[comp],
                                             this->gradients_quad[comp][2]);
            }
          else if (dim == 2)
            {
              // grad x
              apply_gradients<0,true,false> (this->values_dofs[comp],
                                             this->gradients_quad[comp][0]);
              // grad y
              apply_gradients<1,true,false> (this->values_dofs[comp],
                                             this->gradients_quad[comp][1]);
            }
          else if (dim == 1)
            apply_gradients<0,true,false> (this->values_dofs[comp],
                                           this->gradients_quad[comp][0]);
        }
#ifdef DEBUG
  this->gradients_quad_initialized = true;
#endif
    }
  if (evaluate_lapl == true)
    {
      for(unsigned int comp=0;comp<n_components;comp++)
        {
          if (dim == 3)
            {
              // grad x
              this->template apply_hessians<0,true,false> (this->values_dofs[comp],
                                                           this->hessians_quad[comp][0]);
              // grad y
              this->template apply_hessians<1,true,false> (this->values_dofs[comp],
                                                           this->hessians_quad[comp][1]);
              // grad y
              this->template apply_hessians<2,true,false> (this->values_dofs[comp],
                                                           this->hessians_quad[comp][2]);

              vector_t temp1[n_q_points];
              // grad xy
              apply_gradients<0,true,false> (this->values_dofs[comp], temp1);
              apply_gradients<1,true,false> (temp1, this->hessians_quad[comp][3]);
              // grad xz
              apply_gradients<2,true,false> (temp1, this->hessians_quad[comp][4]);
              // grad yz
              apply_gradients<1,true,false> (this->values_dofs[comp], temp1);
              apply_gradients<2,true,false> (temp1, this->hessians_quad[comp][5]);
            }
          else if (dim == 2)
            {
              // grad x
              this->template apply_hessians<0,true,false> (this->values_dofs[comp],
                                                           this->hessians_quad[comp][0]);
              // grad y
              this->template apply_hessians<1,true,false> (this->values_dofs[comp],
                                                           this->hessians_quad[comp][1]);
              vector_t temp1[n_q_points];
              // grad xy
              apply_gradients<0,true,false> (this->values_dofs[comp], temp1);
              apply_gradients<1,true,false> (temp1, this->hessians_quad[comp][2]);
            }
          else if (dim == 1)
            this->template apply_hessians<0,true,false> (this->values_dofs[comp],
                                                         this->hessians_quad[comp][0]);
        }
#ifdef DEBUG
      this->hessians_quad_initialized = true;
#endif
    }
}



template <int dim, int n_points_1d, int n_components, typename Number>
inline
void
FEEvaluationGL<dim,n_points_1d,n_components,Number>::
integrate (bool integrate_val, bool integrate_grad)
{
  Assert (this->cell != numbers::invalid_unsigned_int,
          ExcNotInitialized());
  if (integrate_val == true)
    Assert (this->values_quad_submitted == true,
            internal::ExcAccessToUninitializedField());
  if (integrate_grad == true)
    Assert (this->gradients_quad_submitted == true,
            internal::ExcAccessToUninitializedField());
  if (integrate_val == true)
    std::memcpy (&this->values_dofs[0][0], &this->values_quad[0][0],
                 dofs_per_cell * n_components *
                 sizeof (this->values_dofs[0][0]));
  if (integrate_grad == true)
    {
      for(unsigned int comp=0;comp<n_components;comp++)
        {
          if (dim == 3)
            {
              // grad x: If integrate_val == true we have to add to the previous output
              if (integrate_val == true)
                apply_gradients<0, false, true> (this->gradients_quad[comp][0],
                                                 this->values_dofs[comp]);
              else
                apply_gradients<0, false, false> (this->gradients_quad[comp][0],
                                                  this->values_dofs[comp]);

              // grad y: can sum to temporary x value in temp2
              apply_gradients<1, false, true> (this->gradients_quad[comp][1],
                                               this->values_dofs[comp]);

              // grad z: can sum to temporary x and y value in output
              apply_gradients<2, false, true> (this->gradients_quad[comp][2],
                                               this->values_dofs[comp]);
            }
          else if (dim == 2)
            {
              // grad x: If integrate_val == true we have to add to the previous output
              if (integrate_val == true)
                apply_gradients<0, false, true> (this->gradients_quad[comp][0],
                                                 this->values_dofs[comp]);
              else
                apply_gradients<0, false, false> (this->gradients_quad[comp][0],
                                                  this->values_dofs[comp]);

              // grad y: can sum to temporary x value in temp2
              apply_gradients<1, false, true> (this->gradients_quad[comp][1],
                                               this->values_dofs[comp]);
            }
          else if (dim == 1)
            {
              if (integrate_val == true)
                apply_gradients<0, false, true> (this->gradients_quad[comp][0],
                                                 this->values_dofs[comp]);
              else
                apply_gradients<0, false, false> (this->gradients_quad[comp][0],
                                                  this->values_dofs[comp]);

            }
        }
    }

#ifdef DEBUG
  this->dof_values_initialized = true;
#endif
}



template <int dim, int n_points_1d, int n_components, typename Number>
template <int direction, bool dof_to_quad, bool add>
inline
void
FEEvaluationGL<dim,n_points_1d,n_components,Number>::
apply_gradients (const VectorizedArray<Number> input [],
                 VectorizedArray<Number>       output [])
{
  AssertIndexRange (direction, dim);
  const int mm     = n_points_1d;
  const int nn     = n_points_1d;
  const int n_cols = nn / 2;
  const int mid    = mm / 2;

  const int n_blocks1 = (dim > 1 ? (direction > 0 ? nn : mm) : 1);
  const int n_blocks2 = (dim > 2 ? (direction > 1 ? nn : mm) : 1);
  const int stride    = ((direction > 0 ? nn : 1 ) *
                         (direction > 1 ? nn : 1));

  const vector_t * in = &input[0];
  vector_t * out = &output[0];
  for (int i2=0; i2<n_blocks2; ++i2)
  {
    for (int i1=0; i1<n_blocks1; ++i1)
    {
      for (int col=0; col<n_cols; ++col)
        {
          vector_t val0, val1, res0, res1;
          if (dof_to_quad == true)
            {
              val0 = this->data.shape_gradients[col];
              val1 = this->data.shape_gradients[nn-1-col];
            }
          else
            {
              val0 = this->data.shape_gradients[col*n_points_1d];
              val1 = this->data.shape_gradients[(nn-col-1)*n_points_1d];
            }
          if (mid > 0)
            {
              res0 = val0 * in[0];
              res1 = val1 * in[0];
              res0 -= val1 * in[stride*(mm-1)];
              res1 -= val0 * in[stride*(mm-1)];
              for (int ind=1; ind<mid; ++ind)
                {
                  if (dof_to_quad == true)
                    {
                      val0 = this->data.shape_gradients[ind*n_points_1d+col];
                      val1 = this->data.shape_gradients[ind*n_points_1d+nn-1-col];
                    }
                  else
                    {
                      val0 = this->data.shape_gradients[col*n_points_1d+ind];
                      val1 = this->data.shape_gradients[(nn-col-1)*n_points_1d+ind];
                    }

                                // at inner points, the gradient is zero for
                                // ind==col
                  if (ind == col)
                    {
                      res1 += val1 * in[stride*ind];
                      res0 -= val1 * in[stride*(mm-1-ind)];
                    }
                  else
                    {
                      res0 += val0 * in[stride*ind];
                      res1 += val1 * in[stride*ind];
                      res0 -= val1 * in[stride*(mm-1-ind)];
                      res1 -= val0 * in[stride*(mm-1-ind)];
                    }
                }
            }
          else
            res0 = res1 = vector_t();
          if (mm % 2 == 1)
            {
              if (dof_to_quad == true)
                val0 = this->data.shape_gradients[mid*n_points_1d+col];
              else
                val0 = this->data.shape_gradients[col*n_points_1d+mid];
              val1 = val0 * in[stride*mid];
              res0 += val1;
              res1 -= val1;
            }
          if (add == false)
            {
              out[stride*col]         = res0;
              out[stride*(nn-1-col)]  = res1;
            }
          else
            {
              out[stride*col]        += res0;
              out[stride*(nn-1-col)] += res1;
            }
        }
      if ( nn%2 == 1 )
        {
          vector_t val0, res0;
          if (dof_to_quad == true)
            val0 = this->data.shape_gradients[n_cols];
          else
            val0 = this->data.shape_gradients[n_cols*n_points_1d];
          if (mid > 0)
            {
              res0  = in[0] - in[stride*(mm-1)];
              res0 *= val0;
              for (int ind=1; ind<mid; ++ind)
                {
                  if (dof_to_quad == true)
                    val0 = this->data.shape_gradients[ind*n_points_1d+n_cols];
                  else
                    val0 = this->data.shape_gradients[n_cols*n_points_1d+ind];
                  vector_t val1  = in[stride*ind] - in[stride*(mm-1-ind)];
                  val1 *= val0;
                  res0 += val1;
                }
            }
          else
            res0 = vector_t();
          if (add == false)
            out[stride*n_cols]  = res0;
          else
            out[stride*n_cols] += res0;
        }

                                // increment: in regular case, just go to the
                                // next point in x-direction. for y-part in 3D
                                // and if we are at the end of one chunk in
                                // x-dir, need to jump over to the next layer
                                // in z-direction
      switch (direction)
        {
        case 0:
          in += mm;
          out += nn;
          break;
        case 1:
        case 2:
          ++in;
          ++out;
          break;
        default:
          Assert (false, ExcNotImplemented());
        }
    }

    if (direction == 1)
      {
        in  += nn * (mm-1);
        out += nn * (nn-1);
      }
  }

#endif  // ifndef DOXYGEN

}


DEAL_II_NAMESPACE_CLOSE

#endif
