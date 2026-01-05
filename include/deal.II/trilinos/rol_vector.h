// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_trilinos_rol_vector_h
#define dealii_trilinos_rol_vector_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_TRILINOS_WITH_ROL
#  include <deal.II/base/exceptions.h>
#  include <deal.II/base/index_set.h>
#  include <deal.II/base/mpi.h>
#  include <deal.II/base/types.h>

#  include <deal.II/lac/vector.h>

#  include <ROL_Vector.hpp>

#  include <limits>
#  include <tuple>
#  include <type_traits>


DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{
  /**
   * An adaptor that provides an interface to the
   * <a href="https://trilinos.org/docs/dev/packages/rol/doc/html/index.html">
   * Rapid Optimization Library</a> (ROL), a Trilinos package.
   *
   * This class provides the implementation of the ROL::Vector interface
   * for vectors of type <tt>VectorType</tt>. It supports vectors that satisfy
   * the following requirements:
   *
   * The <tt>VectorType</tt> should contain the following types.
   * ```
   * VectorType::size_type;  // The type for size of the vector.
   * VectorType::value_type; // The type for elements stored in the vector.
   * VectorType::real_type;  // The type for real-valued numbers.
   * ```
   *
   * However, ROL doesn't distinguish ROLVector::value_type from
   * ROLVector::real_type. This is due to ROL's assumption that the
   * ROLVector::value_type itself is a type for real-valued numbers.
   * Therefore, ROLVector supports vectors whose real_type is
   * convertible to value_type in the sense that
   * <code>std::is_convertible_v<real_type, value_type></code> yields
   * <code>true</code>.
   *
   * The <tt>VectorType</tt> should contain the following methods.
   * @code
   * // Reinitialize the current vector using a given vector's
   * // size (and the parallel distribution) without copying
   * // the elements.
   * void VectorType::reinit(const VectorType &, ...);
   *
   * // Copies the data of a given vector to the current.
   * // Resize the current vector if necessary (MPI safe).
   * VectorType& VectorType::operation=(const VectorType &);
   *
   * // Access the value of the ith component.
   * VectorType::value_type&
   * VectorType::operator[](const VectorType::size_type i);
   *
   * // Return whether the vector currently is in a state where
   * // ghost values can be read or not.
   * bool VectorType::has_ghost_elements();
   *
   * // Fills the data field for ghost indices with the values stored
   * // in the respective positions of the owning processor.
   * void VectorType::update_ghost_values();
   *
   * // Compress the vector i.e., flush the buffers of the
   * // vector object if it has any.
   * void VectorType::compress(VectorOperation::insert);
   * @endcode
   *
   * @note The current implementation in ROL doesn't support vector sizes above
   *   the largest value of int type. Some of the vectors in deal.II (see
   *   @ref Vector) may not satisfy the above requirements.
   *
   *
   * <h3>Optimization and finite elements</h3>
   *
   * Numerical optimization can be a useful tool in the context of finite
   * element methods. However, not all degrees of freedom of a finite element
   * discretization should really be part of the optimization space -- typically
   * those that are constrained through boundary conditions and hanging nodes
   * should be excluded from the optimization process.
   *
   * This can be achieved in two different ways: either through constrained
   * optimization; or by reducing the optimization space which allows
   * unconstrained optimization. For the latter, one way to do it is masking
   * those vector elements that correspond to the optimization space. All
   * operations necessary for the optimization process would only \e see these
   * masked entries, which are effectively only a subspace of the finite element
   * discretization.
   *
   * This wrapper class provides a means for masking using an IndexSet upon
   * construction.
   *
   *
   * <h3>Ghost elements on distributed vectors</h3>
   *
   * This wrapper class also manages all distributed vector classes of deal.II.
   *
   * ROL will manipulate the wrapped vector during the optimization progress.
   * For vector classes of the PETScWrappers and TrilinosWrappers namespaces,
   * this manipulation is only allowed on vectors \e without ghost elements.
   * For the LinearAlgebra::distributed::Vector class, manipulation with ghost
   * elements is permitted, although it requires a potentially slow ghost
   * exchange. For more information on ghost vectors, see also the corresponding
   * @ref GlossGhostedVector "glossary entry".
   *
   * To enforce the same behavior of this wrapper on all distributed vectors
   * classes, we only allow manipulating operations for vectors \e without ghost
   * elements.
   */
  template <typename VectorType>
  class ROLVector : public ROL::Vector<typename VectorType::value_type>
  {
    /**
     * An alias for size type of <tt>VectorType</tt>.
     */
    using size_type = typename VectorType::size_type;

    /**
     * An alias for element type stored in the <tt>VectorType</tt>.
     */
    using value_type = typename VectorType::value_type;

    /**
     * An alias for real-valued numbers.
     */
    using real_type = typename VectorType::real_type;

    static_assert(std::is_convertible_v<real_type, value_type>,
                  "The real_type of the current VectorType is not "
                  "convertible to the value_type.");

  private:
    /**
     * ROL pointer to the underlying vector of type <tt>VectorType</tt>.
     */
    ROL::Ptr<VectorType> vector_ptr;

    /**
     * IndexSet that represents the locally owned optimization space.
     */
    IndexSet optimization_space;

    /**
     * Global dimension of the optimization space.
     */
    dealii::types::global_dof_index global_dimension;

    /**
     * Prefix sum of the number of elements in the optimization space of lower
     * rank MPI processes.
     */
    dealii::types::global_dof_index prefix_sum;

  public:
    /**
     * Constructor.
     *
     * The optimization space covers all locally owned degrees of freedom.
     */
    ROLVector(const ROL::Ptr<VectorType> &vector_ptr);

    /**
     * Constructor.
     *
     * The IndexSet @p optimization_space masks entries in the wrapped vector
     * on which the optimization process will be applied only. It needs to
     * contain a subset of locally owned indices of the wrapped vector.
     */
    ROLVector(const ROL::Ptr<VectorType> &vector_ptr,
              const IndexSet             &optimization_space);

    /**
     * Return the ROL pointer to the wrapped vector.
     */
    ROL::Ptr<VectorType>
    getVector();

    /**
     * Return the ROL pointer to the wrapped vector as read-only.
     */
    ROL::Ptr<const VectorType>
    getVector() const;

    /**
     * Return the global dimension of the optimization space.
     */
    int
    dimension() const;

    /**
     * Set the wrapped vector to a given @p rol_vector by overwriting its
     * contents.
     *
     * If the current wrapped vector has ghost elements,
     * then <code> VectorType::operator=(const VectorType&) </code> should still
     * be allowed on it.
     *
     * @note @p rol_vector has to be of type ROLVector.
     */
    void
    set(const ROL::Vector<value_type> &rol_vector);

    /**
     * Add @p rol_vector to the wrapped vector.
     *
     * The operation will only be applied to masked elements. Both vectors need
     * to correspond to the same optimization space.
     *
     * @note @p rol_vector has to be of type ROLVector.
     */
    void
    plus(const ROL::Vector<value_type> &rol_vector);

    /**
     * Scale the wrapped vector by @p alpha and add @p rol_vector to it.
     *
     * The operation will only be applied to masked elements. Both vectors need
     * to correspond to the same optimization space.
     *
     * @note @p rol_vector has to be of type ROLVector.
     */
    void
    axpy(const value_type alpha, const ROL::Vector<value_type> &rol_vector);

    /**
     * Scale the wrapped vector by @p alpha.
     *
     * The operation will only be applied to masked elements.
     */
    void
    scale(const value_type alpha);

    /**
     * Return the dot product with a given @p rol_vector.
     *
     * The operation will only be applied to masked elements. Both vectors need
     * to correspond to the same optimization space.
     *
     * @note @p rol_vector has to be of type ROLVector.
     */
    value_type
    dot(const ROL::Vector<value_type> &rol_vector) const;

    /**
     * Return the $L^{2}$ norm of the wrapped vector.
     *
     * The operation will only be applied to masked elements.
     *
     * The returned type is of ROLVector::value_type so as to maintain
     * consistency with ROL::Vector<ROLVector::value_type> and
     * more importantly to not to create an overloaded version namely,
     * <code> ROLVector::real_type norm() const; </code>
     * if real_type and value_type are not of the same type.
     */
    value_type
    norm() const;

    /**
     * Return a clone of the wrapped vector as a pointer to the parent class.
     */
    ROL::Ptr<ROL::Vector<value_type>>
    clone() const;

    /**
     * Create and return a ROL pointer to the basis vector corresponding to the
     * @p i ${}^{th}$ element of the global optimization space.
     */
    ROL::Ptr<ROL::Vector<value_type>>
    basis(const int i) const;

    /**
     * Apply unary function @p f to all the elements of the wrapped vector.
     *
     * The operation will only be applied to masked elements.
     */
    void
    applyUnary(const ROL::Elementwise::UnaryFunction<value_type> &f);

    /**
     * Apply binary function @p f along with ROL::Vector @p rol_vector to all
     * the elements of the wrapped vector.
     *
     * The operation will only be applied to masked elements. Both vectors need
     * to correspond to the same optimization space.
     *
     * @note @p rol_vector has to be of type ROLVector.
     */
    void
    applyBinary(const ROL::Elementwise::BinaryFunction<value_type> &f,
                const ROL::Vector<value_type>                      &rol_vector);

    /**
     * Return the accumulated value on applying reduction operation @p r on
     * all the elements of the wrapped vector.
     *
     * The operation will only be applied to masked elements.
     */
    value_type
    reduce(const ROL::Elementwise::ReductionOp<value_type> &r) const;

    /**
     * Print the wrapped vector to the output stream @p outStream.
     *
     * This function will print the entire vector, and not just the masked
     * elements.
     */
    void
    print(std::ostream &outStream) const;
  };


  /*------------------------------member definitions--------------------------*/
#  ifndef DOXYGEN


  template <typename VectorType>
  ROLVector<VectorType>::ROLVector(const ROL::Ptr<VectorType> &vector_ptr)
    : ROLVector(vector_ptr, vector_ptr->locally_owned_elements())
  {}



  template <typename VectorType>
  ROLVector<VectorType>::ROLVector(const ROL::Ptr<VectorType> &vector_ptr,
                                   const IndexSet &optimization_space)
    : vector_ptr(vector_ptr)
    , optimization_space(optimization_space)
  {
    Assert(
      optimization_space.is_subset_of(vector_ptr->locally_owned_elements()),
      ExcMessage(
        "Provided IndexSet needs to be a subset of locally owned indices."));

    std::tie(prefix_sum, global_dimension) =
      Utilities::MPI::partial_and_total_sum(optimization_space.n_elements(),
                                            vector_ptr->get_mpi_communicator());

    // The return type of ROL::Vector::dimension() has to be int,
    // so we check this requirement of ROL here.
    Assert(global_dimension < static_cast<dealii::types::global_dof_index>(
                                std::numeric_limits<int>::max()),
           ExcMessage("The number of elements to optimize is greater than the "
                      "largest value of type int."));
  }



  template <typename VectorType>
  ROL::Ptr<VectorType>
  ROLVector<VectorType>::getVector()
  {
    return vector_ptr;
  }



  template <typename VectorType>
  ROL::Ptr<const VectorType>
  ROLVector<VectorType>::getVector() const
  {
    return vector_ptr;
  }



  template <typename VectorType>
  void
  ROLVector<VectorType>::set(const ROL::Vector<value_type> &other_)
  {
    const ROLVector &other = dynamic_cast<const ROLVector &>(other_);

    // Perform a deep copy of the vector pointed to by vector_ptr.
    (*vector_ptr) = *(other.getVector());

    optimization_space = other.optimization_space;
    global_dimension   = other.global_dimension;
    prefix_sum         = other.prefix_sum;
  }



  template <typename VectorType>
  void
  ROLVector<VectorType>::plus(const ROL::Vector<value_type> &other_)
  {
    const ROLVector &other = dynamic_cast<const ROLVector &>(other_);

    Assert(vector_ptr->has_ghost_elements() == false, ExcGhostsPresent());
    Assert(optimization_space.size() == vector_ptr->size(),
           ExcMessage("Optimization space is out-of-sync. "
                      "Please create a new wrapper."));
    Assert(optimization_space == other.optimization_space,
           ExcMessage("Optimization spaces of vectors do not match."));

    for (const auto i : optimization_space)
      (*vector_ptr)[i] += (*other.getVector())[i];

    vector_ptr->compress(VectorOperation::add);
  }



  template <typename VectorType>
  void
  ROLVector<VectorType>::axpy(const value_type               alpha,
                              const ROL::Vector<value_type> &other_)
  {
    const ROLVector &other = dynamic_cast<const ROLVector &>(other_);

    Assert(vector_ptr->has_ghost_elements() == false, ExcGhostsPresent());
    Assert(optimization_space.size() == vector_ptr->size(),
           ExcMessage("Optimization space is out-of-sync. "
                      "Please create a new wrapper."));
    Assert(optimization_space == other.optimization_space,
           ExcMessage("Optimization spaces of vectors do not match."));

    for (const auto i : optimization_space)
      (*vector_ptr)[i] += alpha * (*other.getVector())[i];

    vector_ptr->compress(VectorOperation::add);
  }



  template <typename VectorType>
  int
  ROLVector<VectorType>::dimension() const
  {
    Assert(optimization_space.size() == vector_ptr->size(),
           ExcMessage("Optimization space is out-of-sync. "
                      "Please create a new wrapper."));

    return static_cast<int>(global_dimension);
  }



  template <typename VectorType>
  void
  ROLVector<VectorType>::scale(const value_type alpha)
  {
    Assert(vector_ptr->has_ghost_elements() == false, ExcGhostsPresent());
    Assert(optimization_space.size() == vector_ptr->size(),
           ExcMessage("Optimization space is out-of-sync. "
                      "Please create a new wrapper."));

    for (const auto i : optimization_space)
      (*vector_ptr)[i] *= alpha;

    vector_ptr->compress(VectorOperation::insert);
  }



  template <typename VectorType>
  typename VectorType::value_type
  ROLVector<VectorType>::dot(const ROL::Vector<value_type> &other_) const
  {
    const ROLVector &other = dynamic_cast<const ROLVector &>(other_);

    Assert(optimization_space.size() == vector_ptr->size(),
           ExcMessage("Optimization space is out-of-sync. "
                      "Please create a new wrapper."));
    Assert(optimization_space == other.optimization_space,
           ExcMessage("Optimization spaces of vectors do not match."));

    value_type dot(0);
    for (const auto i : optimization_space)
      dot += (*vector_ptr)[i] * (*other.getVector())[i];

    return Utilities::MPI::sum<value_type>(dot,
                                           vector_ptr->get_mpi_communicator());
  }



  template <typename VectorType>
  typename VectorType::value_type
  ROLVector<VectorType>::norm() const
  {
    return std::sqrt(this->dot(*this));
  }



  template <typename VectorType>
  ROL::Ptr<ROL::Vector<typename VectorType::value_type>>
  ROLVector<VectorType>::clone() const
  {
    ROL::Ptr<VectorType> vec_ptr = ROL::makePtr<VectorType>(*vector_ptr);

    return ROL::makePtr<ROLVector>(vec_ptr, optimization_space);
    // TODO: also somehow pass the partial sums etc?
  }



  template <typename VectorType>
  ROL::Ptr<ROL::Vector<typename VectorType::value_type>>
  ROLVector<VectorType>::basis(const int i) const
  {
    Assert(optimization_space.size() == vector_ptr->size(),
           ExcMessage("Optimization space is out-of-sync. "
                      "Please create a new wrapper."));
    AssertIndexRange(i, global_dimension);

    // clone the internal vector as basis
    ROL::Ptr<VectorType> basis_ptr = ROL::makePtr<VectorType>();
    basis_ptr->reinit(*vector_ptr, false);

    // check whether the i-th basis belongs to us
    if ((i >= prefix_sum) && (i < prefix_sum + optimization_space.n_elements()))
      {
        const IndexSet::size_type local_index = i - prefix_sum;

        // global dof index corresponding to i-th basis in optimization
        const dealii::types::global_dof_index global_dof_index =
          optimization_space.nth_index_in_set(local_index);

        (*basis_ptr)[global_dof_index] = 1.;
      }

    basis_ptr->compress(VectorOperation::insert);

    return ROL::makePtr<ROLVector>(basis_ptr, optimization_space);
    // TODO: also somehow pass the partial sums etc?
  }



  template <typename VectorType>
  void
  ROLVector<VectorType>::applyUnary(
    const ROL::Elementwise::UnaryFunction<value_type> &f)
  {
    Assert(vector_ptr->has_ghost_elements() == false, ExcGhostsPresent());
    Assert(optimization_space.size() == vector_ptr->size(),
           ExcMessage("Optimization space is out-of-sync. "
                      "Please create a new wrapper."));

    for (const auto i : optimization_space)
      (*vector_ptr)[i] = f.apply((*vector_ptr)[i]);

    vector_ptr->compress(VectorOperation::insert);
  }



  template <typename VectorType>
  void
  ROLVector<VectorType>::applyBinary(
    const ROL::Elementwise::BinaryFunction<value_type> &f,
    const ROL::Vector<value_type>                      &other_)
  {
    const ROLVector &other = dynamic_cast<const ROLVector &>(other_);

    Assert(vector_ptr->has_ghost_elements() == false, ExcGhostsPresent());
    Assert(optimization_space.size() == vector_ptr->size(),
           ExcMessage("Optimization space is out-of-sync. "
                      "Please create a new wrapper."));
    Assert(optimization_space == other.optimization_space,
           ExcMessage("Optimization spaces of vectors do not match."));

    for (const auto i : optimization_space)
      (*vector_ptr)[i] = f.apply((*vector_ptr)[i], (*other.getVector())[i]);

    vector_ptr->compress(VectorOperation::insert);
  }



  template <typename VectorType>
  typename VectorType::value_type
  ROLVector<VectorType>::reduce(
    const ROL::Elementwise::ReductionOp<value_type> &r) const
  {
    Assert(optimization_space.size() == vector_ptr->size(),
           ExcMessage("Optimization space is out-of-sync. "
                      "Please create a new wrapper."));

    value_type result = r.initialValue();

    // local reduction
    for (const auto i : optimization_space)
      r.reduce((*vector_ptr)[i], result);

    // global reduction
    const auto combiner = [&r](const value_type a,
                               const value_type b) -> value_type {
      auto out = b;
      r.reduce(a, out);
      return out;
    };

    return Utilities::MPI::all_reduce<value_type>(
      result, vector_ptr->get_mpi_communicator(), combiner);
  }



  template <typename VectorType>
  void
  ROLVector<VectorType>::print(std::ostream &outStream) const
  {
    vector_ptr->print(outStream);
  }


#  endif // DOXYGEN


} // namespace TrilinosWrappers


DEAL_II_NAMESPACE_CLOSE


#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_TRILINOS_WITH_ROL

#endif // dealii_trilinos_rol_vector_h
