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

#ifndef dealii_trilinos_rol_adaptor_h
#define dealii_trilinos_rol_adaptor_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_TRILINOS_WITH_ROL
#  include <deal.II/base/exceptions.h>
#  include <deal.II/base/index_set.h>
#  include <deal.II/base/mpi.h>
#  include <deal.II/base/types.h>

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
   * However, ROL doesn't distinguish ROLAdaptor::value_type from
   * ROLAdaptor::real_type. This is due to ROL's assumption that the
   * ROLAdaptor::value_type itself is a type for real-valued numbers.
   * Therefore, ROLAdaptor supports vectors whose real_type is
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
   * VectorType& VectorType::operator=(const VectorType &);
   *
   * // Return the global size of the current vector.
   * VectorType::size_type VectorType::size();
   *
   * // Access the value of the ith component.
   * VectorType::value_type&
   * VectorType::operator[](const VectorType::size_type i);
   *
   * // Compress the vector i.e., flush the buffers of the
   * // vector object if it has any.
   * void VectorType::compress(VectorOperation);
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
   * unconstrained optimization. For the latter, one way to do it is explicitly
   * defining those vector elements that correspond to the optimization space.
   * All operations necessary for the optimization process would only \e see
   * these specified entries, which are effectively only a subspace of the
   * finite element discretization.
   *
   * This wrapper class provides a means for defining the optimization space
   * using an IndexSet upon construction.
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
  class ROLAdaptor : public ROL::Vector<typename VectorType::value_type>
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
    dealii::types::global_dof_index global_opt_dimension;

    /**
     * Prefix sum of the number of elements in the optimization space of lower
     * rank MPI processes.
     *
     * Denotes the starting index in the optimization space, at which the
     * locally owned elements of the wrapped vector begin.
     */
    dealii::types::global_dof_index local_opt_start_index;

  public:
    /**
     * Constructor.
     *
     * This constructor does not take an index set as argument (for this, see
     * the other constructor) and as a consequence considers all variables
     * represented by the vector as available to optimization.
     */
    ROLAdaptor(const ROL::Ptr<VectorType> &vector_ptr);

    /**
     * Constructor.
     *
     * The IndexSet @p optimization_space specifies entries in the wrapped
     * vector to which the optimization process will be applied. It needs to be
     * equal to or a subset of the locally owned indices of the wrapped vector.
     *
     * Passing in the complete index set (partitioned among processors) is
     * equivalent to the previous constructor.
     */
    ROLAdaptor(const ROL::Ptr<VectorType> &vector_ptr,
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
    dimension() const override;

    /**
     * Set the wrapped vector to a given @p rol_vector by overwriting its
     * contents.
     *
     * If the current wrapped vector has ghost elements,
     * then <code> VectorType::operator=(const VectorType&) </code> should still
     * be allowed on it.
     *
     * @note @p rol_vector has to be of type ROLAdaptor.
     */
    void
    set(const ROL::Vector<value_type> &rol_vector) override;

    /**
     * Add @p rol_vector to the wrapped vector.
     *
     * The operation will only be applied to elements of the specified
     * optimization space. Both vectors need to correspond to the same
     * optimization space.
     *
     * @note @p rol_vector has to be of type ROLAdaptor.
     */
    void
    plus(const ROL::Vector<value_type> &rol_vector) override;

    /**
     * Scale the wrapped vector by @p alpha and add @p rol_vector to it.
     *
     * The operation will only be applied to elements of the specified
     * optimization space. Both vectors need to correspond to the same
     * optimization space.
     *
     * @note @p rol_vector has to be of type ROLAdaptor.
     */
    void
    axpy(const value_type               alpha,
         const ROL::Vector<value_type> &rol_vector) override;

    /**
     * Scale the wrapped vector by @p alpha.
     *
     * The operation will only be applied to elements of the specified
     * optimization space.
     */
    void
    scale(const value_type alpha) override;

    /**
     * Return the dot product with a given @p rol_vector.
     *
     * The operation will only be applied to elements of the specified
     * optimization space. Both vectors need to correspond to the same
     * optimization space.
     *
     * @note @p rol_vector has to be of type ROLAdaptor.
     */
    value_type
    dot(const ROL::Vector<value_type> &rol_vector) const override;

    /**
     * Return the $L^{2}$ norm of the wrapped vector.
     *
     * The operation will only be applied to elements of the specified
     * optimization space.
     *
     * The returned type is of ROLAdaptor::value_type so as to maintain
     * consistency with ROL::Vector<ROLAdaptor::value_type> and
     * more importantly to not to create an overloaded version namely,
     * <code> ROLAdaptor::real_type norm() const; </code>
     * if real_type and value_type are not of the same type.
     */
    value_type
    norm() const override;

    /**
     * Create and return a ROL pointer to a clone of the wrapped vector.
     * The cloned vector has the same size as the wrapped vector, but is
     * initialized with zeros.
     */
    ROL::Ptr<ROL::Vector<value_type>>
    clone() const override;

    /**
     * Create and return a ROL pointer to the basis vector corresponding to the
     * @p i ${}^{th}$ element of the global optimization space.
     */
    ROL::Ptr<ROL::Vector<value_type>>
    basis(const int i) const override;

    /**
     * Apply unary function @p f to all the elements of the wrapped vector.
     *
     * The operation will only be applied to elements of the specified
     * optimization space.
     */
    void
    applyUnary(const ROL::Elementwise::UnaryFunction<value_type> &f) override;

    /**
     * Apply binary function @p f along with ROL::Vector @p rol_vector to all
     * the elements of the wrapped vector.
     *
     * The operation will only be applied to elements of the specified
     * optimization space. Both vectors need to correspond to the same
     * optimization space.
     *
     * @note @p rol_vector has to be of type ROLAdaptor.
     */
    void
    applyBinary(const ROL::Elementwise::BinaryFunction<value_type> &f,
                const ROL::Vector<value_type> &rol_vector) override;

    /**
     * Return the accumulated value on applying reduction operation @p r on
     * all the elements of the wrapped vector.
     *
     * The operation will only be applied to elements of the specified
     * optimization space.
     */
    value_type
    reduce(const ROL::Elementwise::ReductionOp<value_type> &r) const override;

    /**
     * Print the wrapped vector to the output stream @p outStream.
     *
     * This function will print the entire vector, and not just the elements of
     * the specified optimization space.
     */
    void
    print(std::ostream &outStream) const override;
  };


  /*------------------------------member definitions--------------------------*/
#  ifndef DOXYGEN


  template <typename VectorType>
  ROLAdaptor<VectorType>::ROLAdaptor(const ROL::Ptr<VectorType> &vector_ptr)
    : ROLAdaptor(vector_ptr, vector_ptr->locally_owned_elements())
  {}



  template <typename VectorType>
  ROLAdaptor<VectorType>::ROLAdaptor(const ROL::Ptr<VectorType> &vector_ptr,
                                     const IndexSet &optimization_space)
    : vector_ptr(vector_ptr)
    , optimization_space(optimization_space)
  {
    Assert(
      optimization_space.is_subset_of(vector_ptr->locally_owned_elements()),
      ExcMessage(
        "Provided IndexSet needs to be a subset of locally owned indices."));

    std::tie(local_opt_start_index, global_opt_dimension) =
      Utilities::MPI::partial_and_total_sum(optimization_space.n_elements(),
                                            vector_ptr->get_mpi_communicator());

    // The return type of ROL::Vector::dimension() has to be int,
    // so we check this requirement of ROL here.
    Assert(global_opt_dimension < static_cast<dealii::types::global_dof_index>(
                                    std::numeric_limits<int>::max()),
           ExcMessage("The number of elements to optimize is greater than the "
                      "largest value of type int."));
  }



  template <typename VectorType>
  ROL::Ptr<VectorType>
  ROLAdaptor<VectorType>::getVector()
  {
    return vector_ptr;
  }



  template <typename VectorType>
  ROL::Ptr<const VectorType>
  ROLAdaptor<VectorType>::getVector() const
  {
    return vector_ptr;
  }



  template <typename VectorType>
  void
  ROLAdaptor<VectorType>::set(const ROL::Vector<value_type> &other_)
  {
    Assert(dynamic_cast<const ROLAdaptor *>(&other_) != nullptr,
           ExcInternalError());
    const ROLAdaptor &other = dynamic_cast<const ROLAdaptor &>(other_);

    // Perform a deep copy of the vector pointed to by vector_ptr.
    (*vector_ptr) = *(other.getVector());

    optimization_space    = other.optimization_space;
    global_opt_dimension  = other.global_opt_dimension;
    local_opt_start_index = other.local_opt_start_index;
  }



  template <typename VectorType>
  void
  ROLAdaptor<VectorType>::plus(const ROL::Vector<value_type> &other_)
  {
    Assert(dynamic_cast<const ROLAdaptor *>(&other_) != nullptr,
           ExcInternalError());
    const ROLAdaptor &other = dynamic_cast<const ROLAdaptor &>(other_);

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
  ROLAdaptor<VectorType>::axpy(const value_type               alpha,
                               const ROL::Vector<value_type> &other_)
  {
    Assert(dynamic_cast<const ROLAdaptor *>(&other_) != nullptr,
           ExcInternalError());
    const ROLAdaptor &other = dynamic_cast<const ROLAdaptor &>(other_);

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
  ROLAdaptor<VectorType>::dimension() const
  {
    Assert(optimization_space.size() == vector_ptr->size(),
           ExcMessage("Optimization space is out-of-sync. "
                      "Please create a new wrapper."));

    return static_cast<int>(global_opt_dimension);
  }



  template <typename VectorType>
  void
  ROLAdaptor<VectorType>::scale(const value_type alpha)
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
  ROLAdaptor<VectorType>::dot(const ROL::Vector<value_type> &other_) const
  {
    Assert(dynamic_cast<const ROLAdaptor *>(&other_) != nullptr,
           ExcInternalError());
    const ROLAdaptor &other = dynamic_cast<const ROLAdaptor &>(other_);

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
  ROLAdaptor<VectorType>::norm() const
  {
    return std::sqrt(this->dot(*this));
  }



  template <typename VectorType>
  ROL::Ptr<ROL::Vector<typename VectorType::value_type>>
  ROLAdaptor<VectorType>::clone() const
  {
    // create new vector with same size as wrapped one
    ROL::Ptr<VectorType> clone_ptr = ROL::makePtr<VectorType>();
    clone_ptr->reinit(*vector_ptr, false);

    return ROL::makePtr<ROLAdaptor>(clone_ptr, optimization_space);
    // TODO: also somehow pass the partial sums etc?
  }



  template <typename VectorType>
  ROL::Ptr<ROL::Vector<typename VectorType::value_type>>
  ROLAdaptor<VectorType>::basis(const int global_opt_index) const
  {
    Assert(optimization_space.size() == vector_ptr->size(),
           ExcMessage("Optimization space is out-of-sync. "
                      "Please create a new wrapper."));
    AssertIndexRange(global_opt_index, global_opt_dimension);

    // clone the internal vector as basis
    ROL::Ptr<VectorType> basis_ptr = ROL::makePtr<VectorType>();
    basis_ptr->reinit(*vector_ptr, false);

    // check whether the i-th basis corresponds to one of the
    // locally owned elements of the wrapped vector
    if ((global_opt_index >= local_opt_start_index) &&
        (global_opt_index <
         local_opt_start_index + optimization_space.n_elements()))
      {
        const dealii::types::global_dof_index local_opt_index =
          global_opt_index - local_opt_start_index;

        // global dof index corresponding to i-th basis in optimization
        const dealii::types::global_dof_index global_dof_index =
          optimization_space.nth_index_in_set(local_opt_index);

        (*basis_ptr)[global_dof_index] = 1.;
      }

    basis_ptr->compress(VectorOperation::insert);

    return ROL::makePtr<ROLAdaptor>(basis_ptr, optimization_space);
    // TODO: also somehow pass the partial sums etc?
  }



  template <typename VectorType>
  void
  ROLAdaptor<VectorType>::applyUnary(
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
  ROLAdaptor<VectorType>::applyBinary(
    const ROL::Elementwise::BinaryFunction<value_type> &f,
    const ROL::Vector<value_type>                      &other_)
  {
    Assert(dynamic_cast<const ROLAdaptor *>(&other_) != nullptr,
           ExcInternalError());
    const ROLAdaptor &other = dynamic_cast<const ROLAdaptor &>(other_);

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
  ROLAdaptor<VectorType>::reduce(
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
  ROLAdaptor<VectorType>::print(std::ostream &outStream) const
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

#endif // dealii_trilinos_rol_adaptor_h
