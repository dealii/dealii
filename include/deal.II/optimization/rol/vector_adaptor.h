//-----------------------------------------------------------
//
//    Copyright (C) 2017 - 2020 by the deal.II authors
//
//    This file is part of the deal.II library.
//
//    The deal.II library is free software; you can use it, redistribute
//    it, and/or modify it under the terms of the GNU Lesser General
//    Public License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//    The full text of the license can be found in the file LICENSE.md at
//    the top level directory of deal.II.
//
//---------------------------------------------------------------

#ifndef dealii_optimization_rol_vector_adaptor_h
#define dealii_optimization_rol_vector_adaptor_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_TRILINOS_WITH_ROL
#  include <deal.II/base/exceptions.h>

#  include <deal.II/lac/vector.h>

#  include <ROL_Vector.hpp>

#  include <type_traits>


DEAL_II_NAMESPACE_OPEN

/**
 * A namespace that provides an interface to the
 * <a href="https://trilinos.org/docs/dev/packages/rol/doc/html/index.html">
 * Rapid Optimization Library</a> (ROL), a Trilinos package.
 */
namespace Rol
{
  /**
   * An adaptor that provides the implementation of the ROL::Vector interface
   * for vectors of type <tt>VectorType</tt>.
   *
   * This class supports vectors that satisfy the following requirements:
   *
   * The <tt>VectorType</tt> should contain the following types.
   * ```
   * VectorType::size_type;  // The type for size of the vector.
   * VectorType::value_type; // The type for elements stored in the vector.
   * VectorType::real_type;  // The type for real-valued numbers.
   * ```
   *
   * However, ROL doesn't distinguish VectorAdaptor::value_type from
   * VectorAdaptor::real_type. This is due to ROL's assumption that the
   * VectorAdaptor::value_type itself is a type for real-valued numbers.
   * Therefore, VectorAdaptor supports vectors whose real_type is
   * convertible to value_type in the sense that
   * <code>std::is_convertible<real_type, value_type>::value</code> yields
   * <code>true</code>.
   *
   * The <tt>VectorType</tt> should contain the following methods.
   * @code
   * // Reinitialize the current vector using a given vector's
   * // size (and the parallel distribution) without copying
   * // the elements.
   * VectorType::reinit(const VectorType &, ...);
   *
   * // Globally add a given vector to the current.
   * VectorType::operator+=(const VectorType &);
   *
   * // Scale all elements by a given scalar.
   * VectorType::operator*=(const VectorType::value_type &);
   *
   * // Perform dot product with a given vector.
   * VectorType::operator*=(const VectorType &);
   *
   * // Scale all elements of the current vector and globally
   * // add a given vector to it.
   * VectorType::add(const VectorType::value_type, const VectorType &);
   *
   * // Copies the data of a given vector to the current.
   * // Resize the current vector if necessary (MPI safe).
   * VectorType::operation=(const VectorType &);
   *
   * // Return the global size of the current vector.
   * VectorType::size();
   *
   * // Return L^2 norm of the current vector
   * VectorType::l2_norm();
   *
   * // Iterator to the start of the (locally owned) element
   * // of the current vector.
   * VectorType::begin();
   *
   * // Iterator to the one past the last (locally owned)
   * // element of the current vector.
   * VectorType::end();
   *
   * // Compress the vector i.e., flush the buffers of the
   * // vector object if it has any.
   * VectorType::compress(VectorOperation::insert);
   * @endcode
   *
   * @note The current implementation in ROL doesn't support vector sizes above
   * the largest value of int type. Some of the vectors in deal.II (see
   * @ref Vector)
   * may not satisfy the above requirements.
   *
   * @author Vishal Boddu, 2017
   */
  template <typename VectorType>
  class VectorAdaptor : public ROL::Vector<typename VectorType::value_type>
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

    static_assert(std::is_convertible<real_type, value_type>::value,
                  "The real_type of the current VectorType is not "
                  "convertible to the value_type.");

  private:
    /**
     * Teuchos smart reference counting pointer to the underlying vector of type
     * <tt>VectorType</tt>.
     */
    Teuchos::RCP<VectorType> vector_ptr;

  public:
    /**
     * Constructor.
     */
    VectorAdaptor(const Teuchos::RCP<VectorType> &vector_ptr);

    /**
     * Return the Teuchos smart reference counting pointer to
     * the wrapper vector, #vector_ptr.
     */
    Teuchos::RCP<VectorType>
    getVector();

    /**
     * Return the Teuchos smart reference counting pointer to const vector.
     */
    Teuchos::RCP<const VectorType>
    getVector() const;

    /**
     * Return the dimension (global vector size) of the wrapped vector.
     */
    int
    dimension() const;

    /**
     * Set the wrapper vector to a given ROL::Vector @p rol_vector by
     * overwriting its contents.
     *
     * If the current wrapper vector has ghost elements,
     * then <code> VectorType::operator=(const VectorType&) </code> should still
     * be allowed on it.
     */
    void
    set(const ROL::Vector<value_type> &rol_vector);

    /**
     * Perform addition.
     */
    void
    plus(const ROL::Vector<value_type> &rol_vector);

    /**
     * Scale the wrapper vector by @p alpha and add ROL::Vector @p rol_vector
     * to it.
     */
    void
    axpy(const value_type alpha, const ROL::Vector<value_type> &rol_vector);

    /**
     * Scale the wrapper vector.
     */
    void
    scale(const value_type alpha);

    /**
     * Return the dot product with a given ROL::Vector @p rol_vector.
     */
    value_type
    dot(const ROL::Vector<value_type> &rol_vector) const;

    /**
     * Return the $L^{2}$ norm of the wrapped vector.
     *
     * The returned type is of VectorAdaptor::value_type so as to maintain
     * consistency with ROL::Vector<VectorAdaptor::value_type> and
     * more importantly to not to create an overloaded version namely,
     * <code> VectorAdaptor::real_type norm() const; </code>
     * if real_type and value_type are not of the same type.
     */
    value_type
    norm() const;

    /**
     * Return a clone of the wrapped vector.
     */
    Teuchos::RCP<ROL::Vector<value_type>>
    clone() const;

    /**
     * Create and return a Teuchos smart reference counting pointer to the basis
     * vector corresponding to the @p i ${}^{th}$ element of
     * the wrapper vector.
     */
    Teuchos::RCP<ROL::Vector<value_type>>
    basis(const int i) const;

    /**
     * Apply unary function @p f to all the elements of the wrapped vector.
     */
    void
    applyUnary(const ROL::Elementwise::UnaryFunction<value_type> &f);

    /**
     * Apply binary function @p f along with ROL::Vector @p rol_vector to all
     * the elements of the wrapped vector.
     */
    void
    applyBinary(const ROL::Elementwise::BinaryFunction<value_type> &f,
                const ROL::Vector<value_type> &                     rol_vector);

    /**
     * Return the accumulated value on applying reduction operation @p r on
     * all the elements of the wrapped vector.
     */
    value_type
    reduce(const ROL::Elementwise::ReductionOp<value_type> &r) const;

    /**
     * Print the wrapped vector to the output stream @p outStream.
     */
    void
    print(std::ostream &outStream) const;
  };


  /*------------------------------member definitions--------------------------*/
#  ifndef DOXYGEN


  template <typename VectorType>
  VectorAdaptor<VectorType>::VectorAdaptor(
    const Teuchos::RCP<VectorType> &vector_ptr)
    : vector_ptr(vector_ptr)
  {}



  template <typename VectorType>
  Teuchos::RCP<VectorType>
  VectorAdaptor<VectorType>::getVector()
  {
    return vector_ptr;
  }



  template <typename VectorType>
  Teuchos::RCP<const VectorType>
  VectorAdaptor<VectorType>::getVector() const
  {
    return vector_ptr;
  }



  template <typename VectorType>
  void
  VectorAdaptor<VectorType>::set(const ROL::Vector<value_type> &rol_vector)
  {
    const VectorAdaptor &vector_adaptor =
      Teuchos::dyn_cast<const VectorAdaptor>(rol_vector);

    (*vector_ptr) = *(vector_adaptor.getVector());
  }



  template <typename VectorType>
  void
  VectorAdaptor<VectorType>::plus(const ROL::Vector<value_type> &rol_vector)
  {
    Assert(this->dimension() == rol_vector.dimension(),
           ExcDimensionMismatch(this->dimension(), rol_vector.dimension()));

    const VectorAdaptor &vector_adaptor =
      Teuchos::dyn_cast<const VectorAdaptor>(rol_vector);

    *vector_ptr += *(vector_adaptor.getVector());
  }



  template <typename VectorType>
  void
  VectorAdaptor<VectorType>::axpy(const value_type               alpha,
                                  const ROL::Vector<value_type> &rol_vector)
  {
    Assert(this->dimension() == rol_vector.dimension(),
           ExcDimensionMismatch(this->dimension(), rol_vector.dimension()));

    const VectorAdaptor &vector_adaptor =
      Teuchos::dyn_cast<const VectorAdaptor>(rol_vector);

    vector_ptr->add(alpha, *(vector_adaptor.getVector()));
  }



  template <typename VectorType>
  int
  VectorAdaptor<VectorType>::dimension() const
  {
    Assert(vector_ptr->size() < std::numeric_limits<int>::max(),
           ExcMessage("The size of the vector being used is greater than "
                      "largest value of type int."));
    return static_cast<int>(vector_ptr->size());
  }



  template <typename VectorType>
  void
  VectorAdaptor<VectorType>::scale(const value_type alpha)
  {
    (*vector_ptr) *= alpha;
  }



  template <typename VectorType>
  typename VectorType::value_type
  VectorAdaptor<VectorType>::dot(
    const ROL::Vector<value_type> &rol_vector) const
  {
    Assert(this->dimension() == rol_vector.dimension(),
           ExcDimensionMismatch(this->dimension(), rol_vector.dimension()));

    const VectorAdaptor &vector_adaptor =
      Teuchos::dyn_cast<const VectorAdaptor>(rol_vector);

    return (*vector_ptr) * (*vector_adaptor.getVector());
  }



  template <typename VectorType>
  typename VectorType::value_type
  VectorAdaptor<VectorType>::norm() const
  {
    return vector_ptr->l2_norm();
  }



  template <typename VectorType>
  Teuchos::RCP<ROL::Vector<typename VectorType::value_type>>
  VectorAdaptor<VectorType>::clone() const
  {
    Teuchos::RCP<VectorType> vec_ptr = Teuchos::rcp(new VectorType);
    (*vec_ptr)                       = (*vector_ptr);

    return Teuchos::rcp(new VectorAdaptor(vec_ptr));
  }



  template <typename VectorType>
  Teuchos::RCP<ROL::Vector<typename VectorType::value_type>>
  VectorAdaptor<VectorType>::basis(const int i) const
  {
    Teuchos::RCP<VectorType> vec_ptr = Teuchos::rcp(new VectorType);

    // Zero all the entries in dealii vector.
    vec_ptr->reinit(*vector_ptr, false);

    if (vector_ptr->locally_owned_elements().is_element(i))
      vec_ptr->operator[](i) = 1.;

    if (vec_ptr->has_ghost_elements())
      {
        vec_ptr->update_ghost_values();
      }
    else
      {
        vec_ptr->compress(VectorOperation::insert);
      }

    Teuchos::RCP<VectorAdaptor> e = Teuchos::rcp(new VectorAdaptor(vec_ptr));

    return e;
  }



  template <typename VectorType>
  void
  VectorAdaptor<VectorType>::applyUnary(
    const ROL::Elementwise::UnaryFunction<value_type> &f)
  {
    const typename VectorType::iterator vend = vector_ptr->end();

    for (typename VectorType::iterator iterator = vector_ptr->begin();
         iterator != vend;
         iterator++)
      *iterator = f.apply(*iterator);

    if (vector_ptr->has_ghost_elements())
      {
        vector_ptr->update_ghost_values();
      }
    else
      {
        vector_ptr->compress(VectorOperation::insert);
      }
  }



  template <typename VectorType>
  void
  VectorAdaptor<VectorType>::applyBinary(
    const ROL::Elementwise::BinaryFunction<value_type> &f,
    const ROL::Vector<value_type> &                     rol_vector)
  {
    Assert(this->dimension() == rol_vector.dimension(),
           ExcDimensionMismatch(this->dimension(), rol_vector.dimension()));

    const VectorAdaptor &vector_adaptor =
      Teuchos::dyn_cast<const VectorAdaptor>(rol_vector);

    const VectorType &given_rol_vector = *(vector_adaptor.getVector());

    const typename VectorType::iterator       vend   = vector_ptr->end();
    const typename VectorType::const_iterator rolend = given_rol_vector.end();

    typename VectorType::const_iterator r_iterator = given_rol_vector.begin();
    for (typename VectorType::iterator l_iterator = vector_ptr->begin();
         l_iterator != vend && r_iterator != rolend;
         l_iterator++, r_iterator++)
      *l_iterator = f.apply(*l_iterator, *r_iterator);

    if (vector_ptr->has_ghost_elements())
      {
        vector_ptr->update_ghost_values();
      }
    else
      {
        vector_ptr->compress(VectorOperation::insert);
      }
  }



  template <typename VectorType>
  typename VectorType::value_type
  VectorAdaptor<VectorType>::reduce(
    const ROL::Elementwise::ReductionOp<value_type> &r) const
  {
    typename VectorType::value_type result = r.initialValue();

    const typename VectorType::iterator vend = vector_ptr->end();

    for (typename VectorType::iterator iterator = vector_ptr->begin();
         iterator != vend;
         iterator++)
      r.reduce(*iterator, result);
    // Parallel reduction among processes is handled internally.

    return result;
  }



  template <typename VectorType>
  void
  VectorAdaptor<VectorType>::print(std::ostream &outStream) const
  {
    vector_ptr->print(outStream);
  }


#  endif // DOXYGEN


} // namespace Rol


DEAL_II_NAMESPACE_CLOSE


#endif // DEAL_II_TRILINOS_WITH_ROL

#endif // dealii_optimization_rol_vector_adaptor_h
