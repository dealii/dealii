//-----------------------------------------------------------
//
//    Copyright (C) 2020 by the deal.II authors
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
//-----------------------------------------------------------

#ifndef dealii_sundials_n_vector_h
#define dealii_sundials_n_vector_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SUNDIALS
#  if DEAL_II_SUNDIALS_VERSION_GTE(5, 0, 0)

#    include <sundials/sundials_nvector.h>

DEAL_II_NAMESPACE_OPEN

namespace SUNDIALS
{
  namespace internal
  {
    /**
     * Retrieve the underlying vector attached to N_Vector @p v. This call will
     * only succeed if the underlying vector is not const. Use
     * unwrap_nvector_const() for this case.
     *
     * @note Users must ensure that they ask for the correct VectorType when
     *   calling this function and there are no type-safety checks in place.
     *
     * @tparam VectorType type of the vector that is stored in @p v
     * @param v vector to unwrap
     * @return the vector that is stored inside @p v
     */
    template <typename VectorType>
    VectorType *
    unwrap_nvector(N_Vector v);

    /**
     * Retrieve the underlying vector attached to N_Vector @p v as a constant
     * pointer.
     *
     * @note Users must ensure that they ask for the correct VectorType when
     *   calling this function and there are no type-safety checks in place.
     *
     * @tparam VectorType type of the vector that is stored in @p v
     * @param v vector to unwrap
     * @return the vector that is stored inside @p v
     */
    template <typename VectorType>
    const VectorType *
    unwrap_nvector_const(N_Vector v);


    /**
     * A view to an N_Vector which can be used whenever a N_Vector is required.
     * An object of this class will automatically clean up all internal data for
     * the view when it is destroyed. This doesn't mean that the actual vector
     * is deleted though, which completely depends on how the N_Vector has been
     * created.
     *
     * @note N_VDestroy() should not be called on the resulting N_Vector since
     *   this would lead to a double delete. Let the destructor do this work.
     */
    template <typename VectorType>
    class NVectorView
    {
    public:
      /**
       * Constructor. This overload is chosen to view a non-const @p vector.
       */
      NVectorView(VectorType &vector);

      /**
       * Constructor. This overload is chosen to view a const @p vector.
       */
      NVectorView(const VectorType &vector);

      /**
       * Implicit conversion to N_Vector. This makes the NVectorView look like
       * an actual N_Vector and it can be used directly as an argument.
       */
      operator N_Vector() const;

      /**
       * Access the N_Vector that is viewed by this object.
       */
      N_Vector operator->() const;

      /**
       * Destructor. Automatically release all memory that was only allocated
       * for the view.
       */
      ~NVectorView() = default;

    private:
      /**
       * Actual pointer to a vector viewed by this class.
       */
      std::unique_ptr<_generic_N_Vector, std::function<void(N_Vector)>>
        vector_ptr;
    };

  } // namespace internal
} // namespace SUNDIALS

DEAL_II_NAMESPACE_CLOSE

#  endif
#endif
#endif
