// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_sundials_n_vector_h
#define dealii_sundials_n_vector_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SUNDIALS
#  include <sundials/sundials_nvector.h>

#  include <functional>
#  include <memory>

DEAL_II_NAMESPACE_OPEN

#  ifndef DOXYGEN
namespace SUNDIALS
{
  namespace internal
  {
    template <typename VectorType>
    class NVectorView;
  }
} // namespace SUNDIALS
#  endif

namespace SUNDIALS
{
  namespace internal
  {
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
    /**
     * Create a NVectorView of the given @p vector.
     *
     * This call is intended to be used as
     *
     * @code
     *   auto view = make_nvector_view(vector);
     * @endcode
     *
     * The resulting object `view` must be kept around as long as any other
     * object will use the internally viewed N_Vector.
     *
     * @tparam VectorType Type of the viewed vector. This parameter can be
     *   deduced automatically and will respect a potential const-qualifier.
     * @param vector The vector to view.
     * @param nvector_context A SUNDIALS context object to be used for the
     *   operations we want to do on this vector.
     * @return A NVectorView of the @p vector.
     *
     * @related NVectorView
     *
     * @note Versions of this function prior to SUNDIALS 6.0 did not take a
     * SUNContext argument.
     */
    template <typename VectorType>
    NVectorView<VectorType>
    make_nvector_view(VectorType &vector, SUNContext nvector_context);
#  else
    /**
     * Create a NVectorView of the given @p vector.
     *
     * This call is intended to be used as
     *
     * @code
     *   auto view = make_nvector_view(vector);
     * @endcode
     *
     * The resulting object `view` must be kept around as long as any other
     * object will use the internally viewed N_Vector.
     *
     * @tparam VectorType Type of the viewed vector. This parameter can be
     *   deduced automatically and will respect a potential const-qualifier.
     * @param vector The vector to view.
     * @param nvector_context A SUNDIALS context object to be used for the
     *   operations we want to do on this vector.
     * @return A NVectorView of the @p vector.
     *
     * @related NVectorView
     *
     * @note This version of make_nvector_view() is only for versions of
     * SUNDIALS prior to 6.0 which did not need a SUNContext object to set up a
     * vector.
     */
    template <typename VectorType>
    NVectorView<VectorType>
    make_nvector_view(VectorType &vector);
#  endif

    /**
     * Retrieve the underlying vector attached to N_Vector @p v. This call will
     * only succeed if the underlying vector is not const. Use
     * unwrap_nvector_const() for this case.
     *
     * @note Users must ensure that they ask for the correct VectorType when
     *   calling this function and there are no type-safety checks in place.
     *
     * @tparam VectorType Type of the vector that is stored in @p v
     * @param v Vector to unwrap
     * @return The vector that is stored inside @p v
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
     * @tparam VectorType Type of the vector that is stored in @p v
     * @param v Vector to unwrap
     * @return The vector that is stored inside @p v
     */
    template <typename VectorType>
    const VectorType *
    unwrap_nvector_const(N_Vector v);

    /**
     * A view to a vector which can be used whenever a N_Vector is required.
     *
     * Objects of this class should preferably be created by
     * make_nvector_view() as
     *
     * @code
     *   auto view = make_nvector_view(vector);
     * @endcode
     *
     * The resulting N_Vector is a view of the actual vector and not owning
     * memory. Also, N_VDestroy() cannot be called on the resulting N_Vector
     * since this would lead to a double delete in the destructor.
     *
     * @note SUNDIALS will never call N_VDestroy() on a vector it didn't create
     *   itself and thus the above constraint is not limiting the user.
     *
     * @tparam VectorType Type of the vector that is stored in @p v
     */
    template <typename VectorType>
    class NVectorView
    {
    public:
      /**
       * Default constructor.
       *
       * The object is not actually viewing anything and needs to be assigned to
       * with operator=(NVectorView &&).
       */
      NVectorView() = default;

#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
      /**
       * Constructor. Create view of @p vector.
       *
       * @note Versions of this function prior to SUNDIALS 6.0 did not take a
       * SUNContext argument.
       */
      NVectorView(VectorType &vector, SUNContext nvector_context);
#  else
      /**
       * Constructor. Create view of @p vector.
       *
       * @note This constructor is only for versions of SUNDIALS prior to 6.0
       * which do not need a SUNContext object to set up a vector.
       */
      NVectorView(VectorType &vector);
#  endif

      /**
       * Move assignment.
       */
      NVectorView(NVectorView &&) noexcept = default;

      /**
       * Move constructor.
       */
      NVectorView &
      operator=(NVectorView &&) noexcept = default;

      /**
       * Explicitly delete copy ctor. This class is move-only.
       */
      NVectorView(const NVectorView &) = delete;

      /**
       * Explicitly delete copy assignment. This class is move-only.
       */
      NVectorView &
      operator=(const NVectorView &) = delete;

      /**
       * Destructor.
       *
       * @note This will not destroy the viewed vector.
       */
      ~NVectorView() = default;

      /**
       * Implicit conversion to N_Vector. This operator makes the NVectorView
       * look like an actual N_Vector and it can be used directly as an
       * argument in many SUNDIALS functions.
       */
      operator N_Vector() const;

      /**
       * Access the N_Vector that is viewed by this object.
       */
      N_Vector
      operator->() const;

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

#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif
#endif
