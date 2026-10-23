// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#ifndef dealii_sundials_sundials_handle_h
#define dealii_sundials_sundials_handle_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SUNDIALS

#  include <deal.II/base/exceptions.h>

#  include <utility>

#endif

DEAL_II_NAMESPACE_OPEN

#ifdef DEAL_II_WITH_SUNDIALS
namespace SUNDIALS
{
  namespace internal
  {
    /**
     * An owning handle for a C object created by SUNDIALS, for example a
     * SUNLinearSolver or a SUNNonlinearSolver. The object is freed in the
     * destructor by the function provided by the @p Traits class, which has to
     * define
     * - a type `Handle` (the SUNDIALS pointer type), and
     * - a static function `free(Handle)` that releases the object.
     *
     * The handle takes ownership of the object passed to the constructor, so
     * any SUNDIALS creation function can be used:
     * @code
     *   SundialsLinearSolver solver(SUNLinSol_PCG(y, SUN_PREC_NONE, 0, ctx));
     * @endcode
     *
     * The handle can be moved but not copied.
     */
    template <typename Traits>
    class SundialsHandle
    {
    public:
      using Handle = typename Traits::Handle;

      /**
       * Constructor. Takes ownership of @p handle, which must not be null
       * (SUNDIALS returns null if the creation of an object fails).
       */
      explicit SundialsHandle(Handle handle)
        : handle(handle)
      {
        AssertThrow(handle != nullptr,
                    ExcMessage("SUNDIALS failed to create an object."));
      }

      ~SundialsHandle()
      {
        reset();
      }

      SundialsHandle(const SundialsHandle &) = delete;
      SundialsHandle &
      operator=(const SundialsHandle &) = delete;

      SundialsHandle(SundialsHandle &&other) noexcept
        : handle(std::exchange(other.handle, nullptr))
      {}

      SundialsHandle &
      operator=(SundialsHandle &&other) noexcept
      {
        if (this != &other)
          {
            reset();
            handle = std::exchange(other.handle, nullptr);
          }
        return *this;
      }

      /**
       * Implicit conversion to the SUNDIALS pointer. The handle keeps the
       * ownership.
       */
      operator Handle() const
      {
        return handle;
      }

    private:
      void
      reset() noexcept
      {
        if (handle != nullptr)
          {
            Traits::free(handle);
            handle = nullptr;
          }
      }

      Handle handle;
    };
  } // namespace internal
} // namespace SUNDIALS
#endif

DEAL_II_NAMESPACE_CLOSE
#endif
