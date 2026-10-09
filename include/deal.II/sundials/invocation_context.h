// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_sundials_invocation_context_h
#define dealii_sundials_invocation_context_h

#include <deal.II/base/config.h>


#ifdef DEAL_II_WITH_SUNDIALS

#  include <deal.II/base/exceptions.h>

#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
#    include <sundials/sundials_context.h>
#    include <sundials/sundials_types.h>
#  endif

#  include <exception>
#  include <memory>
#  include <type_traits>

#endif // DEAL_II_WITH_SUNDIALS


DEAL_II_NAMESPACE_OPEN

#ifdef DEAL_II_WITH_SUNDIALS
namespace SUNDIALS
{
  namespace internal
  {
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
    /**
     * A SUNContext that is freed when the last copy of this object is
     * destroyed. SUNDIALS objects refer to the context they were created with,
     * even while they are freed, so everything that owns such objects also has
     * to own the context. The object converts to the plain SUNContext that
     * SUNDIALS functions expect.
     */
    class SharedSUNContext
    {
    public:
      /**
       * Construct an object that holds no context.
       */
      SharedSUNContext() = default;

      /**
       * Take ownership of @p context, which is freed with SUNContext_Free().
       */
      explicit SharedSUNContext(SUNContext context)
        : context(context, [](std::remove_pointer_t<SUNContext> *c) {
          SUNContext ctx    = c;
          const int  status = SUNContext_Free(&ctx);
          (void)status;
          AssertNothrow(status == 0, ExcMessage("SUNContext_Free() failed."));
        })
      {}

      /**
       * Access the underlying SUNContext.
       */
      operator SUNContext() const
      {
        return context.get();
      }

    private:
      std::shared_ptr<std::remove_pointer_t<SUNContext>> context;
    };
#  endif

    /**
     * This struct contains a reference to the std::exception_ptr @p
     * pending_exception and, for SUNDIALS 6.0 and newer, the SUNContext @p
     * arkode_ctx. This is required for situations when an object needs to pass
     * the context to the objects it uses. An object that owns such objects
     * keeps a copy of @p arkode_ctx to keep the context alive.
     */
    struct InvocationContext
    {
      std::exception_ptr &pending_exception;
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
      SharedSUNContext arkode_ctx;
#  endif
    };
  } // namespace internal
} // namespace SUNDIALS

#endif // DEAL_II_WITH_SUNDIALS

DEAL_II_NAMESPACE_CLOSE

#endif
