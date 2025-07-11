// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_sundials_utilities_h
#define dealii_sundials_utilities_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>


#ifdef DEAL_II_WITH_SUNDIALS
#  include <exception>


DEAL_II_NAMESPACE_OPEN

namespace SUNDIALS
{
  namespace Utilities
  {
    /**
     * A function that calls the function object given by its first argument
     * with the set of arguments following at the end. If the call returns
     * regularly, the current function returns zero to indicate success. If the
     * call fails with an exception of type RecoverableUserCallbackError, then
     * the current function returns 1 to indicate that the called function
     * object thought the error it encountered is recoverable. If the call fails
     * with any other exception, then the current function returns with an error
     * code of -1. In each of the last two cases, the exception thrown by `f`
     * is captured and `eptr` is set to the exception. In case of success,
     * `eptr` is set to `nullptr`.
     */
    template <typename F, typename... Args>
    int
    call_and_possibly_capture_exception(const F            &f,
                                        std::exception_ptr &eptr,
                                        Args &&...args)
    {
      // See whether there is already something in the exception pointer
      // variable. This can only happen if we had previously had
      // a recoverable exception, and the underlying library actually
      // did recover successfully. In that case, we can abandon the
      // exception previously thrown. If eptr contains anything other,
      // then we really don't know how that could have happened, and
      // should probably bail out:
      if (eptr)
        {
          try
            {
              std::rethrow_exception(eptr);
            }
          catch (const RecoverableUserCallbackError &)
            {
              // ok, ignore, but reset the pointer
              eptr = nullptr;
            }
          catch (...)
            {
              // uh oh:
              AssertThrow(false, ExcInternalError());
            }
        }

      // Call the function and if that succeeds, return zero:
      try
        {
          f(std::forward<Args>(args)...);
          eptr = nullptr;
          return 0;
        }
      // If the call failed with a recoverable error, then
      // ignore the exception for now (but store a pointer to it)
      // and return a positive return value (+1). If the underlying
      // implementation manages to recover
      catch (const RecoverableUserCallbackError &)
        {
          eptr = std::current_exception();
          return 1;
        }
      // For any other exception, capture the exception and
      // return -1:
      catch (const std::exception &)
        {
          eptr = std::current_exception();
          return -1;
        }
    }
  } // namespace Utilities
} // namespace SUNDIALS

DEAL_II_NAMESPACE_CLOSE

#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SUNDIALS

#endif
