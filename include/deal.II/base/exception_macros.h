// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_exception_macros_h
#define dealii_exception_macros_h

#include <deal.II/base/config.h>


/**********************************************************************
 * Preprocessor definitions in support of declaring exception classes.
 *
 * These make reference to classes and functions that are declared in
 * exceptions.h. This is ok as far as concerning preprocessor defines
 * is concerned, but in order to use the material below, you also need
 * to have `#include <deal.II/base/exceptions.h>`.
 */
#ifndef DOXYGEN

/**
 * Declare an exception class derived from ExceptionBase without parameters.
 *
 * @note This and similar macro names are examples of preprocessor definitions
 * in the deal.II library that are not prefixed by a string that likely makes
 * them unique to deal.II. As a consequence, it is possible that other
 * libraries your code interfaces with define the same name, and the result
 * will be name collisions (see
 * https://en.wikipedia.org/wiki/Name_collision). One can <code>\#undef</code>
 * this macro, as well as all other macros defined by deal.II that are not
 * prefixed with either <code>DEAL</code> or <code>deal</code>, by including
 * the header <code>deal.II/base/undefine_macros.h</code> after all other
 * deal.II headers have been included.
 *
 * @ingroup Exceptions
 */
#  define DeclException0(Exception0)                \
    class Exception0 : public dealii::ExceptionBase \
    {}


/**
 * Declare an exception class derived from ExceptionBase that can take one
 * runtime argument, but if none is given in the place where you want to throw
 * the exception, it simply reverts to the default text provided when
 * declaring the exception class through this macro.
 *
 * @note This and similar macro names are examples of preprocessor definitions
 * in the deal.II library that are not prefixed by a string that likely makes
 * them unique to deal.II. As a consequence, it is possible that other
 * libraries your code interfaces with define the same name, and the result
 * will be name collisions (see
 * https://en.wikipedia.org/wiki/Name_collision). One can <code>\#undef</code>
 * this macro, as well as all other macros defined by deal.II that are not
 * prefixed with either <code>DEAL</code> or <code>deal</code>, by including
 * the header <code>deal.II/base/undefine_macros.h</code> after all other
 * deal.II headers have been included.
 *
 * @ingroup Exceptions
 */
#  define DeclExceptionMsg(Exception, defaulttext)    \
    class Exception : public dealii::ExceptionBase    \
    {                                                 \
    public:                                           \
      Exception(const std::string &msg = defaulttext) \
        : arg(msg)                                    \
      {}                                              \
      virtual ~Exception() noexcept                   \
      {}                                              \
      virtual void                                    \
      print_info(std::ostream &out) const override    \
      {                                               \
        out << "    " << arg << std::endl;            \
      }                                               \
                                                      \
    private:                                          \
      const std::string arg;                          \
    }

/**
 * Declare an exception class derived from ExceptionBase with one additional
 * parameter.
 *
 * @note This and similar macro names are examples of preprocessor definitions
 * in the deal.II library that are not prefixed by a string that likely makes
 * them unique to deal.II. As a consequence, it is possible that other
 * libraries your code interfaces with define the same name, and the result
 * will be name collisions (see
 * https://en.wikipedia.org/wiki/Name_collision). One can <code>\#undef</code>
 * this macro, as well as all other macros defined by deal.II that are not
 * prefixed with either <code>DEAL</code> or <code>deal</code>, by including
 * the header <code>deal.II/base/undefine_macros.h</code> after all other
 * deal.II headers have been included.
 *
 * @ingroup Exceptions
 */
#  define DeclException1(Exception1, type1, outsequence) \
    class Exception1 : public dealii::ExceptionBase      \
    {                                                    \
    public:                                              \
      Exception1(type1 const &a1)                        \
        : arg1(a1)                                       \
      {}                                                 \
      virtual ~Exception1() noexcept                     \
      {}                                                 \
      virtual void                                       \
      print_info(std::ostream &out) const override       \
      {                                                  \
        out << "    " outsequence << std::endl;          \
      }                                                  \
                                                         \
    private:                                             \
      type1 const arg1;                                  \
    }


/**
 * Declare an exception class derived from ExceptionBase with two additional
 * parameters.
 *
 * @note This and similar macro names are examples of preprocessor definitions
 * in the deal.II library that are not prefixed by a string that likely makes
 * them unique to deal.II. As a consequence, it is possible that other
 * libraries your code interfaces with define the same name, and the result
 * will be name collisions (see
 * https://en.wikipedia.org/wiki/Name_collision). One can <code>\#undef</code>
 * this macro, as well as all other macros defined by deal.II that are not
 * prefixed with either <code>DEAL</code> or <code>deal</code>, by including
 * the header <code>deal.II/base/undefine_macros.h</code> after all other
 * deal.II headers have been included.
 *
 * @ingroup Exceptions
 */
#  define DeclException2(Exception2, type1, type2, outsequence) \
    class Exception2 : public dealii::ExceptionBase             \
    {                                                           \
    public:                                                     \
      Exception2(type1 const &a1, type2 const &a2)              \
        : arg1(a1)                                              \
        , arg2(a2)                                              \
      {}                                                        \
      virtual ~Exception2() noexcept                            \
      {}                                                        \
      virtual void                                              \
      print_info(std::ostream &out) const override              \
      {                                                         \
        out << "    " outsequence << std::endl;                 \
      }                                                         \
                                                                \
    private:                                                    \
      type1 const arg1;                                         \
      type2 const arg2;                                         \
    }


/**
 * Declare an exception class derived from ExceptionBase with three additional
 * parameters.
 *
 * @note This and similar macro names are examples of preprocessor definitions
 * in the deal.II library that are not prefixed by a string that likely makes
 * them unique to deal.II. As a consequence, it is possible that other
 * libraries your code interfaces with define the same name, and the result
 * will be name collisions (see
 * https://en.wikipedia.org/wiki/Name_collision). One can <code>\#undef</code>
 * this macro, as well as all other macros defined by deal.II that are not
 * prefixed with either <code>DEAL</code> or <code>deal</code>, by including
 * the header <code>deal.II/base/undefine_macros.h</code> after all other
 * deal.II headers have been included.
 *
 * @ingroup Exceptions
 */
#  define DeclException3(Exception3, type1, type2, type3, outsequence) \
    class Exception3 : public dealii::ExceptionBase                    \
    {                                                                  \
    public:                                                            \
      Exception3(type1 const &a1, type2 const &a2, type3 const &a3)    \
        : arg1(a1)                                                     \
        , arg2(a2)                                                     \
        , arg3(a3)                                                     \
      {}                                                               \
      virtual ~Exception3() noexcept                                   \
      {}                                                               \
      virtual void                                                     \
      print_info(std::ostream &out) const override                     \
      {                                                                \
        out << "    " outsequence << std::endl;                        \
      }                                                                \
                                                                       \
    private:                                                           \
      type1 const arg1;                                                \
      type2 const arg2;                                                \
      type3 const arg3;                                                \
    }


/**
 * Declare an exception class derived from ExceptionBase with four additional
 * parameters.
 *
 * @note This and similar macro names are examples of preprocessor definitions
 * in the deal.II library that are not prefixed by a string that likely makes
 * them unique to deal.II. As a consequence, it is possible that other
 * libraries your code interfaces with define the same name, and the result
 * will be name collisions (see
 * https://en.wikipedia.org/wiki/Name_collision). One can <code>\#undef</code>
 * this macro, as well as all other macros defined by deal.II that are not
 * prefixed with either <code>DEAL</code> or <code>deal</code>, by including
 * the header <code>deal.II/base/undefine_macros.h</code> after all other
 * deal.II headers have been included.
 *
 * @ingroup Exceptions
 */
#  define DeclException4(Exception4, type1, type2, type3, type4, outsequence) \
    class Exception4 : public dealii::ExceptionBase                           \
    {                                                                         \
    public:                                                                   \
      Exception4(type1 const &a1,                                             \
                 type2 const &a2,                                             \
                 type3 const &a3,                                             \
                 type4 const &a4)                                             \
        : arg1(a1)                                                            \
        , arg2(a2)                                                            \
        , arg3(a3)                                                            \
        , arg4(a4)                                                            \
      {}                                                                      \
      virtual ~Exception4() noexcept                                          \
      {}                                                                      \
      virtual void                                                            \
      print_info(std::ostream &out) const override                            \
      {                                                                       \
        out << "    " outsequence << std::endl;                               \
      }                                                                       \
                                                                              \
    private:                                                                  \
      type1 const arg1;                                                       \
      type2 const arg2;                                                       \
      type3 const arg3;                                                       \
      type4 const arg4;                                                       \
    }


/**
 * Declare an exception class derived from ExceptionBase with five additional
 * parameters.
 *
 * @note This and similar macro names are examples of preprocessor definitions
 * in the deal.II library that are not prefixed by a string that likely makes
 * them unique to deal.II. As a consequence, it is possible that other
 * libraries your code interfaces with define the same name, and the result
 * will be name collisions (see
 * https://en.wikipedia.org/wiki/Name_collision). One can <code>\#undef</code>
 * this macro, as well as all other macros defined by deal.II that are not
 * prefixed with either <code>DEAL</code> or <code>deal</code>, by including
 * the header <code>deal.II/base/undefine_macros.h</code> after all other
 * deal.II headers have been included.
 *
 * @ingroup Exceptions
 */
#  define DeclException5(                                       \
    Exception5, type1, type2, type3, type4, type5, outsequence) \
    class Exception5 : public dealii::ExceptionBase             \
    {                                                           \
    public:                                                     \
      Exception5(type1 const &a1,                               \
                 type2 const &a2,                               \
                 type3 const &a3,                               \
                 type4 const &a4,                               \
                 type5 const &a5)                               \
        : arg1(a1)                                              \
        , arg2(a2)                                              \
        , arg3(a3)                                              \
        , arg4(a4)                                              \
        , arg5(a5)                                              \
      {}                                                        \
      virtual ~Exception5() noexcept                            \
      {}                                                        \
      virtual void                                              \
      print_info(std::ostream &out) const override              \
      {                                                         \
        out << "    " outsequence << std::endl;                 \
      }                                                         \
                                                                \
    private:                                                    \
      type1 const arg1;                                         \
      type2 const arg2;                                         \
      type3 const arg3;                                         \
      type4 const arg4;                                         \
      type5 const arg5;                                         \
    }

#else /*ifndef DOXYGEN*/

// Dummy definitions for doxygen:

/**
 * Declare an exception class derived from ExceptionBase without parameters.
 *
 * @note This and similar macro names are examples of preprocessor definitions
 * in the deal.II library that are not prefixed by a string that likely makes
 * them unique to deal.II. As a consequence, it is possible that other
 * libraries your code interfaces with define the same name, and the result
 * will be name collisions (see
 * https://en.wikipedia.org/wiki/Name_collision). One can <code>\#undef</code>
 * this macro, as well as all other macros defined by deal.II that are not
 * prefixed with either <code>DEAL</code> or <code>deal</code>, by including
 * the header <code>deal.II/base/undefine_macros.h</code> after all other
 * deal.II headers have been included.
 *
 * @ingroup Exceptions
 */
#  define DeclException0(Exception0) \
    /** @ingroup Exceptions */       \
    static dealii::ExceptionBase &Exception0()

/**
 * Declare an exception class derived from ExceptionBase that can take one
 * runtime argument, but if none is given in the place where you want to throw
 * the exception, it simply reverts to the default text provided when
 * declaring the exception class through this macro.
 *
 * @note This and similar macro names are examples of preprocessor definitions
 * in the deal.II library that are not prefixed by a string that likely makes
 * them unique to deal.II. As a consequence, it is possible that other
 * libraries your code interfaces with define the same name, and the result
 * will be name collisions (see
 * https://en.wikipedia.org/wiki/Name_collision). One can <code>\#undef</code>
 * this macro, as well as all other macros defined by deal.II that are not
 * prefixed with either <code>DEAL</code> or <code>deal</code>, by including
 * the header <code>deal.II/base/undefine_macros.h</code> after all other
 * deal.II headers have been included.
 *
 * @ingroup Exceptions
 */
#  define DeclExceptionMsg(Exception, defaulttext) \
    /** @ingroup Exceptions */                     \
    /** @dealiiExceptionMessage{defaulttext} */    \
    static dealii::ExceptionBase &Exception()

/**
 * Declare an exception class derived from ExceptionBase with one additional
 * parameter.
 *
 * @note This and similar macro names are examples of preprocessor definitions
 * in the deal.II library that are not prefixed by a string that likely makes
 * them unique to deal.II. As a consequence, it is possible that other
 * libraries your code interfaces with define the same name, and the result
 * will be name collisions (see
 * https://en.wikipedia.org/wiki/Name_collision). One can <code>\#undef</code>
 * this macro, as well as all other macros defined by deal.II that are not
 * prefixed with either <code>DEAL</code> or <code>deal</code>, by including
 * the header <code>deal.II/base/undefine_macros.h</code> after all other
 * deal.II headers have been included.
 *
 * @ingroup Exceptions
 */
#  define DeclException1(Exception1, type1, outsequence) \
    /** @ingroup Exceptions */                           \
    /** @dealiiExceptionMessage{outsequence} */          \
    static dealii::ExceptionBase &Exception1(type1 arg1)


/**
 * Declare an exception class derived from ExceptionBase with two additional
 * parameters.
 *
 * @note This and similar macro names are examples of preprocessor definitions
 * in the deal.II library that are not prefixed by a string that likely makes
 * them unique to deal.II. As a consequence, it is possible that other
 * libraries your code interfaces with define the same name, and the result
 * will be name collisions (see
 * https://en.wikipedia.org/wiki/Name_collision). One can <code>\#undef</code>
 * this macro, as well as all other macros defined by deal.II that are not
 * prefixed with either <code>DEAL</code> or <code>deal</code>, by including
 * the header <code>deal.II/base/undefine_macros.h</code> after all other
 * deal.II headers have been included.
 *
 * @ingroup Exceptions
 */
#  define DeclException2(Exception2, type1, type2, outsequence) \
    /** @ingroup Exceptions */                                  \
    /** @dealiiExceptionMessage{outsequence} */                 \
    static dealii::ExceptionBase &Exception2(type1 arg1, type2 arg2)


/**
 * Declare an exception class derived from ExceptionBase with three additional
 * parameters.
 *
 * @note This and similar macro names are examples of preprocessor definitions
 * in the deal.II library that are not prefixed by a string that likely makes
 * them unique to deal.II. As a consequence, it is possible that other
 * libraries your code interfaces with define the same name, and the result
 * will be name collisions (see
 * https://en.wikipedia.org/wiki/Name_collision). One can <code>\#undef</code>
 * this macro, as well as all other macros defined by deal.II that are not
 * prefixed with either <code>DEAL</code> or <code>deal</code>, by including
 * the header <code>deal.II/base/undefine_macros.h</code> after all other
 * deal.II headers have been included.
 *
 * @ingroup Exceptions
 */
#  define DeclException3(Exception3, type1, type2, type3, outsequence) \
    /** @ingroup Exceptions */                                         \
    /** @dealiiExceptionMessage{outsequence} */                        \
    static dealii::ExceptionBase &Exception3(type1 arg1, type2 arg2, type3 arg3)


/**
 * Declare an exception class derived from ExceptionBase with four additional
 * parameters.
 *
 * @note This and similar macro names are examples of preprocessor definitions
 * in the deal.II library that are not prefixed by a string that likely makes
 * them unique to deal.II. As a consequence, it is possible that other
 * libraries your code interfaces with define the same name, and the result
 * will be name collisions (see
 * https://en.wikipedia.org/wiki/Name_collision). One can <code>\#undef</code>
 * this macro, as well as all other macros defined by deal.II that are not
 * prefixed with either <code>DEAL</code> or <code>deal</code>, by including
 * the header <code>deal.II/base/undefine_macros.h</code> after all other
 * deal.II headers have been included.
 *
 * @ingroup Exceptions
 */
#  define DeclException4(Exception4, type1, type2, type3, type4, outsequence) \
    /** @ingroup Exceptions */                                                \
    /** @dealiiExceptionMessage{outsequence} */                               \
    static dealii::ExceptionBase &Exception4(type1 arg1,                      \
                                             type2 arg2,                      \
                                             type3 arg3,                      \
                                             type4 arg4)


/**
 * Declare an exception class derived from ExceptionBase with five additional
 * parameters.
 *
 * @note This and similar macro names are examples of preprocessor definitions
 * in the deal.II library that are not prefixed by a string that likely makes
 * them unique to deal.II. As a consequence, it is possible that other
 * libraries your code interfaces with define the same name, and the result
 * will be name collisions (see
 * https://en.wikipedia.org/wiki/Name_collision). One can <code>\#undef</code>
 * this macro, as well as all other macros defined by deal.II that are not
 * prefixed with either <code>DEAL</code> or <code>deal</code>, by including
 * the header <code>deal.II/base/undefine_macros.h</code> after all other
 * deal.II headers have been included.
 *
 * @ingroup Exceptions
 */
#  define DeclException5(                                       \
    Exception5, type1, type2, type3, type4, type5, outsequence) \
    /** @ingroup Exceptions */                                  \
    /** @dealiiExceptionMessage{outsequence} */                 \
    static dealii::ExceptionBase &Exception5(                   \
      type1 arg1, type2 arg2, type3 arg3, type4 arg4, type5 arg5)

#endif /*ifndef DOXYGEN*/



/**********************************************************************
 * Preprocessor definitions in support of assertions.
 *
 * These make reference to classes and functions that are declared in
 * exceptions.h. This is ok as far as concerning preprocessor defines
 * is concerned, but in order to use the material below, you also need
 * to have `#include <deal.II/base/exceptions.h>`.
 */

/**
 * Wrapper macro around __builtin_expect(). Only used in the assertion macros
 * (Assert(), AssertNothrow(), and AssertThrow()).
 */
#ifdef DEAL_II_HAVE_BUILTIN_EXPECT
#  define DEAL_II_BUILTIN_EXPECT(a, b) __builtin_expect((a), (b))
#else
#  define DEAL_II_BUILTIN_EXPECT(a, b) (a)
#endif


/**
 * A macro that serves as the main routine in the exception mechanism for debug
 * mode error checking. It asserts that a certain condition is fulfilled,
 * otherwise issues an error and aborts the program.
 *
 * A more detailed description can be found in the
 * @ref Exceptions
 * topic. It is first used in step-5 and step-6.
 * See also the <tt>ExceptionBase</tt> class for more information.
 *
 * @note Active in DEBUG mode only
 *
 * @note This and similar macro names are examples of preprocessor definitions
 * in the deal.II library that are not prefixed by a string that likely makes
 * them unique to deal.II. As a consequence, it is possible that other
 * libraries your code interfaces with define the same name, and the result
 * will be name collisions (see
 * https://en.wikipedia.org/wiki/Name_collision). One can <code>\#undef</code>
 * this macro, as well as all other macros defined by deal.II that are not
 * prefixed with either <code>DEAL</code> or <code>deal</code>, by including
 * the header <code>deal.II/base/undefine_macros.h</code> after all other
 * deal.II headers have been included.
 *
 * @ingroup Exceptions
 */
#ifdef DEBUG
#  if DEAL_II_KOKKOS_VERSION_GTE(3, 6, 0)
#    define Assert(cond, exc)                                                \
      do                                                                     \
        {                                                                    \
          KOKKOS_IF_ON_HOST(({                                               \
            if (DEAL_II_BUILTIN_EXPECT(!(cond), false))                      \
              ::dealii::deal_II_exceptions::internals::issue_error_noreturn( \
                ::dealii::deal_II_exceptions::internals::ExceptionHandling:: \
                  abort_or_throw_on_exception,                               \
                __FILE__,                                                    \
                __LINE__,                                                    \
                __PRETTY_FUNCTION__,                                         \
                #cond,                                                       \
                #exc,                                                        \
                exc);                                                        \
          }))                                                                \
          KOKKOS_IF_ON_DEVICE(({                                             \
            if (!(cond))                                                     \
              Kokkos::abort(#cond);                                          \
          }))                                                                \
        }                                                                    \
      while (false)
#  else /*if DEAL_II_KOKKOS_VERSION_GTE(3,6,0)*/
#    ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
#      define Assert(cond, exc)                                              \
        do                                                                   \
          {                                                                  \
            if (DEAL_II_BUILTIN_EXPECT(!(cond), false))                      \
              ::dealii::deal_II_exceptions::internals::issue_error_noreturn( \
                ::dealii::deal_II_exceptions::internals::ExceptionHandling:: \
                  abort_or_throw_on_exception,                               \
                __FILE__,                                                    \
                __LINE__,                                                    \
                __PRETTY_FUNCTION__,                                         \
                #cond,                                                       \
                #exc,                                                        \
                exc);                                                        \
          }                                                                  \
        while (false)
#    else /*#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST*/
#      define Assert(cond, exc)     \
        do                          \
          {                         \
            if (!(cond))            \
              Kokkos::abort(#cond); \
          }                         \
        while (false)
#    endif /*ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST*/
#  endif   /*KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST*/
#else      /*ifdef DEBUG*/
#  define Assert(cond, exc) \
    do                      \
      {                     \
        if (false)          \
          if (!(cond))      \
            {               \
            }               \
      }                     \
    while (false)
#endif /*ifdef DEBUG*/



/**
 * A variant of the <tt>Assert</tt> macro above that exhibits the same runtime
 * behavior as long as disable_abort_on_exception was not called.
 *
 * However, if disable_abort_on_exception was called, this macro merely prints
 * the exception that would be thrown to deallog and continues normally
 * without throwing an exception.
 *
 * A more detailed description can be found in the
 * @ref Exceptions
 * topic, in the discussion about the corner case at the bottom of the page.
 *
 * @note This and similar macro names are examples of preprocessor definitions
 * in the deal.II library that are not prefixed by a string that likely makes
 * them unique to deal.II. As a consequence, it is possible that other
 * libraries your code interfaces with define the same name, and the result
 * will be name collisions (see
 * https://en.wikipedia.org/wiki/Name_collision). One can <code>\#undef</code>
 * this macro, as well as all other macros defined by deal.II that are not
 * prefixed with either <code>DEAL</code> or <code>deal</code>, by including
 * the header <code>deal.II/base/undefine_macros.h</code> after all other
 * deal.II headers have been included.
 *
 * @note Active in DEBUG mode only
 * @ingroup Exceptions
 */
#ifdef DEBUG
#  define AssertNothrow(cond, exc)                                      \
    do                                                                  \
      {                                                                 \
        if (DEAL_II_BUILTIN_EXPECT(!(cond), false))                     \
          ::dealii::deal_II_exceptions::internals::issue_error_nothrow( \
            __FILE__, __LINE__, __PRETTY_FUNCTION__, #cond, #exc, exc); \
      }                                                                 \
    while (false)
#else
#  define AssertNothrow(cond, exc) \
    do                             \
      {                            \
        if (false)                 \
          if (!(cond))             \
            {                      \
            }                      \
      }                            \
    while (false)
#endif

/**
 * A macro that serves as the main routine in the exception mechanism for
 * dynamic error checking. It asserts that a certain condition is fulfilled,
 * otherwise
 * throws an exception via the C++ @p throw mechanism. This exception can
 * be caught via a @p catch clause, as is shown in step-6 and all following
 * tutorial programs.
 *
 * A more detailed description can be found in the
 * @ref Exceptions
 * topic. It is first used in step-9 and step-13.
 * See also the <tt>ExceptionBase</tt> class for more information.
 *
 * @note This and similar macro names are examples of preprocessor definitions
 * in the deal.II library that are not prefixed by a string that likely makes
 * them unique to deal.II. As a consequence, it is possible that other
 * libraries your code interfaces with define the same name, and the result
 * will be name collisions (see
 * https://en.wikipedia.org/wiki/Name_collision). One can <code>\#undef</code>
 * this macro, as well as all other macros defined by deal.II that are not
 * prefixed with either <code>DEAL</code> or <code>deal</code>, by including
 * the header <code>deal.II/base/undefine_macros.h</code> after all other
 * deal.II headers have been included.
 *
 * @note Active in both DEBUG and RELEASE modes
 * @ingroup Exceptions
 */
#define AssertThrow(cond, exc)                                         \
  do                                                                   \
    {                                                                  \
      if (DEAL_II_BUILTIN_EXPECT(!(cond), false))                      \
        ::dealii::deal_II_exceptions::internals::issue_error_noreturn( \
          ::dealii::deal_II_exceptions::internals::ExceptionHandling:: \
            throw_on_exception,                                        \
          __FILE__,                                                    \
          __LINE__,                                                    \
          __PRETTY_FUNCTION__,                                         \
          #cond,                                                       \
          #exc,                                                        \
          exc);                                                        \
    }                                                                  \
  while (false)


/**
 * `DEAL_II_NOT_IMPLEMENTED` is a macro (that looks like a function call
 * when used as in `DEAL_II_NOT_IMPLEMENTED();`) that is used to raise an
 * error in places where a piece of code is not yet implemented. If a code
 * runs into such a place, it will be aborted with an error message that
 * explains the situation, along with a backtrace of how the code ended
 * up in this place. Alternatively, if
 * deal_II_exceptions::internals::ExceptionHandling::abort_or_throw_on_exception
 * is set to ExceptionHandling::throw_on_exception, then the corresponding
 * error is thrown as a C++ exception that can be caught (though in
 * many cases codes will then find it difficult to do what they wanted
 * to do).
 *
 * This macro is first used in step-8 of the tutorial.
 *
 * A typical case where it is used would look as follows: Assume that we want
 * to implement a function that describes the right hand side of an equation
 * that corresponds to a known solution (i.e., we want to use the "Method
 * of manufactured solutions", see step-7). We have computed the right
 * hand side that corresponds to the 1d and 2d solutions, but we have been
 * too lazy so far to do the calculations for the 3d case, perhaps because
 * we first want to test correctness in 1d and 2d before moving on to the 3d
 * case. We could then write this right hand side as follows (the specific
 * formulas in the `return` statements are not important):
 * @code
 *   template <int dim>
 *   double right_hand_side (const Point<dim> &x)
 *   {
 *     if (dim==1)
 *       return x[0]*std::sin(x[0]);
 *     else if (dim==2)
 *       return x[0]*std::sin(x[0])*std::sin(x[1];
 *     else
 *       DEAL_II_NOT_IMPLEMENTED();
 *   }
 * @endcode
 * Here, the call to `DEAL_II_NOT_IMPLEMENTED()` simply indicates that we
 * haven't gotten around to filling in this code block. If someone ends up
 * running the program in 3d, execution will abort in the location with an
 * error message that indicates where this happened and why.
 */
#define DEAL_II_NOT_IMPLEMENTED()                                \
  ::dealii::deal_II_exceptions::internals::issue_error_noreturn( \
    ::dealii::deal_II_exceptions::internals::ExceptionHandling:: \
      abort_or_throw_on_exception,                               \
    __FILE__,                                                    \
    __LINE__,                                                    \
    __PRETTY_FUNCTION__,                                         \
    nullptr,                                                     \
    nullptr,                                                     \
    ::dealii::StandardExceptions::ExcNotImplemented())


/**
 * `DEAL_II_ASSERT_UNREACHABLE` is a macro (that looks like a function call
 * when used as in `DEAL_II_ASSERT_UNREACHABLE();`) that is used to raise an
 * error in places where the programmer believed that execution should
 * never get to. If a code
 * runs into such a place, it will be aborted with an error message that
 * explains the situation, along with a backtrace of how the code ended
 * up in this place. Alternatively, if
 * deal_II_exceptions::internals::ExceptionHandling::abort_or_throw_on_exception
 * is set to ExceptionHandling::throw_on_exception, then the corresponding
 * error is thrown as a C++ exception that can be caught (though in
 * many cases codes will then find it difficult to do what they wanted
 * to do).
 *
 * A typical case where it is used would look as follows. In many cases,
 * one has a finite enumeration of things that can happen, and one runs
 * through those in a sequence of `if`-`else` blocks, or perhaps
 * with a `switch` selection and a number of `case` statements. Of
 * course, if the code is correct, if all possible cases are handled,
 * nothing terrible can happen -- though perhaps it is worth making sure
 * that we have really covered all cases by using `DEAL_II_ASSERT_UNREACHABLE()`
 * as the *last* case. Here is an example:
 * @code
 *   enum OutputFormat { vtk, vtu };
 *
 *   void write_output (const OutputFormat format)
 *   {
 *     if (format == vtk)
 *       {
 *         ... write in VTK format ...
 *       }
 *     else // must not clearly be VTU format
 *       {
 *         ... write in VTU format ...
 *       }
 *   }
 * @endcode
 * The issue here is "Are we really sure it is VTU format if we end up in
 * the `else` block"? There are two reasons that should make us suspicious.
 * First, the authors of the code may later have expanded the list of options
 * in the `OutputFormat` enum, but forgotten to also update the
 * `write_output()` function. We may then end up in the `else` branch even
 * though the argument indicates the now possible third option that was added
 * to `OutputFormat`. The second possibility to consider is that enums are
 * really just fancy ways of using integers; from a language perspective, it
 * is allowed to pass *any* integer to `write_output()`, even values that do
 * not match either `vtk` or `vtu`. This is then clearly a bug in the program,
 * but one that we are better off if we catch it as early as possible.
 *
 * We can guard against both cases by writing the code as follows instead:
 * @code
 *   enum OutputFormat { vtk, vtu };
 *
 *   void write_output (const OutputFormat format)
 *   {
 *     if (format == vtk)
 *       {
 *         ... write in VTK format ...
 *       }
 *     else if (format == vtu)
 *       {
 *         ... write in VTU format ...
 *       }
 *     else // we shouldn't get here, but if we did, abort the program now!
 *       DEAL_II_ASSERT_UNREACHABLE();
 *   }
 * @endcode
 *
 * This macro is first used in step-7, where we show another example of
 * a context where it is frequently used.
 */
#define DEAL_II_ASSERT_UNREACHABLE()                             \
  ::dealii::deal_II_exceptions::internals::issue_error_noreturn( \
    ::dealii::deal_II_exceptions::internals::ExceptionHandling:: \
      abort_or_throw_on_exception,                               \
    __FILE__,                                                    \
    __LINE__,                                                    \
    __PRETTY_FUNCTION__,                                         \
    nullptr,                                                     \
    nullptr,                                                     \
    ::dealii::StandardExceptions::ExcMessage(                    \
      "The program has hit a line of code that the programmer "  \
      "marked with the macro DEAL_II_ASSERT_UNREACHABLE() to "   \
      "indicate that the program should never reach this "       \
      "location. You will have to find out (best done in a "     \
      "debugger) why that happened. Typical reasons include "    \
      "passing invalid arguments to functions (for example, if " \
      "a function takes an 'enum' with two possible values "     \
      "as argument, but you call the function with a third "     \
      "value), or if the programmer of the code that triggered " \
      "the error believed that a variable can only have "        \
      "specific values, but either that assumption is wrong "    \
      "or the computation of that value is buggy."               \
      "\n\n"                                                     \
      "In those latter conditions, where some internal "         \
      "assumption is not satisfied, there may not be very "      \
      "much you can do if you encounter such an exception, "     \
      "since it indicates an error in deal.II, not in your "     \
      "own program. If that is the situation you encounter, "    \
      "try to come up with "                                     \
      "the smallest possible program that still demonstrates "   \
      "the error and contact the deal.II mailing lists with it " \
      "to obtain help."))


/**
 * Special assertion for dimension mismatch.
 *
 * Since this is used very often and always repeats the arguments, we
 * introduce this special assertion for ExcDimensionMismatch in order to keep
 * the user codes shorter.
 *
 * @note This and similar macro names are examples of preprocessor definitions
 * in the deal.II library that are not prefixed by a string that likely makes
 * them unique to deal.II. As a consequence, it is possible that other
 * libraries your code interfaces with define the same name, and the result
 * will be name collisions (see
 * https://en.wikipedia.org/wiki/Name_collision). One can <code>\#undef</code>
 * this macro, as well as all other macros defined by deal.II that are not
 * prefixed with either <code>DEAL</code> or <code>deal</code>, by including
 * the header <code>deal.II/base/undefine_macros.h</code> after all other
 * deal.II headers have been included.
 *
 * @ingroup Exceptions
 */
#define AssertDimension(dim1, dim2)                                           \
  Assert(::dealii::deal_II_exceptions::internals::compare_for_equality(dim1,  \
                                                                       dim2), \
         dealii::ExcDimensionMismatch((dim1), (dim2)))

/**
 * Special assertion for integer conversions.
 *
 * deal.II does not always use the same integer types as its dependencies. For
 * example, PETSc uses signed integers whereas deal.II uses unsigned integers.
 * This assertion checks that we can successfully convert between two index
 * types.
 */
#define AssertIntegerConversion(index1, index2)                         \
  Assert(::dealii::deal_II_exceptions::internals::compare_for_equality( \
           index1, index2),                                             \
         dealii::ExcInvalidIntegerConversion((index1), (index2)))

/**
 * Special assertion for integer conversions which will throw exceptions.
 * Otherwise this is the same as AssertIntegerConversion.
 */
#define AssertThrowIntegerConversion(index1, index2)                         \
  AssertThrow(::dealii::deal_II_exceptions::internals::compare_for_equality( \
                index1, index2),                                             \
              dealii::ExcInvalidIntegerConversion((index1), (index2)))

/**
 * An assertion that tests whether <tt>vec</tt> has size <tt>dim1</tt>, and
 * each entry of the vector is itself an array that has the size <tt>dim2</tt>.
 *
 * @note This and similar macro names are examples of preprocessor definitions
 * in the deal.II library that are not prefixed by a string that likely makes
 * them unique to deal.II. As a consequence, it is possible that other
 * libraries your code interfaces with define the same name, and the result
 * will be name collisions (see
 * https://en.wikipedia.org/wiki/Name_collision). One can <code>\#undef</code>
 * this macro, as well as all other macros defined by deal.II that are not
 * prefixed with either <code>DEAL</code> or <code>deal</code>, by including
 * the header <code>deal.II/base/undefine_macros.h</code> after all other
 * deal.II headers have been included.
 *
 * @ingroup Exceptions
 */
#define AssertVectorVectorDimension(VEC, DIM1, DIM2) \
  AssertDimension(VEC.size(), DIM1);                 \
  for (const auto &subvector : VEC)                  \
    {                                                \
      (void)subvector;                               \
      AssertDimension(subvector.size(), DIM2);       \
    }


/**
 * An assertion that tests that a given index is within the half-open
 * range <code>[0,range)</code>. It throws an exception object
 * <code>ExcIndexRange(index,0,range)</code> if the assertion
 * fails.
 *
 * @note This and similar macro names are examples of preprocessor definitions
 * in the deal.II library that are not prefixed by a string that likely makes
 * them unique to deal.II. As a consequence, it is possible that other
 * libraries your code interfaces with define the same name, and the result
 * will be name collisions (see
 * https://en.wikipedia.org/wiki/Name_collision). One can <code>\#undef</code>
 * this macro, as well as all other macros defined by deal.II that are not
 * prefixed with either <code>DEAL</code> or <code>deal</code>, by including
 * the header <code>deal.II/base/undefine_macros.h</code> after all other
 * deal.II headers have been included.
 *
 * @ingroup Exceptions
 */
#define AssertIndexRange(index, range)                                         \
  Assert(::dealii::deal_II_exceptions::internals::compare_less_than(index,     \
                                                                    range),    \
         ::dealii::ExcIndexRangeType<::dealii::internal::argument_type_t<void( \
           std::common_type_t<decltype(index), decltype(range)>)>>((index),    \
                                                                   0,          \
                                                                   (range)))

/**
 * An assertion that checks whether a number is finite or not. We explicitly
 * cast the number to std::complex to match the signature of the exception
 * (see there for an explanation of why we use std::complex at all) and to
 * satisfy the fact that std::complex has no implicit conversions.
 *
 * @note This and similar macro names are examples of preprocessor definitions
 * in the deal.II library that are not prefixed by a string that likely makes
 * them unique to deal.II. As a consequence, it is possible that other
 * libraries your code interfaces with define the same name, and the result
 * will be name collisions (see
 * https://en.wikipedia.org/wiki/Name_collision). One can <code>\#undef</code>
 * this macro, as well as all other macros defined by deal.II that are not
 * prefixed with either <code>DEAL</code> or <code>deal</code>, by including
 * the header <code>deal.II/base/undefine_macros.h</code> after all other
 * deal.II headers have been included.
 *
 * @ingroup Exceptions
 */
#define AssertIsFinite(number)               \
  Assert(dealii::numbers::is_finite(number), \
         dealii::ExcNumberNotFinite(std::complex<double>(number)))

/**
 * Assert that a geometric object is not used. This assertion is used when
 * constructing triangulations and should normally not be used inside user
 * codes.
 */
#define AssertIsNotUsed(obj) Assert((obj)->used() == false, ExcInternalError())

#ifdef DEAL_II_WITH_MPI
/**
 * An assertion that checks whether or not an error code returned by an MPI
 * function is equal to <code>MPI_SUCCESS</code>. If the check fails then an
 * exception of type ExcMPI is thrown with the given error code as an
 * argument.
 *
 * @note This and similar macro names are examples of preprocessor definitions
 * in the deal.II library that are not prefixed by a string that likely makes
 * them unique to deal.II. As a consequence, it is possible that other
 * libraries your code interfaces with define the same name, and the result
 * will be name collisions (see
 * https://en.wikipedia.org/wiki/Name_collision). One can <code>\#undef</code>
 * this macro, as well as all other macros defined by deal.II that are not
 * prefixed with either <code>DEAL</code> or <code>deal</code>, by including
 * the header <code>deal.II/base/undefine_macros.h</code> after all other
 * deal.II headers have been included.
 *
 * @note Active only if deal.II is compiled with MPI
 * @ingroup Exceptions
 */
#  define AssertThrowMPI(error_code) \
    AssertThrow(error_code == MPI_SUCCESS, dealii::ExcMPI(error_code))
#else
#  define AssertThrowMPI(error_code) \
    do                               \
      {                              \
      }                              \
    while (false)
#endif // DEAL_II_WITH_MPI

#ifdef DEAL_II_TRILINOS_WITH_SEACAS
/**
 * Assertion that checks that the error code produced by calling an ExodusII
 * routine is equal to EX_NOERR (which is zero).
 *
 * @note This and similar macro names are examples of preprocessor definitions
 * in the deal.II library that are not prefixed by a string that likely makes
 * them unique to deal.II. As a consequence, it is possible that other
 * libraries your code interfaces with define the same name, and the result
 * will be name collisions (see
 * https://en.wikipedia.org/wiki/Name_collision). One can <code>\#undef</code>
 * this macro, as well as all other macros defined by deal.II that are not
 * prefixed with either <code>DEAL</code> or <code>deal</code>, by including
 * the header <code>deal.II/base/undefine_macros.h</code> after all other
 * deal.II headers have been included.
 *
 * @ingroup Exceptions
 */
#  define AssertThrowExodusII(error_code) \
    AssertThrow(error_code == 0,          \
                dealii::StandardExceptions::ExcExodusII(error_code));
#endif // DEAL_II_TRILINOS_WITH_SEACAS


#ifdef DEAL_II_WITH_SUNDIALS
/**
 * Assertion that checks that the error code produced by calling one
 * of SUNDIALS' ARKode functions.
 *
 * @note This and similar macro names are examples of preprocessor definitions
 * in the deal.II library that are not prefixed by a string that likely makes
 * them unique to deal.II. As a consequence, it is possible that other
 * libraries your code interfaces with define the same name, and the result
 * will be name collisions (see
 * https://en.wikipedia.org/wiki/Name_collision). One can <code>\#undef</code>
 * this macro, as well as all other macros defined by deal.II that are not
 * prefixed with either <code>DEAL</code> or <code>deal</code>, by including
 * the header <code>deal.II/base/undefine_macros.h</code> after all other
 * deal.II headers have been included.
 *
 * @ingroup Exceptions
 */
#  define AssertARKode(code) Assert(code >= 0, ExcARKodeError(code))

/**
 * Assertion that checks that the error code produced by calling one
 * of SUNDIALS' KINSOL functions.
 *
 * @note This and similar macro names are examples of preprocessor definitions
 * in the deal.II library that are not prefixed by a string that likely makes
 * them unique to deal.II. As a consequence, it is possible that other
 * libraries your code interfaces with define the same name, and the result
 * will be name collisions (see
 * https://en.wikipedia.org/wiki/Name_collision). One can <code>\#undef</code>
 * this macro, as well as all other macros defined by deal.II that are not
 * prefixed with either <code>DEAL</code> or <code>deal</code>, by including
 * the header <code>deal.II/base/undefine_macros.h</code> after all other
 * deal.II headers have been included.
 *
 * @ingroup Exceptions
 */
#  define AssertKINSOL(code) Assert(code >= 0, ExcKINSOLError(code))

/**
 * Assertion that checks that the error code produced by calling one
 * of SUNDIALS' IDA functions.
 *
 * @note This and similar macro names are examples of preprocessor definitions
 * in the deal.II library that are not prefixed by a string that likely makes
 * them unique to deal.II. As a consequence, it is possible that other
 * libraries your code interfaces with define the same name, and the result
 * will be name collisions (see
 * https://en.wikipedia.org/wiki/Name_collision). One can <code>\#undef</code>
 * this macro, as well as all other macros defined by deal.II that are not
 * prefixed with either <code>DEAL</code> or <code>deal</code>, by including
 * the header <code>deal.II/base/undefine_macros.h</code> after all other
 * deal.II headers have been included.
 *
 * @ingroup Exceptions
 */
#  define AssertIDA(code) Assert(code >= 0, ExcIDAError(code))
#endif

#endif
