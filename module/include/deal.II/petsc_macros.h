// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// PETSc defines a couple of function-like macros that are mixed in
// with other stuff in their header files so that in order to get to
// these macros, we would have to include the whole file. This of
// course defeats the purpose of wrapping things into module
// partitions.  Rather -- perhaps imprudently -- we repeat these
// macros here.

#ifndef dealii_petsc_macros_h
#define dealii_petsc_macros_h

#ifdef DEAL_II_WITH_PETSC

#  include <petscmacros.h>

// Taken from petscerror.h:
#  define PetscCall(...)                                        \
    do                                                          \
      {                                                         \
        PetscErrorCode ierr_petsc_call_q_;                      \
        PetscStackUpdateLine;                                   \
        ierr_petsc_call_q_ = __VA_ARGS__;                       \
        if (PetscUnlikely(ierr_petsc_call_q_ != PETSC_SUCCESS)) \
          return PetscError(PETSC_COMM_SELF,                    \
                            __LINE__,                           \
                            PETSC_FUNCTION_NAME,                \
                            __FILE__,                           \
                            ierr_petsc_call_q_,                 \
                            PETSC_ERROR_REPEAT,                 \
                            " ");                               \
      }                                                         \
    while (0)


#  define PetscRegister__FUNCT__()

#  define PetscFunctionBeginUser                                    \
    do                                                              \
      {                                                             \
        PetscStackPushNoCheck(PETSC_FUNCTION_NAME, 2, PETSC_FALSE); \
        PetscRegister__FUNCT__();                                   \
      }                                                             \
    while (0)

#  define PetscFunctionReturn(...)                 \
    do                                             \
      {                                            \
        PetscStackPopNoCheck(PETSC_FUNCTION_NAME); \
        return __VA_ARGS__;                        \
      }                                            \
    while (0)

#  define PetscObjectStateIncrease(obj) ((obj)->state++, PETSC_SUCCESS)

#  if defined(PETSC_CLANG_STATIC_ANALYZER) || defined(__clang_analyzer__)
#    define PetscStackUpdateLine
#    define PetscStackPushNoCheck(funct, petsc_routine, hot)
#    define PetscStackPopNoCheck

#  elif defined(PETSC_USE_DEBUG) && !defined(PETSC_HAVE_THREADSAFETY)

#    define PetscStackUpdateLine                                      \
      do                                                              \
        {                                                             \
          if (petscstack.currentsize > 0 &&                           \
              petscstack.function[petscstack.currentsize - 1] ==      \
                PETSC_FUNCTION_NAME)                                  \
            {                                                         \
              petscstack.line[petscstack.currentsize - 1] = __LINE__; \
            }                                                         \
        }                                                             \
      while (0)

#    define PetscStackPushNoCheck(funct, petsc_routine, hot)            \
      do                                                                \
        {                                                               \
          PetscStackSAWsTakeAccess();                                   \
          PetscStackPush_Private(                                       \
            petscstack, __FILE__, funct, __LINE__, petsc_routine, hot); \
          PetscStackSAWsGrantAccess();                                  \
        }                                                               \
      while (0)
#    define PetscStackPopNoCheck(funct)             \
      do                                            \
        {                                           \
          PetscStackSAWsTakeAccess();               \
          PetscStackPop_Private(petscstack, funct); \
          PetscStackSAWsGrantAccess();              \
        }                                           \
      while (0)

#  else /* PETSC_USE_DEBUG */
#    define PetscStackUpdateLine
#    define PetscStackPushNoCheck(funct, petsc_routine, hot)
#    define PetscStackPopNoCheck(funct)
#  endif

#  define SETERRQ(comm, ierr, ...)                                           \
    do                                                                       \
      {                                                                      \
        PetscErrorCode ierr_seterrq_petsc_ = PetscError(comm,                \
                                                        __LINE__,            \
                                                        PETSC_FUNCTION_NAME, \
                                                        __FILE__,            \
                                                        ierr,                \
                                                        PETSC_ERROR_INITIAL, \
                                                        __VA_ARGS__);        \
        return ierr_seterrq_petsc_ ? ierr_seterrq_petsc_ : PETSC_ERR_RETURN; \
      }                                                                      \
    while (0)



// Taken from petscsys.h:
#  define PetscObjectQueryFunction(obj, name, fptr) \
    PetscObjectQueryFunction_Private((obj), (name), (PetscVoidFunction *)(fptr))
#  define PetscObjectComposeFunction(a, b, ...) \
    PetscObjectComposeFunction_Private((a),     \
                                       (b),     \
                                       (PetscVoidFunction)(__VA_ARGS__))

#  define PetscMalloc(a, b)                \
    ((*PetscTrMalloc)((a),                 \
                      PETSC_FALSE,         \
                      __LINE__,            \
                      PETSC_FUNCTION_NAME, \
                      __FILE__,            \
                      (void **)(b)))
#  define PetscMalloc1(m1, r1)                            \
    PetscMallocA(1,                                       \
                 PETSC_FALSE,                             \
                 __LINE__,                                \
                 PETSC_FUNCTION_NAME,                     \
                 __FILE__,                                \
                 ((size_t)((size_t)m1) * sizeof(**(r1))), \
                 (r1))

#  if !defined(PETSC_HAVE_SAWS)
#    define PetscObjectSAWsGrantAccess(obj) PETSC_SUCCESS
#    define PetscObjectSAWsTakeAccess(obj) PETSC_SUCCESS
#  endif

#endif

#endif
