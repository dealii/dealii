/*
    Copyright 2005-2013 Intel Corporation.  All Rights Reserved.

    This file is part of Threading Building Blocks.

    Threading Building Blocks is free software; you can redistribute it
    and/or modify it under the terms of the GNU General Public License
    version 2 as published by the Free Software Foundation.

    Threading Building Blocks is distributed in the hope that it will be
    useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Threading Building Blocks; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

    As a special exception, you may use this file as part of a free software
    library without restriction.  Specifically, if other files instantiate
    templates or use macros or inline functions from this file, or you compile
    this file and link it with other files to produce an executable, this
    file does not by itself cause the resulting executable to be covered by
    the GNU General Public License.  This exception does not however
    invalidate any other reasons why the executable file might be covered by
    the GNU General Public License.
*/

#ifndef __TBB_tbb_config_H
#define __TBB_tbb_config_H

/** This header is supposed to contain macro definitions and C style comments only.
    The macros defined here are intended to control such aspects of TBB build as
    - presence of compiler features
    - compilation modes
    - feature sets
    - known compiler/platform issues
**/

#define __TBB_GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)

#if __clang__
    /**according to clang documentation version can be vendor specific **/
    #define __TBB_CLANG_VERSION (__clang_major__ * 10000 + __clang_minor__ * 100 + __clang_patchlevel__)
#endif

/** Presence of compiler features **/

#if __INTEL_COMPILER == 9999 && __INTEL_COMPILER_BUILD_DATE == 20110811
/* Intel(R) Composer XE 2011 Update 6 incorrectly sets __INTEL_COMPILER. Fix it. */
    #undef __INTEL_COMPILER
    #define __INTEL_COMPILER 1210
#endif

#if (__TBB_GCC_VERSION >= 40400) && !defined(__INTEL_COMPILER)
    /** warning suppression pragmas available in GCC since 4.4 **/
    #define __TBB_GCC_WARNING_SUPPRESSION_PRESENT 1
#endif

/* Select particular features of C++11 based on compiler version.
   ICC 12.1 (Linux), GCC 4.3 and higher, clang 2.9 and higher
   set __GXX_EXPERIMENTAL_CXX0X__ in c++11 mode.

   Compilers that mimics other compilers (ICC, clang) must be processed before
   compilers they mimic (GCC, MSVC).

   TODO: The following conditions should be extended when new compilers/runtimes
   support added.
 */

#if __INTEL_COMPILER
    /** On Windows environment when using Intel C++ compiler with Visual Studio 2010*,
        the C++0x features supported by Visual C++ 2010 are enabled by default
        TODO: find a way to get know if c++0x mode is specified in command line on windows **/
    #define __TBB_CPP11_VARIADIC_TEMPLATES_PRESENT    ( __GXX_EXPERIMENTAL_CXX0X__ && __VARIADIC_TEMPLATES )
    #define __TBB_CPP11_RVALUE_REF_PRESENT            ( (__GXX_EXPERIMENTAL_CXX0X__ || _MSC_VER >= 1600) && (__INTEL_COMPILER >= 1200) )
    #if  _MSC_VER >= 1600
        #define __TBB_EXCEPTION_PTR_PRESENT           ( __INTEL_COMPILER > 1300                                                \
                                                      /*ICC 12.1 Upd 10 and 13 beta Upd 2 fixed exception_ptr linking  issue*/ \
                                                      || (__INTEL_COMPILER == 1300 && __INTEL_COMPILER_BUILD_DATE >= 20120530) \
                                                      || (__INTEL_COMPILER == 1210 && __INTEL_COMPILER_BUILD_DATE >= 20120410) )
    /** libstc++ that comes with GCC 4.6 use C++11 features not supported by ICC 12.1.
     * Because of that ICC 12.1 does not support C++11 mode with with gcc 4.6. (or higher)
     * , and therefore does not  define __GXX_EXPERIMENTAL_CXX0X__ macro**/
    #elif (__TBB_GCC_VERSION >= 40404) && (__TBB_GCC_VERSION < 40600)
        #define __TBB_EXCEPTION_PTR_PRESENT        ( __GXX_EXPERIMENTAL_CXX0X__ && __INTEL_COMPILER >= 1200 )
    #elif (__TBB_GCC_VERSION >= 40600)
        #define __TBB_EXCEPTION_PTR_PRESENT        ( __GXX_EXPERIMENTAL_CXX0X__ && __INTEL_COMPILER >= 1300 )
    #else
        #define __TBB_EXCEPTION_PTR_PRESENT           0
    #endif
    #define __TBB_MAKE_EXCEPTION_PTR_PRESENT          (_MSC_VER >= 1700 || (__GXX_EXPERIMENTAL_CXX0X__ && __TBB_GCC_VERSION >= 40600))
    #define __TBB_STATIC_ASSERT_PRESENT               ( __GXX_EXPERIMENTAL_CXX0X__ || (_MSC_VER >= 1600) )
    #define __TBB_CPP11_TUPLE_PRESENT                 ( (_MSC_VER >= 1600) || ((__GXX_EXPERIMENTAL_CXX0X__) && (__TBB_GCC_VERSION >= 40300)) )
    /** TODO: re-check for compiler version greater than 12.1 if it supports initializer lists**/
    #define __TBB_INITIALIZER_LISTS_PRESENT           0
    #define __TBB_CONSTEXPR_PRESENT                   0
    #define __TBB_DEFAULTED_AND_DELETED_FUNC_PRESENT  0
#elif __clang__
//TODO: these options need to be rechecked
/** on OS X* the only way to get C++11 is to use clang. For library features (e.g. exception_ptr) libc++ is also
 *  required. So there is no need to check GCC version for clang**/
    #define __TBB_CPP11_VARIADIC_TEMPLATES_PRESENT     __has_feature(__cxx_variadic_templates__)
    #define __TBB_CPP11_RVALUE_REF_PRESENT             __has_feature(__cxx_rvalue_references__)
    #define __TBB_EXCEPTION_PTR_PRESENT               (__GXX_EXPERIMENTAL_CXX0X__ && (__cplusplus >= 201103L))
    #define __TBB_MAKE_EXCEPTION_PTR_PRESENT          (__GXX_EXPERIMENTAL_CXX0X__ && (__cplusplus >= 201103L))
    #define __TBB_STATIC_ASSERT_PRESENT               __has_feature(__cxx_static_assert__)
    /**Clang (preprocessor) has problems with dealing with expression having __has_include in #if's
     * used inside C++ code. (At least version that comes with OS X 10.8) **/
    #if (__GXX_EXPERIMENTAL_CXX0X__ && __has_include(<tuple>))
        #define __TBB_CPP11_TUPLE_PRESENT             1
    #endif
    #if (__has_feature(__cxx_generalized_initializers__) && __has_include(<initializer_list>))
        #define __TBB_INITIALIZER_LISTS_PRESENT       1
    #endif
    #define __TBB_CONSTEXPR_PRESENT                   __has_feature(__cxx_constexpr__)
    #define __TBB_DEFAULTED_AND_DELETED_FUNC_PRESENT  (__has_feature(__cxx_defaulted_functions__) && __has_feature(__cxx_deleted_functions__))
#elif __GNUC__
    #define __TBB_CPP11_VARIADIC_TEMPLATES_PRESENT    __GXX_EXPERIMENTAL_CXX0X__
    #define __TBB_CPP11_RVALUE_REF_PRESENT            __GXX_EXPERIMENTAL_CXX0X__
    /** __GCC_HAVE_SYNC_COMPARE_AND_SWAP_4 here is a substitution for _GLIBCXX_ATOMIC_BUILTINS_4, which is a prerequisite 
        for exception_ptr but cannot be used in this file because it is defined in a header, not by the compiler. 
        If the compiler has no atomic intrinsics, the C++ library should not expect those as well. **/
    #define __TBB_EXCEPTION_PTR_PRESENT               ((__GXX_EXPERIMENTAL_CXX0X__) && (__TBB_GCC_VERSION >= 40404) && __GCC_HAVE_SYNC_COMPARE_AND_SWAP_4)
    #define __TBB_MAKE_EXCEPTION_PTR_PRESENT          ((__GXX_EXPERIMENTAL_CXX0X__) && (__TBB_GCC_VERSION >= 40600))
    #define __TBB_STATIC_ASSERT_PRESENT               ((__GXX_EXPERIMENTAL_CXX0X__) && (__TBB_GCC_VERSION >= 40300))
    #define __TBB_CPP11_TUPLE_PRESENT                 ((__GXX_EXPERIMENTAL_CXX0X__) && (__TBB_GCC_VERSION >= 40300))
    #define __TBB_INITIALIZER_LISTS_PRESENT           ((__GXX_EXPERIMENTAL_CXX0X__) && (__TBB_GCC_VERSION >= 40400))
    /** gcc seems have to support constexpr from 4.4 but tests in (test_atomic) seeming reasonable fail to compile prior 4.6**/
    #define __TBB_CONSTEXPR_PRESENT                   ((__GXX_EXPERIMENTAL_CXX0X__) && (__TBB_GCC_VERSION >= 40400))
    #define __TBB_DEFAULTED_AND_DELETED_FUNC_PRESENT  ((__GXX_EXPERIMENTAL_CXX0X__) && (__TBB_GCC_VERSION >= 40400))
#elif _MSC_VER
    #define __TBB_CPP11_VARIADIC_TEMPLATES_PRESENT    0
    #define __TBB_CPP11_RVALUE_REF_PRESENT            0
    #define __TBB_EXCEPTION_PTR_PRESENT               (_MSC_VER >= 1600)
    #define __TBB_STATIC_ASSERT_PRESENT               (_MSC_VER >= 1600)
    #define __TBB_MAKE_EXCEPTION_PTR_PRESENT          (_MSC_VER >= 1700)
    #define __TBB_CPP11_TUPLE_PRESENT                 (_MSC_VER >= 1600)
    #define __TBB_INITIALIZER_LISTS_PRESENT           0
    #define __TBB_CONSTEXPR_PRESENT                   0
    #define __TBB_DEFAULTED_AND_DELETED_FUNC_PRESENT  0
#else
    #define __TBB_CPP11_VARIADIC_TEMPLATES_PRESENT    0
    #define __TBB_CPP11_RVALUE_REF_PRESENT            0
    #define __TBB_EXCEPTION_PTR_PRESENT               0
    #define __TBB_STATIC_ASSERT_PRESENT               0
    #define __TBB_MAKE_EXCEPTION_PTR_PRESENT          0
    #define __TBB_CPP11_TUPLE_PRESENT                 0
    #define __TBB_INITIALIZER_LISTS_PRESENT           0
    #define __TBB_CONSTEXPR_PRESENT                   0
    #define __TBB_DEFAULTED_AND_DELETED_FUNC_PRESENT  0
#endif

//TODO: not clear how exactly this macro affects exception_ptr - investigate
// On linux ICC fails to find existing std::exception_ptr in libstdc++ without this define
#if __INTEL_COMPILER && __GNUC__ && __TBB_EXCEPTION_PTR_PRESENT && !defined(__GCC_HAVE_SYNC_COMPARE_AND_SWAP_4)
    #define __GCC_HAVE_SYNC_COMPARE_AND_SWAP_4 1
#endif

// Work around a bug in MinGW32
#if __MINGW32__ && __TBB_EXCEPTION_PTR_PRESENT && !defined(_GLIBCXX_ATOMIC_BUILTINS_4)
    #define _GLIBCXX_ATOMIC_BUILTINS_4
#endif

#if __GNUC__ || __SUNPRO_CC || __IBMCPP__
    /* ICC defines __GNUC__ and so is covered */
    #define __TBB_ATTRIBUTE_ALIGNED_PRESENT 1
#elif _MSC_VER && (_MSC_VER >= 1300 || __INTEL_COMPILER)
    #define __TBB_DECLSPEC_ALIGN_PRESENT 1
#endif

/* Actually ICC supports gcc __sync_* intrinsics starting 11.1,
 * but 64 bit support for 32 bit target comes in later ones*/
/* TODO: change the version back to 4.1.2 once macro __TBB_WORD_SIZE become optional */
#if (__TBB_GCC_VERSION >= 40306) || (__INTEL_COMPILER >= 1200)
    /** built-in atomics available in GCC since 4.1.2 **/
    #define __TBB_GCC_BUILTIN_ATOMICS_PRESENT 1
#endif

#if (__INTEL_COMPILER >= 1210)
    /** built-in C++11 style atomics available in compiler since 12.1 **/
    #define __TBB_ICC_BUILTIN_ATOMICS_PRESENT 1
#endif

/** User controlled TBB features & modes **/

#ifndef TBB_USE_DEBUG
#ifdef TBB_DO_ASSERT
#define TBB_USE_DEBUG TBB_DO_ASSERT
#else
#ifdef _DEBUG
#define TBB_USE_DEBUG _DEBUG
#else
#define TBB_USE_DEBUG 0
#endif
#endif /* TBB_DO_ASSERT */
#endif /* TBB_USE_DEBUG */

#ifndef TBB_USE_ASSERT
#ifdef TBB_DO_ASSERT
#define TBB_USE_ASSERT TBB_DO_ASSERT
#else
#define TBB_USE_ASSERT TBB_USE_DEBUG
#endif /* TBB_DO_ASSERT */
#endif /* TBB_USE_ASSERT */

#ifndef TBB_USE_THREADING_TOOLS
#ifdef TBB_DO_THREADING_TOOLS
#define TBB_USE_THREADING_TOOLS TBB_DO_THREADING_TOOLS
#else
#define TBB_USE_THREADING_TOOLS TBB_USE_DEBUG
#endif /* TBB_DO_THREADING_TOOLS */
#endif /* TBB_USE_THREADING_TOOLS */

#ifndef TBB_USE_PERFORMANCE_WARNINGS
#ifdef TBB_PERFORMANCE_WARNINGS
#define TBB_USE_PERFORMANCE_WARNINGS TBB_PERFORMANCE_WARNINGS
#else
#define TBB_USE_PERFORMANCE_WARNINGS TBB_USE_DEBUG
#endif /* TBB_PEFORMANCE_WARNINGS */
#endif /* TBB_USE_PERFORMANCE_WARNINGS */

#if __MIC__ || __MIC2__
#define __TBB_DEFINE_MIC 1
#endif

#if !defined(__EXCEPTIONS) && !defined(_CPPUNWIND) && !defined(__SUNPRO_CC) || defined(_XBOX)
    #if TBB_USE_EXCEPTIONS
        #error Compilation settings do not support exception handling. Please do not set TBB_USE_EXCEPTIONS macro or set it to 0.
    #elif !defined(TBB_USE_EXCEPTIONS)
        #define TBB_USE_EXCEPTIONS 0
    #endif
#elif !defined(TBB_USE_EXCEPTIONS)
    #if __TBB_DEFINE_MIC
    #define TBB_USE_EXCEPTIONS 0
    #else
    #define TBB_USE_EXCEPTIONS 1
    #endif
#elif TBB_USE_EXCEPTIONS && __TBB_DEFINE_MIC
    #error Please do not set TBB_USE_EXCEPTIONS macro or set it to 0.
#endif

#ifndef TBB_IMPLEMENT_CPP0X
    /** By default, use C++0x classes if available **/
    #if __GNUC__==4 && __GNUC_MINOR__>=4 && __GXX_EXPERIMENTAL_CXX0X__
        #define TBB_IMPLEMENT_CPP0X 0
    #elif __clang__ && __cplusplus >= 201103L
        //TODO: consider introducing separate macroses for each file?
        //prevent injection of according tbb names into std:: namespace if native headers are present
        #if __has_include(<thread>) || __has_include(<condition_variable>)
            #define TBB_IMPLEMENT_CPP0X 0
        #else
            #define TBB_IMPLEMENT_CPP0X 1
        #endif
    #else
        #define TBB_IMPLEMENT_CPP0X 1
    #endif
#endif /* TBB_IMPLEMENT_CPP0X */

/* TBB_USE_CAPTURED_EXCEPTION should be explicitly set to either 0 or 1, as it is used as C++ const */
#ifndef TBB_USE_CAPTURED_EXCEPTION
    /**TODO: enable it by default on OS X*, once it is enabled in pre-built binary **/
    /** OS X* and IA64 pre-built TBB binaries do not support exception_ptr. **/
    #if __TBB_EXCEPTION_PTR_PRESENT && !defined(__APPLE__) && !defined(__ia64__)
        #define TBB_USE_CAPTURED_EXCEPTION 0
    #else
        #define TBB_USE_CAPTURED_EXCEPTION 1
    #endif
#else /* defined TBB_USE_CAPTURED_EXCEPTION */
    #if !TBB_USE_CAPTURED_EXCEPTION && !__TBB_EXCEPTION_PTR_PRESENT
        #error Current runtime does not support std::exception_ptr. Set TBB_USE_CAPTURED_EXCEPTION and make sure that your code is ready to catch tbb::captured_exception.
    #endif
#endif /* defined TBB_USE_CAPTURED_EXCEPTION */

/** Check whether the request to use GCC atomics can be satisfied **/
#if (TBB_USE_GCC_BUILTINS && !__TBB_GCC_BUILTIN_ATOMICS_PRESENT)
    #error "GCC atomic built-ins are not supported."
#endif

/** Internal TBB features & modes **/

/** __TBB_WEAK_SYMBOLS_PRESENT denotes that the system supports the weak symbol mechanism **/
#define __TBB_WEAK_SYMBOLS_PRESENT ( !_WIN32 && !__APPLE__ && !__sun && ((__TBB_GCC_VERSION >= 40000) || __INTEL_COMPILER ) )

/** __TBB_DYNAMIC_LOAD_ENABLED describes the system possibility to load shared libraries at run time **/
#ifndef __TBB_DYNAMIC_LOAD_ENABLED
    #define __TBB_DYNAMIC_LOAD_ENABLED 1
#endif

/** __TBB_SOURCE_DIRECTLY_INCLUDED is a mode used in whitebox testing when
    it's necessary to test internal functions not exported from TBB DLLs
**/
#if (_WIN32||_WIN64) && __TBB_SOURCE_DIRECTLY_INCLUDED
    #define __TBB_NO_IMPLICIT_LINKAGE 1
    #define __TBBMALLOC_NO_IMPLICIT_LINKAGE 1
#endif

#ifndef __TBB_COUNT_TASK_NODES
    #define __TBB_COUNT_TASK_NODES TBB_USE_ASSERT
#endif

#ifndef __TBB_TASK_GROUP_CONTEXT
    #define __TBB_TASK_GROUP_CONTEXT 1
#endif /* __TBB_TASK_GROUP_CONTEXT */

#ifndef __TBB_SCHEDULER_OBSERVER
    #define __TBB_SCHEDULER_OBSERVER 1
#endif /* __TBB_SCHEDULER_OBSERVER */

#if !defined(TBB_PREVIEW_TASK_ARENA) && __TBB_BUILD
    #define TBB_PREVIEW_TASK_ARENA __TBB_CPF_BUILD
#endif /* TBB_PREVIEW_TASK_ARENA */
#define __TBB_TASK_ARENA TBB_PREVIEW_TASK_ARENA
#if TBB_PREVIEW_TASK_ARENA
    #define TBB_PREVIEW_LOCAL_OBSERVER 1
    #define __TBB_NO_IMPLICIT_LINKAGE 1
    #define __TBB_RECYCLE_TO_ENQUEUE 1
    #define __TBB_TASK_PRIORITY 0 // TODO: it will be removed in next versions
    #if !__TBB_SCHEDULER_OBSERVER
        #error TBB_PREVIEW_TASK_ARENA requires __TBB_SCHEDULER_OBSERVER to be enabled
    #endif
#endif /* TBB_PREVIEW_TASK_ARENA */

#if !defined(TBB_PREVIEW_LOCAL_OBSERVER) && __TBB_BUILD && __TBB_SCHEDULER_OBSERVER
    #define TBB_PREVIEW_LOCAL_OBSERVER 1
#endif /* TBB_PREVIEW_LOCAL_OBSERVER */

#if TBB_USE_EXCEPTIONS && !__TBB_TASK_GROUP_CONTEXT
    #error TBB_USE_EXCEPTIONS requires __TBB_TASK_GROUP_CONTEXT to be enabled
#endif

#ifndef __TBB_TASK_PRIORITY
    #define __TBB_TASK_PRIORITY __TBB_TASK_GROUP_CONTEXT
#endif /* __TBB_TASK_PRIORITY */

#if __TBB_TASK_PRIORITY && !__TBB_TASK_GROUP_CONTEXT
    #error __TBB_TASK_PRIORITY requires __TBB_TASK_GROUP_CONTEXT to be enabled
#endif

#if TBB_PREVIEW_WAITING_FOR_WORKERS || __TBB_BUILD
    #define __TBB_SUPPORTS_WORKERS_WAITING_IN_TERMINATE 1
#endif

#if !defined(__TBB_SURVIVE_THREAD_SWITCH) && \
          (_WIN32 || _WIN64 || __APPLE__ || (__linux__ && !__ANDROID__))
    #define __TBB_SURVIVE_THREAD_SWITCH 1
#endif /* __TBB_SURVIVE_THREAD_SWITCH */

#ifndef __TBB_DEFAULT_PARTITIONER
#if TBB_DEPRECATED
/** Default partitioner for parallel loop templates in TBB 1.0-2.1 */
#define __TBB_DEFAULT_PARTITIONER tbb::simple_partitioner
#else
/** Default partitioner for parallel loop templates since TBB 2.2 */
#define __TBB_DEFAULT_PARTITIONER tbb::auto_partitioner
#endif /* TBB_DEPRECATED */
#endif /* !defined(__TBB_DEFAULT_PARTITIONER */

#ifdef _VARIADIC_MAX
#define __TBB_VARIADIC_MAX _VARIADIC_MAX
#else
#if _MSC_VER >= 1700
#define __TBB_VARIADIC_MAX 5  /* current VS11 setting, may change. */
#else
#define __TBB_VARIADIC_MAX 10
#endif
#endif

// Define preprocessor symbols used to determine architecture
#if _WIN32||_WIN64
#   if defined(_M_X64)||defined(__x86_64__)  // the latter for MinGW support
#       define __TBB_x86_64 1
#   elif defined(_M_IA64)
#       define __TBB_ipf 1
#   elif defined(_M_IX86)||defined(__i386__) // the latter for MinGW support
#       define __TBB_x86_32 1
#   endif
#else /* Assume generic Unix */
#   if !__linux__ && !__APPLE__
#       define __TBB_generic_os 1
#   endif
#   if __x86_64__
#       define __TBB_x86_64 1
#   elif __ia64__
#       define __TBB_ipf 1
#   elif __i386__||__i386  // __i386 is for Sun OS
#       define __TBB_x86_32 1
#   else
#       define __TBB_generic_arch 1
#   endif
#endif
/** Macros of the form __TBB_XXX_BROKEN denote known issues that are caused by
    the bugs in compilers, standard or OS specific libraries. They should be
    removed as soon as the corresponding bugs are fixed or the buggy OS/compiler
    versions go out of the support list.
**/

#if __ANDROID__ && __TBB_GCC_VERSION <= 40403 && !__GCC_HAVE_SYNC_COMPARE_AND_SWAP_8
    /** Necessary because on Android 8-byte CAS and F&A are not available for some processor architectures,
        but no mandatory warning message appears from GCC 4.4.3. Instead, only a linkage error occurs when
        these atomic operations are used (such as in unit test test_atomic.exe). **/
    #define __TBB_GCC_64BIT_ATOMIC_BUILTINS_BROKEN 1
#endif

#if __GNUC__ && __TBB_x86_64 && __INTEL_COMPILER == 1200
    #define __TBB_ICC_12_0_INL_ASM_FSTCW_BROKEN 1
#endif

#if _MSC_VER && __INTEL_COMPILER && (__INTEL_COMPILER<1110 || __INTEL_COMPILER==1110 && __INTEL_COMPILER_BUILD_DATE < 20091012)
    /** Necessary to avoid ICL error (or warning in non-strict mode):
        "exception specification for implicitly declared virtual destructor is
        incompatible with that of overridden one". **/
    #define __TBB_DEFAULT_DTOR_THROW_SPEC_BROKEN 1
#endif

#if defined(_MSC_VER) && _MSC_VER < 1500 && !defined(__INTEL_COMPILER)
    /** VS2005 and earlier do not allow declaring template class as a friend
        of classes defined in other namespaces. **/
    #define __TBB_TEMPLATE_FRIENDS_BROKEN 1
#endif

//TODO: recheck for different clang versions 
#if __GLIBC__==2 && __GLIBC_MINOR__==3 || __MINGW32__ || (__APPLE__ && (__clang__ || __INTEL_COMPILER==1200 && !TBB_USE_DEBUG))
    /** Macro controlling EH usages in TBB tests.
        Some older versions of glibc crash when exception handling happens concurrently. **/
    #define __TBB_THROW_ACROSS_MODULE_BOUNDARY_BROKEN 1
#else
    #define __TBB_THROW_ACROSS_MODULE_BOUNDARY_BROKEN 0
#endif

#if (_WIN32||_WIN64) && __INTEL_COMPILER == 1110
    /** That's a bug in Intel compiler 11.1.044/IA-32/Windows, that leads to a worker thread crash on the thread's startup. **/
    #define __TBB_ICL_11_1_CODE_GEN_BROKEN 1
#endif

#if __clang__ || (__GNUC__==3 && __GNUC_MINOR__==3 && !defined(__INTEL_COMPILER))
    /** Bugs with access to nested classes declared in protected area */
    #define __TBB_PROTECTED_NESTED_CLASS_BROKEN 1
#endif

#if __MINGW32__ && (__GNUC__<4 || __GNUC__==4 && __GNUC_MINOR__<2)
    /** MinGW has a bug with stack alignment for routines invoked from MS RTLs.
        Since GCC 4.2, the bug can be worked around via a special attribute. **/
    #define __TBB_SSE_STACK_ALIGNMENT_BROKEN 1
#else
    #define __TBB_SSE_STACK_ALIGNMENT_BROKEN 0
#endif

#if __GNUC__==4 && __GNUC_MINOR__==3 && __GNUC_PATCHLEVEL__==0
    /* GCC of this version may rashly ignore control dependencies */
    #define __TBB_GCC_OPTIMIZER_ORDERING_BROKEN 1
#endif

#if __FreeBSD__
    /** A bug in FreeBSD 8.0 results in kernel panic when there is contention
        on a mutex created with this attribute. **/
    #define __TBB_PRIO_INHERIT_BROKEN 1

    /** A bug in FreeBSD 8.0 results in test hanging when an exception occurs
        during (concurrent?) object construction by means of placement new operator. **/
    #define __TBB_PLACEMENT_NEW_EXCEPTION_SAFETY_BROKEN 1
#endif /* __FreeBSD__ */

#if (__linux__ || __APPLE__) && __i386__ && defined(__INTEL_COMPILER)
    /** The Intel compiler for IA-32 (Linux|OS X) crashes or generates
        incorrect code when __asm__ arguments have a cast to volatile. **/
    #define __TBB_ICC_ASM_VOLATILE_BROKEN 1
#endif

#if !__INTEL_COMPILER && (_MSC_VER || __GNUC__==3 && __GNUC_MINOR__<=2)
    /** Bug in GCC 3.2 and MSVC compilers that sometimes return 0 for __alignof(T)
        when T has not yet been instantiated. **/
    #define __TBB_ALIGNOF_NOT_INSTANTIATED_TYPES_BROKEN 1
#endif

/* Actually for Clang it should be name __TBB_CPP11_STD_FORWARD_PRESENT.
 * But in order to check for presence of std:: library feature we need to recognize
 * is standard library actually used stdlibc++ (GNU one) or libc++ (clang one).
 * Unfortunately it is not possible at the moment. So postponing it to later moment.*/
/*TODO: for clang rename it to __TBB_CPP11_STD_FORWARD_PRESENT and re-implement it.*/
#if (__INTEL_COMPILER) || (__clang__ &&  __TBB_GCC_VERSION <= 40300)
    #define __TBB_CPP11_STD_FORWARD_BROKEN 1
#else
    #define __TBB_CPP11_STD_FORWARD_BROKEN 0
#endif

#if __TBB_DEFINE_MIC
    /** Main thread and user's thread have different default thread affinity masks. **/
    #define __TBB_MAIN_THREAD_AFFINITY_BROKEN 1
#endif

/** __TBB_WIN8UI_SUPPORT enables support of New Windows*8 Store Apps and limit a possibility to load
    shared libraries at run time only from application container **/
#if defined(WINAPI_FAMILY) && WINAPI_FAMILY == WINAPI_FAMILY_APP
    #define __TBB_WIN8UI_SUPPORT 1
#else
    #define __TBB_WIN8UI_SUPPORT 0
#endif

#if !defined(__EXCEPTIONS) && __GNUC__==4 && (__GNUC_MINOR__==4 ||__GNUC_MINOR__==5 || (__INTEL_COMPILER==1300 && __TBB_GCC_VERSION>=40600 && __TBB_GCC_VERSION<=40700)) && defined(__GXX_EXPERIMENTAL_CXX0X__)
/* There is an issue for specific GCC toolchain when C++11 is enabled
   and exceptions are disabled:
   exceprion_ptr.h/nested_exception.h are using throw unconditionally.
 */
    #define __TBB_LIBSTDCPP_EXCEPTION_HEADERS_BROKEN 1
#else
    #define __TBB_LIBSTDCPP_EXCEPTION_HEADERS_BROKEN 0
#endif

#if __TBB_x86_32 && (__linux__ || __APPLE__ || _WIN32 || __sun) &&  ((defined(__INTEL_COMPILER) && (__INTEL_COMPILER <= 1300)) || (__GNUC__==3 && __GNUC_MINOR__==3 ) || defined(__SUNPRO_CC))
    // Some compilers for IA-32 fail to provide 8-byte alignment of objects on the stack,
    // even if the object specifies 8-byte alignment.  On such platforms, the IA-32 implementation
    // of 64 bit atomics (e.g. atomic<long long>) use different tactics depending upon
    // whether the object is properly aligned or not.
    #define __TBB_FORCE_64BIT_ALIGNMENT_BROKEN 1
#else
    #define __TBB_FORCE_64BIT_ALIGNMENT_BROKEN 0
#endif

#if (__TBB_DEFAULTED_AND_DELETED_FUNC_PRESENT && (__TBB_GCC_VERSION < 40700) && (!defined(__INTEL_COMPILER) && !defined (__clang__)))
    #define __TBB_ZERO_INIT_WITH_DEFAULTED_CTOR_BROKEN 1
#endif
/** End of __TBB_XXX_BROKEN macro section **/

#define __TBB_ATOMIC_CTORS     (__TBB_CONSTEXPR_PRESENT && __TBB_DEFAULTED_AND_DELETED_FUNC_PRESENT && (!__TBB_ZERO_INIT_WITH_DEFAULTED_CTOR_BROKEN))

#endif /* __TBB_tbb_config_H */
