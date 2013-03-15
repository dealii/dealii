/***************************************************************************\
|* Function Parser for C++ v4.5.1                                          *|
|*-------------------------------------------------------------------------*|
|* Copyright: Juha Nieminen                                                *|
|*                                                                         *|
|* This library is distributed under the terms of the                      *|
|* GNU Lesser General Public License version 3.                            *|
|* (See lgpl.txt and gpl.txt for the license text.)                        *|
\***************************************************************************/

// Configuration file
// ------------------

/* NOTE:
   This file is for the internal use of the function parser only.
   You don't need to include this file in your source files, just
   include "fparser.hh".
*/


/* Uncomment any of these lines or define them in your compiler settings
   to enable the correspondent version of the parser. (These are disabled
   by default because they rely on C99 functions, and non-standard libraries
   in the case pf MPFR and GMP, and they make compiling needlessly slower
   and the resulting binary needlessly larger if they are not used in the
   program.)
*/
//#define FP_SUPPORT_FLOAT_TYPE
//#define FP_SUPPORT_LONG_DOUBLE_TYPE
//#define FP_SUPPORT_LONG_INT_TYPE
//#define FP_SUPPORT_MPFR_FLOAT_TYPE
//#define FP_SUPPORT_GMP_INT_TYPE
//#define FP_SUPPORT_COMPLEX_DOUBLE_TYPE
//#define FP_SUPPORT_COMPLEX_FLOAT_TYPE
//#define FP_SUPPORT_COMPLEX_LONG_DOUBLE_TYPE

/* If you are using FunctionParser_ld or FunctionParser_cld and your compiler
   supports the strtold() function, you should uncomment the following line.
 */
//#define FP_USE_STRTOLD


/* Uncomment this line or define it in your compiler settings if you want
   to disable compiling the basic double version of the library, in case
   one of the above types is used but not the double type. (If the double
   type is not used, then disabling it makes compiling faster and the
   resulting binary smaller.)
 */
//#define FP_DISABLE_DOUBLE_TYPE

/* Uncomment this line or define it in your compiler settings to make the
   parser use C++11 math functions. (Note that these may not be supported
   by all compilers.)
*/
//#define FP_SUPPORT_CPLUSPLUS11_MATH_FUNCS

/*
 Whether to use shortcut evaluation for the & and | operators:
*/
#ifndef FP_DISABLE_SHORTCUT_LOGICAL_EVALUATION
#define FP_ENABLE_SHORTCUT_LOGICAL_EVALUATION
#endif

/*
 Comment out the following lines out if you are not going to use the
 optimizer and want a slightly smaller library. The Optimize() method
 can still be called, but it will not do anything.
 If you are unsure, just leave it. It won't slow down the other parts of
 the library.
*/
#ifndef FP_NO_SUPPORT_OPTIMIZER
#define FP_SUPPORT_OPTIMIZER
#endif

#if defined(FP_SUPPORT_COMPLEX_DOUBLE_TYPE) || defined(FP_SUPPORT_COMPLEX_FLOAT_TYPE) || defined(FP_SUPPORT_COMPLEX_LONG_DOUBLE_TYPE)
#define FP_SUPPORT_COMPLEX_NUMBERS
#endif


/*
 No member function of FunctionParser is thread-safe. Most prominently,
 Eval() is not thread-safe. By uncommenting one of these lines the Eval()
 function can be made thread-safe at the cost of a possible small overhead.
 The second version requires that the compiler supports the alloca() function,
 which is not standard, but is faster.
 */
//#define FP_USE_THREAD_SAFE_EVAL
//#define FP_USE_THREAD_SAFE_EVAL_WITH_ALLOCA
