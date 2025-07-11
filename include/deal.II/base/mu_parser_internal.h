// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_mu_parser_internal_h
#define dealii_mu_parser_internal_h

// This file contains functions used internally by the FunctionParser
// and the TensorFunctionParser class.

#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/thread_local_storage.h>

#include <memory>
#include <string>
#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace FunctionParser
  {
    int
    mu_round(double val);

    double
    mu_if(double condition, double thenvalue, double elsevalue);

    double
    mu_or(double left, double right);

    double
    mu_and(double left, double right);

    double
    mu_int(double value);

    double
    mu_ceil(double value);

    double
    mu_floor(double value);

    double
    mu_cot(double value);

    double
    mu_csc(double value);

    double
    mu_sec(double value);

    double
    mu_log(double value);

    double
    mu_pow(double a, double b);

    double
    mu_erf(double value);

    double
    mu_erfc(double value);

    // returns a random value in the range [0,1] initializing the generator
    // with the given seed
    double
    mu_rand_seed(double seed);

    // returns a random value in the range [0,1]
    double
    mu_rand();

    /**
     * Get the array of all function names.
     */
    std::vector<std::string>
    get_function_names();

    /**
     * @addtogroup Exceptions
     * @{
     */
    DeclException2(ExcParseError,
                   int,
                   std::string,
                   << "Parsing Error at Column " << arg1
                   << ". The parser said: " << arg2);

    /** @} */

    /**
     * deal.II uses muParser as a purely internal dependency. To this end, we do
     * not include any muParser headers in our own headers (and the bundled
     * version of the dependency does not install its headers or compile a
     * separate muparser library). Hence, to interface with muParser, we use the
     * PIMPL idiom here to wrap a pointer to mu::Parser objects.
     */
    class muParserBase
    {
    public:
      virtual ~muParserBase() = default;
    };

    /**
     * Class containing the mutable state required by muParser.
     *
     * @note For performance reasons it is best to put all mutable state in a
     * single object so that, for each function call, we only need to get
     * thread-local data exactly once.
     */
    struct ParserData
    {
      /**
       * Default constructor. Threads::ThreadLocalStorage requires that objects
       * be either default- or copy-constructible: make sure we satisfy the
       * first case by declaring it here.
       */
      ParserData() = default;

      /**
       * std::is_copy_constructible gives the wrong answer for containers with
       * non-copy constructible types (e.g., std::vector<std::unique_ptr<int>>)
       * - for more information, see the documentation of
       * Threads::ThreadLocalStorage. Hence, to avoid compilation failures, just
       * delete the copy constructor completely.
       */
      ParserData(const ParserData &) = delete;

      /**
       * Scratch array used to set independent variables (i.e., x, y, and t)
       * before each muParser call.
       */
      std::vector<double> vars;

      /**
       * The actual muParser parser objects (hidden with PIMPL).
       */
      std::vector<std::unique_ptr<muParserBase>> parsers;
    };

    template <int dim, typename Number>
    class ParserImplementation
    {
    public:
      ParserImplementation();

      virtual ~ParserImplementation() = default;

      /**
       * Initialize the internal state of the object. This is the same as the
       * inheriting class method - see FunctionParser::initialize() for more
       * information.
       */
      virtual void
      initialize(const std::string                   &vars,
                 const std::vector<std::string>      &expressions,
                 const std::map<std::string, double> &constants,
                 const bool                           time_dependent = false);

      /**
       * Set up the internal muParser objects to parse and evaluate mathematical
       * expressions.
       */
      void
      init_muparser() const;

      /**
       * Compute the value of a single component.
       */
      Number
      do_value(const Point<dim> &p,
               const double      time,
               unsigned int      component) const;

      /**
       * Compute the values of all components.
       */
      void
      do_all_values(const Point<dim>  &p,
                    const double       time,
                    ArrayView<Number> &values) const;

      /**
       * An array of function expressions (one per component), required to
       * initialize tfp in each thread.
       */
      std::vector<std::string> expressions;

    private:
      /**
       * The muParser objects (hidden with the PIMPL idiom) for each thread (and
       * one for each component).
       */
      mutable Threads::ThreadLocalStorage<internal::FunctionParser::ParserData>
        parser_data;

      /**
       * An array to keep track of all the constants, required to initialize fp
       * in each thread.
       */
      std::map<std::string, double> constants;

      /**
       * An array for the variable names, required to initialize fp in each
       * thread.
       */
      std::vector<std::string> var_names;

      /**
       * State of usability. This variable is checked every time the function is
       * called for evaluation. It's set to true in the initialize() methods.
       */
      bool initialized;

      /**
       * Number of variables. If this is also a function of time, then the
       * number of variables is dim+1, otherwise it is dim. In the case that
       * this is a time dependent function, the time is supposed to be the last
       * variable. If #n_vars is not identical to the number of the variables
       * parsed by the initialize() method, then an exception is thrown.
       */
      unsigned int n_vars;
    };
  } // namespace FunctionParser
} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
