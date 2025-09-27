// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1999 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_convergence_table_h
#define dealii_convergence_table_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/table_handler.h>
#include <deal.II/base/types.h>

#include <ostream>
#include <string>


DEAL_II_NAMESPACE_OPEN


/**
 * The ConvergenceTable class is an application to the TableHandler class and
 * stores some convergence data, such as residuals of the cg-method, or some
 * evaluated <i>L<sup>2</sup></i>-errors of discrete solutions, etc, and
 * evaluates convergence rates or orders.
 *
 * The already implemented #RateMode's are #reduction_rate, where the
 * convergence rate is the quotient of two following rows, and
 * #reduction_rate_log2, that evaluates the order of convergence. These
 * standard evaluations are useful for global refinement, for local refinement
 * this may not be an appropriate method, as the convergence rates should be
 * set in relation to the number of cells or the number of DoFs. The
 * implementations of these non-standard methods is left to a user.
 *
 * For example, the number of cells may be added to the table by
 * calling `add_value("n cells", n_cells)`. The number of DoFs may be
 * added to the table by calling `add_value("n dofs", n_dofs)`. Of course,
 * one can also add more kinds of information by calling
 * add_value() with other arguments. In any case, before the output of the table
 * the functions evaluate_convergence_rates() and
 * evaluate_all_convergence_rates() may be called.
 *
 * There are two possibilities of how to evaluate the convergence rates of
 * multiple columns in the same RateMode.
 * <ol>
 * <li> call evaluate_convergence_rates() for all wanted columns
 * <li> call omit_column_from_convergence_rate_evaluation() for all columns
 * for which this evaluation is not desired and then
 * evaluate_all_convergence_rates() to evaluate the convergence rates of all
 * columns that have not been flagged for omission.
 * </ol>
 *
 * A detailed discussion of this class can also be found in the step-7 and
 * step-13 example programs. It is also used in step-74.
 *
 * @ingroup textoutput
 */
class ConvergenceTable : public TableHandler
{
public:
  /**
   * Constructor.
   */
  ConvergenceTable() = default;

  /**
   * Rate in relation to the rows.
   */
  enum RateMode
  {
    /**
     * Do not do anything.
     */
    none,
    /**
     * Quotient of values in the previous row and in this row.
     */
    reduction_rate,
    /**
     * Logarithm of #reduction_rate to the base 2 representing the order of
     * convergence when halving the grid size, e.g. from h to h/2.
     */
    reduction_rate_log2
  };

  /**
   * Evaluate the convergence rates of the data column
   * <tt>data_column_key</tt> due to the #RateMode in relation to the
   * reference column <tt>reference_column_key</tt>. Be sure that the value
   * types of the table entries of the data column and the reference data
   * column is a number, i.e. double, float, (unsigned) int, and so on.
   *
   * As this class has no information on the space dimension upon which the
   * reference column vs. the value column is based upon, it needs to be
   * passed as last argument to this method. The <i>default dimension for the
   * reference column</i> is 2, which is appropriate for the number of cells
   * in 2d. If you work in 3d, set the number to 3. If the reference column is
   * $1/h$, remember to set the dimension to 1 also when working in 3d to get
   * correct rates.
   *
   * The new rate column and the data column will be merged to a supercolumn.
   * The tex caption of the supercolumn will be (by default) the same as the
   * one of the data column. This may be changed by using the
   * <tt>set_tex_supercaption (...)</tt> function of the base class
   * TableHandler.
   *
   * This method behaves in the following way:
   *
   * If RateMode is reduction_rate, then the computed output is $
   * \frac{e_{n-1}/k_{n-1}}{e_n/k_n}, $ where $k$ is the reference column (no
   * dimension dependence!).
   *
   * If RateMode is reduction_rate_log2, then the computed output is $ dim
   * \frac{\log |e_{n-1}/e_{n}|}{\log |k_n/k_{n-1}|} $.
   *
   * This is useful, for example, if we use as reference key the number of
   * degrees of freedom or better, the number of cells.  Assuming that the
   * error is proportional to $ C (1/\sqrt{k})^r $ in 2d, then this method
   * will produce the rate $r$ as a result. For general dimension, as
   * described by the last parameter of this function, the formula needs to be
   * $ C (1/\sqrt[dim]{k})^r $.
   *
   * @note Since this function adds columns to the table after several rows
   * have already been filled, it switches off the auto fill mode of the
   * TableHandler base class. If you intend to add further data with auto
   * fill, you will have to re-enable it after calling this function.
   */
  void
  evaluate_convergence_rates(const std::string &data_column_key,
                             const std::string &reference_column_key,
                             const RateMode     rate_mode,
                             const unsigned int dim = 2);


  /**
   * Evaluate the convergence rates of the data column
   * <tt>data_column_key</tt> due to the #RateMode.  Be sure that the value
   * types of the table entries of the data column is a number, i.e. double,
   * float, (unsigned) int, and so on.
   *
   * The new rate column and the data column will be merged to a supercolumn.
   * The tex caption of the supercolumn will be (by default) the same as the
   * one of the data column. This may be changed by using the
   * set_tex_supercaption() function of the base class TableHandler.
   *
   * @note Since this function adds columns to the table after several rows
   * have already been filled, it switches off the auto fill mode of the
   * TableHandler base class. If you intend to add further data with auto
   * fill, you will have to re-enable it after calling this function.
   */
  void
  evaluate_convergence_rates(const std::string &data_column_key,
                             const RateMode     rate_mode);

  /**
   * Omit this column <tt>key</tt> (not supercolumn!) from the evaluation of
   * the convergence rates of `all' columns (see the following two functions).
   *
   * The Column::flag==1 is reserved for omitting the column from convergence
   * rate evaluation.
   */
  void
  omit_column_from_convergence_rate_evaluation(const std::string &key);

  /**
   * Evaluate convergence rates due to the <tt>rate_mode</tt> in relation to
   * the reference column <tt>reference_column_key</tt>. This function
   * evaluates the rates of ALL columns except of the columns that are to be
   * omitted (see previous function) and except of the columns that are
   * previously evaluated rate columns.  This function allows to evaluate the
   * convergence rate for almost all columns of a table without calling
   * evaluate_convergence_rates() for each column separately.
   *
   * Example: Columns like <tt>n cells</tt> or <tt>n dofs</tt> columns may be
   * wanted to be omitted in the evaluation of the convergence rates. Hence
   * they should omitted by calling the
   * omit_column_from_convergence_rate_evaluation().
   */
  void
  evaluate_all_convergence_rates(const std::string &reference_column_key,
                                 const RateMode     rate_mode);

  /**
   * Evaluate convergence rates due to the <tt>rate_mode</tt>. This function
   * evaluates the rates of ALL columns except of the columns that are to be
   * omitted (see previous function) and except of the columns that are
   * previously evaluated rate columns.  This function allows to evaluate the
   * convergence rate for almost all columns of a table without calling
   * evaluate_convergence_rates() for each column separately.
   *
   * Example: Columns like <tt>n cells</tt> or <tt>n dofs</tt> columns may be
   * wanted to be omitted in the evaluation of the convergence rates. Hence
   * they should omitted by calling the
   * omit_column_from_convergence_rate_evaluation().
   */
  void
  evaluate_all_convergence_rates(const RateMode rate_mode);

  /**
   * @addtogroup Exceptions
   * @{
   */

  /**
   * Exception
   */
  DeclException1(ExcRateColumnAlreadyExists,
                 std::string,
                 << "Rate column <" << arg1 << "> does already exist.");
  /** @} */
};


DEAL_II_NAMESPACE_CLOSE

#endif
