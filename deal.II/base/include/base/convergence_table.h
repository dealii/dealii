//----------------------------  convergence_table.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  convergence_table.h  ---------------------------
#ifndef __deal2__convergence_table_h
#define __deal2__convergence_table_h


#include <base/config.h>
#include <base/table_handler.h>


/**
 * The @p{ConvergenceTable} class is an application to the @p{TableHandler} class and
 * stores some convergence data, such as residuals of the cg-method, or some evaluated
 * $L^2$-errors of discrete solutions, etc, and evaluates convergence rates or orders.
 *
 * The already implemented @p{RateMode}s are @p{ConvergenceTable::reduction_rate}, 
 * where the convergence rate is the quotient of two following rows, and 
 * @p{ConvergenceTable::reduction_rate_log2}, that evaluates the order of convergence.
 * These standard evaluations are useful for global refinement, for local refinement
 * this may not be an appropriate method, as the convergence rates should be
 * set in relation to the number of cells or the number of DoFs. The
 * implementations of these non-standard methods is left to a user.
 *
 * @sect3{Usage}
 * The number of cells and the number of DoFs may be added to the table by
 * calling e.g.  @p{add_value("n cells", n_cells)}. The table
 * data is also added by calling @p{add_value(...)}.
 * Before the output of the table the functions @p{evaluate_convergence_rates(...)} and
 * @p{evaluate_all_convergence_rates(...)} may be (also multiply) called.
 *
 * There are two possibilities of how to evaluate the convergence rates of multiple 
 * columns in the same @p{RateMode}.
 * @begin{enumerate}
 * @item call @p{evaluate_convergence_rates(data_column_key, ...)} for all wanted columns
 * @item call @p{omit_column_from_convergence_rate_evaluation(data_column_key)} for all
 * NOT wanted columns and then 
 * @p{evaluate_all_convergence_rates(...)} to evaluate the convergence rates of all columns
 * that are not signed to be omitted.
 * @end{enumerate}
 *
 *
 * @author Ralf Hartmann, 1999
 */
class ConvergenceTable: public TableHandler
{
  public:
				     /**
				      * Constructor.
				      */
    ConvergenceTable();
    
				     /**
				      * Rate in relation to the rows (without a
				      * reference column).
				      * 
				      * @p{reduction_rate}: $value(row-1)/value(row)$.
				      *
				      * @p{reduction_rate_log2}: the order $O(h^r)$ of 
				      * the reduction rate with the assumption
				      * row-1: h, row: h/2. Hence the order 
				      * is evaluated by
				      * $log(value(row-1)/value(row))/log 2$.
				      *
				      * Rate in relation to a reference column:
				      *
				      * @p{reduction_rate}: Not implemented.
				      *
				      * @p{reduction_rate_log2}: Not implemented.
				      */
    enum RateMode {
	  none,
	  reduction_rate,
	  reduction_rate_log2
    };
    
				     /**
				      * Evaluates the convergence rates of the
				      * data column @p{data_column_key} 
				      * due to the @p{RateMode} in relation to
				      * the reference column @p{reference_column_key}.
				      * Be sure that the value types of the
				      * table entries of the
				      * data column and the reference data column
				      * is a number, i.e. double, float,
				      * (unsigned) int, and so on.
				      *
				      * The new rate column and the data column
				      * will be merged to a supercolumn. The
				      * tex caption of the supercolumn will be
				      * (by default) the same as the one of the
				      * data column. This may be changed by using
				      * the @p{set_tex_supercaption (...)} function
				      * of the base class @p{TableHandler}.
				      *
				      * Still not implemented.
				      */
    void evaluate_convergence_rates (const std::string &data_column_key,
				     const std::string &reference_column_key,
				     const RateMode     rate_mode);


				     /**
				      * Evaluates the convergence rates of the
				      * data column @p{data_column_key} 
				      * due to the @p{RateMode}.
				      * Be sure that the value types of the
				      * table entries of the data column
				      * is a number, i.e. double, float,
				      * (unsigned) int, and so on.
				      *
				      * The new rate column and the data column
				      * will be merged to a supercolumn. The
				      * tex caption of the supercolumn will be
				      * (by default) the same as the one of the
				      * data column. This may be changed by using
				      * the @p{set_tex_supercaption (...)} function
				      * of the base class @p{TableHandler}.
				      */
    void evaluate_convergence_rates (const std::string &data_column_key, 
				     const RateMode     rate_mode);

				     /**
				      * Omit this column @p{key} 
				      * (not supercolumn!) from the
				      * evaluation of the convergence rates
				      * of `all' columns (see the following 
				      * two functions).
				      *
				      * The Column::flag==1 is reserved for
				      * omitting the column from convergence
				      * rate evalution.
				      */
    void omit_column_from_convergence_rate_evaluation(const std::string &key);

				     /**
				      * Evaluates convergence rates
				      * due to the @p{rate_mode} in relation
				      * to the reference column 
				      * @p{reference_column_key}. This
				      * function evaluates the rates of
				      * ALL columns except of the 
				      * columns that are to be omitted
				      * (see previous function)
				      * and execpt of the columns that are
				      * previously
				      * evaluated rate columns.
				      * This function allows to evaluate
				      * the convergence rate for almost all
				      * columns of a table without calling
				      * @p{evaluate_convergence_rates(data_column, ...)}
				      * for each column separately.
				      *
				      * Example: 
				      * Columns like @p{n cells} or 
				      * @p{n dofs} columns may be wanted
				      * to be omitted in the evaluation
				      * of the convergence rates. Hence they
				      * should omitted by calling the
				      * @p{omit_column_from_convergence_rate_evaluation(..)}.
				      */
    void evaluate_all_convergence_rates(const std::string &reference_column_key,
					const RateMode     rate_mode);

				     /**
				      * Evaluates convergence rates
				      * due to the @p{rate_mode}. This
				      * function evaluates the rates of
				      * ALL columns except of the 
				      * columns that are to be omitted
				      * (see previous function)
				      * and execpt of the columns that are
				      * previously
				      * evaluated rate columns.
				      * This function allows to evaluate
				      * the convergence rate for almost all
				      * columns of a table without calling
				      * @p{evaluate_convergence_rates(data_column, ...)}
				      * for each column seperately.
				      *
				      * Example: 
				      * Columns like @p{n cells} or 
				      * @p{n dofs} columns may be wanted
				      * to be omitted in the evaluation
				      * of the convergence rates. Hence they
				      * should omitted by calling the
				      * @p{omit_column_from_convergence_rate_evaluation(..)}.
				      */
    void evaluate_all_convergence_rates(const RateMode rate_mode);

    				     /**
				      * Exception
				      */
    DeclException0 (ExcWrongValueType);

    				     /**
				      * Exception
				      */
    DeclException1 (ExcRateColumnAlreadyExists,
		    std::string,
		    << "Rate column <" << arg1 << "> does already exist.");
};


#endif
