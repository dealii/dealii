/*----------------------------   convergence_table.h     ---------------------------*/
/*      $Id$                 */
#ifndef __convergence_table_H
#define __convergence_table_H
/*----------------------------   convergence_table.h     ---------------------------*/


#include <base/table_handler.h>


/**
 * The #ConvergenceTable# class is an application to the #TableHandler# class and
 * stores some convergence data, such as residuals of the cg-method, or some evaluated
 * $L^2$-errors of discrete solutions, etc, and evaluates convergence rates or orders.
 *
 * The already implemented #RateMode#s are #ConvergenceTable::reduction_rate#, 
 * where the convergence rate is the quotient of two following rows, and 
 * #ConvergenceTable::reduction_rate_log2#, that evaluates the order of convergence.
 * These standard evaluations are useful for global refinement, for local refinement
 * this may not be a appropriate method, as the convergence rates should be
 * set in relation to the number of cells or the number of DoFs. The
 * implementations of these non-standard methods is left to a user.
 *
 * \subsection{Usage}
 * The number of cells and the number of DoFs may be added to the table by
 * calling e.g.  #add_value("n cells", n_cells)#. The table
 * data is also added by calling #add_value(...)#.
 * Before the output of the table the functions #evaluate_convergence_rates(...)# and
 * #evaluate_all_convergence_rates(...)# may be (also multiply) called.
 *
 * There are two possibilities of how to evaluate the convergence rates of multiple 
 * columns in the same #RateMode#.
 * \begin{enumeration}
 * \item call #evaluate_convergence_rates(data_column_key, ...)# for all wanted columns
 * \item call #omit_column_from_convergence_rate_evaluation(data_column_key)# for all
 * NOT wanted columns and then 
 * #evaluate_all_convergence_rates(...)# to evaluate the convergence rates of all columns
 * that are not signed to be omitted.
 * \end{enumeration}
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
				      * reference column):
				      * 
				      * #reduction_rate#: $value(row-1)/value(row)$.
				      *
				      * #reduction_rate_log2#: the order $O(h^r) of 
				      * the reduction rate with the assumption
				      * row-1: h, row: h/2. Hence the order 
				      * is evaluated by
				      * $\log(value(row-1)/value(row))/\log 2$.
				      *
				      * Rate in relation to a reference column:
				      *
				      * #reduction_rate#: Not implemented.
				      *
				      * #reduction_rate_log2#: Not implemented.
				      */
    enum RateMode {
	  none,
	  reduction_rate,
	  reduction_rate_log2
    };
    
				     /**
				      * Evaluates the convergence rates of the
				      * data column #data_column_key# 
				      * due to the #RateMode# in relation to
				      * the reference column #reference_column_key#.
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
				      * the #set_tex_supercaption (...)# function
				      * of the base class #TableHandler#.
				      *
				      * Still not implemented.
				      */
    void evaluate_convergence_rates (const string &data_column_key,
				     const string &reference_column_key,
				     const RateMode rate_mode);

   
				     /**
				      * Evaluates the convergence rates of the
				      * data column #data_column_key# 
				      * due to the #RateMode#.
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
				      * the #set_tex_supercaption (...)# function
				      * of the base class #TableHandler#.
				      */
    void evaluate_convergence_rates (const string &data_column_key, 
				     const RateMode rate_mode);

				     /**
				      * Omit this column #key# 
				      * (not supercolumn!) from the
				      * evaluation of the convergence rates
				      * of `all' columns (see the following 
				      * two functions).
				      *
				      * The Column::flag==1 is preserved for
				      * omitting the column from convergence
				      * rate evalution.
				      */
    void omit_column_from_convergence_rate_evaluation(const string &key);

				     /**
				      * Evaluates convergence rates
				      * due to the #rate_mode# in relation
				      * to the reference column 
				      * #reference_column_key#. This
				      * function evaluates the rates of
				      * ALL columns except of the 
				      * #omit_columns# and previously
				      * evaluated rate columns.
				      * This function allows to evaluate
				      * the convergence rate for almost all
				      * columns of a table without calling
				      * #evaluate_convergence_rates(data_column, ...)#
				      * for each column seperately.
				      *
				      * Example: 
				      * Columns like #n cells# or 
				      * #n dofs# columns may be wanted
				      * to be omitted in the evaluation
				      * of the convergence rates. Hence they
				      * should be included into #omit_columns#.
				      */
    void evaluate_all_convergence_rates(const string &reference_column_key,
					const RateMode rate_mode);

				     /**
				      * Evaluates convergence rates
				      * due to the #rate_mode#. This
				      * function evaluates the rates of
				      * ALL columns except of the 
				      * #omit_columns# and previously
				      * evaluated rate columns.
				      * This function allows to evaluate
				      * the convergence rate for almost all
				      * columns of a table without calling
				      * #evaluate_convergence_rates(data_column, ...)#
				      * for each column seperately.
				      *
				      * Example: 
				      * Columns like #n cells# or 
				      * #n dofs# columns may be wanted
				      * to be omitted in the evaluation
				      * of the convergence rates. Hence they
				      * should be included into #omit_columns#.
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
		    string,
		    << "Rate column <" << arg1 << "> does already exist.");
};

    



/*----------------------------   convergence_table.h     ---------------------------*/
/* end of #ifndef __convergence_table_H */
#endif
/*----------------------------   convergence_table.h     ---------------------------*/
