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
 * The already implemented #RateMode#s are #ConvergenceTable::standard#, 
 * where the convergence rate is the quotient of two following rows, and 
 * #ConvergenceTable::standard_order#, that evaluates the order of convergence.
 * These standard evaluations are useful for global refinement, for local refinement
 * this may not be a appropriate method, as the convergence rates should be
 * set in relation to the number of cells or the number of DoFs. The
 * implementations of these non-standard methods is left to the fantasy of a
 * fanatic user.
 *
 * \subsection{Usage}
 * The number of cells and the number of DoFs are added to the table by
 * calling #add_run(unsigned int ncells, unsigned int ndofs)#. The
 * data is added by #add_value(...)# of the base class #TableHandler#.
 * Before the output of the table the function #evaluate_convergence_rates#
 * may be (also multiply) called. 
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
    
    enum RateMode {
	  none, standard, standard_order, n_cells, n_dofs
    };

				     /**
				      * Adds the basic information
				      * (number of cells, number of DoFs)
				      * of a new run.
				      */
    void add_run (unsigned int ncells, 
		  unsigned int ndofs);

				     /**
				      * Evaluates the convergence rates of the
				      * column #key# due to a #RateMode# and
				      * merges the rate column and this column
				      * to a supercolumn.
				      */
    void evaluate_convergence_rates (const string &key, 
				     const RateMode conv_rate);

				     /**
				      * Evaluate the convergence rates of
				      * all (not "n cells" or "n dofs") 
				      * columns due to a #RateMode#.
				      */
    void evaluate_convergence_rates(const RateMode conv_rate);

    				     /**
				      * Exception
				      */
    DeclException0 (ExcWrongValueType)

  private:
				     /**
				      * Preset string "n cells".
				      */
    string n_cells_string;

				     /**
				      * Preset string "n dofs".
				      */
    string n_dofs_string;
};

    



/*----------------------------   convergence_table.h     ---------------------------*/
/* end of #ifndef __convergence_table_H */
#endif
/*----------------------------   convergence_table.h     ---------------------------*/
