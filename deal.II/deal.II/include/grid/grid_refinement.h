//----------------------------  grid_refinement.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  grid_refinement.h  ---------------------------
#ifndef __deal2__grid_refinement_h
#define __deal2__grid_refinement_h


#include <lac/vector.h>
#include <vector>

// forward declarations
template <int dim> class Triangulation;
template <class T> class Vector;



/**
 *   This class provides several function that flag certain cells for
 *   coarsening or refinement based on a vector of ``error''
 *   indicators and some selection algorithm.  The central function is
 *   @p{refine (const Vector<float> &criterion, const double threshold)}:
 *   it takes a vector of values, one per active cell,
 *   which denote the criterion according to which the triangulation
 *   is to be refined. It marks all cells for which the criterion is
 *   greater than the threshold being given as the second
 *   argument. Analogously,
 *   @p{coarsen (const Vector<float> &criterion, const double threshold)}
 *   flags those cells for
 *   coarsening for which the criterion is less than the threshold.
 *
 *   There are two variations of these functions, which rely on @p{refine} and
 *   @p{coarsen} by computing the thresholds from other information:
 *   @begin{itemize}
 *   @item @p{refine_and_coarsen_fixed_number}: this function takes a vector as
 *     above and two values between zero and one denoting the fractions of cells to
 *     be refined and coarsened. For this purpose, it sorts the criteria per cell
 *     and takes the threshold to be the one belonging to the cell with the
 *     @p{fraction times n_active_cells} highest criterion. For example, if
 *     the fraction is $0.3$, the threshold is computed to a value such that
 *     30 per cent of cells have a criterion higher than the threshold and are
 *     thus flagged for refinement. The flagging for refinement is done through
 *     the central @p{refine} function. For coarsening, the same holds.
 *
 *     The sorting of criteria is not done actually, since we only need one
 *     value, in the example above the criterion of the cell which is at
 *     30 per cent in the sorted list of cells. The order of cells with higher
 *     and of those with lower criteria is irrelevant. Getting this value is
 *     accomplished by the @p{nth_element} function of the @p{C++} standard
 *     library, which takes only linear time in the number of elements, rather
 *     than @p{N log N} for sorting all values.
 *
 *     A typical value for the fraction of cells to be refined is 0.3.
 *     However, for singular functions or singular error functionals, you may
 *     want to chose a smaller value to avoid overrefinement in regions which
 *     do not contribute much to the error.
 *
 *   @item @p{refine_and_coarsen_fixed_fraction}: this function computes the
 *     threshold such that the number of cells getting flagged for refinement
 *     makes up for a certain fraction of the total error. If this fraction is 50
 *     per cent, for example, the threshold is computed such that the cells with
 *     a criterion greater than the threshold together account for half of the
 *     total error. The definition of the fraction is a bit counterintuitive, since
 *     the total error is the sum over all cells of the local contribution
 *     squared. We define that the fraction $\alpha$ be such that those
 *     elements with the greatest error are refined for which the condition
 *     $\sum \eta_K^2 \le \alpha\eta^2$ holds. Note that $\alpha$ is not
 *     squared. The sum runs over the mentioned
 *     cells, $\eta_K$ are the local error indicators and $\eta$ is the global
 *     indicator with $\eta^2 = \sum \eta_K^2$, with here the sum running over
 *     all cells.
 *
 *     For the bottom fraction the same holds: the threshold for coarsening is
 *     computed such that the cells with criterion less than the threshold
 *     together make up for the fraction of the total error specified.
 *
 *     This strategy is more suited for singular functions and error
 *     functionals, but may lead to very slow convergence of the grid
 *     if only few cells are refined in each step.
 *
 *     From the point of view of implementation, this time we really need to
 *     sort the array of criteria.
 *     Just like the other strategy described above, this function only
 *     computes the threshold values and then passes over to @p{refine} and
 *     @p{coarsen}.
 *
 *     A typical value for the fraction of the total error is 0.5.
 *   @end{itemize}
 *
 *   There are other functions relying on different methods to flag
 *   cells for refinement or coarsening. See their documentation for
 *   further information.
 *
 *   For a more thorough discussion of advantages and disadvantages of the
 *   different strategies for refinement, see the paper of R. Becker and
 *   R. Rannacher titled "A Feed-Back Approach to Error Control in Finite
 *   Element Methods: Basic Analysis and Examples".
 *
 *   It is assumed that the criterion is a value in a certain norm over each
 *   element, such that the square of the total error is the sum over the
 *   squares of the criteria on the cells. The criteria shall be positive.
 *
 *   You can suppress coarsening or refining by giving zero as the fraction
 *   for one of the operations.
 *
 * @author Wolfgang Bangerth, 1998, 2000; Thomas Richter, 2000
 */
class GridRefinement
{
  public:
				     /**
				      * Refine the triangulation
				      * according to the given
				      * criteria. The criterion is a
				      * @p{double} value for each cell
				      * which determines which cells
				      * are to be refined by
				      * comparison with the threshold:
				      * if the value for a cell is
				      * larger than the threshold, the
				      * cell is flagged for
				      * refinement. It is your duty to
				      * guarantee that the threshold
				      * value is in a resonable
				      * range. Please note that the
				      * @p{criteria} array may contain
				      * negative values (sometimes,
				      * error estimators are evaluated
				      * in a way which produces
				      * positive and negative values),
				      * but the comparison with
				      * @p{threshold} is done only on
				      * the absolute values of the
				      * criteria.
				      *
				      * The cells are only flagged for
				      * refinement, they are not
				      * actually refined. To do so,
				      * you have to call the
				      * @p{execute_coarsening_and_refinement}
				      * function.
				      *
				      * There are more sophisticated
				      * strategies for mesh
				      * refinement; refer to the
				      * following functions and to the
				      * general doc for this class for
				      * more information.
				      *
				      * Note that this function takes
				      * a vector of @p{float}s, rather
				      * than the usual @p{double}s,
				      * since accuracy is not so much
				      * needed here and saving memory
				      * may be a good goal when using
				      * many cells.
				      */
    template <int dim, typename number>
    static void refine (Triangulation<dim>   &tria,
			const Vector<number> &criteria,
			const double         threshold);

				     /**
				      * Analogue to the @p{refine}
				      * function: flag all cells for
				      * coarsening for which the
				      * absolute value of the
				      * criterion is less than the
				      * given threshold.
				      *
				      * Note that this function takes
				      * a vector of @p{float}s, rather
				      * than the usual @p{double}s,
				      * since accuracy is not so much
				      * needed here and saving memory
				      * may be a good goal when using
				      * many cells.
				      */
    template <int dim, typename number>
    static void coarsen (Triangulation<dim>   &tria,
			 const Vector<number> &criteria,
			 const double         threshold);
    
				     /**
				      * Refine the triangulation by
				      * refining a certain fraction
				      * @p{top_fraction_of_cells} with
				      * the highest error. Likewise
				      * coarsen the fraction
				      * @p{bottom_fraction_of_cells}
				      * with the least error. To
				      * actually perform the
				      * refinement, call
				      * @p{execute_coarsening_and_refinement}.
				      *
				      * @p{fraction_of_cells} shall be
				      * a value between zero and one.
				      *
				      * Refer to the general doc of
				      * this class for more
				      * information.
				      *
				      * Note that this function takes
				      * a vector of @p{float}s, rather
				      * than the usual @p{double}s,
				      * since accuracy is not so much
				      * needed here and saving memory
				      * may be a good goal when using
				      * many cells.
				      */
    template <int dim, typename number>
    static void refine_and_coarsen_fixed_number (Triangulation<dim>   &tria,
						 const Vector<number> &criteria,
						 const double         top_fraction_of_cells,
						 const double         bottom_fraction_of_cells);
    
				     /**
				      * Refine the triangulation by
				      * flagging those cells which
				      * make up a certain
				      * @p{top_fraction} of the total
				      * error.  Likewise, coarsen all
				      * cells which make up only
				      * @p{bottom_fraction}.  To
				      * actually perform the
				      * refinement, call
				      * @p{execute_coarsening_and_refinement}.
				      *
				      * @p{*_fraction} shall be a
				      * values between zero and one.
				      *
				      * Refer to the general doc of
				      * this class for more
				      * information.
				      *
				      * Note that this function takes
				      * a vector of @p{float}s, rather
				      * than the usual @p{double}s,
				      * since accuracy is not so much
				      * needed here and saving memory
				      * may be a good goal when using
				      * many cells.
				      */
    template<int dim, typename number>
    static void refine_and_coarsen_fixed_fraction (Triangulation<dim>   &tria,
						   const Vector<number> &criteria,
						   const double         top_fraction,
						   const double         bottom_fraction);



				     /**
				      * Refine the triangulation by
				      * flagging certain cells to reach
				      * an optimal grid:
				      * We try to minimize the error
				      * multiplied with the number of
				      * cells in the new grid. All cells
				      * with large error indicator are
				      * refined to generate an optimal
				      * grid in the above sense.
				      * We assume that the error in one
				      * cell is reduced to a quarter
				      * after refinement.
				      * The new triangulation has three
				      * new cells for every flagged cell.
				      *
				      * Refer to the general doc of
				      * this class for more
				      * information.
				      *
				      * Note that this function takes
				      * a vector of @p{float}s, rather
				      * than the usual @p{double}s,
				      * since accuracy is not so much
				      * needed here and saving memory
				      * may be a good goal when using
				      * many cells.
				      */
    
    template<int dim, typename number>
    static void refine_and_coarsen_optimize (Triangulation<dim>   &tria,
					     const Vector<number> &criteria);

				     /**
				      * Exception
				      */
    DeclException2 (ExcInvalidVectorSize,
		    int, int,
		    << "The given vector has " << arg1
		    << " elements, but " << arg2 << " were expected.");
				     /**
				      * Exception
				      */
    DeclException0 (ExcInvalidParameterValue);


  private:
    
				     /**
				      * Sorts the vector @p{ind} as an
				      * index vector of @p{a} in
				      * increasing order.  This
				      * implementation of quicksort
				      * seems to be faster than the
				      * STL version and is needed in
				      * @p{refine_and_coarsen_optimize}
				      */
    template<typename number>
    static void qsort_index(const Vector<number>       &a,
			    std::vector<unsigned int>  &ind,
			    int                         l,
			    int                         r);
};



#endif //__deal2__grid_refinement_h

