//----------------------------  sparse_vanka.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  sparse_vanka.h  ---------------------------
#ifndef __deal2__sparse_vanka_h
#define __deal2__sparse_vanka_h


/* Copyright Guido Kanschat, Wolfgang Bangerth 1999, 2000 */


#include <base/smartpointer.h>
#include <base/multithread_info.h>

#include <vector>
#include <map>

template <typename number> class FullMatrix;
template <typename number> class SparseMatrix;
template <typename number> class Vector;

/**
 * Point-wise Vanka preconditioning.
 * This class does Vanka preconditioning  on a point-wise base.
 * Vanka preconditioners are used for saddle point problems like Stoke's
 * problem or problems arising in optimization where Lagrange multiplier
 * occur and let Netwon's matrix have a zero block. With these matrices the
 * application of Jacobi or Gauss-Seidel methods is impossible, because
 * some diagonal elements are zero in the rows of the Lagrange multiplier.
 * The approach of Vanka is to solve a small (usually indefinite) system
 * of equations for each Langrange multiplie variable (we will also call
 * the pressure in Stoke's equation a Langrange multiplier since it
 * can be interpreted as such).
 *
 * Objects of this class are constructed by passing a vector of indices
 * of the degrees of freedom of the Lagrange multiplier. In the actual
 * preconditioning method, these rows are traversed in the order in which
 * the appear in the matrix. Since this is a Gauﬂ-Seidel like procedure,
 * remember to have a good ordering in advance (for transport dominated
 * problems, Cuthill-McKee algorithms are a good means for this, if points
 * on the inflow boundary are chosen as starting points for the renumbering).
 *
 * For each selected degree of freedom, a local system of equations is built
 * by the degree of freedom itself and all other values coupling immediately,
 * i.e. the set of degrees of freedom considered for the local system of
 * equations for degree of freedom #i# is #i# itself and all #j# such that
 * the element #(i,j)# is a nonzero entry in the sparse matrix under 
 * consideration. The elements #(j,i)# are not considered. We now pick all
 * matrix entries from rows and columns out of the set of degrees of freedom
 * just described out of the global matrix and put it into a local matrix,
 * which is subsequently inverted. This system may be of different size for
 * each degree of freedom, depending for example on the local neighborhood of
 * the respective node on a computational grid.
 *
 * The right hand side is built up in the same way, i.e. by copying
 * all entries that coupled with the one under present consideration,
 * but it is augmented by all degrees of freedom coupling with the
 * degrees from the set described above (i.e. the DoFs coupling second
 * order to the present one). The reason for this is, that the local
 * problems to be solved shall have Dirichlet boundary conditions on
 * the second order coupling DoFs, so we have to take them into
 * account but eliminate them before actually solving; this
 * elimination is done by the modification of the right hand side, and
 * in the end these degrees of freedom do not occur in the matrix and
 * solution vector any more at all.
 *
 * This local system is solved and the values are updated into the
 * destination vector.
 *
 * Remark: the Vanka method is a non-symmetric preconditioning method.
 *
 *
 * \subsection{Example of Use}
 * This little example is taken from a program doing parameter optimization.
 * The Lagrange multiplier is the third component of the finite element
 * used. The system is solved by the GMRES method.
 * \begin{verbatim}
 *                        // tag the Lagrange multiplier variable
 *    vector<bool> signature(3);
 *    signature[0] = signature[1] = false;
 *    signature[2] = true;
 *
 *                        // tag all dofs belonging to the
 *                        // Lagrange multiplier
 *    vector<bool> selected_dofs (dof.n_dofs(), false);
 *    DoFTools::extract_dofs(dof, signature, p_select);
 *                        // create the Vanka object
 *    SparseVanka<double> vanka (global_matrix, selected_dofs);
 *
 *                        // create the solver
 *    SolverGMRES<PreconditionedSparseMatrix<double>,
 *                Vector<double> >    gmres(control,memory,504);
 *    
 *			  // solve
 *    gmres.solve (global_matrix, solution, right_hand_side,
 *	           vanka);
 * \end{verbatim}
 *
 *
 * \subsubsection{Implementor's remark}
 * At present, the local matrices are built up such that the degree of freedom
 * associated with the local Lagrange multiplier is the first one. Thus, usually
 * the upper left entry in the local matrix is zero. It is not clear to me (W.B.)
 * whether this might pose some problems in the inversion of the local matrices.
 * Maybe someone would like to check this.
 *
 *
 * \section{On template instantiations}
 *
 * Member functions of this class are either implemented in this file
 * or in a file of the same name with suffix ``.templates.h''. For the
 * most common combinations of the template parameters, instantiations
 * of this class are provided in a file with suffix ``.cc'' in the
 * ``source'' directory. If you need an instantiation that is not
 * listed there, you have to include this file along with the
 * corresponding ``.templates.h'' file and instantiate the respective
 * class yourself.
 *
 * @author Guido Kanschat, Wolfgang Bangerth; 1999, 2000
 */
template<typename number>
class SparseVanka
{
  public:
				     /**
				      * Constructor. Gets the matrix
				      * for preconditioning and a bit
				      * vector with entries #true# for
				      * all rows to be updated. A
				      * reference to this vector will
				      * be stored, so it must persist
				      * longer than the Vanka
				      * object. The same is true for
				      * the matrix.
				      *
				      * The matrix #M# which is passed
				      * here may or may not be the
				      * same matrix for which this
				      * object shall act as
				      * preconditioner. In particular,
				      * it is conceivable that the
				      * preconditioner is build up for
				      * one matrix once, but is used
				      * for subsequent steps in a
				      * nonlinear process as well,
				      * where the matrix changes in
				      * each step slightly.
				      *
				      * If #conserve_mem# is #false#,
				      * then the inverses of the local
				      * systems are computed now, if
				      * the flag is #true#, then they
				      * are computed every time the
				      * preconditioner is
				      * applied. This saves some
				      * memory, but makes
				      * preconditioning very
				      * slow. Note also, that if the
				      * flag is #false#, the the
				      * contents of the matrix #M# at
				      * the time of calling this
				      * constructor are used, while if
				      * the flag is #true#, then the
				      * values in #M# at the time of
				      * preconditioning are used. This
				      * may lead to different results,
				      * obviously, of #M# changes.
				      *
				      * The parameter #n_threads#
				      * determines how many threads
				      * shall be used in parallel when
				      * building the inverses of the
				      * diagonal blocks. This
				      * parameter is ignored if not in
				      * multithreaded mode. By
				      * default, this variable is set
				      * to the value of
				      * #multithread_info.n_default_threads#.
				      */
    SparseVanka(const SparseMatrix<number> &M,
		const vector<bool>         &selected,
		const bool                  conserve_memory = false,
		const unsigned int          n_threads       = multithread_info.n_default_threads);
    
				     /**
				      * Destructor.
				      * Delete all allocated matrices.
				      */
    ~SparseVanka();

				     /**
				      * Do the preconditioning.
				      * This function takes the residual
				      * in #src# and returns the resulting
				      * update vector in #dst#.
				      */
    template<typename number2>
    void vmult (Vector<number2>       &dst,
		const Vector<number2> &src) const;

				     /**
				      * Exception
				      */
    DeclException0 (ExcMatrixNotSquare);
				     /**
				      * Exception
				      */
    DeclException2 (ExcInvalidVectorSize,
		    unsigned int, unsigned int,
		    << "The dimensions of vectors and matrices, "
		    << arg1 << " and " << arg2 << " do not match.");

  protected:
				     /**
				      * Apply the inverses
				      * corresponding to those degrees
				      * of freedom that have a #true#
				      * value in #dof_mask#, to the
				      * #src# vector and move the
				      * result into #dst#. Actually,
				      * only values for allowed
				      * indices are written to #dst#,
				      * so the application of this
				      * function only does what is
				      * announced in the general
				      * documentation if the given
				      * mask sets all values to zero
				      *
				      * The reason for providing the
				      * mask anyway is that in derived
				      * classes we may want to apply
				      * the preconditioner to parts of
				      * the matrix only, in order to
				      * parallelize the
				      * application. Then, it is
				      * important to only write to
				      * some slices of #dst#, in order
				      * to eliminate the dependencies
				      * of threads of each other.
				      *
				      * If a null pointer is passed
				      * instead of a pointer to the
				      * #dof_mask# (as is the default
				      * value), then it is assumed
				      * that we shall work on all
				      * degrees of freedom. This is
				      * then equivalent to calling the
				      * function with a
				      * #vector<bool>(n_dofs,true)#.
				      *
				      * The #vmult# of this class
				      * of course calls this function
				      * with a null pointer
				      */
    template<typename number2>
    void apply_preconditioner (Vector<number2>       &dst,
			       const Vector<number2> &src,
			       const vector<bool>    *dof_mask = 0) const;    
    
  private:
				     /**
				      * Pointer to the matrix.
				      */
    SmartPointer<const SparseMatrix<number> > matrix;
    
				     /**
				      * Conserve memory flag.
				      */
    const bool conserve_mem;

				     /**
				      * Indices of those degrees of
				      * freedom that we shall work on.
				      */
    const vector<bool> &selected;

				     /**
				      * Number of threads to be used
				      * when building the
				      * inverses. Only relevant in
				      * multithreaded mode.
				      */
    const unsigned int n_threads;
    
				     /**
				      * Array of inverse matrices,
				      * one for each degree of freedom.
				      * Only those elements will be used
				      * that are tagged in #selected#.
				      */
    mutable vector<SmartPointer<FullMatrix<float> > > inverses;

				     /**
				      * Compute the inverses of all
				      * selected diagonal elements.
				      */
    void compute_inverses ();

				     /**
				      * Compute the inverses at
				      * positions in the range
				      * #[begin,end)#. In
				      * non-multithreaded mode,
				      * #compute_inverses()# calls
				      * this function with the whole
				      * range, but in multithreaded
				      * mode, several copies of this
				      * function are spawned.
				      */
    void compute_inverses (const unsigned int begin,
			   const unsigned int end);
    
				     /**
				      * Compute the inverse of the
				      * block located at position
				      * #row#. Since the vector is
				      * used quite often, it is
				      * generated only once in the
				      * caller of this function and
				      * passed to this function which
				      * first clears it. Reusing the
				      * vector makes the process
				      * significantly faster than in
				      * the case where this function
				      * re-creates it each time.
				      */
    void compute_inverse (const unsigned int    row,
			  vector<unsigned int> &local_indices);

};



/**
 * Block version of the sparse Vanka preconditioner. This class
 * divides the matrix into blocks and works on the diagonal blocks
 * only, which of course reduces the efficiency as preconditioner, but
 * is perfectly parallelizable. The constructor takes a parameter into
 * how many blocks the matrix shall be subdivided and then lets the
 * underlying class do the work. Division of the matrix is done in
 * several ways which are described in detail below.
 *
 * This class is probably useless if you don't have a multiprocessor
 * system, since then the amount of work per preconditioning step is
 * the same as for the #SparseVanka# class, but preconditioning
 * properties are worse. On the other hand, if you have a
 * multiprocessor system, the worse preconditioning quality (leading
 * to more iterations of the linear solver) usually is well balanced
 * by the increased speed of application due to the parallelization,
 * leading to an overall decrease in elapsed wall-time for solving
 * your linear system. It should be noted that the quality as
 * preconditioner reduces with growing number of blocks, so there may
 * be an optimal value (in terms of wall-time per linear solve) for
 * the number of blocks.
 *
 * To facilitate writing portable code, if the number of blocks into
 * which the matrix is to be subdivided, is set to one, then this
 * class acts just like the #SparseVanka# class. You may therefore
 * want to set the number of blocks equal to the number of processors
 * you have.
 *
 * Note that the parallelization is done if #deal.II# was configured
 * for multithread use and that the number of threads which is spawned
 * equals the number of blocks. This is reasonable since you will not
 * want to set the number of blocks unnecessarily large, since, as
 * mentioned, this reduces the preconditioning properties.
 *
 *
 * \subsection{Splitting the matrix into blocks}
 *
 * Splitting the matrix into blocks is always done in a way such that
 * the blocks are not necessarily of equal size, but such that the
 * number of selected degrees of freedom for which a local system is
 * to be solved is equal between blocks. The reason for this strategy
 * to subdivision is load-balancing for multithreading. There are
 * several possibilities to actually split the matrix into blocks,
 * which are selected by the flag #blocking_strategy# that is passed
 * to the constructor. By a block, we will in the sequel denote a list
 * of indices of degrees of freedom; the algorithm will work on each
 * block separately, i.e. the solutions of the local systems
 * corresponding to a degree of freedom of one block will only be used
 * to update the degrees of freedom belonging to the same block, but
 * never to update degrees of freedoms of other blocks. A block can be
 * a consecutive list of indices, as in the first alternative below,
 * or a nonconsecutive list of indices. Of course, we assume that the
 * intersection of each two blocks is empty and that the union of all
 * blocks equals the interval #[0,N)#, where #N# is the number of
 * degrees of freedom of the system of equations.
 *
 * \begin{itemize}
 * \item #index_intervals#:
 *    Here, we chose the blocks to be intervals #[a_i,a_{i+1})#,
 *    i.e. consecutive degrees of freedom are usually also within the
 *    same block. This is a reasonable strategy, if the degrees of
 *    freedom have, for example, be re-numbered using the
 *    Cuthill-McKee algorithm, in which spatially neighboring degrees
 *    of freedom have neighboring indices. In that case, coupling in
 *    the matrix is usually restricted to the vicinity of the diagonal
 *    as well, and we can simply cut the matrix into blocks.
 *
 *    The bounds of the intervals, i.e. the #a_i# above, are chosen
 *    such that the number of degrees of freedom on which we shall
 *    work (i.e. usually the degrees of freedom corresponding to
 *    Lagrange multipliers) is about the same in each block; this does
 *    not mean, however, that the sizes of the blocks are equal, since
 *    the blocks also comprise the other degrees of freedom for which
 *    no local system is solved. In the extreme case, consider that
 *    all Lagrange multipliers are sorted to the end of the range of
 *    DoF indices, then the first block would be very large, since it
 *    comprises all other DoFs and some Lagrange multipliers, while
 *    all other blocks are rather small and comprise only Langrange
 *    multipliers. This strategy therefore does not only depend on the
 *    order in which the Lagrange DoFs are sorted, but also on the
 *    order in which the other DoFs are sorted. It is therefore
 *    necessary to note that this almost renders the capability as
 *    preconditioner useless if the degrees of freedom are numbered by
 *    component, i.e. all Lagrange multipliers en bloc.
 *
 * \item #adaptive#: This strategy is a bit more clever in cases where
 *    the Langrange DoFs are clustered, as in the example above. It
 *    works as follows: it first groups the Lagrange DoFs into blocks,
 *    using the same strategy as above. However, instead of grouping
 *    the other DoFs into the blocks of Lagrange DoFs with nearest DoF
 *    index, it decides for each non-Lagrange DoF to put it into the
 *    block of Lagrange DoFs which write to this non-Lagrange DoF most
 *    often. This makes it possible to even sort the Lagrange DoFs to
 *    the end and still associate spatially neighboring non-Lagrange
 *    DoFs to the same blocks where the respective Lagrange DoFs are,
 *    since they couple to each other while spatially distant DoFs
 *    don't couple.
 *
 *    The additional computational effort to sorting the non-Lagrange
 *    DoFs is not very large compared with the inversion of the local
 *    systems and applying the preconditioner, so this strategy might
 *    be reasonable if you want to sort your degrees of freedom by
 *    component. If the degrees of freedom are not sorted by
 *    component, the results of the both strategies outlined above
 *    does not differ much. However, unlike the first strategy, the
 *    performance of the second strategy does not deteriorate if the
 *    DoFs are renumbered by component.
 * \end{itemize}
 *
 *
 * \subsection{Typical results}
 *
 * As a prototypical test case, we use a nonlinear problem from
 * optimization, which leads to a series of saddle point problems,
 * each of which is solved using GMRES with Vanka as
 * preconditioner. The equation had approx. 850 degrees of
 * freedom. With the non-blocked version #SparseVanka# (or
 * #SparseBlockVanka# with #n_blocks==1#), the following numbers of
 * iterations is needed to solver the linear system in each nonlinear
 * step:
 * \begin{verbatim}
 *   101 68 64 53 35 21
 * \end{verbatim}
 *
 * With four blocks, we need the following numbers of iterations
 * \begin{verbatim}
 *   124 88 83 66 44 28
 * \end{verbatim}
 * As can be seen, more iterations are needed. However, in terms of
 * computing time, the first version needs 72 seconds wall time (and
 * 79 seconds CPU time, which is more than wall time since some other
 * parts of the program were parallelized as well), while the second
 * version needed 53 second wall time (and 110 seconds CPU time) on a
 * four processor machine. The total time is in both cases dominated
 * by the linear solvers. In this case, it is therefore worth while
 * using the blocked version of the preconditioner if wall time is
 * more important than CPU time.
 *
 * The results with the block version above were obtained with the
 * first blocking strategy and the degrees of freedom were not
 * numbered by component. Using the second strategy does not much
 * change the numbers of iterations (at most by one in each step) and
 * they also do not change when the degrees of freedom are sorted
 * by component, while the first strategy significantly deteriorated.
 *
 * @author Wolfgang Bangerth, 2000
 */
template<typename number>
class SparseBlockVanka : public SparseVanka<number>
{
  public:
				     /**
				      * Enumeration of the different
				      * methods by which the DoFs are
				      * distributed to the blocks on
				      * which we are to work.
				      */
    enum BlockingStrategy {
	  index_intervals, adaptive
    };
    
				     /**
				      * Constructor. Pass all
				      * arguments except for
				      * #n_blocks# to the base class.
				      */
    SparseBlockVanka (const SparseMatrix<number> &M,
		      const vector<bool>         &selected,
		      const unsigned int          n_blocks,
		      const BlockingStrategy      blocking_strategy,
		      const bool                  conserve_memory = false,
		      const unsigned int          n_threads       = multithread_info.n_default_threads);

				     /**
				      * Apply the preconditioner.
				      */
    template<typename number2>
    void vmult (Vector<number2>       &dst,
		     const Vector<number2> &src) const;
    
  private:
				     /**
				      * Store the number of blocks.
				      */
    const unsigned int n_blocks;

				     /**
				      * In this field, we precompute
				      * for each block which degrees
				      * of freedom belong to it. Thus,
				      * if #dof_masks[i][j]==true#,
				      * then DoF #j# belongs to block
				      * #i#. Of course, no other
				      * #dof_masks[l][j]# may be
				      * #true# for #l!=i#. This
				      * computation is done in the
				      * constructor, to avoid
				      * recomputing each time the
				      * preconditioner is called.
				      */
    vector<vector<bool> > dof_masks;

				     /**
				      * Compute the contents of the
				      * field #dof_masks#. This
				      * function is called from the
				      * constructor.
				      */
    void compute_dof_masks (const SparseMatrix<number> &M,
			    const vector<bool>         &selected,
			    const BlockingStrategy      blocking_strategy);
};


#endif
