//----------------------------  sparse_direct.h  ---------------------------
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
//----------------------------  sparse_direct.h  ---------------------------
#ifndef __deal2__sparse_direct_h
#define __deal2__sparse_direct_h



#include <base/config.h>
#include <base/exceptions.h>
#include <base/subscriptor.h>
#include <base/thread_management.h>
#include <lac/vector.h>
#include <lac/sparse_matrix.h>



/**
 * This class provides an interface to the sparse direct solver MA27
 * from the Harwell Subroutine Library. MA27 is a direct solver
 * specialized for sparse symmetric indefinite systems of linear equations and
 * uses a modified form of Gauss elimination. It is included in the
 * Harwell Subroutine Library (see
 * @url{http://www.cse.clrc.ac.uk/Activity/HSL}) and is written in
 * Fortran. The present class only transforms the data stored in
 * @ref{SparseMatrix} objects into the form which is required by the
 * functions resembling MA27, calls these Fortran functions, and
 * interprets some of the returned values indicating error codes,
 * etc. It also manages allocation of the right amount of temporary
 * storage required by these functions.
 *
 * For a description of the steps necessary for the installation of
 * HSL subroutines, read the section on external libraries in the
 * deal.II ReadMe file.
 *
 * @sect3{Interface and Method}
 *
 * For the meaning of the three functions @p{initialize},
 * @p{factorize}, and @p{solve}, as well as for the method used in
 * MA27, please see the documentation of these functions, which can be
 * obtained from @url{http://www.cse.clrc.ac.uk/Activity/HSL}. In
 * practice, one will most often call the second @p{solve} function,
 * which solves the linear system for a given right hand sidem but one
 * can as well call the three functions separately if, for example,
 * one would like to solve the same matrix for several right hand side
 * vectors; the MA27 solver can do this efficiently, as it computes a
 * decomposition of the matrix, so that subsequent solves only amount
 * to a forward-backward substitution which is significantly less
 * costly than the decomposition process.
 *
 *
 * @sect3{Parameters to the constructor}
 *
 * The constructor of this class takes several arguments. The meaning
 * is the following: the MA27 functions require the user to allocate
 * and pass a certain amount of memory for temporary variables or for
 * data to be passed to subsequent functions. The sizes of these
 * arrays are denoted by the variables @p{LIW1}, @p{LIW2}, and @p{LA},
 * where @p{LIW1} denotes the size of the @p{IW} array in the call to
 * @p{MA27A}, while @p{LIW2} is the array size in the call to
 * @p{MA27B}. The documentation of the MA27 functions gives ways to
 * obtain estimates for their values, e.g. by evaluating values
 * returned by functions called before. However, the documentation
 * only states that the values have to be @em{at least as large} as
 * the estimates, a hint that is not very useful oftentimes (in my
 * humble opinion, the lack of dynamic memory allocation mechanism is
 * a good reason not to program in Fortran 77...).
 *
 * In our experience, it is often necessary to go beyond the proposed
 * values (most often for @p{LA}, but also for @p{LIW1}). The first
 * three parameters of the constructor denote by which factor the
 * initial estimates shall be increased. The default values are 1.2
 * (the documentation recommends this value, 1, and 1.5, values which
 * have often worked for us. Note that the value of @p{LIW} is only
 * changed in the second call if the recommended value times
 * @p{LIW_factor_2} is larger than the array size already is from the
 * call to @p{MA27A}; otherwise, @p{LIW_factor_2} is ignored.
 *
 * If the values thus constructed fail to work, we try to restart the
 * called function with larger values until the calls succeed. The
 * second triple of values passed to the constructor denotes by which
 * factor we shall increase the array sizes. If the increment factors
 * are less than or equal to one, then we only try to call the
 * respective calls to the functions once and abort by throwing an
 * error. Note that the @p{MA27C} function writes out an error message
 * if the value of @p{LA} is too small and gives an indication to
 * which size it should be increased. However, most often the
 * indicated value is far too small and can not be relied upon.
 *
 *
 * @sect3{Note on parallelization}
 *
 * Due to the use of global variables through COMMON blocks, the calls
 * to the sparse direct solver routines is not multithreading-capable,
 * i.e. at each time there may only be one call to these functions
 * active. You have to synchronise your calls to the functions
 * provided by this class using mutexes (see the @ref{Threads}
 * namespace for such classes) to avoid multiple active calls at the
 * same time if you use multithreading. Since you may use this class
 * in different parts of your program, and may not want to use a
 * global variable for locking, this class has a lock as static member
 * variable, which may be accessed using the
 * @p{get_synchronisation_lock} function. Note however, that this
 * class does not perform the synchronisation for you within its
 * member functions. The reason is that you will usually want to
 * synchronise over the calls to @p{initialize} and @p{factorize},
 * since there should probably not be a call to one of these function
 * with another matrix between the calls for one matrix. (The author
 * does not really know whether this is true, but it is probably safe
 * to assume that.) Since such cross-function synchronisation can only
 * be performed from outside, it is left to the user of this class to
 * do so.
 *
 * @author Wolfgang Bangerth, 2000, 2001
 */
class SparseDirectMA27 : public Subscriptor
{
  public:
				     /**
				      * Constructor. See the
				      * documentation of this class
				      * for the meaning of the
				      * parameters to this function.
				      */
    SparseDirectMA27 (const double LIW_factor_1          = 1.2,
		      const double LIW_factor_2          = 1,
		      const double LA_factor             = 1.5,
		      const double LIW_increase_factor_1 = 1.2,
		      const double LIW_increase_factor_2 = 1.2,
		      const double LA_increase_factor    = 1.2,
		      const bool   suppress_output       = true);

				     /**
				      * Initialize some data
				      * structures. This function
				      * computes symbolically some
				      * information based on the
				      * sparsity pattern, but does not
				      * actually use the values of the
				      * matrix, so only the sparsity
				      * pattern has to be passed as
				      * argument.
				      */
    void initialize (const SparsityPattern &sparsity_pattern);

				     /**
				      * Actually factorize the
				      * matrix. This function may be
				      * called multiple times for
				      * different matrices, after the
				      * object of this class has been
				      * initialized for a certain
				      * sparsity pattern. You may
				      * therefore save some computing
				      * time if you want to invert
				      * several matrices with the same
				      * sparsity pattern. However,
				      * note that the bulk of the
				      * computing time is actually
				      * spent in the factorization, so
				      * this functionality may not
				      * always be of large benefit.
				      *
				      * If the initialization step has
				      * not been performed yet, then
				      * the @p{initialize} function is
				      * called at the beginning of
				      * this function.
				      */
    void factorize (const SparseMatrix<double> &matrix);

				     /**
				      * Solve for a certain right hand
				      * side vector. This function may
				      * be called multiple times for
				      * different right hand side
				      * vectors after the matrix has
				      * been factorized. This yields a
				      * big saving in computing time,
				      * since the actual solution is
				      * fast, compared to the
				      * factorization of the matrix.
				      *
				      * The solution will be returned
				      * in place of the right hand
				      * side vector.
				      *
				      * If the factorization has not
				      * happened before, strange
				      * things will happen. Note that
				      * we can't actually call the
				      * @p{factorize} function from
				      * here if it has not yet been
				      * called, since we have no
				      * access to the actual matrix.
				      */
    void solve (Vector<double> &rhs_and_solution) const;

				     /**
				      * Call the three functions above
				      * in that order, i.e. perform
				      * the whole solution process for
				      * the given right hand side
				      * vector.
				      *
				      * The solution will be returned
				      * in place of the right hand
				      * side vector.
				      */
    void solve (const SparseMatrix<double> &matrix,
		Vector<double>             &rhs_and_solution);

				     /**
				      * Return an estimate of the
				      * memory used by this class.
				      */
    unsigned int memory_consumption () const;
    
				     /**
				      * Get a reference to the
				      * synchronisation lock which can
				      * be used for this class. See
				      * the general description of
				      * this class for more
				      * information.
				      */
    Threads::ThreadMutex & get_synchronisation_lock () const;

				     /**
				      * Exception.
				      */
    DeclException1 (ExcMA27AFailed,
		    int,
		    << "The function MA27A failed with an exit code of " << arg1);
				     /**
				      * Exception.
				      */
    DeclException1 (ExcMA27BFailed,
		    int,
		    << "The function MA27B failed with an exit code of " << arg1);
				     /**
				      * Exception.
				      */
    DeclException1 (ExcMA27CFailed,
		    int,
		    << "The function MA27C failed with an exit code of " << arg1);
				     /**
				      * Exception
				      */
    DeclException0 (ExcInitializeAlreadyCalled);
				     /**
				      * Exception
				      */
    DeclException0 (ExcFactorizeNotCalled);
				     /**
				      * Exception
				      */
    DeclException0 (ExcDifferentSparsityPatterns);
    
  private:
				     /**
				      * Store the three values passed
				      * to the cinstructor. See the
				      * documentation of this class
				      * for the meaning of these
				      * variables.
				      */
    const double LIW_factor_1;
    const double LIW_factor_2;
    const double LA_factor;

				     /**
				      * Increase factors in case a
				      * call to a function fails.
				      */
    const double LIW_increase_factor_1;
    const double LIW_increase_factor_2;
    const double LA_increase_factor;
    
				     /**
				      * Flags storing whether the
				      * first two functions have
				      * already been called.
				      */
    bool initialize_called;
    bool factorize_called;

				     /**
				      * Store a pointer to the
				      * sparsity pattern, to make sure
				      * that we use the same thing for
				      * all calls.
				      */
    SmartPointer<const SparsityPattern> sparsity_pattern;
    
				     /**
				      * Number of nonzero elements in
				      * the sparsity pattern on and
				      * above the diagonal.
				      */
    unsigned int n_nonzero_elements;
    
				     /**
				      * Arrays holding row and column
				      * indices.
				      */
    std::vector<unsigned int> row_numbers;
    std::vector<unsigned int> column_numbers;

				     /**
				      * Array to hold the matrix
				      * elements, and later the
				      * elements of the factors.
				      */
    std::vector<double>       A;

				     /**
				      * Length of the @p{A} array.
				      */
    unsigned int         LA;
    
				     /**
				      * Scratch arrays and variables
				      * used by the MA27 functions. We
				      * keep to the names introduced
				      * in the documentation of these
				      * functions, in all uppercase
				      * letters as is usual in
				      * Fortran.
				      */
    unsigned int LIW;
    std::vector<unsigned int> IW;
    std::vector<unsigned int> IKEEP;
    std::vector<unsigned int> IW1;
    
    unsigned int NSTEPS;
    unsigned int MAXFRT;

				     /**
				      * Two values that live inside a
				      * COMMON block of the Fortran
				      * code and are mirrored at these
				      * locations. They are used to
				      * transport information about
				      * the required length of arrays
				      * from the Fortran functions to
				      * the outside world.
				      */
    unsigned int NRLNEC;
    unsigned int NIRNEC;
    
				     /**
				      * Flag indicating the level of
				      * output desired and returning
				      * error values if error occured.
				      */
    int IFLAG;

				     /**
				      * Mutex for synchronising access
				      * to this class.
				      */
    static Threads::ThreadMutex synchronisation_lock;
    
				     /**
				      * Fill the @p{A} array from the
				      * symmetric part of the given
				      * matrix.
				      */
    void fill_A (const SparseMatrix<double> &matrix);
};






/**
 * This class provides an interface to the sparse direct solver MA47
 * from the Harwell Subroutine Library. MA47 is a direct solver
 * specialized for sparse symmetric indefinite systems of linear equations and
 * uses a frontal elimination method. It is included in the Harwell
 * Subroutine Library (see
 * @url{http://www.cse.clrc.ac.uk/Activity/HSL}) and is written in
 * Fortran. The present class only transforms the data stored in
 * @ref{SparseMatrix} objects into the form which is required by the
 * functions resembling MA47, calls these Fortran functions, and
 * interprets some of the returned values indicating error codes,
 * etc. It also manages allocation of the right amount of temporary
 * storage required by these functions.
 *
 *
 * @sect3{Interface and Method}
 *
 * For the meaning of the three functions @p{initialize},
 * @p{factorize}, and @p{solve}, as well as for the method used in
 * MA47, please see the documentation of these functions, which can be
 * obtained from @url{http://www.cse.clrc.ac.uk/Activity/HSL}. In
 * practice, one will most often call the second @p{solve} function,
 * which solves the linear system for a given right hand sidem but one
 * can as well call the three functions separately if, for example,
 * one would like to solve the same matrix for several right hand side
 * vectors; the MA47 solver can do this efficiently, as it computes a
 * decomposition of the matrix, so that subsequent solves only amount
 * to a forward-backward substitution which is significantly less
 * costly than the decomposition process.
 *
 *
 * @sect3{Parameters to the constructor}
 *
 * The constructor of this class takes several arguments. Their
 * meaning is equivalent to those of the constructor of the
 * @ref{SparseDirectMA27} class; see there for more information.
 *
 *
 * @sect3{Note on parallelization}
 *
 * Due to the use of global variables through COMMON blocks, the calls
 * to the sparse direct solver routines is not multithreading-capable,
 * i.e. at each time there may only be one call to these functions
 * active. You have to synchronise your calls to the functions
 * provided by this class using mutexes (see the @ref{Threads}
 * namespace for such classes) to avoid multiple active calls at the
 * same time if you use multithreading. Since you may use this class
 * in different parts of your program, and may not want to use a
 * global variable for locking, this class has a lock as static member
 * variable, which may be accessed using the
 * @p{get_synchronisation_lock} function. Note however, that this
 * class does not perform the synchronisation for you within its
 * member functions. The reason is that you will usually want to
 * synchronise over the calls to @p{initialize} and @p{factorize},
 * since there should probably not be a call to one of these function
 * with another matrix between the calls for one matrix. (The author
 * does not really know whether this is true, but it is probably safe
 * to assume that.) Since such cross-function synchronisation can only
 * be performed from outside, it is left to the user of this class to
 * do so.
 *
 * @author Wolfgang Bangerth, 2000, 2001
 */
class SparseDirectMA47 : public Subscriptor
{
  public:
				     /**
				      * Constructor. See the
				      * documentation of this class
				      * for the meaning of the
				      * parameters to this function.
				      *
				      * This function already calls
				      * the initialization function
				      * @p{MA47ID} to set up some
				      * values.
				      */
    SparseDirectMA47 (const double LIW_factor_1 = 1.4,
		      const double LIW_factor_2 = 1,
		      const double LA_factor    = 3,
		      const double LIW_increase_factor_1 = 1.2,
		      const double LIW_increase_factor_2 = 1.2,
		      const double LA_increase_factor    = 1.2,
		      const bool   suppress_output       = true);

				     /**
				      * Initialize some data
				      * structures. This function
				      * computes symbolically some
				      * information based on the
				      * sparsity pattern, but does not
				      * actually use the values of the
				      * matrix, so only the sparsity
				      * pattern has to be passed as
				      * argument.
				      *
				      * Since the MA47 solver requires
				      * us to omit zero-entries in the
				      * matrix (even if they are in
				      * the sparsity pattern), we have
				      * to actually use the matrix
				      * here, as opposed to the MA27
				      * solver that only required the
				      * sparsity pattern.
				      */
    void initialize (const SparseMatrix<double> &matrix);

				     /**
				      * Actually factorize the
				      * matrix. Unlike for the MA27
				      * solver, this function may not
				      * be called multiple times for
				      * different matrices, since we
				      * have eliminated entries from
				      * the sparsity pattern where
				      * matrix entries happen to be
				      * zero. Since this is likely to
				      * change between matrices
				      * although they have the same
				      * sparsity pattern.
				      *
				      * If the initialization step has
				      * not been performed yet, then
				      * the @p{initialize} function is
				      * called at the beginning of
				      * this function.
				      */
    void factorize (const SparseMatrix<double> &matrix);

				     /**
				      * Solve for a certain right hand
				      * side vector. This function may
				      * be called multiple times for
				      * different right hand side
				      * vectors after the matrix has
				      * been factorized. This yields a
				      * big saving in computing time,
				      * since the actual solution is
				      * fast, compared to the
				      * factorization of the matrix.
				      *
				      * The solution will be returned
				      * in place of the right hand
				      * side vector.
				      *
				      * If the factorization has not
				      * happened before, strange
				      * things will happen. Note that
				      * we can't actually call the
				      * @p{factorize} function from
				      * here if it has not yet been
				      * called, since we have no
				      * access to the actual matrix.
				      */
    void solve (Vector<double> &rhs_and_solution);

				     /**
				      * Call the three functions above
				      * in that order, i.e. perform
				      * the whole solution process for
				      * the given right hand side
				      * vector.
				      *
				      * The solution will be returned
				      * in place of the right hand
				      * side vector.
				      */
    void solve (const SparseMatrix<double> &matrix,
		Vector<double>             &rhs_and_solution);

				     /**
				      * Return an estimate of the
				      * memory used by this class.
				      */
    unsigned int memory_consumption () const;

				     /**
				      * Get a reference to the
				      * synchronisation lock which can
				      * be used for this class. See
				      * the general description of
				      * this class for more
				      * information.
				      */
    Threads::ThreadMutex & get_synchronisation_lock () const;
    
				     /**
				      * Exception.
				      */
    DeclException1 (ExcMA47AFailed,
		    int,
		    << "The function MA47A failed with an exit code of " << arg1);
				     /**
				      * Exception.
				      */
    DeclException1 (ExcMA47BFailed,
		    int,
		    << "The function MA47B failed with an exit code of " << arg1);
				     /**
				      * Exception.
				      */
    DeclException1 (ExcMA47CFailed,
		    int,
		    << "The function MA47C failed with an exit code of " << arg1);
				     /**
				      * Exception
				      */
    DeclException0 (ExcInitializeAlreadyCalled);
				     /**
				      * Exception
				      */
    DeclException0 (ExcFactorizeNotCalled);
				     /**
				      * Exception
				      */
    DeclException0 (ExcCantFactorizeAgain);
				     /**
				      * Exception
				      */
    DeclException0 (ExcDifferentMatrices);
    
  private:
				     /**
				      * Store the three values passed
				      * to the cinstructor. See the
				      * documentation of this class
				      * for the meaning of these
				      * variables.
				      */
    const double LIW_factor_1;
    const double LIW_factor_2;
    const double LA_factor;

				     /**
				      * Increase factors in case a
				      * call to a function fails.
				      */
    const double LIW_increase_factor_1;
    const double LIW_increase_factor_2;
    const double LA_increase_factor;

				     /**
				      * Flags storing whether the
				      * first two functions have
				      * already been called.
				      */
    bool initialize_called;
    bool factorize_called;

				     /**
				      * Store a pointer to the matrix,
				      * to make sure that we use the
				      * same thing for all calls.
				      */
    SmartPointer<const SparseMatrix<double> > matrix;
    
				     /**
				      * Number of nonzero elements in
				      * the sparsity pattern on and
				      * above the diagonal.
				      */
    unsigned int n_nonzero_elements;

				     /**
				      * Control values set by @p{MA47ID}.
				      */
    double       CNTL[2];
    unsigned int ICNTL[7];

				     /**
				      * Info field filled by the MA47
				      * functions and (partially) used
				      * for subsequent MA47 calls.
				      */
    int          INFO[24];

				     /**
				      * Arrays holding row and column
				      * indices.
				      */
    std::vector<unsigned int> row_numbers;
    std::vector<unsigned int> column_numbers;

				     /**
				      * Array to hold the matrix
				      * elements, and later the
				      * elements of the factors.
				      */
    std::vector<double>       A;

				     /**
				      * Length of the @p{A} array.
				      */
    unsigned int         LA;
    
				     /**
				      * Scratch arrays and variables
				      * used by the MA47 functions. We
				      * keep to the names introduced
				      * in the documentation of these
				      * functions, in all uppercase
				      * letters as is usual in
				      * Fortran.
				      */
    unsigned int LIW;
    std::vector<unsigned int> IW;
    std::vector<unsigned int> KEEP;
    std::vector<unsigned int> IW1;

				     /**
				      * Mutex for synchronising access
				      * to this class.
				      */
    static Threads::ThreadMutex synchronisation_lock;
    
				     /**
				      * Fill the @p{A} array from the
				      * symmetric part of the given
				      * matrix.
				      */
    void fill_A (const SparseMatrix<double> &matrix);
};




#endif
