/*----------------------------   solver_selector.h     ---------------------------*/
/*      $Id$                 */
/*                Ralf Hartmann, University of Heidelberg                         */
#ifndef __solver_selector_H
#define __solver_selector_H
/*----------------------------   solver_selector.h     ---------------------------*/

#include <base/smartpointer.h>
#include <lac/solver.h>
#include <lac/vector.h>
#include <lac/vector_memory.h>
#include <lac/solver_control.h>
#include <lac/solver_cg.h>
#include <lac/solver_bicgstab.h>
#include <lac/solver_gmres.h>
#include <lac/vector_memory.h>
#include <lac/solver_richardson.h>
#include <lac/precondition.h>
#include <lac/forward_declarations.h>


/**
 * By calling the #solve# function of this #SolverSelector#, it selects
 * the #solve# function of that #Solver# that was specified in the constructor
 * of this class.
 *
 * \subsection{Usage}
 * The simplest use of this class is the following:
 * \begin{verbatim}
 *                                  // generate a #SolverControl# and
 *                                  // a #VectorMemory#
 * SolverControl control;
 * VectorMemory<Vector<double> > memory;
 *                                  // Line 3:
 *                                  //
 *                                  // generate a #SolverSelector# that
 *                                  // calls the #SolverCG#
 * SolverSelector<SparseMatrix<double>, Vector<double> > 
 *   solver_selector("cg", control, memory);
 *                                  // generate e.g. a #PreconditionRelaxation#
 * PreconditionRelaxation<SparseMatrix<double>, Vector<double> >
 *   preconditioning(A, &SparseMatrix<double>
 *                   ::template precondition_SSOR<double>,0.8);
 *                                  // call the #solve# function with this
 *                                  // preconditioning as last argument
 * solver_selector.solve(A,x,b,preconditioning);
 * \end{verbatim}
 * But the full usefulness of the #SolverSelector# class is not clear
 * until the presentation of the following
 * example that assumes the user using the #ParameterHandler# class and having
 * declared a "solver" entry, e.g. with
 * \begin{verbatim}
 * Parameter_Handler prm;
 * prm.declare_entry ("solver", "none",
 *                    Patterns::Sequence(SolverSelector::get_solver_names()));
 * ...
 * \end{verbatim}
 * Assuming that in the users parameter file there exists the line
 * \begin{verbatim}
 * set solver = cg
 * \end{verbatim}
 * then `Line 3' of the above example reads
 * \begin{verbatim}
 * SolverSelector<SparseMatrix<double>, Vector<double> > 
 *   solver_selector(prm.get("solver"), control, memory);
 * \end{verbatim}
 *
 *
 * If at some time there exists a new solver "xyz" then the user does not need
 * to change his program. Only in the implementation of the #SolverSelector#
 * the calling of this solver has to be added and each user with program lines
 * quoted above only needs to 'set solver = xyz' in his parameter file to get
 * access to that new solver.  :-)
 *
 * (By the way, thanks to Wolfgang for implementing the #ParameterHandler#.)
 * 
 * @author Ralf Hartmann, 1999
 */
template <class Matrix = SparseMatrix<double>,
          class Vector = Vector<double> >
class SolverSelector
{
  public:

				     /**
				      * Constructor. Takes the #SolverName#,
				      * the #SolverControl# and the 
				      * #VectorMemory# as argument.
				      */
    SolverSelector(const string solvername,
		   SolverControl &control,
		   VectorMemory<Vector> &vectorm);
  
				     /**
				      * Destructor
				      */
    ~SolverSelector();

				     /**
				      * Solver procedure. Calls the #solve#
				      * function
				      * of the #solver# whose #SolverName#
				      * was specified in the constructor.
				      * 
				      */
    template<class Preconditioner>
    typename Solver<Matrix,Vector>::ReturnState solve(const Matrix &A,
						      Vector &x,
						      const Vector &b,
						      const Preconditioner &precond) const;
    
				     /**
				      * Set the additional data. For more
				      * info see the #Solver# class.
				      */
    void set_data(const typename SolverRichardson<Matrix,Vector>
		  ::AdditionalData &data);

				     /**
				      * Set the additional data. For more
				      * info see the #Solver# class.
				      */
    void set_data(const typename SolverCG<Matrix,Vector>
		  ::AdditionalData &data);

				     /**
				      * Set the additional data. For more
				      * info see the #Solver# class.
				      */
    void set_data(const typename SolverBicgstab<Matrix,Vector>
		  ::AdditionalData &data); 

				     /**
				      * Set the additional data. For more
				      * info see the #Solver# class.
				      */
    void set_data(const typename SolverGMRES<Matrix,Vector>
		  ::AdditionalData &data);

				     /**
				      * Get the names of all implemented
				      * solvers.
				      */
    static string get_solver_names ();

				     /**
				      * Exception.
				      */
    DeclException1 (ExcSolverDoesNotExist,
		    string, << "Solver " << arg1 << " does not exist. Use one of "
		    << endl << get_solver_names());

  protected:
				     /**
				      * Stores the Name of the solver.
				      */
    string solver_name;
    
				     /**
				      * Stores the #SolverControl# that
				      * is needed in the constructor of
				      * each #Solver# class.
				      */
    SmartPointer<SolverControl> control;

				     /**
				      * Stores the #VectorMemory# that
				      * is needed in the constructor of
				      * each #Solver# class.
				      */
    SmartPointer<VectorMemory<Vector> > vector_memory;

  private:
				     /**
				      * Stores the additional data.
				      */
    typename SolverRichardson<Matrix,Vector>::AdditionalData richardson_data;

				     /**
				      * Stores the additional data.
				      */
    typename SolverCG<Matrix,Vector>::AdditionalData cg_data;

				     /**
				      * Stores the additional data.
				      */
    typename SolverBicgstab<Matrix,Vector>::AdditionalData bicgstab_data;

				     /**
				      * Stores the additional data.
				      */
    typename SolverGMRES<Matrix,Vector>::AdditionalData gmres_data;
};



/* --------------------- Inline and template functions ------------------- */


template <class Matrix, class Vector>
SolverSelector<Matrix, Vector>::SolverSelector(string solver_name,
					       SolverControl &control,
					       VectorMemory<Vector> &vector_memory) :
		solver_name(solver_name),
		control(&control),
		vector_memory(&vector_memory)
{};



template <class Matrix, class Vector>
SolverSelector<Matrix, Vector>::~SolverSelector()
{};



template <class Matrix, class Vector>
template<class Preconditioner>
typename Solver<Matrix,Vector>::ReturnState 
SolverSelector<Matrix, Vector>::solve(const Matrix &A,
				      Vector &x,
				      const Vector &b,
				      const Preconditioner &precond) const
{
  if (solver_name=="richardson")
    {
      SolverRichardson<Matrix,Vector> solver(*control,*vector_memory,
					     richardson_data);
      return solver.solve(A,x,b,precond);
    }       
  else if (solver_name=="cg")
    {
      SolverCG<Matrix,Vector> solver(*control,*vector_memory,
				     cg_data);
      return solver.solve(A,x,b,precond);
    }
  else if (solver_name=="bicgstab")
    {
      SolverBicgstab<Matrix,Vector> solver(*control,*vector_memory,
					   bicgstab_data);
      return solver.solve(A,x,b,precond);
    }
  else if (solver_name=="gmres")
    {
      SolverGMRES<Matrix,Vector> solver(*control,*vector_memory,
					gmres_data);
      return solver.solve(A,x,b,precond);
    }
  else
    Assert(false,ExcSolverDoesNotExist(solver_name));

  return Solver<Matrix,Vector>::breakdown;
};



template <class Matrix, class Vector>
string SolverSelector<Matrix, Vector>::get_solver_names()
{
  return "richardson|cg|bicgstab|gmres";
};



template <class Matrix, class Vector>
void SolverSelector<Matrix, Vector>::set_data(
  const typename SolverGMRES<Matrix,Vector>::AdditionalData &data)
{ 
  gmres_data=data; 
}


template <class Matrix, class Vector>
void SolverSelector<Matrix, Vector>::set_data(
  const typename SolverRichardson<Matrix,Vector>::AdditionalData &data)
{ 
  richardson_data=data; 
}


template <class Matrix, class Vector>
void SolverSelector<Matrix, Vector>::set_data(
  const typename SolverCG<Matrix,Vector>::AdditionalData &data) 
{ 
  cg_data=data; 
}


template <class Matrix, class Vector>
void SolverSelector<Matrix, Vector>::set_data(
  const typename SolverBicgstab<Matrix,Vector>::AdditionalData &data) 
{ 
  bicgstab_data=data; 
};



/*----------------------------   solver_selector.h     ---------------------------*/
/* end of #ifndef __solver_selector_H */
#endif
/*----------------------------   solver_selector.h     ---------------------------*/
