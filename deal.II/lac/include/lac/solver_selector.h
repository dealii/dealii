//----------------------------  solver_selector.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  solver_selector.h  ---------------------------
#ifndef __deal2__solver_selector_h
#define __deal2__solver_selector_h


#include <base/config.h>
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


/**
 * By calling the @p{solve} function of this @p{SolverSelector}, it selects
 * the @p{solve} function of that @p{Solver} that was specified in the constructor
 * of this class.
 *
 * @sect3{Usage}
 * The simplest use of this class is the following:
 * @begin{verbatim}
 *                                  // generate a @p{SolverControl} and
 *                                  // a @p{VectorMemory}
 * SolverControl control;
 * VectorMemory<Vector<double> > memory;
 *                                  // Line 3:
 *                                  //
 *                                  // generate a @p{SolverSelector} that
 *                                  // calls the @p{SolverCG}
 * SolverSelector<SparseMatrix<double>, Vector<double> > 
 *   solver_selector("cg", control, memory);
 *                                  // generate e.g. a @p{PreconditionRelaxation}
 * PreconditionRelaxation<SparseMatrix<double>, Vector<double> >
 *   preconditioning(A, &SparseMatrix<double>
 *                   ::template precondition_SSOR<double>,0.8);
 *                                  // call the @p{solve} function with this
 *                                  // preconditioning as last argument
 * solver_selector.solve(A,x,b,preconditioning);
 * @end{verbatim}
 * But the full usefulness of the @p{SolverSelector} class is not
 * clear until the presentation of the following example that assumes
 * the user using the @p{ParameterHandler} class and having declared a
 * "solver" entry, e.g. with
 * @begin{verbatim}
 * Parameter_Handler prm;
 * prm.declare_entry ("solver", "none",
 *                    Patterns::Sequence(SolverSelector::get_solver_names()));
 * ...
 * @end{verbatim}
 * Assuming that in the users parameter file there exists the line
 * @begin{verbatim}
 * set solver = cg
 * @end{verbatim}
 * then `Line 3' of the above example reads
 * @begin{verbatim}
 * SolverSelector<SparseMatrix<double>, Vector<double> > 
 *   solver_selector(prm.get("solver"), control, memory);
 * @end{verbatim}
 *
 *
 * If at some time there exists a new solver "xyz" then the user does not need
 * to change his program. Only in the implementation of the @p{SolverSelector}
 * the calling of this solver has to be added and each user with program lines
 * quoted above only needs to 'set solver = xyz' in his parameter file to get
 * access to that new solver.  :-)
 *
 * (By the way, thanks to Wolfgang for implementing the @p{ParameterHandler}.)
 * 
 * @author Ralf Hartmann, 1999
 */
template <class Vector = Vector<double> >
class SolverSelector : public Subscriptor
{
  public:

				     /**
				      * Constructor. Takes the @p{SolverName},
				      * the @p{SolverControl} and the 
				      * @p{VectorMemory} as argument.
				      */
    SolverSelector (const std::string    &solvername,
		    SolverControl        &control,
		    VectorMemory<Vector> &vectorm);
  
				     /**
				      * Destructor
				      */
    ~SolverSelector();

				     /**
				      * Solver procedure. Calls the @p{solve}
				      * function
				      * of the @p{solver} whose @p{SolverName}
				      * was specified in the constructor.
				      * 
				      */
    template<class Matrix, class Preconditioner>
    void solve(const Matrix &A,
	       Vector &x,
	       const Vector &b,
	       const Preconditioner &precond) const;
    
				     /**
				      * Set the additional data. For more
				      * info see the @p{Solver} class.
				      */
    void set_data(const typename SolverRichardson<Vector>
		  ::AdditionalData &data);

				     /**
				      * Set the additional data. For more
				      * info see the @p{Solver} class.
				      */
    void set_data(const typename SolverCG<Vector>
		  ::AdditionalData &data);

				     /**
				      * Set the additional data. For more
				      * info see the @p{Solver} class.
				      */
    void set_data(const typename SolverBicgstab<Vector>
		  ::AdditionalData &data); 

				     /**
				      * Set the additional data. For more
				      * info see the @p{Solver} class.
				      */
    void set_data(const typename SolverGMRES<Vector>
		  ::AdditionalData &data);

				     /**
				      * Get the names of all implemented
				      * solvers.
				      */
    static std::string get_solver_names ();

				     /**
				      * Exception.
				      */
    DeclException1 (ExcSolverDoesNotExist,
		    std::string, << "Solver " << arg1 << " does not exist. Use one of "
		    << std::endl << get_solver_names());

  protected:
				     /**
				      * Stores the Name of the solver.
				      */
    std::string solver_name;
    
				     /**
				      * Stores the @p{SolverControl} that
				      * is needed in the constructor of
				      * each @p{Solver} class.
				      */
    SmartPointer<SolverControl> control;

				     /**
				      * Stores the @p{VectorMemory} that
				      * is needed in the constructor of
				      * each @p{Solver} class.
				      */
    SmartPointer<VectorMemory<Vector> > vector_memory;

  private:
				     /**
				      * Stores the additional data.
				      */
    typename SolverRichardson<Vector>::AdditionalData richardson_data;

				     /**
				      * Stores the additional data.
				      */
    typename SolverCG<Vector>::AdditionalData cg_data;

				     /**
				      * Stores the additional data.
				      */
    typename SolverBicgstab<Vector>::AdditionalData bicgstab_data;

				     /**
				      * Stores the additional data.
				      */
    typename SolverGMRES<Vector>::AdditionalData gmres_data;
};


/* --------------------- Inline and template functions ------------------- */


template <class Vector>
SolverSelector<Vector>::SolverSelector(const std::string    &solver_name,
				       SolverControl        &control,
				       VectorMemory<Vector> &vector_memory) :
		solver_name(solver_name),
		control(&control),
		vector_memory(&vector_memory)
{};


template <class Vector>
SolverSelector<Vector>::~SolverSelector()
{};


template <class Vector>
template<class Matrix, class Preconditioner>
void
SolverSelector<Vector>::solve(const Matrix &A,
			      Vector &x,
			      const Vector &b,
			      const Preconditioner &precond) const
{
  if (solver_name=="richardson")
    {
      SolverRichardson<Vector> solver(*control,*vector_memory,
					     richardson_data);
      solver.solve(A,x,b,precond);
    }       
  else if (solver_name=="cg")
    {
      SolverCG<Vector> solver(*control,*vector_memory,
				     cg_data);
      solver.solve(A,x,b,precond);
    }
  else if (solver_name=="bicgstab")
    {
      SolverBicgstab<Vector> solver(*control,*vector_memory,
					   bicgstab_data);
      solver.solve(A,x,b,precond);
    }
  else if (solver_name=="gmres")
    {
      SolverGMRES<Vector> solver(*control,*vector_memory,
					gmres_data);
      solver.solve(A,x,b,precond);
    }
  else
    Assert(false,ExcSolverDoesNotExist(solver_name));
};


template <class Vector>
std::string SolverSelector<Vector>::get_solver_names()
{
  return "richardson|cg|bicgstab|gmres";
};


template <class Vector>
void SolverSelector<Vector>::set_data(
  const typename SolverGMRES<Vector>::AdditionalData &data)
{ 
  gmres_data=data; 
}


template <class Vector>
void SolverSelector<Vector>::set_data(
  const typename SolverRichardson<Vector>::AdditionalData &data)
{ 
  richardson_data=data; 
}


template <class Vector>
void SolverSelector<Vector>::set_data(
  const typename SolverCG<Vector>::AdditionalData &data) 
{ 
  cg_data=data; 
}


template <class Vector>
void SolverSelector<Vector>::set_data(
  const typename SolverBicgstab<Vector>::AdditionalData &data) 
{ 
  bicgstab_data=data; 
};


#endif
