/*----------------------------   problem_base.h     ---------------------------*/
/*      $Id$                 */
#ifndef __problem_base_H
#define __problem_base_H
/*----------------------------   problem_base.h     ---------------------------*/


#include <lac/dsmatrix.h>
#include <base/exceptions.h>
#include <grid/dof_constraints.h>


// forward declaration
template <int dim> class Triangulation;
template <int dim> class DoFHandler;
template <int dim> class FiniteElement;
template <int dim> class Quadrature;
template <int dim> class DataOut;




/**
  Base class for user problems. This class stores the system matrix and right
  hand side vectors as well as a solution vector. It initiates the assemblage
  process of matrix and vectors and so on.

  This class is not extremely versatile as could certainly be. For example
  it presently only supports sparse matrices and has no multigrid features.
  However, all these things depend strongly on the problem and it seems
  best to implement many of these things yourself. Thus, this class is more
  a display of concept haw to work with deal.II.


  {\bf Assemblage}

  The #assemble# member function does the assemblage of the system matrix and
  the given number of right hand sides. It does the following steps:
  \begin{itemize}
    \item Create sparsity pattern of the system matrix and condense it with
      the constraints induced by hanging nodes.
    \item Initialize an assembler object.
    \item Loop over all cells and assemble matrix and vectors using the given
      quadrature formula and the equation object which contains the weak
      formulation of the equation.
    \item Condense the system matrix and right hand side with the constraints
      induced by hanging nodes.
  \end{itemize}
  */
template <int dim>
class ProblemBase {
  public:
				     /**
				      * Constructor. Use this triangulation and
				      * degree of freedom object during the
				      * lifetime of this object. The dof
				      * object must refer to the given
				      * triangulation.
				      */
    ProblemBase (Triangulation<dim> *tria,
		 DoFHandler<dim>    *dof_handler);

				     /**
				      * Initiate the process of assemblage of
				      * vectors and system matrix. Use the
				      * given equation object and the given
				      * quadrature formula.
				      *
				      * For what exactly happens here, refer to
				      * the general doc of this class.
				      */
    virtual void assemble (const Equation<dim>      &equation,
			   const Quadrature<dim>    &q,
			   const FiniteElement<dim> &fe);


    void solve ();
    void fill_data (DataOut<dim> &) const;
    
				     /**
				      * Exception
				      */
    DeclException0 (ExcDofAndTriaDontMatch);
				     /**
				      * Exception
				      */
    DeclException0 (ExcNoMemory);
    
  protected:
				     /**
				      * Pointer to the triangulation to work on.
				      */
    Triangulation<dim> *tria;

				     /**
				      * Pointer to the degree of freedom handler
				      * to work with.
				      */
    DoFHandler<dim>    *dof_handler;

				     /**
				      * Sparsity pattern of the system matrix.
				      */
    dSMatrixStruct      system_sparsity;

				     /**
				      * System matrix.
				      */
    dSMatrix            system_matrix;

				     /**
				      * Vector storing the right hand side.
				      */
    dVector             right_hand_side;

				     /**
				      * Solution vector.
				      */
    dVector             solution;

				     /**
				      * List of constraints introduced by
				      * hanging nodes.
				      */
    ConstraintMatrix    constraints;
    
  friend class Assembler<dim>;
};

    


/*----------------------------   problem_base.h     ---------------------------*/
/* end of #ifndef __problem_base_H */
#endif
/*----------------------------   problem_base.h     ---------------------------*/
