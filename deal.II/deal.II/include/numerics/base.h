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
template <int dim> class Function;



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


  {\bf Solving}

  Calling the #solve# function with a solver object, the system of equations
  which results after having called the #assemble# function is solved. After
  this, the solution vector is distributed again, i.e. the constrained nodes
  are given their correct values.
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
				      * Destructor. Declare this only to have
				      * a virtual destructor, which is safer
				      * as we have virtual functions.
				      * It actually does nothing spectacular.
				      */
    virtual ~ProblemBase ();

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


				     /**
				      * Solve the system of equations.
				      */
    virtual void solve ();

				     /**
				      * Integrate the difference between
				      * the solution computed before and
				      * the reference solution, which
				      * is given as a continuous function
				      * object. The integration is
				      * performed using the given quadrature
				      * rule and assumes that the given
				      * finite element objects equals
				      * that used for the computation of
				      * the solution.
				      *
				      * The result ist stored in the
				      * #difference# object, where each
				      * entry equals the L_1 norm of
				      * the difference on one cell. The
				      * order of entries is the same as
				      * a #cell_iterator# takes when
				      * started with #begin_active# and
				      * promoted with the #++# operator.
				      *
				      * You can use the
				      * #distribute_cell_to_dof_vector#
				      * function to convert cell based
				      * data to a data vector with values
				      * on the degrees of freedom, which
				      * can then be attached to a
				      * #DataOut# object.
				      *
				      * To get the global L_1 error,
				      * you have to sum up the entries
				      * in #difference#, e.g. using the STL
				      * function
				      * #accumulate (d.begin(), d.end(), 0)#,
				      * if #d# is the difference vector.
				      */
    void integrate_L1_difference (const Function<dim>      &exact_solution,
				  vector<double>           &difference,
				  const Quadrature<dim>    &q,
				  const FiniteElement<dim> &fe) const;
    
				     /**
				      * Initialize the #DataOut# object with
				      * the grid and DoF handler used for this
				      * computation, as well as with the
				      * solution. Overload this function if
				      * you have multiple data sets to be
				      * written, or alternativelt call this
				      * function and attach the additional
				      * vectors directly to the #DataOut#
				      * object.
				      *
				      * Solution name and physical units are
				      * derived by calling the virtual
				      * function #get_solution_name#.
				      */
    virtual void fill_data (DataOut<dim> &) const;


				     /**
				      * Return solution name and
				      * physical units as a pair of
				      * #char*#. The default implementation
				      * returns #make_pair ("solution","")#,
				      * which results in "<dimensionless>"
				      * upon output.
				      * Overload this function, if you
				      * want anything else.
				      */
    virtual pair<char*,char*> get_solution_name () const;
    
				     /**
				      * Exception
				      */
    DeclException0 (ExcDofAndTriaDontMatch);
				     /**
				      * Exception
				      */
    DeclException0 (ExcNoMemory);
				     /**
				      * Exception
				      */
    DeclException0 (ExcInvalidFE);
    
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
