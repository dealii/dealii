//----------------------------  base.h  ---------------------------
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
//----------------------------  base.h  ---------------------------
#ifndef __deal2__base_h
#define __deal2__base_h


/*----------------------------   problem_base.h     ---------------------------*/


#include <lac/sparse_matrix.h>
#include <base/exceptions.h>
#include <dofs/dof_constraints.h>
#include <fe/fe_update_flags.h>
#include <map>
#include <string>


/**
 * The use of this class is now deprecated!
 *
 * Base class for user problems. This class stores the system matrix and right
 * hand side vectors as well as a solution vector. It initiates the assemblage
 * process of matrix and vectors and so on.
 *
 * This class is not extremely versatile as could certainly be. For example
 * it presently only supports sparse matrices and has no multigrid features.
 * However, all these things depend strongly on the problem and it seems
 * best to implement many of these things yourself. Thus, this class is more
 * a display of concept how to work with deal.II.
 *
 *
 * \subsection{Assemblage}
 *
 * The #assemble# member function does the assemblage of the system matrix and
 * the given number of right hand sides. It does the following steps:
 * \begin{itemize}
 *   \item Initialize solution vector with zero entries.
 *   \item Create sparsity pattern of the system matrix and condense it with
 *     the constraints induced by hanging nodes.
 *   \item Initialize an assembler object.
 *   \item Loop over all cells and assemble matrix and vectors using the given
 *     quadrature formula and the equation object which contains the weak
 *     formulation of the equation.
 *   \item Apply Dirichlet boundary conditions. See the section on boundary
 *     conditions for more details.
 *   \item Condense the system matrix and right hand side with the constraints
 *     induced by hanging nodes.
 * \end{itemize}
 *
 * The #assemble# function needs an object describing the boundary of the domain,
 * since for higher order finite elements, we may be tempted to use curved faces
 * of cells for better approximation of the boundary. In this case, the
 * transformation from the unit cell to the real cell requires knowledge of
 * the exact boundary of the domain.
 * 
 *
 * \subsection{Solving}
 *
 * Calling the #solve# function with a solver object, the system of equations
 * which results after having called the #assemble# function is solved. After
 * this, the solution vector is distributed again, i.e. the constrained nodes
 * are given their correct values.
 *
 *
 * \subsection{Boundary conditions}
 *
 * During assemblage of matrices and right hand side, use is made of dirichlet
 * boundary conditions (in short: bc) specified to the #assemble# function. You
 * can specify a list of pairs of boundary indicators (of type #unsigned char#;
 * see the section in the documentation of the \Ref{Triangulation} class for more
 * details) and the according functions denoting the dirichlet boundary values
 * of the nodes on boundary faces with this boundary indicator.
 *
 * To actually apply the boundary conditions, use is made of the
 * #MatrixTools::apply_boundary_values# function and by interpolation of
 * the #boundary_values# using the #MatrixTool::interpolate_boundary_values#
 * function. See there for more information.
 *
 *
 * @author Wolfgang Bangerth, 1998
 */
template <int dim>
class ProblemBase
{
  public:
				     /**
				      *	Declare a data type which denotes a
				      *	mapping between a boundary indicator
				      *	and the function denoting the boundary
				      *	values on this part of the boundary.
				      *	Only one boundary function may be given
				      *	for each boundary indicator, which is
				      *	guaranteed by the #map# data type.
				      *	
				      *	See the general documentation of this
				      *	class for more detail.
				      */
    typedef map<unsigned char,const Function<dim>*> FunctionMap;
				     /**
				      * Typdedef an iterator which assembles
				      * matrices and vectors.
				      */
    typedef TriaActiveIterator<dim, Assembler<dim> > active_assemble_iterator;
    
				     /**
				      * Constructor.
				      */
    ProblemBase ();

				     /**
				      * Destructor. Declare this only to have
				      * a virtual destructor, which is safer
				      * as we have virtual functions.
				      * It actually does nothing spectacular.
				      */
    virtual ~ProblemBase ();

				     /**
				      * Use this triangulation and
				      * degree of freedom object during the
				      * lifetime of this object. The dof
				      * object must refer to the given
				      * triangulation.
				      */
    void set_tria_and_dof (Triangulation<dim> *tria,
			   DoFHandler<dim>    *dof_handler);

				     /**
				      * Reset all fields to a state as if we
				      * were right after calling the
				      * constructor. This is useful if you
				      * want to use an object derived from
				      * this base class for multiple
				      * successive calculations. In special,
				      * all aquired memory should be freed
				      * until it is needed again.
				      */
    void clear ();
    
				     /**
				      * Initiate the process of assemblage of
				      * vectors and system matrix. Use the
				      * given equation object and the given
				      * quadrature formula. Also use the list
				      * of dirichlet boundary value functions
				      * (by default, no dirichlet bc are assumed
				      * which means that all bc are included
				      * into the weak formulation).
				      *
				      * For what exactly happens here, refer to
				      * the general doc of this class.
				      */
    virtual void assemble (const Equation<dim>      &equation,
			   const Quadrature<dim>    &q,
			   const UpdateFlags         update_flags,
			   const FunctionMap        &dirichlet_bc = FunctionMap());
    
				     /**
				      * Solve the system of equations. This uses
				      * a simple CG method.
				      */
    virtual void solve ();

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
				      * The solution name are
				      * derived by calling the virtual
				      * function #get_solution_name#.
				      */
    virtual void fill_data (DataOut<dim> &) const;


				     /**
				      * Return the name of the solution as a
				      * #string#. The default implementation
				      * returns #"solution"#.
				      * Overload this function, if you
				      * want anything else.
				      */
    virtual string get_solution_name () const;

				     /**
				      * Exception
				      */
    DeclException0 (ExcDofAndTriaDontMatch);
				     /**
				      * Exception
				      */
    DeclException0 (ExcNoTriaSelected);
    
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
    SparsityPattern      system_sparsity;

				     /**
				      * System matrix.
				      */
    SparseMatrix<double> system_matrix;

				     /**
				      * Vector storing the right hand side.
				      */
    Vector<double>       right_hand_side;

				     /**
				      * Solution vector.
				      */
    Vector<double>       solution;

				     /**
				      * List of constraints introduced by
				      * hanging nodes.
				      */
    ConstraintMatrix    constraints;

    friend class Assembler<dim>;
};


/*----------------------------   problem_base.h     ---------------------------*/

#endif
/*----------------------------   problem_base.h     ---------------------------*/
