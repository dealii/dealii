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

template <class Key, class T, class Compare = less<Key>, class Alloc = alloc> class map;




/**
   Denote which norm/integral is to be computed. The following possibilities
   are implemented:
   \begin{itemize}
   \item #mean#: the function or difference of functions is integrated
     on each cell.
   \item #L1_norm#: the absolute value of the function is integrated.
   \item #L2_norm#: the square of the function is integrated on each
     cell; afterwards the root is taken of this value.
   \end{itemize}
*/
enum NormType {
      mean,
      L1_norm,
      L2_norm,
      Linfty_norm
};




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
    \item Initialize solution vector with zero entries.
    \item Create sparsity pattern of the system matrix and condense it with
      the constraints induced by hanging nodes.
    \item Initialize an assembler object.
    \item Loop over all cells and assemble matrix and vectors using the given
      quadrature formula and the equation object which contains the weak
      formulation of the equation.
    \item Apply Dirichlet boundary conditions. See the section on boundary
      conditions for more details.
    \item Condense the system matrix and right hand side with the constraints
      induced by hanging nodes.
  \end{itemize}


  {\bf Solving}

  Calling the #solve# function with a solver object, the system of equations
  which results after having called the #assemble# function is solved. After
  this, the solution vector is distributed again, i.e. the constrained nodes
  are given their correct values.


  {\bf Boundary conditions}

  During assemblage of matrices and right hand side, use is made of dirichlet
  boundary conditions (in short: bc) specified to the #assemble# function. You
  can specify a list of pairs of boundary indicators (of type #unsigned char#;
  see the section in the documentation of the \Ref{Triangulation} class for more
  details) and the according functions denoting the dirichlet boundary values
  of the nodes on boundary faces with this boundary indicator.

  Usually, all other boundary conditions, such as inhomogeneous Neumann values
  or mixed boundary conditions are handled in the weak formulation. No attempt
  is made to include this into the process of assemblage therefore.

  The inclusion into the assemblage process is as follows: when the matrix and
  vectors are set up, a list of nodes subject to dirichlet bc is made and
  matrix and vectors are changed accordingly. This is done by deleting all
  entries in the matrix in the line of this degree of freedom, setting the
  main diagonal entry to one and the right hand side element to the
  boundary value at this node. This forces this node's value to be as specified.
  To decouple the remaining linear system of equations and to make the system
  symmetric again (at least if it was before), one Gauss elimination
  step is performed with this line, by adding this (now almost empty) line to
  all other lines which couple with the given degree of freedom and thus
  eliminating all coupling between this degree of freedom and others. Now
  also the column consists only of zeroes, apart from the main diagonal entry.

  To make solving faster, we preset the solution vector with the right boundary
  values. Since boundary nodes can never be hanging nodes, and since all other
  entries of the solution vector are zero, we need not condense the solution
  vector if the condensation process is done in-place. If done by copying
  matrix and vectors to smaller ones, it would also be necessary to condense
  the solution vector to preserve the preset boundary values.
  
  It it not clear whether the deletion of coupling between the boundary degree
  of freedom and other dofs really forces the corresponding entry in the
  solution vector to have the right value when using iterative solvers,
  since their search directions may contains components in the direction
  of the boundary node. For this reason, we perform a very simple line
  balancing by not setting the main diagonal entry to unity, but rather
  to the value it had before deleting this line. Of course we have to change
  the right hand side appropriately. This is not a very good
  strategy, but it at least should give the main diagonal entry a value
  in the right order of dimension, which makes the solving process a bit
  more stable. A refined algorithm would set the entry to the mean of the
  other diagonal entries, but this seems to be too expensive.

  Because of the mentioned question, whether or not a preset solution value
  which does not couple with other degrees of freedom remains its value or
  not during solving iteratively, it may or may not be necessary to set
  the correct value after solving again. This question is an open one as of
  now and may be answered by future experience.
  

  {\bf Computing errors}

  The function #integrate_difference# performs the calculation of the error
  between the finite element solution and a given (continuous) reference
  function in different norms. The integration is performed using a given
  quadrature formulae and assumes that the given finite element objects equals
  that used for the computation of the solution.
  
  The result ist stored in a vector (named #difference#), where each entry
  equals the given norm of the difference on one cell. The order of entries
  is the same as a #cell_iterator# takes when started with #begin_active# and
  promoted with the #++# operator.
  
  You can use the #distribute_cell_to_dof_vector# function to convert cell
  based data to a data vector with values on the degrees of freedom, which
  can then be attached to a #DataOut# object to be printed.
  
  Presently, there is the possibility to compute the following values from the
  difference, on each cell: #mean#, #L1_norm#, #L2_norm#, #Linfty_norm#.
  For the mean difference value, the reference function minus the numerical
  solution is computed, not the other way round.

  The infinity norm of the difference on a given cell returns the maximum
  absolute value of the difference at the quadrature points given by the
  quadrature formula parameter. This will in some cases not be too good
  an approximation, since for example the Gauss quadrature formulae do
  not evaluate the difference at the end or corner points of the cells.
  You may want to chose a quadrature formula with more quadrature points
  or one with another distribution of the quadrature points in this case.
  
  To get the {\it global} L_1 error, you have to sum up the entries in
  #difference#, e.g. using the STL function
  #accumulate (d.begin(), d.end(), 0)#, if #d# is the difference vector.
  For the global L_2 difference, you have to sum up the squares of the
  entries and take the root of the sum. These two operations represent the
  l_1 and l_2 norms of the vectors, but you need not take the absolute
  value of each entry, since the cellwise norms are already positive.
  
  To get the global mean difference, simply sum up the elements as above.
  To get the L_\infty norm, take the maximum of the vector elements, e.g.
  using the STL function #max_element (d.begin(), d.end())#.
  */
template <int dim>
class ProblemBase {
  public:
				     /**
				      * Declare a data type which denotes a
				      * mapping between a boundary indicator
				      * and the function denoting the boundary
				      * values on this part of the boundary.
				      * Only one boundary function may be given
				      * for each boundary indicator, which is
				      * guaranteed by the #map# data type.
				      *
				      * See the general documentation of this
				      * class for more detail.
				      */
    typedef map<unsigned char,Function<dim> > DirichletBC;
    
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
			   const FiniteElement<dim> &fe,
			   const DirichletBC        &dirichlet_bc = DirichletBC());
    
				     /**
				      * Solve the system of equations.
				      */
    virtual void solve ();

				     /**
				      * Integrate the difference between
				      * the solution computed before and
				      * the reference solution, which
				      * is given as a continuous function
				      * object.
				      *
				      * See the general documentation of this
				      * class for more information.
				      */
    void integrate_difference (const Function<dim>      &exact_solution,
			       vector<double>           &difference,
			       const Quadrature<dim>    &q,
			       const FiniteElement<dim> &fe,
			       const NormType           &norm) const;
    
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
				     /**
				      * Exception
				      */
    DeclException0 (ExcNotImplemented);
    
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

				     /**
				      * Apply dirichlet boundary conditions
				      * to the system matrix and vectors
				      * as described in the general
				      * documentation.
				      */
    void apply_dirichlet_bc (dSMatrix &matrix,
			     dVector  &solution,
			     dVector  &right_hand_side,
			     const DirichletBC &dirichlet_bc);
    
    friend class Assembler<dim>;
};

    


/*----------------------------   problem_base.h     ---------------------------*/
/* end of #ifndef __problem_base_H */
#endif
/*----------------------------   problem_base.h     ---------------------------*/
