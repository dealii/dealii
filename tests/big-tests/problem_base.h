/*----------------------------   problem_base.h     ---------------------------*/
/* Previously deal.II/include/numerics/base.h and deal.II/source/numerics/base.cc,
   but deprecated before deal.II version 3.0 and removed after version 3.1

   Kept here to ensure that the example programs in big-tests still
   compile. It is put into this directory to allow simple access from
   all the subdirectories.
*/

#ifndef __deal2__base_h
#define __deal2__base_h




#include <lac/sparse_matrix.h>
#include <base/exceptions.h>
#include <dofs/dof_constraints.h>
#include <fe/fe_update_flags.h>
#include <map>
#include <string>

template <int dim> class DataOut;
template <int dim> class Assembler;
template <int dim> class Equation;

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
 * @sect3{Assemblage}
 *
 * The @p{assemble} member function does the assemblage of the system matrix and
 * the given number of right hand sides. It does the following steps:
 * @begin{itemize}
 *   @item Initialize solution vector with zero entries.
 *   @item Create sparsity pattern of the system matrix and condense it with
 *     the constraints induced by hanging nodes.
 *   @item Initialize an assembler object.
 *   @item Loop over all cells and assemble matrix and vectors using the given
 *     quadrature formula and the equation object which contains the weak
 *     formulation of the equation.
 *   @item Apply Dirichlet boundary conditions. See the section on boundary
 *     conditions for more details.
 *   @item Condense the system matrix and right hand side with the constraints
 *     induced by hanging nodes.
 * @end{itemize}
 *
 * The @p{assemble} function needs an object describing the boundary of the domain,
 * since for higher order finite elements, we may be tempted to use curved faces
 * of cells for better approximation of the boundary. In this case, the
 * transformation from the unit cell to the real cell requires knowledge of
 * the exact boundary of the domain.
 * 
 *
 * @sect3{Solving}
 *
 * Calling the @p{solve} function with a solver object, the system of equations
 * which results after having called the @p{assemble} function is solved. After
 * this, the solution vector is distributed again, i.e. the constrained nodes
 * are given their correct values.
 *
 *
 * @sect3{Boundary conditions}
 *
 * During assemblage of matrices and right hand side, use is made of dirichlet
 * boundary conditions (in short: bc) specified to the @p{assemble} function. You
 * can specify a list of pairs of boundary indicators (of type @p{unsigned char};
 * see the section in the documentation of the @ref{Triangulation} class for more
 * details) and the according functions denoting the dirichlet boundary values
 * of the nodes on boundary faces with this boundary indicator.
 *
 * To actually apply the boundary conditions, use is made of the
 * @p{MatrixTools::apply_boundary_values} function and by interpolation of
 * the @p{boundary_values} using the @p{MatrixTool::interpolate_boundary_values}
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
				      *	guaranteed by the @p{map} data type.
				      *	
				      *	See the general documentation of this
				      *	class for more detail.
				      */
    typedef typename std::map<unsigned char,const Function<dim>*> FunctionMap;
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
				      * Initialize the @p{DataOut} object with
				      * the grid and DoF handler used for this
				      * computation, as well as with the
				      * solution. Overload this function if
				      * you have multiple data sets to be
				      * written, or alternativelt call this
				      * function and attach the additional
				      * vectors directly to the @p{DataOut}
				      * object.
				      *
				      * The solution name are
				      * derived by calling the virtual
				      * function @p{get_solution_name}.
				      */
    virtual void fill_data (DataOut<dim> &) const;


				     /**
				      * Return the name of the solution as a
				      * @p{string}. The default implementation
				      * returns @p{"solution"}.
				      * Overload this function, if you
				      * want anything else.
				      */
    virtual std::string get_solution_name () const;

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




#include <numerics/assembler.h>
#include <numerics/matrices.h>
#include <numerics/vectors.h>
#include <dofs/dof_constraints.h>
#include <grid/tria_iterator.h>
#include <numerics/data_out.h>
#include <dofs/dof_tools.h>
#include <base/function.h>
#include <fe/fe.h>
#include <base/quadrature.h>
#include <lac/vector.h>
#include <lac/precondition.h>
#include <lac/solver_cg.h>
#include <lac/vector_memory.h>

#include <map>
#include <numeric>
#include <algorithm>
#include <cmath>


template <int dim>
ProblemBase<dim>::ProblemBase () :
		tria(0),
		dof_handler(0),
		system_sparsity(),        // dummy initialisation, is later reinit'd
		system_matrix()           // dummy initialisation, is later reinit'd
{};


template <int dim>
void ProblemBase<dim>::set_tria_and_dof (Triangulation<dim> *t,
					 DoFHandler<dim>    *d)
{
  tria        = t;
  dof_handler = d;

				   // allow a user to reset both tria and
				   // dof to null pointers, but don't allow
				   // something other
  if ((tria != 0) || (dof_handler != 0))
    {
      Assert ((tria!=0) && (dof_handler!=0), ExcNoTriaSelected());
      Assert (tria == &dof_handler->get_tria(), ExcDofAndTriaDontMatch());
    };
};


template <int dim>
void ProblemBase<dim>::clear () {
  if (tria)        { delete tria;         tria        = 0; };
  if (dof_handler) { delete dof_handler;  dof_handler = 0; };
  system_sparsity.reinit (0,0,1);
  system_matrix.clear ();
  right_hand_side.reinit (0);
  solution.reinit (0);
  constraints.clear ();
};


template <int dim>
ProblemBase<dim>::~ProblemBase () {};


template <int dim>
void ProblemBase<dim>::assemble (const Equation<dim>      &equation,
				 const Quadrature<dim>    &quadrature,
				 const UpdateFlags         update_flags,
				 const FunctionMap        &dirichlet_bc)
{
  Assert ((tria!=0) && (dof_handler!=0), ExcNoTriaSelected());
  
  system_sparsity.reinit (dof_handler->n_dofs(),
			  dof_handler->n_dofs(),
			  dof_handler->max_couplings_between_dofs());
  right_hand_side.reinit (dof_handler->n_dofs());
  
				   // make up sparsity pattern and
				   // compress with constraints
  constraints.clear ();
  DoFTools::make_hanging_node_constraints (*dof_handler, constraints);
  constraints.close ();
  DoFTools::make_sparsity_pattern (*dof_handler, system_sparsity);
  constraints.condense (system_sparsity);

				   // reinite system matrix
  system_matrix.reinit (system_sparsity);
				   // reinit solution vector, preset
				   // with zeroes.
  solution.reinit (dof_handler->n_dofs());
  
				   // create assembler
  Assembler<dim>::AssemblerData data (*dof_handler,
				      true, true, //assemble matrix and rhs
				      system_matrix,
				      right_hand_side,
				      quadrature,
				      update_flags);
  active_assemble_iterator assembler (tria,
				      tria->begin_active()->level(),
				      tria->begin_active()->index(),
				      &data);
				   // loop over all cells, fill matrix and rhs
  do 
    {
      assembler->assemble (equation);
    }
  while ((++assembler).state() == valid);

				   // condense system matrix in-place
  constraints.condense (system_matrix);

				   // condense right hand side in-place
  constraints.condense (right_hand_side);

				   // apply Dirichlet bc as described
				   // in the docs
  std::map<unsigned int, double> boundary_value_list;

  for (typename FunctionMap::const_iterator dirichlet = dirichlet_bc.begin() ;
       dirichlet != dirichlet_bc.end() ; ++dirichlet)
    VectorTools::interpolate_boundary_values (*dof_handler,
					      dirichlet->first,
					      *(dirichlet->second), 
					      boundary_value_list);
  MatrixTools<dim>::apply_boundary_values (boundary_value_list,
					   system_matrix, solution,
					   right_hand_side,
					   true);
};


template <int dim>
void ProblemBase<dim>::solve ()
{
  Assert ((tria!=0) && (dof_handler!=0), ExcNoTriaSelected());
  
  SolverControl           control(4000, 1e-16);
  PrimitiveVectorMemory<> memory;
  SolverCG<>              cg(control,memory);

  PreconditionSSOR<> prec;
  prec.initialize (system_matrix, 1.2);

				   // solve
  cg.solve (system_matrix, solution, right_hand_side, prec);
				   // distribute solution
  constraints.distribute (solution);
};


template <int dim>
void ProblemBase<dim>::fill_data (DataOut<dim> &out) const
{
  Assert ((tria!=0) && (dof_handler!=0), ExcNoTriaSelected());
  
  out.clear ();
  out.attach_dof_handler (*dof_handler);

  out.add_data_vector (solution, get_solution_name());
};


template <int dim>
std::string ProblemBase<dim>::get_solution_name () const
{
  return "solution";
};




#endif

