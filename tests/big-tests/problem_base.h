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



/**
 * Equation class to be passed to the @ref{Assembler} if you want to make up the
 * mass matrix for your problem. The mass matrix is the matrix with
 * $m_{ij} = \int_\Omega \phi_i(x) \phi_j(x) dx$.
 *
 * You may pass a coefficient function to the constructor. If you do so, the
 * assemble routines compute the matrix
 * $m_{ij} = \int_\Omega a(x) \phi_i(x) \phi_j(x) dx$
 * instead. The coefficient will in many cases be a strictly positive function.
 *
 * The class also has functions to create a right hand side
 * $f_i = \int_\Omega f(x) \phi_i(x) dx$. The function $f(x)$ has to be
 * given to the constructor; if none is given, an error is issued if you
 * try to create a right hand side vector. The function to create right
 * hand side vectors is the same for all the matrix class in this file,
 * since it does not depend on the operator.
 *
 * The defaults for both right hand side and coefficient function is a
 * @p{NULL} pointer. If you need a coefficient but no right hand side object,
 * simply pass a @p{NULL} pointer to the constructor for its first argument.
 *
 *
 * @sect3{Other possibilities}
 *
 * You will usually want to use this object only if you have coefficients
 * which vary over each cell. If you have coefficients which are constant
 * on each cell or even on the whole domain, you can get the local mass
 * matrix easier by calling the @ref{FiniteElement}@p{::get_local_mass_matrix} and
 * then scaling this one on each cell. This has the additional benefit that
 * the mass matrix is evaluated exactly, i.e. not using a quadrature formula
 * and is normally much faster since it can be precomputed and needs only
 * be scaled appropriately.
 *
 * The useful use of this object is therefore probable one of the following
 * cases:
 * @begin{itemize}
 * @item Mass lumping: use an @ref{Assembler} object and a special quadrature
 *   formula to voluntarily evaluate the mass matrix incorrect. For example
 *   by using the trapezoidal formula, the mass matrix will become a
 *   diagonal (at least if no hanging nodes are considered). However, there
 *   may be easier ways to set up the resulting matrix, for example by
 *   scaling the diagonal elements of the unit matrix by the area element
 *   of the respective cell.
 *
 * @item Nonconstant coefficient: if the coefficient varies considerably over
 *   each element, there is no way around this class. However, there are many
 *   cases where it is sufficient to assume that the function be constant on
 *   each cell (taking on its mean value throughout the cell for example, or
 *   more easily computed, its value at the center of mass of the element).
 *   A proper analysis of the error introduced by an assumed constant
 *   coefficient may be worth the effort.
 *
 *   Nonconstant coefficients to the mass matrix occur in mechanical problems
 *   if the density or other mechanical properties vary with the space
 *   coordinate.
 *
 * @item Simple plugging together of system matrices: if the system matrix has
 *    the form $s_{ij} = m_{ij} + \alpha a_{ij}$, for example, with $M$ and
 *    $A$ being the mass and laplace matrix, respectively (this matrix $S$
 *    occurs in the discretization of the heat and the wave equation, amoung
 *    others), once could conceive an equation object in which the @p{assemble}
 *    functions do nothing but sum up the contributions delivered by the
 *    @p{assemble} functions of the @ref{MassMatrix} and @ref{LaplaceMatrix} classes.
 *    Since numerical quadrature is necessary here anyway, this way is
 *    justifyable to quickly try something out. In the further process it
 *    may be useful to replace this behaviour by more sophisticated methods,
 *    however.
 * @end{itemize}
 */
template <int dim>
class MassMatrix :  public Equation<dim> {
  public:
				     /**
				      * Constructor. Pass a function object if
				      * you want to create a right hand side
				      * vector, pass a function pointer (default
				      * is a NULL pointer). It is your duty to
				      * guarantee that the function object for
				      * the right hand side lives at least as
				      * long as this object does.
				      *
				      * You may also pass a function describing
				      * the weight to the integral (see the
				      * general docs for more information). The
				      * same applies for this object as said
				      * above.
				      */
    MassMatrix (const Function<dim> * const rhs = 0,
		const Function<dim> * const a = 0);

				     /**
				      * Assemble the cell matrix and right hand
				      * side vector for this cell. You need to
				      * give a right hand side object to the
				      * constructor to use this function. If
				      * a coefficient was given to the
				      * constructor, it is used.
				      *
				      * This function assumes the cell matrix
				      * and right hand side to have the right
				      * size and to be empty.
				      */
    virtual void assemble (FullMatrix<double>  &cell_matrix,
			   Vector<double>      &rhs,
			   const FEValues<dim> &fe_values,
			   const typename DoFHandler<dim>::cell_iterator &) const;

				     /**
				      * Construct the cell matrix for this cell.
				      * If a coefficient was given to the
				      * constructor, it is used.
				      */
    virtual void assemble (FullMatrix<double>  &cell_matrix,
			   const FEValues<dim> &fe_values,
			   const typename DoFHandler<dim>::cell_iterator &) const;

				     /**
				      * Only construct the right hand side
				      * vector for this cell. You need to give
				      * a right hand side function to the
				      * constructor in order to call this
				      * function.
				      */
    virtual void assemble (Vector<double>      &rhs,
			   const FEValues<dim> &fe_values,
			   const typename DoFHandler<dim>::cell_iterator &) const;
    
				     /**
				      * Exception
				      */
    DeclException0 (ExcNoRHSSelected);
    
  protected:
				     /**
				      * Pointer to a function describing the
				      * right hand side of the problem. Should
				      * be zero if not given to the constructor
				      * and should then not be used.
				      */
    const Function<dim> * const right_hand_side;

				     /**
				      * Pointer to a function describing the
				      * coefficient to the integral for the
				      * matrix entries. Should
				      * be zero if not given to the constructor
				      * and should then not be used.
				      */
    const Function<dim> * const coefficient;
};


/**
 * Equation class to be passed to the @ref{Assembler} if you want to make up the
 * laplace matrix for your problem. The laplace matrix is the matrix with
 * $a_{ij} = \int_\Omega \nabla\phi_i(x) \cdot \nabla\phi_j(x) dx$.
 *
 * You may pass a coefficient function to the constructor. If you do so, the
 * assemble routines compute the matrix
 * $m_{ij} = \int_\Omega a(x) \nabla\phi_i(x) \cdot \nabla\phi_j(x) dx$
 * instead. The coefficient will in many cases be a strictly positive function.
 *
 * The class also has functions to create a right hand side
 * $f_i = \int_\Omega f(x) \phi_i(x) dx$. The function $f(x)$ has to be
 * given to the constructor; if none is given, an error is issued if you
 * try to create a right hand side vector. The function to create right
 * hand side vectors is the same for all the matrix class in this file,
 * since it does not depend on the operator.
 *
 * The defaults for both right hand side and coefficient function is a
 * @p{NULL} pointer. If you need a coefficient but no right hand side object,
 * simply pass a @p{NULL} pointer to the constructor for its first argument.
 */
template <int dim>
class LaplaceMatrix :  public Equation<dim> {
  public:
				     /**
				      * Constructor. Pass a function object if
				      * you want to create a right hand side
				      * vector, pass a function pointer (default
				      * is a NULL pointer). It is your duty to
				      * guarantee that the function object for
				      * the right hand side lives at least as
				      * long as this object does.
				      *
				      * You may also pass a function describing
				      * the weight to the integral (see the
				      * general docs for more information). The
				      * same applies for this object as said
				      * above.
				      */
    LaplaceMatrix (const Function<dim> * const rhs = 0,
		   const Function<dim> * const a = 0);

				     /**
				      * Assemble the cell matrix and right hand
				      * side vector for this cell. You need to
				      * give a right hand side object to the
				      * constructor to use this function. If
				      * a coefficient was given to the
				      * constructor, it is used.
				      */
    virtual void assemble (FullMatrix<double>  &cell_matrix,
			   Vector<double>      &rhs,
			   const FEValues<dim> &fe_values,
			   const typename DoFHandler<dim>::cell_iterator &) const;

				     /**
				      * Construct the cell matrix for this cell.
				      * If a coefficient was given to the
				      * constructor, it is used.
				      */
    virtual void assemble (FullMatrix<double>  &cell_matrix,
			   const FEValues<dim> &fe_values,
			   const typename DoFHandler<dim>::cell_iterator &) const;

				     /**
				      * Only construct the right hand side
				      * vector for this cell. You need to give
				      * a right hand side function to the
				      * constructor in order to call this
				      * function.
				      */
    virtual void assemble (Vector<double>      &rhs,
			   const FEValues<dim> &fe_values,
			   const typename DoFHandler<dim>::cell_iterator &) const;

				     /**
				      * Exception
				      */
    DeclException0 (ExcNoRHSSelected);
    
  protected:
				     /**
				      * Pointer to a function describing the
				      * right hand side of the problem. Should
				      * be zero if not given to the constructor
				      * and should then not be used.
				      */
    const Function<dim> * const right_hand_side;

    				     /**
				      * Pointer to a function describing the
				      * coefficient to the integral for the
				      * matrix entries. Should
				      * be zero if not given to the constructor
				      * and should then not be used.
				      */
    const Function<dim> * const coefficient;
};





template <int dim>
MassMatrix<dim>::MassMatrix (const Function<dim> * const rhs,
			     const Function<dim> * const a) :
		Equation<dim> (1),
		right_hand_side (rhs),
		coefficient (a)
{};



template <int dim>
void MassMatrix<dim>::assemble (FullMatrix<double>      &cell_matrix,
				const FEValues<dim>     &fe_values,
				const typename DoFHandler<dim>::cell_iterator &) const
{
  const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
		     n_q_points    = fe_values.n_quadrature_points;
  const FiniteElement<dim>    &fe  = fe_values.get_fe();
  const unsigned int n_components  = fe.n_components();

  Assert (cell_matrix.n() == dofs_per_cell,
	  Equation<dim>::ExcWrongSize(cell_matrix.n(), dofs_per_cell));
  Assert (cell_matrix.m() == dofs_per_cell,
	  Equation<dim>::ExcWrongSize(cell_matrix.m(), dofs_per_cell));
  Assert (cell_matrix.all_zero(),
	  Equation<dim>::ExcObjectNotEmpty());
  
  const FullMatrix<double> &values    = fe_values.get_shape_values ();
  const std::vector<double>     &weights   = fe_values.get_JxW_values ();


  if (coefficient != 0)
    {
      if (coefficient->n_components == 1)
					 // scalar coefficient given
	{
	  std::vector<double> coefficient_values (fe_values.n_quadrature_points);
	  coefficient->value_list (fe_values.get_quadrature_points(),
				   coefficient_values);
	  for (unsigned int i=0; i<dofs_per_cell; ++i) 
	    for (unsigned int j=0; j<dofs_per_cell; ++j)
	      if ((n_components == 1)
		  ||
		  (fe.system_to_component_index(i).first ==
		   fe.system_to_component_index(j).first))
		{
		  for (unsigned int point=0; point<n_q_points; ++point)
		    cell_matrix(i,j) += (values(i,point) *
					 values(j,point) *
					 weights[point] *
					 coefficient_values[point]);
		};
	}
      else
					 // vectorial coefficient
					 // given
	{
	  std::vector<Vector<double> > coefficient_values (fe_values.n_quadrature_points,
						      Vector<double>(n_components));
	  coefficient->vector_value_list (fe_values.get_quadrature_points(),
					  coefficient_values);
	  for (unsigned int i=0; i<dofs_per_cell; ++i) 
	    for (unsigned int j=0; j<dofs_per_cell; ++j)
	      if ((n_components == 1)
		  ||
		  (fe.system_to_component_index(i).first ==
		   fe.system_to_component_index(j).first))
		{
		  for (unsigned int point=0; point<n_q_points; ++point)
		    cell_matrix(i,j) += (values(i,point) *
					 values(j,point) *
					 weights[point] *
					 coefficient_values[point](
					   fe.system_to_component_index(i).first));
		};
	};
      
    }
  else
				     // no coefficient given
    for (unsigned int i=0; i<dofs_per_cell; ++i) 
      for (unsigned int j=0; j<dofs_per_cell; ++j)
	if ((n_components == 1)
	    ||
	    (fe.system_to_component_index(i).first ==
	     fe.system_to_component_index(j).first))
	  {
	    for (unsigned int point=0; point<n_q_points; ++point)
	      cell_matrix(i,j) += (values(i,point) *
				   values(j,point) *
				   weights[point]);
	  };
};



template <int dim>
void MassMatrix<dim>::assemble (FullMatrix<double>  &cell_matrix,
				Vector<double>      &rhs,
				const FEValues<dim> &fe_values,
				const typename DoFHandler<dim>::cell_iterator &) const
{
  Assert (right_hand_side != 0, ExcNoRHSSelected());

  const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
		     n_q_points    = fe_values.n_quadrature_points;
  const FiniteElement<dim>    &fe  = fe_values.get_fe();
  const unsigned int n_components  = fe.n_components();

				   // for system elements: not
				   // implemented at present
  Assert (n_components==1, ExcNotImplemented());
  
  Assert (cell_matrix.n() == dofs_per_cell,
	  Equation<dim>::ExcWrongSize(cell_matrix.n(), dofs_per_cell));
  Assert (cell_matrix.m() == dofs_per_cell,
	  Equation<dim>::ExcWrongSize(cell_matrix.m(), dofs_per_cell));
  Assert (rhs.size() == dofs_per_cell,
	  Equation<dim>::ExcWrongSize(rhs.size(), dofs_per_cell));
  Assert (cell_matrix.all_zero(),
	  Equation<dim>::ExcObjectNotEmpty());
  Assert (rhs.all_zero(),
	  Equation<dim>::ExcObjectNotEmpty());

  const FullMatrix<double> &values    = fe_values.get_shape_values ();
  const std::vector<double>     &weights   = fe_values.get_JxW_values ();
  std::vector<double>            rhs_values (fe_values.n_quadrature_points);
  right_hand_side->value_list (fe_values.get_quadrature_points(), rhs_values);

  if (coefficient != 0)
    {
      std::vector<double> coefficient_values (n_q_points);
      coefficient->value_list (fe_values.get_quadrature_points(),
			       coefficient_values);
      for (unsigned int point=0; point<n_q_points; ++point)
	for (unsigned int i=0; i<dofs_per_cell; ++i) 
	  {
	    for (unsigned int j=0; j<dofs_per_cell; ++j)
	      cell_matrix(i,j) += (values(i,point) *
				   values(j,point) *
				   weights[point] *
				   coefficient_values[point]);
	    rhs(i) += values(i,point) *
		      rhs_values[point] *
		      weights[point];
	  };
    }
  else
    for (unsigned int point=0; point<n_q_points; ++point)
      for (unsigned int i=0; i<dofs_per_cell; ++i) 
	{
	  for (unsigned int j=0; j<dofs_per_cell; ++j)
	    cell_matrix(i,j) += (values(i,point) *
				 values(j,point) *
				 weights[point]);
	  rhs(i) += values(i,point) *
		    rhs_values[point] *
		    weights[point];
	};
};



template <int dim>
void MassMatrix<dim>::assemble (Vector<double>      &rhs,
				const FEValues<dim> &fe_values,
				const typename DoFHandler<dim>::cell_iterator &) const
{
  Assert (right_hand_side != 0, ExcNoRHSSelected());

  const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
		     n_q_points    = fe_values.n_quadrature_points;
  const FiniteElement<dim>    &fe  = fe_values.get_fe();
  const unsigned int n_components  = fe.n_components();

				   // for system elements: not
				   // implemented at present
  Assert (n_components==1, ExcNotImplemented());

  Assert (rhs.size() == dofs_per_cell,
	  Equation<dim>::ExcWrongSize(rhs.size(), dofs_per_cell));
  Assert (rhs.all_zero(),
	  Equation<dim>::ExcObjectNotEmpty());

  const FullMatrix<double> &values    = fe_values.get_shape_values ();
  const std::vector<double>     &weights   = fe_values.get_JxW_values ();
  std::vector<double>            rhs_values(fe_values.n_quadrature_points);
  right_hand_side->value_list (fe_values.get_quadrature_points(), rhs_values);

  for (unsigned int point=0; point<n_q_points; ++point)
    for (unsigned int i=0; i<dofs_per_cell; ++i) 
      rhs(i) += values(i,point) *
		rhs_values[point] *
		weights[point];
};





template <int dim>
LaplaceMatrix<dim>::LaplaceMatrix (const Function<dim> * const rhs,
				   const Function<dim> * const a) :
		Equation<dim> (1),
		right_hand_side (rhs),
		coefficient (a) {};



template <int dim>
void LaplaceMatrix<dim>::assemble (FullMatrix<double>         &cell_matrix,
				   Vector<double>             &rhs,
				   const FEValues<dim>        &fe_values,
				   const typename DoFHandler<dim>::cell_iterator &) const
{
  Assert (right_hand_side != 0, ExcNoRHSSelected());
  
  const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
		     n_q_points    = fe_values.n_quadrature_points;
  const FiniteElement<dim>    &fe  = fe_values.get_fe();
  const unsigned int n_components  = fe.n_components();

				   // for system elements: might be
				   // not so useful, not implemented
				   // at present
  Assert (n_components==1, ExcNotImplemented());

  Assert (cell_matrix.n() == dofs_per_cell,
	  Equation<dim>::ExcWrongSize(cell_matrix.n(), dofs_per_cell));
  Assert (cell_matrix.m() == dofs_per_cell,
	  Equation<dim>::ExcWrongSize(cell_matrix.m(), dofs_per_cell));
  Assert (rhs.size() == dofs_per_cell,
	  Equation<dim>::ExcWrongSize(rhs.size(), dofs_per_cell));
  Assert (cell_matrix.all_zero(),
	  Equation<dim>::ExcObjectNotEmpty());
  Assert (rhs.all_zero(),
	  Equation<dim>::ExcObjectNotEmpty());

  const std::vector<std::vector<Tensor<1,dim> > >&gradients = fe_values.get_shape_grads ();
  const FullMatrix<double>             &values    = fe_values.get_shape_values ();
  std::vector<double>        rhs_values(fe_values.n_quadrature_points);
  const std::vector<double> &weights   = fe_values.get_JxW_values ();
  right_hand_side->value_list (fe_values.get_quadrature_points(), rhs_values);

  if (coefficient != 0)
    {
      std::vector<double> coefficient_values(n_q_points);
      coefficient->value_list (fe_values.get_quadrature_points(),
			       coefficient_values);
      for (unsigned int point=0; point<n_q_points; ++point)
	for (unsigned int i=0; i<dofs_per_cell; ++i) 
	  {
	    for (unsigned int j=0; j<dofs_per_cell; ++j)
	      cell_matrix(i,j) += (gradients[i][point] *
				   gradients[j][point]) *
				  weights[point] *
				  coefficient_values[point];
	    rhs(i) += values(i,point) *
		      rhs_values[point] *
		      weights[point];
	  };
    }
  else
    for (unsigned int point=0; point<n_q_points; ++point)
      for (unsigned int i=0; i<dofs_per_cell; ++i) 
	{
	  for (unsigned int j=0; j<dofs_per_cell; ++j)
	    cell_matrix(i,j) += (gradients[i][point] *
				 gradients[j][point]) *
				weights[point];
	  rhs(i) += values(i,point) *
		    rhs_values[point] *
		    weights[point];
	};

};



template <int dim>
void LaplaceMatrix<dim>::assemble (FullMatrix<double>  &cell_matrix,
				   const FEValues<dim> &fe_values,
				   const typename DoFHandler<dim>::cell_iterator &) const
{
  const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
		     n_q_points    = fe_values.n_quadrature_points;

  const FiniteElement<dim>    &fe  = fe_values.get_fe();
  const unsigned int n_components  = fe.n_components();

				   // for system elements: might be
				   // not so useful, not implemented
				   // at present
  Assert ((n_components==1) || (coefficient==0), ExcNotImplemented());

  Assert (cell_matrix.n() == dofs_per_cell,
	  Equation<dim>::ExcWrongSize(cell_matrix.n(), dofs_per_cell));
  Assert (cell_matrix.m() == dofs_per_cell,
	  Equation<dim>::ExcWrongSize(cell_matrix.m(), dofs_per_cell));
  Assert (cell_matrix.all_zero(),
	  Equation<dim>::ExcObjectNotEmpty());
  
  const std::vector<std::vector<Tensor<1,dim> > >&gradients = fe_values.get_shape_grads ();
  const std::vector<double> &weights   = fe_values.get_JxW_values ();
   
  if (coefficient != 0)
    {
      std::vector<double> coefficient_values(n_q_points);
      coefficient->value_list (fe_values.get_quadrature_points(),
			       coefficient_values);
      for (unsigned int point=0; point<n_q_points; ++point)
	for (unsigned int i=0; i<dofs_per_cell; ++i) 
	  for (unsigned int j=0; j<dofs_per_cell; ++j)
	    cell_matrix(i,j) += (gradients[i][point] *
				 gradients[j][point]) *
				weights[point] *
				coefficient_values[point];
    }
  else
    for (unsigned int i=0; i<dofs_per_cell; ++i) 
      for (unsigned int j=0; j<dofs_per_cell; ++j)
	if ((n_components==1)
	    ||
	    (fe.system_to_component_index(i).first ==
	     fe.system_to_component_index(j).first))
	  {
	    for (unsigned int point=0; point<n_q_points; ++point)
	      cell_matrix(i,j) += (gradients[i][point] *
				   gradients[j][point]) *
				  weights[point];
	  };
};



template <int dim>
void LaplaceMatrix<dim>::assemble (Vector<double>      &rhs,
				   const FEValues<dim> &fe_values,
				   const typename DoFHandler<dim>::cell_iterator &) const
{
  Assert (right_hand_side != 0, ExcNoRHSSelected());

  const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
		     n_q_points    = fe_values.n_quadrature_points;
  const FiniteElement<dim>    &fe  = fe_values.get_fe();
  const unsigned int n_components  = fe.n_components();

				   // for system elements: might be
				   // not so useful, not implemented
				   // at present
  Assert (n_components==1, ExcNotImplemented());

  Assert (rhs.size() == dofs_per_cell,
	  Equation<dim>::ExcWrongSize(rhs.size(), dofs_per_cell));
  Assert (rhs.all_zero(),
	  Equation<dim>::ExcObjectNotEmpty());

  const FullMatrix<double> &values    = fe_values.get_shape_values ();
  const std::vector<double>     &weights   = fe_values.get_JxW_values ();
  std::vector<double>        rhs_values(fe_values.n_quadrature_points);
  right_hand_side->value_list (fe_values.get_quadrature_points(), rhs_values);
   
  for (unsigned int point=0; point<n_q_points; ++point)
    for (unsigned int i=0; i<dofs_per_cell; ++i) 
      rhs(i) += values(i,point) *
		rhs_values[point] *
		weights[point];
};



#include <grid/tria_iterator.h>
#include <grid/tria_iterator.templates.h>
#include <fe/fe.h>
#include <lac/full_matrix.h>
#include <lac/vector.h>
#include <lac/sparse_matrix.h>
#include <base/quadrature.h>
#include <fe/mapping_q1.h>


// if necessary try to work around a bug in the IBM xlC compiler
#ifdef XLC_WORK_AROUND_STD_BUG
using namespace std;
#endif


//TODO: purge this variable
static const MappingQ1<deal_II_dimension> mapping;

template <int dim>
Assembler<dim>::AssemblerData::AssemblerData (const DoFHandler<dim>    &dof,
					      const bool                assemble_matrix,
					      const bool                assemble_rhs,
					      SparseMatrix<double>     &matrix,
					      Vector<double>           &rhs_vector,
					      const Quadrature<dim>    &quadrature,
					      const UpdateFlags        &update_flags) :
		dof(dof),
		assemble_matrix(assemble_matrix),
		assemble_rhs(assemble_rhs),
		matrix(matrix),
		rhs_vector(rhs_vector),
		quadrature(quadrature),
		update_flags(update_flags)
{};


template <int dim>
Assembler<dim>::Assembler (Triangulation<dim>  *tria,
			   const int            level,
			   const int            index,
			   const AssemblerData *local_data) :
		DoFCellAccessor<dim> (tria,level,index, &local_data->dof),
		cell_matrix (dof_handler->get_fe().dofs_per_cell),
		cell_vector (Vector<double>(dof_handler->get_fe().dofs_per_cell)),
		assemble_matrix (local_data->assemble_matrix),
		assemble_rhs (local_data->assemble_rhs),
		matrix(local_data->matrix),
		rhs_vector(local_data->rhs_vector),
		fe_values (mapping, dof_handler->get_fe(),
			   local_data->quadrature,
			   local_data->update_flags)
{
  Assert (!assemble_matrix || (matrix.m() == dof_handler->n_dofs()),
	  ExcInvalidData());
  Assert (!assemble_matrix || (matrix.n() == dof_handler->n_dofs()),
	  ExcInvalidData());
  Assert (!assemble_rhs || (rhs_vector.size()==dof_handler->n_dofs()),
	  ExcInvalidData());
};


template <int dim>
void Assembler<dim>::assemble (const Equation<dim> &equation) {
				   // re-init fe values for this cell
  fe_values.reinit (DoFHandler<dim>::cell_iterator (*this));
  const unsigned int n_dofs = dof_handler->get_fe().dofs_per_cell;

  if (assemble_matrix)
    cell_matrix.clear ();
  if (assemble_rhs)
    cell_vector.clear ();


// fill cell matrix and vector if required
  DoFHandler<dim>::cell_iterator this_cell (*this);
  if (assemble_matrix && assemble_rhs) 
    equation.assemble (cell_matrix, cell_vector, fe_values, this_cell);
  else
    if (assemble_matrix)
      equation.assemble (cell_matrix, fe_values, this_cell);
    else
      if (assemble_rhs)
	equation.assemble (cell_vector, fe_values, this_cell);
      else
	Assert (false, ExcNoAssemblingRequired());


// get indices of dofs
  std::vector<unsigned int> dofs (n_dofs);
  get_dof_indices (dofs);

				   // one could use the
				   // @p{distribute_local_to_global} functions
				   // here, but they would require getting the
				   // dof indices twice, so we leave it the
				   // way it was originally programmed.
  
				   // distribute cell matrix
  if (assemble_matrix)
    for (unsigned int i=0; i<n_dofs; ++i)
      for (unsigned int j=0; j<n_dofs; ++j)
	matrix.add(dofs[i], dofs[j], cell_matrix(i,j));

				   // distribute cell vector
  if (assemble_rhs)
    for (unsigned int j=0; j<n_dofs; ++j)
      rhs_vector(dofs[j]) += cell_vector(j);
};


// explicit instantiations
template class Assembler<deal_II_dimension>;

template class TriaRawIterator<deal_II_dimension,Assembler<deal_II_dimension> >;
template class TriaIterator<deal_II_dimension,Assembler<deal_II_dimension> >;
template class TriaActiveIterator<deal_II_dimension,Assembler<deal_II_dimension> >;



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

