/* $Id$ */
/* Author: Wolfgang Bangerth, University of Heidelberg, 2000 */

/*    $Id$       */
/*    Version: $Name$                                          */
/*                                                                */
/*    Copyright (C) 2000, 2001, 2002 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */

				 // As usual, the first few include
				 // files are already known, so we
				 // will not comment on them further.
#include <base/quadrature_lib.h>
#include <base/function.h>
#include <base/logstream.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/vector_memory.h>
#include <lac/precondition.h>
#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <grid/grid_generator.h>
#include <grid/grid_refinement.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary_lib.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <fe/fe_values.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <numerics/data_out.h>
#include <dofs/dof_constraints.h>
#include <numerics/error_estimator.h>

				 // In this example, we need
				 // vector-valued finite elements. The
				 // support for these can be found in
				 // the following include file:
#include <fe/fe_system.h>
				 // We will compose the vector-valued
				 // finite elements from regular Q1
				 // elements which can be found here,
				 // as usual:
#include <fe/fe_q.h>

				 // This again is C++:
#include <fstream>


				 // The main class is, except for its
				 // name, almost unchanged with
				 // respect to the step-6 example. The
				 // only change is the use of a
				 // different class for the ``fe''
				 // variable.
template <int dim>
class ElasticProblem 
{
  public:
    ElasticProblem ();
    ~ElasticProblem ();
    void run ();
    
  private:
    void setup_system ();
    void assemble_system ();
    void solve ();
    void refine_grid ();
    void output_results (const unsigned int cycle) const;

    Triangulation<dim>   triangulation;
    DoFHandler<dim>      dof_handler;

				     // Instead of a concrete finite
				     // element class such as
				     // ``FE_Q'', we now use a more
				     // generic one, ``FESystem''. In
				     // fact, it is not a finite
				     // element itself, but rather a
				     // class that can be used to
				     // stack several usual elements
				     // together to form one
				     // vector-valued finite
				     // element. In our case, we will
				     // compose the vector-valued
				     // element of ``FE_Q(1)'' objects,
				     // as shown below in the
				     // constructor of this class.
    FESystem<dim>        fe;

    ConstraintMatrix     hanging_node_constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double>       solution;
    Vector<double>       system_rhs;
};


				 // Before going over to the
				 // implementation of the main class,
				 // we declare and define the class
				 // which describes the right hand
				 // side. This time, the right hand
				 // side is vector-valued, as is the
				 // solution, so we will describe the
				 // new elements in some more detail.
template <int dim>
class RightHandSide :  public Function<dim> 
{
  public:
				     // The first thing is that
				     // vector-valued functions have a
				     // constructor, since they need
				     // to pass down to the base class
				     // of how many components the
				     // function consists. The default
				     // value in the constructor of
				     // the base class is one, so we
				     // need not define a constructor
				     // for the usual scalar function.
    RightHandSide ();
    
				     // The next function is a
				     // replacement for the ``value''
				     // function of the previous
				     // examples. There, a second
				     // parameter ``component'' was
				     // given, which denoted which
				     // component was requested. Here,
				     // we implement a function that
				     // returns the whole vector of
				     // values at the given place at
				     // once.
    virtual void vector_value (const Point<dim> &p,
			       Vector<double>   &values) const;

				     // Then, in analogy to the
				     // ``value_list'' function, there
				     // is a function
				     // ``vector_value_list'', which
				     // returns the values of the
				     // vector-valued function at
				     // several points at once:
    virtual void vector_value_list (const std::vector<Point<dim> > &points,
				    std::vector<Vector<double> >   &value_list) const;
};


				 // This is the constructor of the
				 // right hand side class. As said
				 // above, it only passes down to the
				 // base class the number of
				 // components, which is ``dim'' in
				 // the present case. Note that
				 // although the implementation is
				 // very short here, we do not move it
				 // into the class declaration, since
				 // our style guides require that
				 // inside the class declaration only
				 // declarations have to happen and
				 // that definitions are always to be
				 // found outside.
template <int dim>
RightHandSide<dim>::RightHandSide () :
		Function<dim> (dim)
{};


				 // This is the function that returns
				 // the whole vector of values at the
				 // point ``p'' at once:
template <int dim>
inline
void RightHandSide<dim>::vector_value (const Point<dim> &p,
				       Vector<double>   &values) const 
{
				   // To prevent cases where the
				   // return value has not previously
				   // been set to the right size
				   // (which is kind of a convention
				   // in the deal.II library), we test
				   // for this case and otherwise
				   // throw an exception:
  Assert (values.size() == dim, 
	  ExcDimensionMismatch (values.size(), dim));
				   // Likewise, if by some accident
				   // someone tried to compile and run
				   // the program in only one space
				   // dimension (in which the elastic
				   // equations do not make much sense
				   // since they reduce to the
				   // ordinary Laplace equation), we
				   // terminate the program if the
				   // dimension is not as expected.
  Assert (dim >= 2, ExcInternalError());
  
				   // The rest of the function is as
				   // would probably be expected given
				   // the form of the right hand side
				   // function. First we define the
				   // centers of the two points around
				   // which are the sources of
				   // x-displacement, i.e. (0.5,0) and
				   // (-0.5,0). Note that upon
				   // construction of the ``Point''
				   // objects, all components are set
				   // to zero.
  Point<dim> point_1, point_2;
  point_1(0) = 0.5;
  point_2(0) = -0.5;
  
				   // If now the point ``p'' is in the
				   // circle of radius 0.2 around one
				   // of these points, then set the
				   // force in x-direction to one,
				   // otherwise to zero:
  if (((p-point_1).square() < 0.2*0.2) ||
      ((p-point_2).square() < 0.2*0.2))
    values(0) = 1;
  else
    values(0) = 0;
  
				   // Likewise, if ``p'' is in the
				   // vicinity of the origin, then set
				   // the y-force to 1, otherwise to
				   // zero:
  if (p.square() < 0.2*0.2)
    values(1) = 1;
  else
    values(1) = 0;    
};



				 // Now, this is the function of the
				 // right hand side class that returns
				 // the values at several points at
				 // once.
template <int dim>
void RightHandSide<dim>::vector_value_list (const std::vector<Point<dim> > &points,
					    std::vector<Vector<double> >   &value_list) const 
{
				   // First we define an abbreviation
				   // for the number of points which
				   // we shall work on:
  const unsigned int n_points = points.size();

				   // Then we check whether the number
				   // of output slots has been set
				   // correctly, i.e. to the number of
				   // input points:
  Assert (value_list.size() == n_points, 
	  ExcDimensionMismatch (value_list.size(), n_points));

				   // Finally we treat each of the
				   // points. In one of the previous
				   // examples, we have explained why
				   // the
				   // ``value_list''/``vector_value_list''
				   // function had been introduced: to
				   // prevent us from calling virtual
				   // functions too frequently. On the
				   // other hand, we now need to
				   // implement the same function
				   // twice, which can lead to
				   // confusion if one function is
				   // changed but the other is
				   // not. However, we can prevent
				   // this situation using the
				   // following construct:
  for (unsigned int p=0; p<n_points; ++p)
    RightHandSide<dim>::vector_value (points[p],
				      value_list[p]);
				   // It calls the ``vector_value''
				   // function defined above for each
				   // point, and thus preempts all
				   // chances for inconsistency. It is
				   // important to note how the
				   // function was called: using the
				   // full class qualification using
				   // ``RightHandSide::'', since this
				   // calls the function directly and
				   // not using the virtual function
				   // table. The call is thus as fast
				   // as a call to any non-virtual
				   // function. In addition, we have
				   // declared the ``vector_value''
				   // function ``inline'', i.e. the
				   // compiler can remove the function
				   // call altogether and the
				   // resulting code can in principle
				   // be as fast as if we had
				   // duplicated the code.
};




template <int dim>
ElasticProblem<dim>::ElasticProblem () :
		dof_handler (triangulation),
						 // As said before, we
						 // would like to
						 // construct one
						 // vector-valued
						 // finite element as
						 // outer product of
						 // several scalar
						 // finite
						 // elements. Of
						 // course, the number
						 // of scalar finite
						 // element we would
						 // like to stack
						 // together equals
						 // the number of
						 // components the
						 // solution function
						 // has, which is
						 // ``dim'' since we
						 // consider
						 // displacement in
						 // each space
						 // direction. The
						 // ``FESystem'' class
						 // can handle this:
						 // we pass it the
						 // finite element of
						 // which we would
						 // like to compose
						 // the system of, and
						 // how often it shall
						 // be repeated:
		fe (FE_Q<dim>(1), dim)
				 // In fact, the ``FESystem'' class
				 // has several more constructors
				 // which can perform more complex
				 // operations than just stacking
				 // together several scalar finite
				 // elements of the same type into
				 // one; we will get to know these
				 // possibilities in later examples.
				 //
				 // It should be noted that the
				 // ``FESystem'' object thus created
				 // does not actually use the finite
				 // element which we have passed to it
				 // as first parameter. We could thus
				 // use an anonymous object created
				 // in-place. The ``FESystem''
				 // constructor only needs the
				 // parameter to generate a copy of
				 // the finite element from this.
{};



template <int dim>
ElasticProblem<dim>::~ElasticProblem () 
{
  dof_handler.clear ();
};


				 // Setting up the system of equations
				 // is equal to the function used in
				 // the step-6 example. The
				 // ``DoFHandler'' class and all other
				 // classes used take care of the
				 // vector-valuedness of the finite
				 // element themselves (in fact, the
				 // do not do so, since they only take
				 // care how many degrees of freedom
				 // there are per vertex, line and
				 // cell, and they do not ask what they
				 // represent, i.e. whether the finite
				 // element under consideration is
				 // vector-valued or whether it is,
				 // for example, a scalar Hermite
				 // element with several degrees of
				 // freedom on each vertex).
template <int dim>
void ElasticProblem<dim>::setup_system ()
{
  dof_handler.distribute_dofs (fe);
  hanging_node_constraints.clear ();
  DoFTools::make_hanging_node_constraints (dof_handler,
					   hanging_node_constraints);
  hanging_node_constraints.close ();
  sparsity_pattern.reinit (dof_handler.n_dofs(),
			   dof_handler.n_dofs(),
			   dof_handler.max_couplings_between_dofs());
				   // When making the sparsity
				   // pattern, there is some potential
				   // for optimization if not all
				   // components couple to all
				   // others. However, this is not the
				   // case for the elastic equations,
				   // so we use the standard call:
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);

  hanging_node_constraints.condense (sparsity_pattern);

  sparsity_pattern.compress();

  system_matrix.reinit (sparsity_pattern);

  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());
};


				 // The big changes in this program
				 // are in the creation of matrix and
				 // right hand side, since they are
				 // problem-dependent. We will go
				 // through that process step-by-step,
				 // since it is a bit more complicated
				 // than in previous examples.
template <int dim>
void ElasticProblem<dim>::assemble_system () 
{  
				   // First thing: the quadrature
				   // formula does not need
				   // modification since we still deal
				   // with bilinear functions.
  QGauss2<dim>  quadrature_formula;
				   // Also, the ``FEValues'' objects
				   // takes care of everything for us
				   // (or better: it does not really
				   // so; as in the comment in the
				   // function setting up the system,
				   // here as well the ``FEValues''
				   // object computes the same data on
				   // each cell, but it has some
				   // functionality to access data
				   // stored inside the finite element
				   // where they are precomputed upon
				   // construction).
  FEValues<dim> fe_values (fe, quadrature_formula, 
			   UpdateFlags(update_values    |
				       update_gradients |
				       update_q_points  |
				       update_JxW_values));

				   // The number of degrees of freedom
				   // per cell we now obviously ask
				   // from the composed finite element
				   // rather than from the underlying
				   // scalar Q1 element. Here, it is
				   // ``dim'' times the number of
				   // degrees of freedom per cell of
				   // the Q1 element, but this is not
				   // something we need to care about.
  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.n_quadrature_points;

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

				   // As was shown in previous
				   // examples as well, we need a
				   // place where to store the values
				   // of the coefficients at all the
				   // quadrature points on a cell. In
				   // the present situation, we have
				   // two coefficients, lambda and mu.
  std::vector<double>     lambda_values (n_q_points);
  std::vector<double>     mu_values (n_q_points);

				   // Well, we could as well have
				   // omitted the above two arrays
				   // since we will use constant
				   // coefficients for both lambda and
				   // mu, which can be declared like
				   // this. They both represent
				   // functions always returning the
				   // constant value 1.0. Although we
				   // could omit the respective
				   // factors in the assemblage of the
				   // matrix, we use them here for
				   // purpose of demonstration.
  ConstantFunction<dim> lambda(1.), mu(1.);

				   // Then again, we need to have the
				   // same for the right hand
				   // side. This is exactly as before
				   // in previous examples. However,
				   // we now have a vector-valued
				   // right hand side, which is why
				   // the data type of the
				   // ``rhs_values'' array is
				   // changed. We initialize it by
				   // ``n_q_points'' elements, each of
				   // which is a ``Vector<double>''
				   // with ``dim'' elements.
  RightHandSide<dim>      right_hand_side;
  std::vector<Vector<double> > rhs_values (n_q_points,
					   Vector<double>(dim));


				   // Now we can begin with the loop
				   // over all cells:
  typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
						 endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      cell_matrix.clear ();
      cell_rhs.clear ();

      fe_values.reinit (cell);

				       // As in previous examples, we
				       // define some abbreviations
				       // for the various data that
				       // the ``FEValues'' class
				       // offers:
      const FullMatrix<double> 
	& shape_values = fe_values.get_shape_values();
      const std::vector<std::vector<Tensor<1,dim> > >
	& shape_grads  = fe_values.get_shape_grads();
      const std::vector<double>
	& JxW_values   = fe_values.get_JxW_values();
      const std::vector<Point<dim> >
	& q_points     = fe_values.get_quadrature_points();
      
				       // Next we get the values of
				       // the coefficients at the
				       // quadrature points:
      lambda.value_list (q_points, lambda_values);
      mu.value_list     (q_points, mu_values);

				       // Then assemble the entries of
				       // the local stiffness matrix
				       // and right hand side
				       // vector. This follows almost
				       // one-to-one the pattern
				       // described in the
				       // introduction of this example
				       // and will not comment much on
				       // this.
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	{
					   // One of the few comments
					   // in place is how we access
					   // the function ``comp(i)''
					   // used in the
					   // introduction. This is
					   // possible as follows:
	  const unsigned int 
	    component_i = fe.system_to_component_index(i).first;
					   // By accessing the
					   // ``first'' variable of
					   // the return value of the
					   // ``system_to_component_index''
					   // function, you might
					   // already have guessed
					   // that there is more in
					   // it. In fact, the
					   // function returns a
					   // ``std::pair<unsigned int,
					   // unsigned int>'', of
					   // which the first element
					   // is ``comp(i)'' and the
					   // second is the value
					   // ``base(i)'' also noted
					   // in the text. You will
					   // rather seldom need to
					   // access this second
					   // value, but the first is
					   // important when using
					   // vector valued elements.
	  
	  for (unsigned int j=0; j<dofs_per_cell; ++j) 
	    {
	      const unsigned int 
		component_j = fe.system_to_component_index(j).first;
	      
	      for (unsigned int q_point=0; q_point<n_q_points;
		   ++q_point)
		{
						   // Now add up the
						   // contribution of
						   // this cell to the
						   // local matrix:
		  cell_matrix(i,j) 
		    += 
						     // This first term is
						     // ((lambda+mu) d_i u_i, d_j v_j).
						     // Note that
						     // ``shape_grads[i][q_point]''
						     // returns the
						     // gradient of
						     // the i-th shape
						     // function at
						     // quadrature
						     // point
						     // q_point. The
						     // component
						     // ``comp(i)'',
						     // which is the
						     // derivative of
						     // the i-th shape
						     // function with
						     // respect to the
						     // comp(i)th
						     // coordinate is
						     // accessed by
						     // the appended
						     // brackets.
		    (
		      (shape_grads[i][q_point][component_i] *
		       shape_grads[j][q_point][component_j] *
		       (lambda_values[q_point] +
			mu_values[q_point]))
		      +
						       // The second term is
						       // (mu nabla u_i, nabla v_j).
						       // We need not
						       // access a
						       // specific
						       // component of
						       // the
						       // gradient,
						       // since we
						       // only have to
						       // compute the
						       // scalar
						       // product of
						       // the two
						       // gradients,
						       // of which an
						       // overloaded
						       // version of
						       // the
						       // operator*
						       // takes care,
						       // as in
						       // previous
						       // examples.
						       //
						       // Note that by
						       // using the ?:
						       // operator, we
						       // only do this
						       // if comp(i)
						       // equals
						       // comp(j),
						       // otherwise a
						       // zero is
						       // added (which
						       // will be
						       // optimized
						       // away by the
						       // compiler).
		      ((component_i == component_j) ?
		       (shape_grads[i][q_point] *
			shape_grads[j][q_point] *
			mu_values[q_point])  :
		       0)
		    )
		    *
		    JxW_values[q_point];
		};
	    };
	};

				       // Assembling the right hand
				       // side is also just as
				       // discussed in the
				       // introduction. We will
				       // therefore not discuss it
				       // further.
      right_hand_side.vector_value_list (q_points, rhs_values);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	{
	  const unsigned int 
	    component_i = fe.system_to_component_index(i).first;
	  
	  for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	    cell_rhs(i) += shape_values(i,q_point) *
			   rhs_values[q_point](component_i) *
			   JxW_values[q_point];
	};

				       // The transfer from local
				       // degrees of freedom into the
				       // global matrix and right hand
				       // side vector does not depend
				       // on the equation under
				       // consideration, and is thus
				       // the same as in all previous
				       // examples.
      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	{
	  for (unsigned int j=0; j<dofs_per_cell; ++j)
	    system_matrix.add (local_dof_indices[i],
			       local_dof_indices[j],
			       cell_matrix(i,j));
	  
	  system_rhs(local_dof_indices[i]) += cell_rhs(i);
	};
    };

  hanging_node_constraints.condense (system_matrix);
  hanging_node_constraints.condense (system_rhs);

				   // The interpolation of the
				   // boundary values needs a small
				   // modification: since the solution
				   // function is vector-valued, so
				   // needs to be the boundary
				   // values. The ``ZeroFunction''
				   // constructor accepts a parameter
				   // that tells it that it shall
				   // represent a vector valued,
				   // constant zero function with that
				   // many components. By default,
				   // this parameter is equal to one,
				   // in which case the
				   // ``ZeroFunction'' object would
				   // represent a scalar
				   // function. Since the solution
				   // vector has ``dim'' components,
				   // we need to pass ``dim'' as
				   // number of components to the zero
				   // function as well.
  std::map<unsigned int,double> boundary_values;
  VectorTools::interpolate_boundary_values (dof_handler,
					    0,
					    ZeroFunction<dim>(dim),
					    boundary_values);
  MatrixTools::apply_boundary_values (boundary_values,
				      system_matrix,
				      solution,
				      system_rhs);
};



				 // The solver does not care about
				 // where the system of equations
				 // comes, as long as it stays
				 // positive definite and symmetric
				 // (which are the requirements for
				 // the use of the CG solver), which
				 // the system is. Therefore, we need
				 // not change anything.
template <int dim>
void ElasticProblem<dim>::solve () 
{
  SolverControl           solver_control (1000, 1e-12);
  PrimitiveVectorMemory<> vector_memory;
  SolverCG<>              cg (solver_control, vector_memory);

  PreconditionSSOR<> preconditioner;
  preconditioner.initialize(system_matrix, 1.2);

  cg.solve (system_matrix, solution, system_rhs,
	    preconditioner);

  hanging_node_constraints.distribute (solution);
};



				 // The function that does the
				 // refinement of the grid is the same
				 // as in the step-6 example. The
				 // quadrature formula is adapted to
				 // the linear elements again. Note
				 // that the error estimator by
				 // default adds up the estimated
				 // obtained from all components of
				 // the finite element solution, that
				 // is it uses the displacement in all
				 // directions with the same
				 // weight. If we would like the grid
				 // to be adapted to the
				 // x-displacement only, we could pass
				 // the function an additional
				 // parameter which tells it to do so
				 // and do not consider the
				 // displacements in all other
				 // directions for the error
				 // indicators.
template <int dim>
void ElasticProblem<dim>::refine_grid ()
{
  Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

  FunctionMap<dim>::type neumann_boundary;
  KellyErrorEstimator<dim>::estimate (dof_handler,
				      QGauss2<dim-1>(),
				      neumann_boundary,
				      solution,
				      estimated_error_per_cell);

  GridRefinement::refine_and_coarsen_fixed_number (triangulation,
						   estimated_error_per_cell,
						   0.3, 0.03);

  triangulation.execute_coarsening_and_refinement ();
};


				 // The output happens mostly as has
				 // been shown in previous examples
				 // already. The only difference is
				 // that the solution function is
				 // vector valued. The ``DataOut''
				 // class takes care of this
				 // automatically, but we have to give
				 // each component of the solution
				 // vector a different name.
template <int dim>
void ElasticProblem<dim>::output_results (const unsigned int cycle) const
{
  std::string filename = "solution-";
  filename += ('0' + cycle);
  Assert (cycle < 10, ExcInternalError());
  
  filename += ".gmv";
  std::ofstream output (filename.c_str());

  DataOut<dim> data_out;
  data_out.attach_dof_handler (dof_handler);

 

				   // As said above, we need a
				   // different name for each
				   // component of the solution
				   // function. To pass one name for
				   // each component, a vector of
				   // strings is used. Since the
				   // number of components is the same
				   // as the number of dimensions we
				   // are working in, the following
				   // ``switch'' statement is used.
				   //
				   // We note that some graphics
				   // programs have restriction as to
				   // what characters are allowed in
				   // the names of variables. The
				   // library therefore supports only
				   // the minimal subset of these
				   // characters that is supported by
				   // all programs. Basically, these
				   // are letters, numbers,
				   // underscores, and some other
				   // characters, but in particular no
				   // whitespace and minus/hyphen. The
				   // library will throw an exception
				   // otherwise, at least if in debug
				   // mode.
  std::vector<std::string> solution_names;
  switch (dim)
    {
      case 1:
	    solution_names.push_back ("displacement");
	    break;
      case 2:
	    solution_names.push_back ("x_displacement");	    
	    solution_names.push_back ("y_displacement");
	    break;
      case 3:
	    solution_names.push_back ("x_displacement");	    
	    solution_names.push_back ("y_displacement");
	    solution_names.push_back ("z_displacement");
	    break;
					     // It is good style to
					     // let the program die if
					     // we run upon a case
					     // which we did not
					     // consider. Remember
					     // that the ``Assert''
					     // macro throws an
					     // exception if the
					     // condition in the first
					     // parameter is not
					     // satisfied. Of course,
					     // the condition
					     // ``false'' can never be
					     // satisfied, so the
					     // program will always
					     // abort whenever it gets
					     // to this statement:
      default:
	    Assert (false, ExcInternalError());
    };
	     
				   // After setting up the names for
				   // the different components of the
				   // solution vector, we can add the
				   // solution vector to the list of
				   // data vectors scheduled for
				   // output. Note that the following
				   // function takes a vector of
				   // strings as second argument,
				   // whereas the one which we have
				   // used in all previous examples
				   // accepted a string there. In
				   // fact, the latter function is
				   // only a shortcut for the function
				   // which we call here: it puts the
				   // single string that is passed to
				   // it into a vector of strings with
				   // only one element and forwards
				   // that to the other function.
  data_out.add_data_vector (solution, solution_names);
  data_out.build_patches ();
  data_out.write_gmv (output);
};



template <int dim>
void ElasticProblem<dim>::run () 
{
  for (unsigned int cycle=0; cycle<8; ++cycle)
    {
      std::cout << "Cycle " << cycle << ':' << std::endl;

      if (cycle == 0)
	{
					   // As in previous examples,
					   // we use the unit square
					   // (or cube) as domain.
	  GridGenerator::hyper_cube (triangulation, -1, 1);
					   // This time, we have to
					   // refine the coarse grid
					   // twice before we first
					   // solve on it. The reason
					   // is the following: we use
					   // the ``Gauss2''
					   // quadrature formula for
					   // integration of the right
					   // hand side; that means
					   // that there are four
					   // quadrature points on
					   // each cell (in 2D). If we
					   // only refine the initial
					   // grid once globally, then
					   // there will be only four
					   // quadrature points in
					   // each direction on the
					   // domain. However, the
					   // right hand side function
					   // was chosen to be rather
					   // localized and in that
					   // case all quadrature
					   // points lie outside the
					   // support of the right
					   // hand side function. The
					   // right hand side vector
					   // will then contain only
					   // zeroes and the solution
					   // of the system of
					   // equations is the zero
					   // vector, i.e. a finite
					   // element function that it
					   // zero everywhere. We
					   // should not be surprised
					   // about such things
					   // happening, since we have
					   // chosen an initial grid
					   // that is totally
					   // unsuitable for the
					   // problem at hand.
					   //
					   // The unfortunate thing is
					   // that if the discrete
					   // solution is constant,
					   // then the error
					   // indicators computed by
					   // the
					   // ``KellyErrorEstimator''
					   // class are zero for each
					   // cell as well, and the
					   // call to
					   // ``refine_and_coarsen_fixed_number''
					   // on the ``triangulation''
					   // object will not flag any
					   // cells for refinement
					   // (why should it if the
					   // indicated error is zero
					   // for each cell?). The
					   // grid in the next
					   // iteration will therefore
					   // consist of four cells
					   // only as well, and the
					   // same problem occurs
					   // again.
					   //
					   // The conclusion needs to
					   // be: while of course we
					   // will not choose the
					   // initial grid to be
					   // well-suited for the
					   // accurate solution of the
					   // problem, we must at
					   // least choose it such
					   // that it has the chance
					   // to capture the most
					   // striking features of the
					   // solution. In this case,
					   // it needs to be able to
					   // see the right hand
					   // side. Thus, we refine
					   // twice globally. (Note
					   // that the
					   // ``refine_global''
					   // function is not part of
					   // the ``GridRefinement''
					   // class in which
					   // ``refine_and_coarsen_fixed_number''
					   // is declared, for
					   // example. The reason is
					   // first that it is not an
					   // algorithm that computed
					   // refinement flags from
					   // indicators, but more
					   // importantly that it
					   // actually performs the
					   // refinement, in contrast
					   // to the functions in
					   // ``GridRefinement'' that
					   // only flag cells without
					   // actually refining the
					   // grid.)
	  triangulation.refine_global (2);
	}
      else
	refine_grid ();

      std::cout << "   Number of active cells:       "
		<< triangulation.n_active_cells()
		<< std::endl;

      setup_system ();

      std::cout << "   Number of degrees of freedom: "
		<< dof_handler.n_dofs()
		<< std::endl;
      
      assemble_system ();
      solve ();
      output_results (cycle);
    };
};


				 // The main function is again exactly
				 // like in step-6 (apart from the
				 // changed class names, of course).
int main () 
{
  try
    {
      deallog.depth_console (0);

      ElasticProblem<2> elastic_problem_2d;
      elastic_problem_2d.run ();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      std::cerr << "Exception on processing: " << std::endl
		<< exc.what() << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      
      return 1;
    }
  catch (...) 
    {
      std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      std::cerr << "Unknown exception!" << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      return 1;
    };

  return 0;
};
