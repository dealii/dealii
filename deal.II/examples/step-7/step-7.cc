/* $Id$ */
/* Author: Wolfgang Bangerth and Ralf Hartmann, University of Heidelberg, 2000 */

				 // These first include files have all
				 // been treated in previous examples,
				 // so we won't explain what is in
				 // them again.
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
#include <grid/grid_generator.h>
#include <grid/grid_refinement.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary_lib.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_constraints.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <fe/fe_q.h>
#include <numerics/matrices.h>
#include <numerics/error_estimator.h>
#include <numerics/data_out.h>

				 // In this example, we will not use
				 // the numeration scheme which is
				 // used per default by the
				 // ``DoFHandler'' class, but will
				 // renumber them using the
				 // Cuthill-McKee algorithm. The
				 // necessary functions are declared
				 // in the following file:
#include <numerics/dof_renumbering.h>
				 // Then we will show a little trick
				 // how we can make sure that objects
				 // are not deleted while they are
				 // still in use. For this purpose,
				 // there is the ``SmartPointer''
				 // helper class, which is declared in
				 // this file:
#include <base/smartpointer.h>
				 // Then we will want to use the
				 // ``integrate_difference'' function
				 // mentioned in the introduction. It
				 // comes from this file:
#include <numerics/vectors.h>
				 // We are going to use a
				 // ``ConvergenceTable'' that collects
				 // all important data during a run
				 // and prints it at the end as a
				 // table.
#include <base/convergence_table.h>
				 // And finally, we need to use the
				 // ``FEFaceValues'' class, which is
				 // declare in the same file as the
				 // ``FEValues'' class:
#include <fe/fe_values.h>

				 // We need one more include from
				 // standard C++, which is necessary
				 // when we try to find out the actual
				 // type behind a pointer to a base
				 // class. We will explain this in
				 // slightly more detail below.
#include <typeinfo>
#include <fstream>



				 // Since we want to compare the
				 // exactly known continuous solution
				 // to the computed one, we need a
				 // function object which represents
				 // the continuous solution. On the
				 // other hand, we need the right hand
				 // side function, and that one of
				 // course shares some characteristics
				 // with the solution. In order to
				 // reduce dependencies which arise if
				 // we have to change something in
				 // both classes at the same time, we
				 // exclude the common characteristics
				 // of both functions into a base
				 // class.
				 //
				 // The common characteristics for the
				 // given solution, which as explained
				 // in the introduction is a sum of
				 // three exponentials, are here: the
				 // number of exponentials, their
				 // centers, and their half width. We
				 // declare them in the following
				 // class. Since the number of
				 // exponentials is a constant scalar
				 // integral quantity, C++ allows its
				 // definition (i.e. assigning a
				 // value) right at the place of
				 // declaration (i.e. where we declare
				 // that such a variable exists).
template <int dim>
class SolutionBase 
{
  protected:
    static const unsigned int n_source_centers = 3;    
    static const Point<dim>   source_centers[n_source_centers];
    static const double       width;
};


				 // The variables which denote the
				 // centers and the width of the
				 // exponentials have just been
				 // declared, now we still need to
				 // assign values to them. Here, we
				 // can show another small piece of
				 // template sorcery, namely how we
				 // can assign different values to
				 // these variables depending on the
				 // dimension. We will only use the 2d
				 // case in the program, but we show
				 // the 1d case for exposition of a
				 // useful technique.
				 //
				 // First we assign values to the
				 // centers for the 1d case, where we
				 // place the centers equidistantly at
				 // -1/3, 0, and 1/3:
template <>
const Point<1>
SolutionBase<1>::source_centers[SolutionBase<1>::n_source_centers]
= { Point<1>(-1.0 / 3.0), 
    Point<1>(0.0), 
    Point<1>(+1.0 / 3.0)   };

				 // Then we place the centers for the
				 // 2d case as follows:
template <>
const Point<2>
SolutionBase<2>::source_centers[SolutionBase<2>::n_source_centers]
= { Point<2>(-0.5, +0.5), 
    Point<2>(-0.5, -0.5), 
    Point<2>(+0.5, -0.5)   };

				 // There remains to assign a value to
				 // the half-width of the
				 // exponentials. We would like to use
				 // the same value for all dimensions,
				 // so here is how that works:
template <int dim>
const double SolutionBase<dim>::width = 1./3.;



				 // After declaring and defining the
				 // characteristics of solution and
				 // right hand side, we can declare
				 // the classes representing these
				 // two. They both represent
				 // continuous functions, so they are
				 // derived from the ``Function<dim>''
				 // base class, and they also inherit
				 // the characteristics defined in the
				 // ``SolutionBase'' class.
				 //
				 // The actual classes are declared in
				 // the following. Note that in order
				 // to compute the error of the
				 // numerical solution against the
				 // continuous one in the L2 and H1
				 // norms, we have to export value and
				 // gradient of the exact solution,
				 // which is done by overloading the
				 // respective virtual member
				 // functions in the ``Function'' base
				 // class.
template <int dim>
class Solution : public Function<dim>,
		 protected SolutionBase<dim>
{
  public:
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;
    virtual Tensor<1,dim> gradient (const Point<dim>   &p,
				    const unsigned int  component = 0) const;
};


				 // The actual definition of the
				 // values and gradients of the exact
				 // solution class is according to
				 // their mathematical definition and
				 // probably needs not much
				 // explanation.
template <int dim>
double Solution<dim>::value (const Point<dim>   &p,
			     const unsigned int) const
{
  double return_value = 0;
  for (unsigned int i=0; i<n_source_centers; ++i)
    {
				       // One of the few things worth
				       // mentioning is the following
				       // variables, which represents
				       // the vector (x-x_i). It is
				       // computed in the way that one
				       // would intuitively expect:
      const Point<dim> shifted_point = p-source_centers[i];
      
				       // The ``Point<dim>'' class
				       // offers a member function
				       // ``square'' that does what
				       // it's name suggests.
      return_value += exp(-shifted_point.square() / (width*width));
    };
  
  return return_value;
};



template <int dim>
Tensor<1,dim> Solution<dim>::gradient (const Point<dim>   &p,
				       const unsigned int) const
{
				   // In order to accumulate the
				   // gradient from the contributions
				   // of the exponentials, we allocate
				   // an object which denotes the
				   // mathematical quantity of a
				   // tensor of rank ``1'' and
				   // dimension ``dim''. Its default
				   // constructor sets it to the
				   // vector containing only zeroes,
				   // so we need not explicitly care
				   // for its initialization.
  Tensor<1,dim> return_value;
				   // Note that we could as well have
				   // taken the type of the object to
				   // be ``Point<dim>''. Tensors of
				   // rank 1 and points are almost
				   // exchangeable, and have only very
				   // slightly different mathematical
				   // meanings. In fact, the
				   // ``Point<dim>'' class is derived
				   // from the ``Tensor<1,dim>''
				   // class, which makes up for their
				   // mutual exchange ability.

  for (unsigned int i=0; i<n_source_centers; ++i)
    {
      const Point<dim> shifted_point = p-source_centers[i];
      
				       // For the gradient, note that
				       // it's direction is along
				       // (x-x_i), so we add up
				       // multiples of this distance
				       // vector, where the factor is
				       // given by the exponentials.
      return_value += (-2 / (width*width) *
		       exp(-shifted_point.square() / (width*width)) *
		       shifted_point);
    };
  
  return return_value;
};



				 // Besides the function that
				 // represents the exact solution, we
				 // also need a function which we can
				 // use as right hand side when
				 // assembling the linear system of
				 // discretized equations. This is
				 // accomplished using the following
				 // class and the following definition
				 // of its function. Note that here we
				 // only need the value of the
				 // function, not its gradients or
				 // higher derivatives.
template <int dim>
class RightHandSide : public Function<dim>,
		      protected SolutionBase<dim>
{
  public:
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;
};


				 // The value of the right hand side
				 // is given by the negative Laplacian
				 // of the solution plus the solution
				 // itself, since we wanted to solve
				 // Helmholtz's equation:
template <int dim>
double RightHandSide<dim>::value (const Point<dim>   &p,
				  const unsigned int) const
{
  double return_value = 0;
  for (unsigned int i=0; i<n_source_centers; ++i)
    {
      const Point<dim> shifted_point = p-source_centers[i];
      
				       // The first contribution is
				       // the Laplacian:
      return_value += ((2*dim - 4*shifted_point.square()/(width*width)) / 
		       (width*width) *
		       exp(-shifted_point.square() / (width*width)));
				       // And the second is the
				       // solution itself:
      return_value += exp(-shifted_point.square() / (width*width));
    };
  
  return return_value;
};



				 // Then we need the class that does
				 // all the work. It is mostly the
				 // same as in previous examples, and
				 // we will discuss the differences
				 // only when we declare the
				 // respective functions or variables
				 // below.
template <int dim>
class LaplaceProblem 
{
  public:
				     // We will use this class in
				     // several modes: for different
				     // finite elements, as well as
				     // for adaptive and global
				     // refinement. The decision
				     // whether global or adaptive
				     // refinement shall be used is
				     // communicated to the
				     // constructor of this class
				     // through an enumeration type,
				     // which we declare here:
    enum RefinementMode {
	  global_refinement, adaptive_refinement
    };
    
				     // This is the constructor of the
				     // class, it takes the finite
				     // element and the refinement
				     // mode as parameter and stores
				     // them in local variables.
    LaplaceProblem (const FiniteElement<dim> &fe,
		    const RefinementMode      refinement_mode);

				     // The following two functions
				     // are the same as in previous
				     // examples.
    ~LaplaceProblem ();

    void run ();
    
  private:
				     // As are these:
    void setup_system ();
    void assemble_system ();
    void solve ();
    void refine_grid ();

				     // After the solution has been
				     // computed, we perform some
				     // analysis on it, such as
				     // computing the error in various
				     // norms. This is done in the
				     // following function. To enable
				     // some output, we pass it the
				     // number of the refinement
				     // cycle.
    void process_solution (const unsigned int cycle);

				     // Now for the data elements of
				     // this class:
    Triangulation<dim>                      triangulation;
    DoFHandler<dim>                         dof_handler;

				     // The finite elements which the
				     // objects of this class operate
				     // on are passed to the
				     // constructor of this class. It
				     // has to store a pointer to the
				     // finite element for the member
				     // functions to use. Now, for the
				     // present class there is no big
				     // deal in that, but since we
				     // want to show techniques rather
				     // than solutions in these
				     // programs, we will here point
				     // out a problem that often
				     // occurs -- and of course the
				     // right solution as well.
				     //
				     // Consider the following
				     // situation that occurs in all
				     // the example programs: we have
				     // a triangulation object, and we
				     // have a finite element object,
				     // and we also have an object of
				     // type ``DoFHandler'' that uses
				     // both of the first two. These
				     // three objects all have a
				     // lifetime that is rather long
				     // compared to most other
				     // objects: they are basically
				     // set at the beginning of the
				     // program or an outer loop, and
				     // they are destroyed at the very
				     // end. The question is: can we
				     // guarantee that the two objects
				     // which the ``DoFHandler'' uses,
				     // live at least as long as they
				     // are in use? This means that
				     // the ``DoFHandler'' must have a
				     // kind of lock on the
				     // destruction of the other
				     // objects, and it can only
				     // release this lock once it has
				     // cleared all active references
				     // to these objects. We have seen
				     // what happens if we violate
				     // this order of destruction in
				     // the previous example program:
				     // an exception is thrown that
				     // terminates the program in
				     // order to notify the programmer
				     // of this potentially dangerous
				     // state where an object is
				     // pointed to that no longer
				     // persists.
				     //
				     // We will show here how the
				     // library managed to find out
				     // that there are still active
				     // references to an
				     // object. Basically, the method
				     // is along the following line:
				     // all objects that are subject
				     // to such potentially dangerous
				     // pointers are derived from a
				     // class called
				     // ``Subscriptor''. For example,
				     // the ``Triangulation'',
				     // ``DoFHandler'', and a base
				     // class of the ``FiniteElement''
				     // class are derived from
				     // ``Subscriptor``. This latter
				     // class does not offer much
				     // functionality, but it has a
				     // built-in counter which we can
				     // subscribe to, thus the name of
				     // the class. Whenever we
				     // initialize a pointer to that
				     // object, we can increase it use
				     // counter, and when we move away
				     // our pointer or do not need it
				     // any more, we decrease the
				     // counter again. This way, we
				     // can always check how many
				     // objects still use that
				     // object. If an object of a
				     // class that is derived from the
				     // ``Subscriptor'' class is
				     // destroyed, it also has to call
				     // the destructor of the
				     // ``Subscriptor'' class; this
				     // will then check whether the
				     // counter is really zero. If
				     // yes, then there are no active
				     // references to this object any
				     // more, and we can safely
				     // destroy it. If the counter is
				     // non-zero, however, then the
				     // destruction would result in
				     // stale and thus potentially
				     // dangerous pointers, and we
				     // rather throw an exception to
				     // alert the programmer that she
				     // is doing something dangerous
				     // and better had her program
				     // fixed.
				     //
				     // While this certainly all
				     // sounds very well, it has some
				     // problems in terms of
				     // usability: what happens if I
				     // forget to increase the counter
				     // when I let a pointer point to
				     // such an object? And what
				     // happens if I forget to
				     // decrease it again? Note that
				     // this may lead to extremely
				     // difficult to find bugs, since
				     // the place where we have
				     // forgotten something may be
				     // very far away from the place
				     // where the check for zeroness
				     // of the counter upon
				     // destruction actually
				     // fails. This kind of bug is
				     // very annoying and usually very
				     // hard to fix.
				     //
				     // The solution to this problem
				     // is to again use some C++
				     // trickery: we create a class
				     // that acts just like a pointer,
				     // i.e. can be dereferenced, can
				     // be assigned to and from other
				     // pointers, and so on. This can
				     // be done by overloading the
				     // several dereferencing
				     // operators of that
				     // class. Within the
				     // constructors, destructors, and
				     // assignment operators of that
				     // class, we can however also
				     // manage increasing or
				     // decreasing the use counters of
				     // the objects we point
				     // to. Objects of that class
				     // therefore can be used just
				     // like ordinary pointers to
				     // objects, but they also serve
				     // to change the use counters of
				     // those objects without the need
				     // for the programmer to do so
				     // herself. The class that
				     // actually does all this is
				     // called ``SmartPointer'' and
				     // takes as template parameter
				     // the data type of the object
				     // which it shall point to. The
				     // latter type may be any class,
				     // as long as it is derived from
				     // the ``Subscriptor'' class.
				     //
				     // In the present example
				     // program, we protect object
				     // using the pointer to the
				     // finite element, i.e. the
				     // following member variable,
				     // from the situation that for
				     // some reason the finite element
				     // pointed to is destroyed while
				     // still in use. Note that the
				     // pointer is assigned at
				     // construction time of this
				     // object, and destroyed upon
				     // destruction of this object, so
				     // the lock on the destruction of
				     // the finite element object is
				     // basically all through the
				     // lifetime of this object.
    SmartPointer<const FiniteElement<dim> > fe;

				     // The next few member variables
				     // are unspectacular, since they
				     // have already been discussed in
				     // detail:
    ConstraintMatrix                        hanging_node_constraints;

    SparsityPattern                         sparsity_pattern;
    SparseMatrix<double>                    system_matrix;

    Vector<double>                          solution;
    Vector<double>                          system_rhs;

				     // The second last variable
				     // stores the refinement mode
				     // passed to the
				     // constructor. Since it is only
				     // set in the constructor, we can
				     // declare this variable
				     // constant, to avoid that
				     // someone sets it involuntarily
				     // (e.g. in an `if'-statement
				     // where == was written as = by
				     // chance).
    const RefinementMode                    refinement_mode;

				     // For each refinement level some
				     // important data (like the
				     // number of cells, or the L2
				     // error of the numerical
				     // solution) is printed out. The
				     // ``TableHandler'' can be used
				     // to collect all this data and
				     // to output it at the end of the
				     // run as a table in a simple
				     // text format or in Tex
				     // format. Here we don't only use
				     // the ``TableHandler'' but we
				     // use the derived class
				     // ``ConvergenceTable'' that
				     // additionally evaluates rates
				     // of convergence.
    ConvergenceTable                        convergence_table;
};



				 // In the constructor of this class,
				 // we only set the variables passed
				 // to this object, and associate the
				 // DoF handler object with the
				 // triangulation (which is empty at
				 // present, however).
template <int dim>
LaplaceProblem<dim>::LaplaceProblem (const FiniteElement<dim> &fe,
				     const RefinementMode refinement_mode) :
		dof_handler (triangulation),
		fe (&fe),
		refinement_mode (refinement_mode)
{};



template <int dim>
LaplaceProblem<dim>::~LaplaceProblem () 
{
  dof_handler.clear ();
};


				 // The following function sets up the
				 // degrees of freedom, sizes of
				 // matrices and vectors, etc. Most of
				 // its functionality has been showed
				 // in previous examples, the only
				 // difference being the renumbering
				 // step.
template <int dim>
void LaplaceProblem<dim>::setup_system ()
{
  dof_handler.distribute_dofs (*fe);
				   // Renumbering the degrees of
				   // freedom is not overly difficult,
				   // as long as you use one of the
				   // algorithms included in the
				   // library. It requires just one
				   // line of code, namely the
				   // following:
  DoFRenumbering::Cuthill_McKee (dof_handler);
				   // Note, however, that when you
				   // renumber the degrees of freedom,
				   // you must do so immediately after
				   // distributing them, since such
				   // things as hanging nodes, the
				   // sparsity pattern etc. depend on
				   // the absolute numbers which are
				   // altered by renumbering.
				   //
				   // Renumbering does not serve any
				   // specific purpose in this
				   // example, it is done only for
				   // exposition of the technique. To
				   // see the effect of renumbering on
				   // the sparsity pattern of the
				   // matrix, refer to the second
				   // example program.

				   // The rest of the function is
				   // almost identically taken over
				   // from previous examples:
  hanging_node_constraints.clear ();
  DoFTools::make_hanging_node_constraints (dof_handler,
					   hanging_node_constraints);
  hanging_node_constraints.close ();

  sparsity_pattern.reinit (dof_handler.n_dofs(),
			   dof_handler.n_dofs(),
			   dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
  hanging_node_constraints.condense (sparsity_pattern);
  sparsity_pattern.compress();

  system_matrix.reinit (sparsity_pattern);

  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());
};



				 // Assembling the system of equations
				 // for the problem at hand is mostly
				 // as for the example programs
				 // before. However, some things have
				 // changed anyway, so we comment on
				 // this function fairly extensively.
template <int dim>
void LaplaceProblem<dim>::assemble_system () 
{  
				   // First we need to define objects
				   // which will be used as quadrature
				   // formula for domain and face
				   // integrals.
				   //
				   // Note the way in which we define
				   // a quadrature rule for the faces:
				   // it is simply a quadrature rule
				   // for one dimension less!
  QGauss3<dim>   quadrature_formula;
  QGauss3<dim-1> face_quadrature_formula;
				   // For simpler use later on, we
				   // alias the number of quadrature
				   // points to local variables:
  const unsigned int n_q_points    = quadrature_formula.n_quadrature_points;
  const unsigned int n_face_q_points = face_quadrature_formula.n_quadrature_points;
  
				   // Then we need objects which can
				   // evaluate the values, gradients,
				   // etc of the shape functions at
				   // the quadrature points. While it
				   // seems that it should be feasible
				   // to do it with one object for
				   // both domain and face integrals,
				   // there is a subtle difference
				   // since the weights in the domain
				   // integrals include the measure of
				   // the cell in the domain, while
				   // the face integral quadrature
				   // requires the measure of the face
				   // in a lower-dimensional
				   // manifold. Internally these two
				   // classes are rooted on a common
				   // base class which does most of
				   // the work; that, however, is
				   // something that you need not
				   // worry about.
				   //
				   // For the domain integrals in the
				   // bilinear form for Helmholtz's
				   // equation, we need to compute the
				   // values and gradients, as well as
				   // the weights at the quadrature
				   // points. Furthermore, we need the
				   // quadrature points on the real
				   // cell (rather than on the unit
				   // cell) to evaluate the right hand
				   // side function.
  FEValues<dim>  fe_values (*fe, quadrature_formula, 
			    UpdateFlags(update_values    |
					update_gradients |
					update_q_points  |
					update_JxW_values));

				   // For the face integrals, we only
				   // need the values of the shape
				   // functions, as well as the
				   // weights. We also need the normal
				   // vectors and quadrature points on
				   // the real cell since we want to
				   // determine the Neumann values
				   // from the exact solution object
				   // (see below).
  FEFaceValues<dim> fe_face_values (*fe, face_quadrature_formula, 
				    UpdateFlags(update_values         |
						update_q_points       |
						update_normal_vectors |
						update_JxW_values));

				   // In order to make programming
				   // more readable below, we alias
				   // the number of degrees of freedom
				   // per cell to a local variable, as
				   // already done for the number of
				   // quadrature points above:
  const unsigned int   dofs_per_cell = fe->dofs_per_cell;

				   // Then we need some objects
				   // already known from previous
				   // examples: An object denoting the
				   // right hand side function, its
				   // values at the quadrature points
				   // on a cell, the cell matrix and
				   // right hand side, and the indices
				   // of the degrees of freedom on a
				   // cell.
  RightHandSide<dim>   right_hand_side;
  std::vector<double>       rhs_values (n_q_points);

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

				   // Then we define an object
				   // denoting the exact solution
				   // function. We will use it to
				   // compute the Neumann values at
				   // the boundary from it. Usually,
				   // one would of course do so using
				   // a separate object, in particular
				   // since the exact solution is not
				   // known while the Neumann values
				   // are prescribed. We will,
				   // however, be a little bit lazy
				   // and use what we already have in
				   // information. Real-life programs
				   // would to go other ways here, of
				   // course.
  Solution<dim> exact_solution;

				   // Now for the main loop over all
				   // cells. This is mostly unchanged
				   // from previous examples, so we
				   // only comment on the things that
				   // have changed.
  typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
						 endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      cell_matrix.clear ();
      cell_rhs.clear ();

      fe_values.reinit (cell);
      const FullMatrix<double> 
	& shape_values = fe_values.get_shape_values();
      const std::vector<std::vector<Tensor<1,dim> > >
	& shape_grads  = fe_values.get_shape_grads();
      const std::vector<double>
	& JxW_values   = fe_values.get_JxW_values();
      const std::vector<Point<dim> >
	& q_points     = fe_values.get_quadrature_points();

      right_hand_side.value_list (q_points, rhs_values);
      
      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  {
	    for (unsigned int j=0; j<dofs_per_cell; ++j)
					       // The first thing that
					       // has changed is the
					       // bilinear form. It
					       // now contains the
					       // additional term from
					       // the Helmholtz
					       // equation, namely the
					       // scalar products of
					       // the two function
					       // values, rather than
					       // their gradients,
					       // which is the second
					       // term below:
	      cell_matrix(i,j) += ((shape_grads[i][q_point] *
				    shape_grads[j][q_point] *
				    JxW_values[q_point])
				   +
				   (shape_values(i,q_point) *
				    shape_values(j,q_point) *
				    JxW_values[q_point]));

	    cell_rhs(i) += (shape_values (i,q_point) *
			    rhs_values [q_point] *
			    fe_values.JxW (q_point));
	  };

				       // Then there is that second
				       // term on the right hand side,
				       // the contour integral. First
				       // we have to find out whether
				       // the intersection of the face
				       // of this cell with the
				       // boundary part Gamma2 is
				       // nonzero. To this end, we
				       // loop over all faces and
				       // check whether its boundary
				       // indicator equals ``1'',
				       // which is the value that we
				       // have assigned to that
				       // portions of the boundary
				       // composing Gamma2 in a
				       // function further below. The
				       // default value of boundary
				       // indicators is ``0'' for
				       // external faces, and ``255''
				       // for internal faces (the
				       // latter value should never be
				       // changed, and there is also
				       // no need to do so), so faces
				       // can only have an indicator
				       // equal to ``1'' if we have
				       // explicitly set it.
      for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
	if (cell->face(face)->boundary_indicator() == 1)
	  {
					     // If we came into here,
					     // then we have found an
					     // external face
					     // belonging to
					     // Gamma2. Next, we have
					     // to compute the values
					     // of the shape functions
					     // and the other
					     // quantities which we
					     // will need for the
					     // computation of the
					     // contour integral. This
					     // is done using the
					     // ``reinit'' function
					     // which we already know
					     // from the ``FEValue''
					     // class:
	    fe_face_values.reinit (cell, face);

					     // Then, for simpler
					     // access, we alias the
					     // various quantities to
					     // local variables:
	    const FullMatrix<double> 
	      & face_shape_values   = fe_face_values.get_shape_values();
	    const std::vector<double>
	      & face_JxW_values     = fe_face_values.get_JxW_values();
	    const std::vector<Point<dim> >
	      & face_q_points       = fe_face_values.get_quadrature_points();
	    const std::vector<Point<dim> >
	      & face_normal_vectors = fe_face_values.get_normal_vectors ();

					     // And we can then
					     // perform the
					     // integration by using a
					     // loop over all
					     // quadrature points.
	    for (unsigned int q_point=0; q_point<n_face_q_points; ++q_point)
	      {
						 // On each quadrature
						 // point, we first
						 // compute the value
						 // of the normal
						 // derivative. We do
						 // so using the
						 // gradient of the
						 // exact solution and
						 // the normal vector
						 // to the face at the
						 // present quadrature
						 // point:
		const double neumann_value
		  = (exact_solution.gradient (face_q_points[q_point]) *
		     face_normal_vectors[q_point]);

						 // Using this, we can
						 // compute the
						 // contribution of
						 // this face for each
						 // shape function:
		for (unsigned int i=0; i<dofs_per_cell; ++i)
		  cell_rhs(i) += (neumann_value *
				  face_shape_values(i,q_point) *
				  face_JxW_values[q_point]);
	      };
	  };

				       // Now that we have the
				       // contributions of the present
				       // cell, we can transfer it to
				       // the global matrix and right
				       // hand side vector, as in the
				       // examples before.
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

				   // The rest of the function has
				   // also been shown previously:
  hanging_node_constraints.condense (system_matrix);
  hanging_node_constraints.condense (system_rhs);

				   // Only with the interpolation of
				   // boundary values, there is one
				   // notable thing, namely that now
				   // the boundary indicator for which
				   // we interpolate boundary values
				   // (denoted by the second parameter
				   // to
				   // ``interpolate_boundary_values'')
				   // does not represent the whole
				   // boundary an more. Rather, it is
				   // that portion of the boundary
				   // which we have not assigned
				   // another indicator (see
				   // below). The degrees of freedom
				   // at the boundary that do not
				   // belong to Gamma1 are therefore
				   // excluded from the interpolation
				   // of boundary values.
  std::map<unsigned int,double> boundary_values;
  VectorTools::interpolate_boundary_values (dof_handler,
					    0,
					    Solution<dim>(),
					    boundary_values);
  MatrixTools<dim>::apply_boundary_values (boundary_values,
					   system_matrix,
					   solution,
					   system_rhs);
};


				 // Solving the system of equations is
				 // done in the same way as before.
template <int dim>
void LaplaceProblem<dim>::solve () 
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


				 // Now for the function doing grid
				 // refinement. Depending on the
				 // refinement mode passed to the
				 // constructor, we do global or
				 // adaptive refinement.
template <int dim>
void LaplaceProblem<dim>::refine_grid ()
{
  switch (refinement_mode) 
    {
				       // If global refinement is
				       // required, this is simple:
      case global_refinement:
      {
	triangulation.refine_global (1);
	break;
      };

					// In case of adaptive
					// refinement, we use the same
					// functions and classes as in
					// the previous example
					// program. Note that one
					// could treat Neumann
					// boundaries differently than
					// Dirichlet boundaries, and
					// one should in fact do so
					// here since we have Neumann
					// boundary conditions on part
					// of the boundaries, but
					// since we don't have a
					// function here that
					// describes the Neumann
					// values (we only construct
					// these values from the exact
					// solution when assembling
					// the matrix), we omit this
					// detail here.
      case adaptive_refinement:
      {
	Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

	FunctionMap<dim>::type neumann_boundary;
	KellyErrorEstimator<dim>::estimate (dof_handler,
					    QGauss3<dim-1>(),
					    neumann_boundary,
					    solution,
					    estimated_error_per_cell);

	GridRefinement::refine_and_coarsen_fixed_number (triangulation,
							 estimated_error_per_cell,
							 0.3, 0.03);
	
	triangulation.execute_coarsening_and_refinement ();

	break;
      };
    };
};



				 // Finally process the solution after
				 // it has been computed. For this, we
				 // integrate the error in various
				 // norms, and we generate tables that
				 // will be later used to display the
				 // convergence against the continuous
				 // solution in a nice format.
template <int dim>
void LaplaceProblem<dim>::process_solution (const unsigned int cycle)
{
				   // In order to integrate the
				   // difference between computed
				   // numerical solution and the
				   // continuous solution (described
				   // by the ``Solution'' class
				   // defined at the top of this
				   // file), we first need a vector
				   // that will hold the norm of the
				   // error on each cell. Since
				   // accuracy with 16 digits is not
				   // so important for these
				   // quantities, we save some memory
				   // by using ``float'' instead of
				   // ``double'' values.
  Vector<float> difference_per_cell (triangulation.n_active_cells());

				   // Next we use a function from the
				   // library which computes the error
				   // in the L2 norm on each cell. We
				   // have to pass it the DoF handler
				   // object, the vector holding the
				   // nodal values of the numerical
				   // solution, the continuous
				   // solution as a function object,
				   // the vector into which it shall
				   // place the norm of the error on
				   // each cell, a quadrature rule by
				   // which this norm shall be
				   // computed, and the type of norm
				   // to be used. Here, we use a Gauss
				   // formula with three points in
				   // each space direction, and
				   // compute the L2 norm.
  VectorTools::integrate_difference (dof_handler,
				     solution,
				     Solution<dim>(),
				     difference_per_cell,
				     QGauss3<dim>(),
				     L2_norm);
				   // Finally, we want to get the
				   // global L2 norm. This can of
				   // course be obtained by summing
				   // the squares of the norms on each
				   // cell, and taking the square root
				   // of that value. This is
				   // equivalent to taking the l2
				   // (lower case ``l'') norm of the
				   // vector of norms on each cell:
  const double L2_error = difference_per_cell.l2_norm();

				   // The same procedure is done to
				   // get the H1 semi-norm:
  VectorTools::integrate_difference (dof_handler,
				     solution,
				     Solution<dim>(),
				     difference_per_cell,
				     QGauss3<dim>(),
				     H1_seminorm);
  const double H1_error = difference_per_cell.l2_norm();

				   // Finally, we compute the maximum
				   // norm. Of course, we can't
				   // actually use the true maximum,
				   // but only the maximum at the
				   // quadrature points. Since this
				   // quite sensitively depends on the
				   // quadrature rule being used, and
				   // since we would like to avoid
				   // false results due to
				   // super-convergence effects at
				   // some points, we use a special
				   // quadrature rule that is obtained
				   // by iterating the trapezoidal
				   // rule five times in each space
				   // direction. Note that the
				   // constructor of the ``QIterated''
				   // class takes a one-dimensional
				   // quadrature rule and a number
				   // that tells it how often it shall
				   // use this rule in each space
				   // direction.
  QTrapez<1>     q_trapez;
  QIterated<dim> q_iterated (q_trapez, 5);

				   // Using this special quadrature
				   // rule, we can now try to find the
				   // maximal error on each cell:
  VectorTools::integrate_difference (dof_handler,
				     solution,
				     Solution<dim>(),
				     difference_per_cell,
				     q_iterated,
				     Linfty_norm);
				   // Obviously, the maximal error
				   // globally is the maximum over the
				   // maximal errors on each cell:
  const double Linfty_error = difference_per_cell.linfty_norm();

				   // After all these errors have been
				   // computed, we finally write some
				   // output and put all the data into
				   // a table.
  const unsigned int n_active_cells=triangulation.n_active_cells();
  const unsigned int n_dofs=dof_handler.n_dofs();
  
  std::cout << "Cycle " << cycle << ':' 
	    << std::endl
	    << "   Number of active cells:       "
	    << n_active_cells
	    << std::endl
	    << "   Number of degrees of freedom: "
	    << n_dofs
	    << std::endl;

				   // Add the important data to the
				   // ``TableHandler'' by giving the key
				   // of the column and the value.
				   // You don't need to define the keys
				   // beforehand, just add the values,
				   // and the column will be
				   // introduced into the table in the
				   // order the values are added the
				   // first time.
  convergence_table.add_value("cycle", cycle);
  convergence_table.add_value("cells", n_active_cells);
  convergence_table.add_value("dofs", n_dofs);
  convergence_table.add_value("L2", L2_error);
  convergence_table.add_value("H1", H1_error);
  convergence_table.add_value("Linfty", Linfty_error);
				   // You may set the precision with
				   // which the values will be written
				   // upon output.
  convergence_table.set_precision("L2", 3);
  convergence_table.set_precision("H1", 3);
  convergence_table.set_precision("Linfty", 3);
				   // The default notation is fixed
				   // point. For the columns you'd
				   // like to see in scientific notation
				   // set the `scientific_flag' `true'
				   // by the following lines:
  convergence_table.set_scientific("L2", true);
  convergence_table.set_scientific("H1", true);
  convergence_table.set_scientific("Linfty", true);
				   // For the output of a table into a
				   // LaTeX file, the default captions of
				   // the columns are the keys given
				   // as argument to the ``add_value''
				   // functions.  If you'd like to
				   // have TeX captions that differ
				   // from the default ones you can
				   // specify them by the following.
  convergence_table.set_tex_caption("cells", "\\# cells");
  convergence_table.set_tex_caption("dofs", "\\# dofs");
  convergence_table.set_tex_caption("L2", "$L^2$-error");
  convergence_table.set_tex_caption("H1", "$H^1$-error");
  convergence_table.set_tex_caption("Linfty", "$L^\\infty$-error");
				   // Note, that `\\' is reduced to
				   // `\' by the compiler such that the
				   // real TeX caption is e.g.
				   // `$L^\infty$-error'.
				   //
				   // The default TeX format of each
				   // column of the table is `c'
				   // (centered). To specify a
				   // different (e.g. `right') one,
				   // the following function may be
				   // used:
  convergence_table.set_tex_format("cells", "r");
  convergence_table.set_tex_format("dofs", "r");
};



				 // The following function is the main
				 // one which controls the flow of
				 // execution. The basic layout is as
				 // in previous examples: an outer
				 // loop over successively refined
				 // grids, and in this loop first
				 // problem setup, assemblage of the
				 // linear system, solution, and
				 // post-processing.
template <int dim>
void LaplaceProblem<dim>::run () 
{
  for (unsigned int cycle=0; cycle<7; ++cycle)
    {
				       // The first action in each
				       // iteration of the outer loop
				       // is setting up the grid on
				       // which we will solve in this
				       // iteration. In the first
				       // iteration, the coarsest grid
				       // is generated, in later
				       // iterations it is refined,
				       // for which we call the
				       // ``refine_grid'' function.
      if (cycle == 0)
	{
					   // Setting up the coarse
					   // grid is done as in
					   // previous examples: we
					   // first create an initial
					   // grid, which is the unit
					   // square [-1,1]x[-1,1] in
					   // the present case. Then
					   // we refine it globally a
					   // specific number of
					   // times.
	  GridGenerator::hyper_cube (triangulation, -1, 1);
	  triangulation.refine_global (1);

					   // However, here we have to
					   // do something else in
					   // addition: mark those
					   // faces that belong to the
					   // different components of
					   // the boundary, Gamma1 and
					   // Gamma2. We will use the
					   // following convention:
					   // Faces belonging to
					   // Gamma1 will have the
					   // boundary indicator ``0''
					   // (which is the default,
					   // so we don't have to set
					   // it explicitely), and
					   // faces belonging to
					   // Gamma2 will use ``1'' as
					   // boundary indicator.
					   //
					   // To set these values, we
					   // loop over all cells,
					   // then over all faces of a
					   // given cell, check
					   // whether it belongs to
					   // the boundary Gamma2, and
					   // if so set its boundary
					   // indicator to ``1''.
					   //
					   // It is worth noting that
					   // we have to loop over all
					   // cells here, not only the
					   // active ones. The reason
					   // is that upon refinement,
					   // newly created faces
					   // inherit the boundary
					   // indicator of their
					   // parent face. If we now
					   // only set the boundary
					   // indicator for active
					   // faces, coarsen some
					   // cells and refine them
					   // later on, they will
					   // again have the boundary
					   // indicator of the parent
					   // cell which we have not
					   // modified, instead of the
					   // one we
					   // intended. Therefore, we
					   // have to change the
					   // boundary indicators of
					   // all faces on Gamma2,
					   // irrespective whether
					   // they are active or not.
	  typename Triangulation<dim>::cell_iterator cell = triangulation.begin (),
						     endc = triangulation.end();
	  for (; cell!=endc; ++cell)
	    for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
	      if ((cell->face(face)->center()(0) == -1)
		  ||
		  (cell->face(face)->center()(1) == -1))
		cell->face(face)->set_boundary_indicator (1);
	}
      else
	{
					   // If this is not the first
					   // step, the we call
					   // ``refine_grid'' to
					   // actually refine the grid
					   // according to the
					   // refinement mode passed to
					   // the constructor.
	  refine_grid ();
	};
      

				       // The next steps you already
				       // know from previous
				       // examples. This is mostly the
				       // basic set-up of every finite
				       // element program:
      setup_system ();
      
      assemble_system ();
      solve ();

				       // The last step in this chain
				       // of function calls is usually
				       // evaluation of the computed
				       // solution for the quantities
				       // one is interested in. This
				       // is done in the following
				       // function. We pass the number
				       // of the loop iteration since
				       // that might be of interest to
				       // see in the logs which this
				       // function produces.
      process_solution (cycle);
    };
  
				   // After the last iteration we
				   // output the solution on the
				   // finest grid. This is done using
				   // the following sequence of
				   // statements which you have
				   // already seen in previous
				   // examples:
  std::string filename;
  switch (refinement_mode)
    {
      case global_refinement:
	    filename = "solution-global";
	    break;
      case adaptive_refinement:
	    filename = "solution-adaptive";
	    break;
      default:
	    Assert (false, ExcInternalError());
    };

				   // We augment the filename by a
				   // postfix denoting the finite
				   // element which we have used in
				   // the computation. Finding out
				   // which finite element we are
				   // actually using is not that
				   // simple here, since we only have
				   // a pointer to the common base
				   // class of all finite elements,
				   // which does not know anything
				   // about polynomial
				   // degrees. However, we actually
				   // know that we have generated a
				   // finite element of class
				   // ``FE_Q'', so we can use some C++
				   // feature to actually get a
				   // reference to the ``FE_Q''
				   // element pointed to by the
				   // reference and ask it for the
				   // polynomial degree. Note that if
				   // for whatever reason the object
				   // referenced behind the pointer to
				   // the base class should not be of
				   // type ``FE_Q'', then the C++
				   // language lets the
				   // ``dynamic_cast'' operator
				   // applied to a reference type
				   // throw an exception (if it were a
				   // pointer type, then a null
				   // pointer would be returned, which
				   // would then yield a segmentation
				   // fault when dereferenced in the
				   // subsequent call to
				   // ``get_order'').
  switch (dynamic_cast<const FE_Q<dim>&>(*fe).get_degree())
    {
      case 1:
	    filename += "-q1";
	    break;
      case 2:
	    filename += "-q2";
	    break;

      default:
				       // The finite element is
				       // neither Q1 nor Q2. This
				       // should not have happened,
				       // but maybe someone has tried
				       // to change this in ``main'',
				       // so it might happen. We catch
				       // this case and throw an
				       // exception, since we don't
				       // know how to name the
				       // respective output file
	    Assert (false, ExcInternalError());
    };
  
    
  filename += ".gmv";
	    
  std::ofstream output (filename.c_str());


  DataOut<dim> data_out;
  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution, "solution");

				   // Now building the intermediate
				   // format as before is the next
				   // step. We introduce one more
				   // feature of deal.II here. The
				   // background is the following: in
				   // some of the runs of this
				   // function, we have used
				   // biquadratic finite
				   // elements. However, since almost
				   // all output formats only support
				   // bilinear data, the data is
				   // written only bilinear, and
				   // information is lost
				   // therefore. Of course, we can't
				   // change the format in which
				   // graphic programs accept their
				   // inputs, but we can write the
				   // data differently such that we
				   // more closely resemble the
				   // information available in the
				   // quadratic approximation. We can,
				   // for example, write each cell as
				   // four sub-cells with bilinear data
				   // each, such that we have nine
				   // data points for each cell in the
				   // triangulation. The graphic
				   // programs will, of course,
				   // display this data still only
				   // bilinear, but at least we have
				   // given some more of the
				   // information we have.
				   //
				   // In order to allow writing more
				   // than one sub-cell per actual
				   // cell, the ``build_patches''
				   // function accepts a parameter
				   // (the default is ``1'', which is
				   // why you haven't seen this
				   // parameter in previous
				   // examples). This parameter
				   // denotes into how many sub-cells
				   // per space direction each cell
				   // shall be subdivided for
				   // output. For example, if you give
				   // ``2'', this leads to 4 cells in
				   // 2D and 8 cells in 3D. For
				   // quadratic elements, two
				   // sub-cells per space direction is
				   // obviously the right choice, so
				   // this is what we choose. In
				   // general, for elements of
				   // polynomial order ``q'', we use
				   // ``q'' subdivisions, and the
				   // order of the elements is
				   // determined in the same way as
				   // above in the ``run'' function:
  const unsigned int
    n_subcells = dynamic_cast<const FE_Q<dim>&>(*fe).get_degree();
  data_out.build_patches (n_subcells);

				   // Finally write out the data in
				   // GMV format.
  data_out.write_gmv (output);

				   // In each cycle values were added
				   // to the TableHandler. Now write
				   // the table to the standard output
				   // stream `std::cout'. Note, that the
				   // output in text format is a quite
				   // simple one and the captions may
				   // not be printed directly above
				   // the specific columns.
  convergence_table.write_text(std::cout);
				   // The table can also be written
				   // into a Tex file.  The (nicely)
				   // formatted table can be viewed at
				   // after calling `latex filename'
				   // and e.g. `xdvi filename', where
				   // filename is the name of the file
				   // which we construct from the name
				   // of the finite element and the
				   // refinement mode, as above
  if (true)
    {
      std::string filename = "error";
      switch (refinement_mode)
	{
	  case global_refinement:
		filename += "-global";
		break;
	  case adaptive_refinement:
		filename += "-adaptive";
		break;
	  default:
		Assert (false, ExcInternalError());
	};

      switch (dynamic_cast<const FE_Q<dim>&>(*fe).get_degree())
	{
	  case 1:
		filename += "-q1";
		break;
	  case 2:
		filename += "-q2";
		break;
	  default:
		Assert (false, ExcInternalError());
	};
  
      filename += ".tex";
      
      std::ofstream table_file(filename.c_str());
      convergence_table.write_tex(table_file);
      table_file.close();
    }
				   // In case you want the same
				   // caption for several columns, you
				   // can merge some columns to a
				   // super column by
  convergence_table.add_column_to_supercolumn("cycle", "n cells");
  convergence_table.add_column_to_supercolumn("cells", "n cells");
				   // You don't always need to output
				   // all columns. Also you don't need
				   // to restrict the order of the
				   // columns in the table to the
				   // order the columns were
				   // originally added during the run.
				   // Select and re-order the columns
				   // by adding the columns or the
				   // supercolumns to a new string
				   // vector.
  std::vector<std::string> new_order;
  new_order.push_back("n cells");
  new_order.push_back("H1");
  new_order.push_back("L2");
				   // and call
  convergence_table.set_column_order (new_order);

				   // In case of global refinement, it
				   // might be of interest to also
				   // output the convergence
				   // rates. This may be done by the
				   // functionality the
				   // ``ConvergenceTable'' offers over
				   // the regular
				   // ``TableHandler''. However, we do
				   // it only for global refinement,
				   // since for adaptive refinement
				   // the determination of something
				   // like an order of convergence is
				   // somewhat more involved.
  if (refinement_mode==global_refinement)
    {
				       // For everything that happened to
				       // the `ConvergenceTable' until
				       // this point, it would have been
				       // sufficient to use a simple
				       // `TableHandler'. Indeed, the
				       // `ConvergenceTable' is derived
				       // from the `TableHandler' but it
				       // offers the additional
				       // functionality of automatically
				       // evaluating convergence rates
      convergence_table.evaluate_convergence_rates(
	"L2", ConvergenceTable::reduction_rate);
				       // and/or the order of convergence.
      convergence_table.evaluate_convergence_rates(
	"L2", ConvergenceTable::reduction_rate_log2);
      convergence_table.evaluate_convergence_rates(
	"H1", ConvergenceTable::reduction_rate_log2);
				       // Each of the last three
				       // function calls produces an
				       // additional column that is
				       // merged with the original
				       // column (in our example the
				       // `L2' and the `H1' column) to
				       // a supercolumn.
    }

				   // Finally, the convergence chart
				   // is written. The filename is
				   // again constructed as above.
  convergence_table.write_text(std::cout);

  if (true)
    {
      std::string filename = "convergence";
      switch (refinement_mode)
	{
	  case global_refinement:
		filename += "-global";
		break;
	  case adaptive_refinement:
		filename += "-adaptive";
		break;
	  default:
		Assert (false, ExcInternalError());
	};
      switch (dynamic_cast<const FE_Q<dim>&>(*fe).get_degree())
	{
	  case 1:
		filename += "-q1";
		break;
	  case 2:
		filename += "-q2";
		break;
	  default:
		Assert (false, ExcInternalError());
	};
      filename += ".tex";

      std::ofstream table_file(filename.c_str());
      convergence_table.write_tex(table_file);
      table_file.close();
    }
};


				 // The main function is mostly as
				 // before. The only difference is
				 // that we solve three times, once
				 // for Q1 and adaptive refinement,
				 // once for Q1 elements and global
				 // refinement, and once for Q2
				 // elements and global refinement.
int main () 
{
  try
    {
      deallog.depth_console (0);

				       // Now for the three calls to
				       // the main class. Each call is
				       // blocked into curly braces in
				       // order to destroy the
				       // respective objects (i.e. the
				       // finite element and the
				       // LaplaceProblem object) at
				       // the end of the block and
				       // before we go to the next
				       // run.
      {
	std::cout << "Solving with Q1 elements, adaptive refinement" << std::endl
		  << "=============================================" << std::endl
		  << std::endl;
	
	FE_Q<2> fe(1);
	LaplaceProblem<2> laplace_problem_2d (fe, LaplaceProblem<2>::adaptive_refinement);
	laplace_problem_2d.run ();

	std::cout << std::endl;
      };
	
      {
	std::cout << "Solving with Q1 elements, global refinement" << std::endl
		  << "===========================================" << std::endl
		  << std::endl;
	
	FE_Q<2> fe(1);
	LaplaceProblem<2> laplace_problem_2d (fe, LaplaceProblem<2>::global_refinement);
	laplace_problem_2d.run ();

	std::cout << std::endl;
      };
       
      {
	std::cout << "Solving with Q2 elements, global refinement" << std::endl
		  << "===========================================" << std::endl
		  << std::endl;
	
	FE_Q<2> fe(2);
	LaplaceProblem<2> laplace_problem_2d (fe, LaplaceProblem<2>::global_refinement);
	laplace_problem_2d.run ();

	std::cout << std::endl;
      };
       
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
