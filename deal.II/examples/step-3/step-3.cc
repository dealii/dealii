/* $Id$ */

				 // These include files are already
				 // known to you. They declare the
				 // classes which handle
				 // triangulations and enumerate the
				 // degrees of freedom.
#include <grid/tria.h>
#include <dofs/dof_handler.h>
				 // And this is the file in which the
				 // functions are declared which
				 // create grids.
#include <grid/grid_generator.h>

				 // The next three files contain
				 // classes which are needed for loops
				 // over all cells and to get the
				 // information from the cell objects.
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>

				 // In this file are the finite
				 // element descriptions.
#include <fe/fe_lib.lagrange.h>

				 // And this file is needed for the
				 // creation of sparsity patterns of
				 // sparse matrices, as shown in
				 // previous examples:
#include <dofs/dof_tools.h>

				 // The next two file are needed for
				 // assembling the matrix using
				 // quadrature on each cell. The
				 // classes declared in them will be
				 // explained below.
#include <fe/fe_values.h>
#include <base/quadrature_lib.h>

				 // The following three include files
				 // we need for the treatment of
				 // boundary values:
#include <base/function.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>

				 // These include files are for the
				 // linear algebra which we employ to
				 // solve the system of equations
				 // arising from the finite element
				 // discretization of the Laplace
				 // equation. We will use vectors and
				 // full matrices for assembling the
				 // system of equations locally on
				 // each cell, and transfer the
				 // results into a sparse matrix. We
				 // will then use a Conjugate Gradient
				 // solver to solve the problem, for
				 // which we need a preconditioner (in
				 // this program, we use the identity
				 // preconditioner which does nothing,
				 // but we need to include the file
				 // anyway), and a class which
				 // provides the solver with some
				 // memory for temporary vectors.
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/vector_memory.h>
#include <lac/precondition.h>

				 // Finally, this is for output to a
				 // file.
#include <numerics/data_out.h>
#include <fstream>


				 // Instead of the procedural
				 // programming of previous examples,
				 // we encapsulate everything into a
				 // class for this program. The class
				 // consists of functions which do
				 // certain aspects of a finite
				 // element program, a `main' function
				 // which controls what is done first
				 // and what is done next, and a list
				 // of member variables.
class LaplaceProblem 
{
  public:
				     // This is the constructor:
    LaplaceProblem ();

				     // And the top-level function,
				     // which is called from the
				     // outside to start the whole
				     // program (see the `main'
				     // function at the bottom of this
				     // file):
    void run ();
    
  private:
				     // Then there are some member
				     // functions that mostly do what
				     // their names suggest. Since
				     // they do not need to be called
				     // from outside, they are made
				     // private to this class.
  private:
    void make_grid_and_dofs ();
    void assemble_system ();
    void solve ();
    void output_results ();

				     // And then we have the member
				     // variables. There are variables
				     // describing the triangulation
				     // and the numbering of the
				     // degrees of freedom...
    Triangulation<2>     triangulation;
    FEQ1<2>              fe;
    DoFHandler<2>        dof_handler;

				     // ...variables for the sparsity
				     // pattern and values of the
				     // system matrix resulting from
				     // the discretization of the
				     // Laplace equation...
    SparseMatrixStruct   sparsity_pattern;
    SparseMatrix<double> system_matrix;

				     // ...and variables which will
				     // hold the right hand side and
				     // solution vectors.
    Vector<double>       solution;
    Vector<double>       system_rhs;
};


				 // Here comes the constructor. It
				 // does not much more than associate
				 // the dof_handler variable to the
				 // triangulation we use. All the
				 // other member variables of the
				 // LaplaceProblem class have a
				 // default constructor which does all
				 // we want.
LaplaceProblem::LaplaceProblem () :
		dof_handler (triangulation)
{};


				 // Now, the first thing we've got to
				 // do is to generate the
				 // triangulation on which we would
				 // like to do our computation and
				 // number each vertex with a degree
				 // of freedom. We have seen this in
				 // the previous examples before. Then
				 // we have to set up space for the
				 // system matrix and right hand side
				 // of the discretized problem. This
				 // is what this function does:
void LaplaceProblem::make_grid_and_dofs ()
{
				   // First create the grid and refine
				   // all cells five times. Since the
				   // initial grid (which is the
				   // square [-1,1]x[-1,1]) consists
				   // of only one cell, the final grid
				   // has 32 times 32 cells, for a
				   // total of 1024.
  GridGenerator::hyper_cube (triangulation, -1, 1);
  triangulation.refine_global (5);
				   // Unsure that 1024 is the correct
				   // number? Let's see:
				   // n_active_cells return the number
				   // of terminal cells. By terminal
				   // we mean the cells on the finest
				   // grid.
  cout << "Number of active cells: "
       << triangulation.n_active_cells()
       << endl;
				   // We stress the adjective
				   // `terminal' or `active', since
				   // there are more cells, namely the
				   // parent cells of the finest
				   // cells, their parents, etc, up to
				   // the one cell which made up the
				   // initial grid. Of course, on the
				   // next coarser level, the number
				   // of cells is one quarter that of
				   // the cells on the finest level,
				   // i.e. 256, then 64, 16, 4, and
				   // 1. We can get the total number
				   // of cells like this:
  cout << "Total number of cells: "
       << triangulation.n_cells()
       << endl;
				   // Note the distinction between
				   // n_active_cells() and n_cells().
  
				   // Next we enumerate all the
				   // degrees of freedom. This is done
				   // by using the distribute_dofs
				   // function, as we have seen in
				   // previous examples. Since we use
				   // the FEQ1 class, i.e. bilinear
				   // elements, this associates one
				   // degree of freedom with each
				   // vertex.
  dof_handler.distribute_dofs (fe);

				   // Now that we have the degrees of
				   // freedom, we can take a look at
				   // how many there are:
  cout << "Number of degrees of freedom: "
       << dof_handler.n_dofs()
       << endl;
				   // There should be one DoF for each
				   // vertex. Since we have a 32 times
				   // 32 grid, the number of DoFs
				   // should be 33 times 33, or 1089.

				   // As we have seen in the previous
				   // example, we set up a sparse
				   // matrix for the system matrix and
				   // tag those entries that might be
				   // nonzero. Since that has already
				   // been done, we won't discuss the
				   // next few lines:
  sparsity_pattern.reinit (dof_handler.n_dofs(),
			   dof_handler.n_dofs(),
			   dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
  sparsity_pattern.compress();

				   // Now the sparsity pattern is
				   // built and fixed (after
				   // `compress' has been called, you
				   // can't add nonzero entries
				   // anymore; the sparsity pattern is
				   // `sealed', so to say), and we can
				   // initialize the matrix itself
				   // with it. Note that the
				   // SparseMatrixStruct object does
				   // not hold the values of the
				   // matrix, it only stores the
				   // places where entries are. The
				   // entries are themselves stored in
				   // objects of type SparseMatrix, of
				   // which our variable system_matrix
				   // is one.
				   //
				   // The distinction between sparsity
				   // pattern and matrix was made to
				   // allow several matrices to use
				   // the same sparsity pattern. This
				   // may not seem relevant, but when
				   // you consider the size which
				   // matrices can have, and that it
				   // may take some time to build the
				   // sparsity pattern, this becomes
				   // important in large-scale
				   // problems.
  system_matrix.reinit (sparsity_pattern);

				   // The last thing to do in this
				   // function is to set the sizes of
				   // the right hand side vector and
				   // the solution vector to the right
				   // values:
  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());
};


				 // Now comes the difficult part:
				 // assembling matrices and
				 // vectors. In fact, this is not
				 // overly difficult, but it is
				 // something that the library can't
				 // do for you as for most of the
				 // other things in the functions
				 // above and below.
				 //
				 // The general way to assemble
				 // matrices and vectors is to loop
				 // over all cells, and on each cell
				 // compute the contribution of that
				 // cell to the global matrix and
				 // right hand side by quadrature. The
				 // idea now is that since we only
				 // need the finite element shape
				 // functions on the quadrature points
				 // of each cell, we don't need the
				 // shape functions of the finite
				 // element themselves any
				 // more. Therefore, we won't deal
				 // with the finite element object
				 // `fe' (which was of type FEQ1), but
				 // with another object which only
				 // provides us with the values,
				 // gradients, etc of the shape
				 // functions at the quadrature
				 // points. The objects which do this
				 // are of type FEValues.
void LaplaceProblem::assemble_system () 
{
				   // Ok, let's start: we need a
				   // quadrature formula for the
				   // evaluation of the integrals on
				   // each cell. Let's take a Gauss
				   // formula with three quadrature
				   // points in each direction, i.e. a
				   // total of nine points since we
				   // are in 2D:
  QGauss3<2>  quadrature_formula;
				   // And we initialize the object
				   // which we have briefly talked
				   // about above. It needs to be told
				   // which the finite element is that
				   // we want to use, the quadrature
				   // points and their
				   // weights. Finally, we have to
				   // tell it what we want it to
				   // compute on each cell: we need
				   // the values of the shape
				   // functions at the quadrature
				   // points, their gradients, and
				   // also the weights of the
				   // quadrature points and the
				   // determinants of the Jacobian
				   // transformations from the unit
				   // cell to the real cells. The
				   // values of the shape functions
				   // computed by specifying
				   // update_values; the gradients are
				   // done alike, using
				   // update_gradients. The
				   // determinants of the Jacobians
				   // and the weights are always used
				   // together, so only the products
				   // (Jacobians times weights, or
				   // short JxW) are computed; since
				   // we also need them, we have to
				   // list them as well:
  FEValues<2> fe_values (fe, quadrature_formula, 
			 UpdateFlags(update_values    |
				     update_gradients |
				     update_JxW_values));

				   // For use further down below, we
				   // define two short cuts for the
				   // number of degrees of freedom on
				   // each cell (since we are in 2D
				   // and degrees of freedom are
				   // associated with vertices only,
				   // this number is four). We also
				   // define an abbreviation for the
				   // number of quadrature points
				   // (here that should be nine). In
				   // general, it is a good idea to
				   // use their symbolic names instead
				   // of hard-coding these number even
				   // if you know them, since you may
				   // want to change the quadrature
				   // formula and/or finite element at
				   // some time; the program will just
				   // work with these changes, without
				   // the need to change the matrix
				   // assemblage.
				   //
				   // The shortcuts, finally, are only
				   // defined to make the following
				   // loops a bit more readable. You
				   // will see them in many places in
				   // larger programs, and
				   // `dofs_per_cell' and `n_q_points'
				   // are more or less standard names
				   // for these purposes.
  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q_points    = quadrature_formula.n_quadrature_points;

				   // Now, we said that we wanted to
				   // assemble the global matrix and
				   // vector cell-by-cell. We could
				   // write the results directly into
				   // the global matrix, but this is
				   // not very efficient since access
				   // to the elements of a sparse
				   // matrix is slow. Rather, we first
				   // compute the contribution of each
				   // ell in a small matrix with the
				   // degrees of freedom on the
				   // present cell, and only transfer
				   // them to the global matrix when
				   // the copmutations are finished
				   // for this cell. We do the same
				   // for the right hand side vector,
				   // although access times are not so
				   // problematic for them.
  FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs (dofs_per_cell);

				   // When assembling the
				   // contributions of each cell, we
				   // do this with the local numbering
				   // of the degrees of freedom
				   // (i.e. the number running from
				   // zero through
				   // dofs_per_cell-1). However, when
				   // we transfer the result into the
				   // global matrix, we have to know
				   // the global numbers of the
				   // degrees of freedom. When we get
				   // them, we need a scratch array
				   // for these numbers:
  vector<int>        local_dof_indices (dofs_per_cell);

				   // Now for th loop over all
				   // cells. You have seen before how
				   // this works, so this should be
				   // familiar to you:
  DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active(),
				      endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
				       // We are on one cell, and we
				       // would like the values and
				       // gradients of the shape
				       // functions be computed, as
				       // well as the determinants of
				       // the Jacobian matrices of the
				       // mapping between unit cell
				       // and true cell, at the
				       // quadrature points. Since all
				       // these values depend on the
				       // geometry of the cell, we
				       // have to have the FEValues
				       // object re-compute them on
				       // each cell:
      fe_values.reinit (cell);

				       // Reset the values of the
				       // contributions of this cell
				       // to global matrix and global
				       // right hand side to zero,
				       // before we fill them.
      cell_matrix.clear ();
      cell_rhs.clear ();

				       // Assemble the matrix: For the
				       // Laplace problem, the matrix
				       // on each cell is the integral
				       // over the gradients of shape
				       // function i and j. Since we
				       // do not integrate, but rather
				       // use quadrature, this is the
				       // sum over all quadrature
				       // points of the integrands
				       // times the determinant of the
				       // Jacobian matrix at the
				       // quadrature point times the
				       // weight of this quadrature
				       // point. You can get the
				       // gradient of shape function i
				       // at quadrature point q_point
				       // by using
				       // fe_values.shape_grad(i,q_point);
				       // this gradient is a
				       // 2-dimensional vector (in
				       // fact it is of type
				       // Tensor<1,dim>, with here
				       // dim=2) and the product of
				       // two such vectors is the
				       // scalar product, i.e. the
				       // product of the two
				       // shape_grad function calls is
				       // the dot product.
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=0; j<dofs_per_cell; ++j)
	  for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	    cell_matrix(i,j) += (fe_values.shape_grad (i, q_point) *
				 fe_values.shape_grad (j, q_point) *
				 fe_values.JxW (q_point));

				       // We then do the same thing
				       // for the right hand
				       // side. Here, the integral is
				       // over the shape function i
				       // times the right hand side
				       // function, which we choose to
				       // be the function with
				       // constant value one (more
				       // interesting examples will be
				       // considered in the following
				       // programs). Again, we compute
				       // the integral by quadrature,
				       // which transforms the
				       // integral to a sum over all
				       // quadrature points of the
				       // value of the shape function
				       // at that point times the
				       // right hand side function
				       // (i.e. 1) times the Jacobian
				       // determinant times the weight
				       // of that quadrature point:
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	  cell_rhs(i) += (fe_values.shape_value (i, q_point) *
			  1 *
			  fe_values.JxW (q_point));

				       // Now that we have the
				       // contribution of this cell,
				       // we have to transfer it to
				       // the global matrix and right
				       // hand side. To this end, we
				       // first have to find out which
				       // global numbers the degrees
				       // of freedom on this cell
				       // have. Let's simply ask the
				       // cell for that information:
      cell->get_dof_indices (local_dof_indices);

				       // Then again loop over all
				       // shape functions i and j and
				       // transfer the local elements
				       // to the global matrix. The
				       // global numbers can be
				       // obtained using
				       // local_dof_indices[i]:
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=0; j<dofs_per_cell; ++j)
	  system_matrix.add (local_dof_indices[i],
			     local_dof_indices[j],
			     cell_matrix(i,j));

				       // And again, we do the same
				       // thing for the right hand
				       // side vector.
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	system_rhs(local_dof_indices[i]) += cell_rhs(i);
    };


				   // Now almost everything is set up
				   // for the solution of the discrete
				   // system. However, we have not yet
				   // taken care of boundary values
				   // (in fact, Laplace's equation
				   // without Dirichlet boundary
				   // values is not even uniquely
				   // solvable, since you can add an
				   // arbitrary constant to the
				   // discrete solution). We therefore
				   // have to take into account
				   // boundary values.
				   //
				   // For this, we first obtain a list
				   // of the degrees of freedom on the
				   // boundary and the value the shape
				   // function shall have there. For
				   // simplicity, we only interpolate
				   // the boundary value function,
				   // rather than projecting them onto
				   // the boundary. There is a
				   // function in the library which
				   // does exactly this:
				   // interpolate_boundary_values. Its
				   // parameters are (omitting
				   // parameters for which default
				   // values exist which are
				   // sufficient here): the DoFHandler
				   // object to get the global numbers
				   // of the degrees of freedom on the
				   // boundary; the component of the
				   // boundary where the boundary
				   // values shall be interpolated;
				   // the boundary value function
				   // itself; and the output object.
				   //
				   // The component of the boundary is
				   // meant as follows: in many cases,
				   // you may want to impose certain
				   // boundary values only on parts of
				   // the boundary. For example, you
				   // may have inflow and outflow
				   // boundaries in fluid dynamics,
				   // are clamped and free parts of
				   // bodies in deformation
				   // computations of bodies. Then you
				   // will want to denote these
				   // different parts of the boundary
				   // by different numbers and tell
				   // the interpolate_boundary_values
				   // function to only compute the
				   // boundary values on a certain
				   // part of the boundary (e.g. the
				   // clamped part, or the inflow
				   // boundary). By default, all
				   // boundaries have the number `0',
				   // and since we have not changed
				   // that, this is still so;
				   // therefore, if we give `0' as the
				   // desired portion of the boundary,
				   // this means we get the whole
				   // boundary.
				   //
				   // The function describing the
				   // boundary values is an object of
				   // type `Function' or of a derived
				   // class. One of the derived
				   // classes is ZeroFunction, which
				   // described a function which is
				   // zero everywhere. We create such
				   // an object in-place and pass it
				   // to the
				   // interpolate_boundary_values
				   // function.
				   //
				   // Finally, the output object is a
				   // list of pairs of global degree
				   // of freedom numbers (i.e. the
				   // number of the degrees of freedom
				   // on the boundary) and their
				   // boundary values (which are zero
				   // here for all entries). This
				   // mapping of DoF numbers to
				   // boundary values is done by the
				   // `map' class.
  map<int,double> boundary_values;
  VectorTools::interpolate_boundary_values (dof_handler,
					    0,
					    ZeroFunction<2>(),
					    boundary_values);
				   // Now that we got the list of
				   // boundary DoFs and their
				   // respective boundary values,
				   // let's use them to modify the
				   // system of equations
				   // accordingly. This is done by the
				   // following function call:
  MatrixTools<2>::apply_boundary_values (boundary_values,
					 system_matrix,
					 solution,
					 system_rhs);
};


				 // The following function simply
				 // solves the discretized
				 // equation. As the system is quite a
				 // large one for direct solvers such
				 // as Gauss elimination or LU
				 // decomposition, we use a Conjugate
				 // Gradient algorithm. You should
				 // remember that the number of
				 // variables here (only 1089) is a
				 // very small number for finite
				 // element computations, where
				 // 100.000 is a more usual number;
				 // for this number of variables,
				 // direct methods are no longer
				 // usable and you are forced to use
				 // methods like CG.
void LaplaceProblem::solve () 
{
				   // We need to tell the algorithm
				   // where to stop. This is done by
				   // using a SolverControl object,
				   // and as stopping criterion we
				   // say: maximally 1000 iterations
				   // (which is far more than is
				   // needed for 1089 variables; see
				   // the results section to find out
				   // how many were really used), and
				   // stop if the norm of the residual
				   // is below 1e-12. In practice, the
				   // latter criterion will be the one
				   // which stops the iteration.
  SolverControl           solver_control (1000, 1e-12);
				   // Furthermore, the CG algorithm
				   // needs some space for temporary
				   // vectors. Rather than allocating
				   // it on the stack or heap itself,
				   // it relies on helper objects,
				   // which can sometimes do a better
				   // job at this. The
				   // PrimitiveVectorMemory class is
				   // such a helper class which the
				   // solver can ask for memory. The
				   // angle brackets indicate that
				   // this class really takes a
				   // template parameter (here the
				   // data type of the vectors we
				   // use), which however has a
				   // default value, which is
				   // appropriate here.
  PrimitiveVectorMemory<> vector_memory;
				   // Then we need the solver
				   // itself. The template parameters
				   // here are the matrix type and the
				   // type of the vectors. They
				   // default to the ones we use here.
  SolverCG<>              cg (solver_control, vector_memory);

				   // Now solve the system of
				   // equations. The CG solver takes a
				   // preconditioner, but we don't
				   // want to use one, so we tell it
				   // to use the identity operation as
				   // preconditioner.
  cg.solve (system_matrix, solution, system_rhs,
	    PreconditionIdentity());
				   // Now that the solver has done its
				   // job, the solution variable
				   // contains the nodal values of the
				   // solution function.
};


				 // The last part of a typical finite
				 // element program is to output the
				 // results and maybe do some
				 // postprocessing (for example
				 // compute the maximal stress values
				 // at the boundary, or the average
				 // flux across the outflow, etc). We
				 // have no such postprocessing here,
				 // but we would like to write the
				 // solution to a file.
void LaplaceProblem::output_results () 
{
				   // To write the output to a file,
				   // we need an object which knows
				   // about output formats and the
				   // like. This is the DataOut class,
				   // and we need an object of that
				   // type:
  DataOut<2> data_out;
				   // Now we have to tell it where to
				   // take the values from which it
				   // shall write. We tell it which
				   // DoFHandler object to use, and we
				   // add the solution vector (and the
				   // name by which it shall be
				   // written to disk) to the list of
				   // data that is to be written. If
				   // we had more than one vector
				   // which we would like to look at
				   // in the output (for example right
				   // hand sides, errors per cell,
				   // etc) we would add them as well:
  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution, "solution");
				   // After the DataOut object knows
				   // which data it is to work on, we
				   // have to tell it to process them
				   // into something the backends can
				   // handle. The reason is that we
				   // have separated the frontend
				   // (which knows about how to treat
				   // DoFHandler objects and data
				   // vectors) from the backend (which
				   // knows several output formats)
				   // and use an intermediate data
				   // format to transfer data from the
				   // front- to the backend. The data
				   // is transformed into this
				   // intermediate format by the
				   // following function:
  data_out.build_patches ();

				   // Now we have everything in place
				   // for the actual output. Just open
				   // a file and write the data into
				   // it, using GNUPLOT format (there
				   // are other functions which write
				   // their data in postscript, AVS,
				   // GMV, or some other format):
  ofstream output ("solution.gpl");
  data_out.write_gnuplot (output);
};


				 // The following function is the main
				 // function which calls all the other
				 // functions of the LaplaceProblem
				 // class. The order in which this is
				 // done resembles the order in which
				 // most finite element programs
				 // work. Since the names are mostly
				 // self-explanatory, there is not
				 // much to comment about:
void LaplaceProblem::run () 
{
  make_grid_and_dofs ();
  assemble_system ();
  solve ();
  output_results ();
};

    

				 // This is the main function of the
				 // program. Since the concept of a
				 // main function is mostly a remnant
				 // from the pre-object era in C/C++
				 // programming, it often does not
				 // much more than creating an object
				 // of the top-level class and calling
				 // it principle function. This is
				 // what is done here as well.
int main () 
{
  LaplaceProblem laplace_problem;
  laplace_problem.run ();
  return 0;
};
