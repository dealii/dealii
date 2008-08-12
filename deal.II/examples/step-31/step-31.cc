/* $Id$ */
/* Author: Wolfgang Bangerth, Texas A&M University, 2007 */

/*    $Id$       */
/*                                                                */
/*    Copyright (C) 2007, 2008 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */

				 // @sect3{Include files}

				 // We include the functionality
				 // of these well-known deal.II
				 // library files and some C++
				 // header files.
#include <base/quadrature_lib.h>
#include <base/logstream.h>
#include <base/function.h>
#include <base/utilities.h>

#include <lac/block_vector.h>
#include <lac/full_matrix.h>
#include <lac/block_sparse_matrix.h>
#include <lac/solver_gmres.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>
#include <lac/sparse_direct.h>
#include <lac/sparse_ilu.h>
#include <lac/block_matrix_array.h>

#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary_lib.h>
#include <grid/grid_tools.h>
#include <grid/grid_refinement.h>

#include <dofs/dof_handler.h>
#include <dofs/dof_renumbering.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <dofs/dof_constraints.h>

#include <fe/fe_q.h>
#include <fe/fe_dgq.h>
#include <fe/fe_system.h>
#include <fe/fe_values.h>
#include <fe/mapping_q1.h>

#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <numerics/data_out.h>
#include <numerics/derivative_approximation.h>
#include <numerics/solution_transfer.h>

#include <fstream>
#include <sstream>

                                 // This is Trilinos
#include <Epetra_SerialComm.h>
#include <Epetra_Map.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include <Amesos.h>
#include <ml_include.h>
#include <ml_MultiLevelPreconditioner.h>

				 // Next, we import all deal.II
				 // names into global namespace
using namespace dealii;


				 // @sect3{Defining the inner preconditioner type}

				 // This class creates a local typedef that
				 // specifies the preconditioner we're
				 // going to use in the code below, depending
				 // on the space dimension. This
				 // is in complete analogy to step-22.
template <int dim>
struct InnerPreconditioner;

template <>
struct InnerPreconditioner<2>
{
    typedef SparseDirectUMFPACK type;
};

template <>
struct InnerPreconditioner<3>
{
    typedef SparseILU<double> type;
};



				 // This implements an AMG
				 // preconditioner based on the
				 // Trilinos ML implementation.
				 // What this class does is twofold.
				 // When the constructor of the class
				 // is invoked, a ML preconditioner
				 // object is created based on the
				 // DoFHandler and matrix
				 // that we want the preconditioner to
				 // be based on. A call of
				 // the respective <code>vmult</code>
				 // function does call the respective
				 // operation in the Trilinos package,
				 // where it is called 
				 // <code>ApplyInverse</code>.
 
				 // There are a few pecularities in
				 // the constructor. Since the
				 // Trilinos objects we want to use are
				 // heavily dependent on Epetra objects,
				 // the fundamental construction
				 // routines for vectors and 
				 // matrices in Trilinos, we do a 
				 // copy of our deal.II preconditioner
				 // matrix to a Epetra matrix. This 
				 // is of course not optimal, but for
				 // the time being there is no direct
				 // support for our data interface.
				 // When doing this time-consuming 
				 // operation, we can still profit 
				 // from the fact that some of the
				 // entries in the preconditioner matrix
				 // are zero and hence can be 
				 // neglected.
template <int dim>
class TrilinosAmgPreconditioner : public Subscriptor
{
  public:
    TrilinosAmgPreconditioner (
	  const DoFHandler<dim>      &DofHandler,
	  const SparseMatrix<double> &PreconditionerMatrix,
	  const bool                 VectorValuedProblem);

  void vmult (Vector<double>       &dst,
              const Vector<double> &src) const;

  private:
    const DoFHandler<dim>      *dof_handler;
    const SparseMatrix<double> *preconditioner_matrix;
    const SparsityPattern      *sparsity_pattern;

    const unsigned int n_u;
    
    const bool vector_valued_problem;
    
    ML_Epetra::MultiLevelPreconditioner* ml_precond;
    
    Epetra_SerialComm  communicator;
    std::auto_ptr<Epetra_Map>       Map;
    std::auto_ptr<Epetra_CrsMatrix> Matrix;
};


template <int dim>
TrilinosAmgPreconditioner<dim>::TrilinosAmgPreconditioner(
		const DoFHandler<dim>       &DofHandler,
		const SparseMatrix<double>  &PreconditionerMatrix,
		const bool                  VectorValuedProblem
		)
		:
		dof_handler           (&DofHandler),
		preconditioner_matrix (&PreconditionerMatrix),
		sparsity_pattern      (&preconditioner_matrix->get_sparsity_pattern()),
		n_u                   (preconditioner_matrix->m()),
		vector_valued_problem (VectorValuedProblem)
{
  
				 // Init Epetra Matrix, avoid 
				 // storing the nonzero elements.
  {
    Map.reset (new Epetra_Map(n_u, 0, communicator));
    
    std::vector<int> row_lengths (n_u);
    for (unsigned int row=0; row<n_u; ++row)
      {
	const unsigned int temporary_row_length = 
	    sparsity_pattern->row_length (row);
	unsigned int local_length = 0;
	for (unsigned int col=0; col<temporary_row_length; ++col)
	  {
	    unsigned int col_index = sparsity_pattern->column_number (row, col);
	    if (std::abs((*preconditioner_matrix) (row, col_index)) > 1e-15)
	      local_length += 1;
	  }
	row_lengths[row] = local_length;
      }
  
    Matrix.reset (new Epetra_CrsMatrix(Copy, *Map, &row_lengths[0], true));
  
    const unsigned int max_nonzero_entries
      = *std::max_element (row_lengths.begin(), row_lengths.end());
  
    std::vector<double> values(max_nonzero_entries, 0);
    std::vector<int> row_indices(max_nonzero_entries);
  
    for (unsigned int row=0; row<n_u; ++row)
      {
	const unsigned int temporary_row_length = 
	    sparsity_pattern->row_length (row);
	
	row_indices.resize (row_lengths[row], 0);
	values.resize (row_lengths[row], 0.);
  
	int col_counter = 0;
	for (unsigned int col=0; col<temporary_row_length; ++col)
	  {
	    unsigned int col_index = sparsity_pattern->column_number (row, col);
	    if (std::abs((*preconditioner_matrix) (row, col_index)) > 1e-15)
	      {
		row_indices[col_counter] = 
		    sparsity_pattern->column_number (row, col);
		values[col_counter] = 
		    (*preconditioner_matrix) (row, row_indices[col_counter]);
		++col_counter;
	      }
	    Assert (col_counter == row_lengths[row],
		    ExcMessage("Filtering out zeros could not "
			       "be successfully finished!"));
	  }
  
	Matrix->InsertGlobalValues(row, row_lengths[row],
				   &values[0], &row_indices[0]);
      }
      
    Matrix->FillComplete();
  }
  
				 // And now build the AMG
				 // preconditioner.
  const bool output_amg_info = false;
  Teuchos::ParameterList MLList;
  
				 // set default values for 
				 // smoothed aggregation in MLList,
				 // the standard choice for elliptic
				 // (Laplace-type) problems
  ML_Epetra::SetDefaults("SA",MLList);
  
  if (output_amg_info)
    MLList.set("ML output", 10);
  else
    MLList.set("ML output", 0);
  
				 // modify some AMG parameters from default to
				 // get better performance on FE induced Laplace
				 // type matrices
  MLList.set("max levels",10);
  MLList.set("increasing or decreasing", "increasing");
  MLList.set("aggregation: type", "Uncoupled");
  MLList.set("smoother: type", "symmetric Gauss-Seidel");
  MLList.set("smoother: sweeps", 1);
  MLList.set("smoother: damping factor", 4./3.);
  MLList.set("smoother: pre or post", "both");
  MLList.set("coarse: type","Amesos-KLU");
  
				 // Build null space, i.e. build dim vectors
				 // of ones in each velocity component.
  if (vector_valued_problem)
    {
      std::vector<double> null_vectors (dim * n_u, 0.);
      {
	unsigned int n_ud = n_u/dim;
	Assert (n_ud * dim == n_u,
	    ExcMessage("Cannot find portions of single velocity components!"));
	
	std::vector<bool> velocity_d_dofs (dof_handler->n_dofs(), false);
	std::vector<bool> velocity_mask (dim + 2, false);
	for (unsigned int d=0; d<dim; d++)
	  {
	    velocity_mask[d] = true;
	    DoFTools::extract_dofs(*dof_handler, velocity_mask, 
				    velocity_d_dofs);
	    velocity_mask[d] = false;
	    
	    unsigned int counter = 0;
	    for (unsigned int i=0; i<dof_handler->n_dofs(); ++i)
	      {
		if (velocity_d_dofs[i])
		  {
		    Assert(i < n_u,
			    ExcMessage("Could not correctly locate velocity "
				      "dofs in velocity system!"));
		    null_vectors [d* n_u + i] = 1.;
		    ++counter;
		  }
	      }
	    Assert (counter == n_ud,
		    ExcMessage("Failed to extract correct components "
				"that should consitute null space!"));
	  }
	MLList.set("null space: dimension", dim);
	MLList.set("null space: vectors", &null_vectors[0]);
	MLList.set("null space: type", "pre-computed");
      }
    }
	
  ml_precond = new ML_Epetra::MultiLevelPreconditioner(*Matrix, MLList, true);

  if (output_amg_info)
    ml_precond->PrintUnused(0);
}

				 // For the implementation of the
				 // <code>vmult</code> function we
				 // note that invoking a call of 
				 // the Trilinos preconditioner 
				 // requires us to use Epetra vectors
				 // as well. Luckily, it is sufficient
				 // to provide a view, i.e., feed 
				 // Trilinos with a pointer to the
				 // data, so we avoid copying the
				 // content of the vectors during
				 // the iteration. In the declaration
				 // of the right hand side, we need
				 // to cast the source vector (that
				 // is <code>const</code> in all deal.II 
				 // calls) to non-constant value, as
				 // this is the way Trilinos wants to
				 // have them.
template <int dim>
void TrilinosAmgPreconditioner<dim>::vmult (Vector<double>       &dst,
					    const Vector<double> &src) const
{
  Epetra_Vector LHS (View, *Map, dst.begin());
  Epetra_Vector RHS (View, *Map, const_cast<double*>(src.begin()));
  
  int res = ml_precond->ApplyInverse (RHS, LHS);
  
  Assert (res == 0,
	  ExcMessage ("Trilinos AMG MultiLevel preconditioner returned "
		      "errorneously!"));
}



				 // @sect3{The <code>BoussinesqFlowProblem</code> class template}

				 // The definition of this class is
				 // mainly based on the step-22 tutorial
				 // program. Most of the data types are
				 // the same as there. However, we
				 // deal with a time-dependent system now,
				 // and there is temperature to take care
				 // of as well, so we need some additional
				 // function and variable declarations.
				 // Furthermore, we have a slightly more
				 // sophisticated solver we are going to
				 // use, so there is a second pointer
				 // to a sparse ILU for a pressure
				 // mass matrix as well.
template <int dim>
class BoussinesqFlowProblem
{
  public:
    BoussinesqFlowProblem (const unsigned int degree);
    void run ();

  private:
    void setup_dofs (const bool setup_matrices);
    void assemble_preconditioner ();
    void assemble_system ();
    void assemble_rhs_T ();
    double get_maximal_velocity () const;
    void solve ();
    void output_results () const;
    void refine_mesh ();

    const unsigned int   degree;

    Triangulation<dim>   triangulation;
    FESystem<dim>        fe;
    DoFHandler<dim>      dof_handler;

    ConstraintMatrix     hanging_node_constraints;

    BlockSparsityPattern      sparsity_pattern;
    BlockSparseMatrix<double> system_matrix;
    BlockSparseMatrix<double> preconditioner_matrix;

    double time_step;
    unsigned int timestep_number;

    BlockVector<double> solution;
    BlockVector<double> old_solution;
    BlockVector<double> system_rhs;

    //boost::shared_ptr<typename InnerPreconditioner<dim>::type> A_preconditioner;
    boost::shared_ptr<TrilinosAmgPreconditioner<dim> >  Amg_preconditioner;
    boost::shared_ptr<SparseILU<double> > Mp_preconditioner;

    bool rebuild_matrices;
    bool rebuild_preconditioner;
};




				 // @sect3{Equation data}

				 // Again, the next stage in the program
				 // is the definition of the equation 
				 // data, that is, the various
				 // boundary conditions, the right hand
				 // side and the initial condition (remember
				 // that we're about to solve a time-
				 // dependent system). The basic strategy
				 // for this definition is the same as in
				 // step-22. Regarding the details, though,
				 // there are some differences.

				 // The first
				 // thing is that we don't set any boundary
				 // conditions on the velocity, as is
				 // explained in the introduction. So
				 // what is left are two conditions for
				 // pressure <i>p</i> and temperature
				 // <i>T</i>.

				 // Secondly, we set an initial
				 // condition for all problem variables,
				 // i.e., for <b>u</b>, <i>p</i> and <i>T</i>,
				 // so the function has <i>dim+2</i>
				 // components.
				 // In this case, we choose a very simple
				 // test case, where everything is zero.

				 // @sect4{Boundary values}
template <int dim>
class PressureBoundaryValues : public Function<dim>
{
  public:
    PressureBoundaryValues () : Function<dim>(1) {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
};


template <int dim>
double
PressureBoundaryValues<dim>::value (const Point<dim>  &/*p*/,
                                    const unsigned int /*component*/) const
{
  return 0;
}



template <int dim>
class TemperatureBoundaryValues : public Function<dim>
{
  public:
    TemperatureBoundaryValues () : Function<dim>(1) {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
};



template <int dim>
double
TemperatureBoundaryValues<dim>::value (const Point<dim> &p,
                                      const unsigned int /*component*/) const
{
//TODO: leftover from olden times. replace by something sensible once we have
//diffusion in the temperature field
  if (p[0] == 0)
    return 1;
  else
    return 0;
}



				 // @sect4{Initial values}
template <int dim>
class InitialValues : public Function<dim>
{
  public:
    InitialValues () : Function<dim>(dim+2) {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;

    virtual void vector_value (const Point<dim> &p,
                               Vector<double>   &value) const;
};


template <int dim>
double
InitialValues<dim>::value (const Point<dim>  &,
                           const unsigned int) const
{
  return 0;
}


template <int dim>
void
InitialValues<dim>::vector_value (const Point<dim> &p,
                                  Vector<double>   &values) const
{
  for (unsigned int c=0; c<this->n_components; ++c)
    values(c) = InitialValues<dim>::value (p, c);
}



				 // @sect4{Right hand side}
				 // 
				 // The last definition of this kind
				 // is the one for the right hand
				 // side function. Again, the content
				 // of the function is very
				 // basic and zero in most of the
				 // components, except for a source
				 // of temperature in some isolated
				 // regions near the bottom of the
				 // computational domain, as is explained
				 // in the problem description in the
				 // introduction.
template <int dim>
class RightHandSide : public Function<dim>
{
  public:
    RightHandSide () : Function<dim>(dim+2) {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;

    virtual void vector_value (const Point<dim> &p,
                               Vector<double>   &value) const;
};


template <int dim>
double
RightHandSide<dim>::value (const Point<dim>  &p,
                           const unsigned int component) const
{
  if (component == dim+1)
    return ((p.distance (Point<dim>(.3,.1)) < 1./32)
	    ||
	    (p.distance (Point<dim>(.45,.1)) < 1./32)
	    ||
	    (p.distance (Point<dim>(.75,.1)) < 1./32)
	    ?
	    1
	    :
	    0);
  else
    return 0;
}


template <int dim>
void
RightHandSide<dim>::vector_value (const Point<dim> &p,
                                  Vector<double>   &values) const
{
  for (unsigned int c=0; c<this->n_components; ++c)
    values(c) = RightHandSide<dim>::value (p, c);
}




				 // @sect3{Linear solvers and preconditioners}

				 // This section introduces some
				 // objects that are used for the
				 // solution of the linear equations of
				 // Stokes system that we need to
				 // solve in each time step. The basic
				 // structure is still the same as
				 // in step-20, where Schur complement
				 // based preconditioners and solvers
				 // have been introduced, with the 
				 // actual interface taken from step-22.

				 // @sect4{The <code>InverseMatrix</code> class template}

				 // This class is an interface to
				 // calculate the action of an
				 // "inverted" matrix on a vector
				 // (using the <code>vmult</code>
				 // operation)
				 // in the same way as the corresponding
				 // function in step-22: when the
				 // product of an object of this class
				 // is requested, we solve a linear
				 // equation system with that matrix
				 // using the CG method, accelerated
				 // by a preconditioner of (templated) class
				 // <code>Preconditioner</code>.
template <class Matrix, class Preconditioner>
class InverseMatrix : public Subscriptor
{
  public:
    InverseMatrix (const Matrix         &m,
		   const Preconditioner &preconditioner);

    void vmult (Vector<double>       &dst,
                const Vector<double> &src) const;

  private:
    const SmartPointer<const Matrix> matrix;
    const Preconditioner &preconditioner;
};


template <class Matrix, class Preconditioner>
InverseMatrix<Matrix,Preconditioner>::InverseMatrix (const Matrix &m,
						     const Preconditioner &preconditioner)
                :
                matrix (&m),
		preconditioner (preconditioner)
{}



template <class Matrix, class Preconditioner>
void InverseMatrix<Matrix,Preconditioner>::vmult (Vector<double>       &dst,
						  const Vector<double> &src) const
{
  SolverControl solver_control (src.size(), 1e-6*src.l2_norm());
  SolverCG<> cg (solver_control);

  dst = 0;

  try
    {
      cg.solve (*matrix, dst, src, preconditioner);
    }
  catch (std::exception &e)
    {
      Assert (false, ExcMessage(e.what()));
    }
}

				 // @sect4{Schur complement preconditioner}

				 // This is the implementation
				 // of the Schur complement
				 // preconditioner as described
				 // in the section on improved
				 // solvers in step-22.
				 // 
				 // The basic 
				 // concept of the preconditioner is 
				 // different to the solution 
				 // strategy used in step-20 and 
				 // step-22. There, the Schur
				 // complement was used for a 
				 // two-stage solution of the linear
				 // system. Recall that the process
				 // in the Schur complement solver is
				 // a Gaussian elimination of
				 // a 2x2 block matrix, where each
				 // block is solved iteratively. 
				 // Here, the idea is to let 
				 // an iterative solver act on the
				 // whole system, and to use 
				 // a Schur complement for 
				 // preconditioning. As usual when
				 // dealing with preconditioners, we
				 // don't intend to exacly set up a 
				 // Schur complement, but rather use
				 // a good approximation to the
				 // Schur complement for the purpose of
				 // preconditioning.
				 // 
				 // So the question is how we can
				 // obtain a good preconditioner.
				 // Let's have a look at the 
				 // preconditioner matrix <i>P</i>
				 // acting on the block system, built
				 // as
				 // @f{eqnarray*}
				 //   P^{-1}
				 //   = 
				 //   \left(\begin{array}{cc}
				 //     A^{-1} & 0 \\ S^{-1} B A^{-1} & -S^{-1}
				 //   \end{array}\right)
				 // @f}
				 // using the Schur complement 
				 // $S = B A^{-1} B^T$. If we apply
				 // this matrix in the solution of 
				 // a linear system, convergence of
				 // an iterative Krylov-based solver
				 // will be governed by the matrix
				 // @f{eqnarray*}
				 //   P^{-1}\left(\begin{array}{cc}
				 //     A & B^T \\ B & 0
				 //   \end{array}\right) 
 				 //  = 
				 //   \left(\begin{array}{cc}
				 //     I & A^{-1} B^T \\ 0 & 0
				 //   \end{array}\right),
				 // @f}
				 // which turns out to be very simple.
				 // A GMRES solver based on exact
				 // matrices would converge in two
				 // iterations, since there are
				 // only two distinct eigenvalues.
				 // Such a preconditioner for the
				 // blocked Stokes system has been 
				 // proposed by Silvester and Wathen,
				 // Fast iterative solution of 
				 // stabilised Stokes systems part II. 
				 // Using general block preconditioners.
				 // (SIAM J. Numer. Anal., 31 (1994),
				 // pp. 1352-1367).
				 // 
				 // The deal.II users who have already
				 // gone through the step-20 and step-22 
				 // tutorials can certainly imagine
				 // how we're going to implement this.
				 // We replace the inverse matrices
				 // in $P^{-1}$ using the InverseMatrix
				 // class, and the inverse Schur 
				 // complement will be approximated
				 // by the pressure mass matrix $M_p$.
				 // Having this in mind, we define a
				 // preconditioner class with a 
				 // <code>vmult</code> functionality,
				 // which is all we need for the
				 // interaction with the usual solver
				 // functions further below in the
				 // program code.
				 // 
				 // First the declarations. These
				 // are similar to the definition of
				 // the Schur complement in step-20,
				 // with the difference that we need
				 // some more preconditioners in
				 // the constructor.
template <class PreconditionerA, class PreconditionerMp>
class BlockSchurPreconditioner : public Subscriptor
{
  public:
    BlockSchurPreconditioner (const BlockSparseMatrix<double>         &S,
          const InverseMatrix<SparseMatrix<double>,PreconditionerMp>  &Mpinv,
          const PreconditionerA &Apreconditioner);

  void vmult (BlockVector<double>       &dst,
              const BlockVector<double> &src) const;

  private:
    const SmartPointer<const BlockSparseMatrix<double> > system_matrix;
    const SmartPointer<const InverseMatrix<SparseMatrix<double>,
                       PreconditionerMp > > m_inverse;
    const PreconditionerA &a_preconditioner;

    mutable Vector<double> tmp;

};

template <class PreconditionerA, class PreconditionerMp>
BlockSchurPreconditioner<PreconditionerA, PreconditionerMp>::BlockSchurPreconditioner(
          const BlockSparseMatrix<double>                            &S,
          const InverseMatrix<SparseMatrix<double>,PreconditionerMp> &Mpinv,
          const PreconditionerA &Apreconditioner
          )
                :
                system_matrix           (&S),
                m_inverse               (&Mpinv),
                a_preconditioner        (Apreconditioner),
                tmp                     (S.block(1,1).m())
{
}


				 // This is the <code>vmult</code>
				 // function. We implement
				 // the action of $P^{-1}$ as described
				 // above in three successive steps.
				 // The first step multiplies
				 // the velocity vector by a 
				 // preconditioner of the matrix <i>A</i>.
				 // The resuling velocity vector
				 // is then multiplied by $B$ and
				 // subtracted from the pressure.
				 // This second step only acts on 
				 // the pressure vector and is 
				 // accomplished by the command
				 // SparseMatrix::residual. Next, 
				 // we change the sign in the 
				 // temporary pressure vector and
				 // finally multiply by the pressure
				 // mass matrix to get the final
				 // pressure vector.
template <class PreconditionerA, class PreconditionerMp>
void BlockSchurPreconditioner<PreconditionerA, PreconditionerMp>::vmult (
                                     BlockVector<double>       &dst,
                                     const BlockVector<double> &src) const
{
  a_preconditioner.vmult (dst.block(0), src.block(0));
  system_matrix->block(1,0).residual(tmp, dst.block(0), src.block(1));
  tmp *= -1;
  m_inverse->vmult (dst.block(1), tmp);
}



				 // @sect3{BoussinesqFlowProblem class implementation}

				 // @sect4{BoussinesqFlowProblem::BoussinesqFlowProblem}
				 // 
				 // The constructor of this class is
				 // an extension of the constructor
				 // in step-22. We need to include 
				 // the temperature in the definition
				 // of the finite element. As discussed
				 // in the introduction, we are going 
				 // to use discontinuous elements 
				 // of one degree less than for pressure
				 // there. Moreover, we initialize
				 // the time stepping as well as the
				 // options for the matrix assembly 
				 // and preconditioning.
template <int dim>
BoussinesqFlowProblem<dim>::BoussinesqFlowProblem (const unsigned int degree)
                :
                degree (degree),
                fe (FE_Q<dim>(degree+1), dim,
                    FE_Q<dim>(degree), 1,
                    FE_DGQ<dim>(degree-1), 1),
                dof_handler (triangulation),
                time_step (0),
		rebuild_matrices (true),
		rebuild_preconditioner (true)
{}



				 // @sect4{BoussinesqFlowProblem::setup_dofs}
				 // 
				 // This function does the same as
				 // in most other tutorial programs. 
				 // As a slight difference, the 
				 // program is called with a 
				 // parameter <code>setup_matrices</code>
				 // that decides whether to 
				 // recreate the sparsity pattern
				 // and the associated stiffness
				 // matrix.
				 // 
				 // The body starts by assigning dofs on 
				 // basis of the chosen finite element,
				 // and then renumbers the dofs 
				 // first using the Cuthill_McKee
				 // algorithm (to generate a good
				 // quality ILU during the linear
				 // solution process) and then group
				 // components of velocity, pressure
				 // and temperature together. This 
				 // happens in complete analogy to
				 // step-22.
				 // 
				 // We then proceed with the generation
				 // of the hanging node constraints
				 // that arise from adaptive grid
				 // refinement. Next we impose
				 // the no-flux boundary conditions
				 // $\vec{u}\cdot \vec{n}=0$ by adding
				 // a respective constraint to the
				 // hanging node constraints
				 // matrix. The second parameter in 
				 // the function describes the first 
				 // of the velocity components
				 // in the total dof vector, which is 
				 // zero here. The parameter 
				 // <code>no_normal_flux_boundaries</code>
				 // sets the no flux b.c. to those
				 // boundaries with boundary indicator
				 // zero.
template <int dim>
void BoussinesqFlowProblem<dim>::setup_dofs (const bool setup_matrices)
{
  dof_handler.distribute_dofs (fe);
  DoFRenumbering::Cuthill_McKee (dof_handler);
  std::vector<unsigned int> block_component (dim+2,0);
  block_component[dim] = 1;
  block_component[dim+1] = 2;
  DoFRenumbering::component_wise (dof_handler, block_component);

  hanging_node_constraints.clear ();
  DoFTools::make_hanging_node_constraints (dof_handler,
					   hanging_node_constraints);
  std::set<unsigned char> no_normal_flux_boundaries;
  no_normal_flux_boundaries.insert (0);
  VectorTools::compute_no_normal_flux_constraints (dof_handler, 0,
						   no_normal_flux_boundaries,
						   hanging_node_constraints);
  hanging_node_constraints.close ();

				 // The next step is, as usual, 
				 // to write some information
				 // to the screen. The information
				 // that is most interesting during
				 // the calculations is the
				 // number of degrees of freedom
				 // in the individual components,
				 // so we count them. The function 
				 // to do this is the same as the
				 // one used in step-22, which 
				 // uses the grouping of all
				 // velocity components into
				 // one block as introduced
				 // above.
  std::vector<unsigned int> dofs_per_block (3);
  DoFTools::count_dofs_per_block (dof_handler, dofs_per_block, block_component);  
  const unsigned int n_u = dofs_per_block[0],
                     n_p = dofs_per_block[1],
		     n_T = dofs_per_block[2];

  std::cout << "Number of active cells: "
            << triangulation.n_active_cells()
            << std::endl
            << "Number of degrees of freedom: "
            << dof_handler.n_dofs()
            << " (" << n_u << '+' << n_p << '+'<< n_T <<')'
            << std::endl
            << std::endl;

				 // The next step is to 
				 // create the sparsity 
				 // pattern for the system matrix 
				 // based on the Boussinesq 
				 // system. As in step-22, 
				 // we choose to create the
				 // pattern not as in the
				 // first tutorial programs,
				 // but by using the blocked
				 // version of 
				 // CompressedSetSparsityPattern.
				 // The reason for doing this 
				 // is mainly a memory issue,
				 // that is, the basic procedures
				 // consume too much memory
				 // when used in three spatial
				 // dimensions as we intend
				 // to do for this program.
				 // 
				 // So, in case we need
				 // to recreate the matrices,
				 // we first release the
				 // stiffness matrix from the
				 // sparsity pattern and then
				 // set up an object of the 
				 // BlockCompressedSetSparsityPattern
				 // consisting of three blocks. 
				 // Each of these blocks is
				 // initialized with the
				 // respective number of 
				 // degrees of freedom. 
				 // Once the blocks are 
				 // created, the overall size
				 // of the sparsity pattern
				 // is initiated by invoking 
				 // the <code>collect_sizes()</code>
				 // command, and then the
				 // sparsity pattern can be
				 // filled with information.
				 // Then, the hanging
				 // node constraints are applied
				 // to the temporary sparsity
				 // pattern, which is finally
				 // then completed and copied
				 // into the general sparsity
				 // pattern structure.
  
				 // Observe that we use a 
				 // coupling argument for 
				 // telling the function
				 // <code>make_sparsity_pattern</code>
				 // which components actually
				 // will hold data and which 
				 // we're going to neglect.
				 // 
				 // After these actions, we 
				 // need to reassign the 
				 // system matrix structure to
				 // the sparsity pattern.
  if (setup_matrices == true)
    {
      system_matrix.clear ();
      preconditioner_matrix.clear ();

      BlockCompressedSetSparsityPattern csp (3,3);
 
      csp.block(0,0).reinit (n_u, n_u);
      csp.block(0,1).reinit (n_u, n_p);
      csp.block(0,2).reinit (n_u, n_T);
      csp.block(1,0).reinit (n_p, n_u);
      csp.block(1,1).reinit (n_p, n_p);
      csp.block(1,2).reinit (n_p, n_T);
      csp.block(2,0).reinit (n_T, n_u);
      csp.block(2,1).reinit (n_T, n_p);
      csp.block(2,2).reinit (n_T, n_T);
      
      csp.collect_sizes ();

      Table<2,DoFTools::Coupling> coupling (dim+2, dim+2);
      
      for (unsigned int component = 0; component<dim+2; ++component)
	for (unsigned int component2 = 0; component2<dim+2; ++component2)
	  coupling[component][component2] = DoFTools::always;
      
      for (unsigned int component = 0; component<dim+1; ++component)
	{
	  coupling[dim+1][component] = DoFTools::none;
	  coupling[component][dim+1] = DoFTools::none;
	}
      
      DoFTools::make_sparsity_pattern (dof_handler, coupling, csp);
      hanging_node_constraints.condense (csp);
      sparsity_pattern.copy_from (csp);

      system_matrix.reinit (sparsity_pattern);
    }

				 // As last action in this function,
				 // we need to set the vectors
				 // for the solution, the old 
				 // solution (required for 
				 // time stepping) and the system
				 // right hand side to the 
				 // three-block structure given
				 // by velocity, pressure and
				 // temperature.
  solution.reinit (3);
  solution.block(0).reinit (n_u);
  solution.block(1).reinit (n_p);
  solution.block(2).reinit (n_T);
  solution.collect_sizes ();

  old_solution.reinit (3);
  old_solution.block(0).reinit (n_u);
  old_solution.block(1).reinit (n_p);
  old_solution.block(2).reinit (n_T);
  old_solution.collect_sizes ();

  system_rhs.reinit (3);
  system_rhs.block(0).reinit (n_u);
  system_rhs.block(1).reinit (n_p);
  system_rhs.block(2).reinit (n_T);
  system_rhs.collect_sizes ();
}



template <int dim>
double scalar_product (const Tensor<2,dim> &t1,
		       const Tensor<2,dim> &t2)
{
  double s = 0;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      s += t1[i][j] * t2[i][j];
  return s;
}


template <int dim>
void
BoussinesqFlowProblem<dim>::assemble_preconditioner ()
{
  preconditioner_matrix = 0;

  QGauss<dim>   quadrature_formula(degree+2);
  QGauss<dim-1> face_quadrature_formula(degree+2);

  FEValues<dim> fe_values (fe, quadrature_formula,
			   update_JxW_values |
			   update_gradients);
  const unsigned int   dofs_per_cell   = fe.dofs_per_cell;

  const unsigned int   n_q_points      = quadrature_formula.size();

  FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

  std::vector<Tensor<2,dim> > phi_grad_u (dofs_per_cell);

  const FEValuesExtractors::Vector velocities (0);

  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      local_matrix = 0;

      for (unsigned int q=0; q<n_q_points; ++q)
	{
	  for (unsigned int k=0; k<dofs_per_cell; ++k)
	    phi_grad_u[k] = fe_values[velocities].gradient(k,q);

	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    for (unsigned int j=0; j<dofs_per_cell; ++j)
	      local_matrix(i,j) += scalar_product (phi_grad_u[i], phi_grad_u[j])
				   * fe_values.JxW(q);
	}

      cell->get_dof_indices (local_dof_indices);

      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=0; j<dofs_per_cell; ++j)
	  preconditioner_matrix.add (local_dof_indices[i],
				     local_dof_indices[j],
				     local_matrix(i,j));
    }
  
  hanging_node_constraints.condense (preconditioner_matrix);
}



				 // @sect4{BoussinesqFlowProblem::assemble_system}
				 // 
				 // The assembly of the Boussinesq 
				 // system is acutally a two-step
				 // procedure. One is to create
				 // the Stokes system matrix and
				 // right hand side for the 
				 // velocity-pressure system as
				 // well as the mass matrix for
				 // temperature, and
				 // the second is to create the
				 // rhight hand side for the temperature
				 // dofs. The reason for doing this
				 // in two steps is simply that 
				 // the time stepping we have chosen
				 // needs the result from the Stokes
				 // system at the current time step
				 // for building the right hand
				 // side of the temperature equation.
				 // 
				 // This function does the 
				 // first of these two tasks.
				 // There are two different situations
				 // for calling this function. The
				 // first one is when we reset the
				 // mesh, and both the matrix and
				 // the right hand side have to
				 // be generated. The second situation
				 // only sets up the right hand
				 // side. The reason for having 
				 // two different accesses is that
				 // the matrix of the Stokes system
				 // does not change in time unless
				 // the mesh is changed, so we can
				 // save a considerable amount of
				 // work by doing the full assembly
				 // only when it is needed.
				 // 
				 // Regarding the technical details
				 // of implementation, not much has
				 // changed from step-22. We reset
				 // matrix and vector, create 
				 // a quadrature formula on the 
				 // cells and one on cell faces
				 // (for implementing Neumann 
				 // boundary conditions). Then,
				 // we create a respective
				 // FEValues object for both the 
				 // cell and the face integration.
				 // For the the update flags of
				 // the first, we perform the
				 // calculations of basis function
				 // derivatives only in
				 // case of a full assembly, since
				 // they are not needed otherwise,
				 // which makes the call of
				 // the FEValues::reinit function
				 // further down in the program 
				 // more efficient.
				 // 
				 // The declarations proceed 
				 // with some shortcuts for 
				 // array sizes, the creation of
				 // the local matrix and right 
				 // hand side as well as the
				 // vector for the indices of
				 // the local dofs compared to
				 // the global system.
template <int dim>
void BoussinesqFlowProblem<dim>::assemble_system ()
{
  if (rebuild_matrices == true)
    system_matrix=0;

  system_rhs=0;

  QGauss<dim>   quadrature_formula(degree+2);
  QGauss<dim-1> face_quadrature_formula(degree+2);

  FEValues<dim> fe_values (fe, quadrature_formula,
			   update_values    |
			   update_quadrature_points  |
			   update_JxW_values |
			   (rebuild_matrices == true
			    ?
			    update_gradients
			    :
			    UpdateFlags(0)));
  FEFaceValues<dim> fe_face_values (fe, face_quadrature_formula,
				    update_values    | 
				    update_normal_vectors |
				    update_quadrature_points  | 
				    update_JxW_values);

  const unsigned int   dofs_per_cell   = fe.dofs_per_cell;

  const unsigned int   n_q_points      = quadrature_formula.size();
  const unsigned int   n_face_q_points = face_quadrature_formula.size();

  FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       local_rhs (dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

				 // These few declarations provide
				 // the structures for the evaluation
				 // of inhomogeneous Neumann boundary
				 // conditions from the function
				 // declaration made above.
				 // The vector <code>old_solution_values</code>
				 // evaluates the solution 
				 // at the old time level, since
				 // the temperature from the
				 // old time level enters the 
				 // Stokes system as a source
				 // term in the momentum equation.
				 // 
				 // Then, we create a variable
				 // to hold the Rayleigh number,
				 // the measure of buoyancy.
				 // 
				 // The set of vectors we create
				 // next hold the evaluations of
				 // the basis functions that will
				 // be used for creating the
				 // matrices. This gives faster
				 // access to that data, which
				 // increases the performance
				 // of the assembly. See step-22 
				 // for details.
				 // 
				 // The last few declarations 
				 // are used to extract the 
				 // individual blocks (velocity,
				 // pressure, temperature) from
				 // the total FE system.
  const PressureBoundaryValues<dim> pressure_boundary_values;
  std::vector<double>               boundary_values (n_face_q_points);

  std::vector<Vector<double> >      old_solution_values(n_q_points,
							Vector<double>(dim+2));

  const double Rayleigh_number = 10;

  std::vector<Tensor<1,dim> >          phi_u       (dofs_per_cell);
  std::vector<SymmetricTensor<2,dim> > phi_grads_u (dofs_per_cell);
  std::vector<double>                  div_phi_u   (dofs_per_cell);
  std::vector<double>                  phi_p       (dofs_per_cell);
  std::vector<double>                  phi_T       (dofs_per_cell);
  std::vector<Tensor<1,dim> >          grad_phi_T  (dofs_per_cell);

  const FEValuesExtractors::Vector velocities (0);
  const FEValuesExtractors::Scalar pressure (dim);
  const FEValuesExtractors::Scalar temperature (dim+1);

				 // Now starts the loop over
				 // all cells in the problem.
				 // The first commands are all
				 // very familiar, doing the
				 // evaluations of the element
				 // basis functions, resetting
				 // the local arrays and 
				 // getting the values of the
				 // old solution at the
				 // quadrature point. Then we
				 // are ready to loop over
				 // the quadrature points 
				 // on the cell.
  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      local_matrix = 0;
      local_rhs = 0;

      fe_values.get_function_values (old_solution, old_solution_values);

      for (unsigned int q=0; q<n_q_points; ++q)
	{
	  const double old_temperature = old_solution_values[q](dim+1);

				// Extract the basis relevant
				// terms in the inner products
				// once in advance as shown
				// in step-22 in order to 
				// accelerate assembly.
				// 
				// Once this is done, we 
				// start the loop over the
				// rows and columns of the
				// local matrix and feed
				// the matrix with the relevant
				// products. The right hand
				// side is filled with the 
				// forcing term driven by
				// temperature in direction
				// of gravity (which is 
				// vertical in our example).
				// Note that the right hand 
				// side term is always generated,
				// whereas the matrix 
				// contributions are only
				// updated when it is 
				// requested by the
				// <code>rebuild_matrices</code>
				// flag.
	  for (unsigned int k=0; k<dofs_per_cell; ++k)
	    {
	      phi_u[k] = fe_values[velocities].value (k,q);
	      if (rebuild_matrices)
	        {
		  phi_grads_u[k] = fe_values[velocities].symmetric_gradient(k,q);
		  div_phi_u[k]   = fe_values[velocities].divergence (k, q);
		  phi_p[k]       = fe_values[pressure].value (k, q);
		  phi_T[k]       = fe_values[temperature].value (k, q);
		  grad_phi_T[k]  = fe_values[temperature].gradient (k, q);
		}
	    }

	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    {

	      const Tensor<1,dim> phi_i_u = fe_values[velocities].value (i, q);

	      if (rebuild_matrices)
		for (unsigned int j=0; j<dofs_per_cell; ++j)
		  local_matrix(i,j) += (phi_grads_u[i] * phi_grads_u[j]
					- div_phi_u[i] * phi_p[j]
					- phi_p[i] * div_phi_u[j]
					+ phi_p[i] * phi_p[j]
					+ phi_T[i] * phi_T[j])
				       * fe_values.JxW(q);

	      const Point<dim> gravity = ( (dim == 2) ? (Point<dim> (0,1)) : 
						        (Point<dim> (0,1,0)) );

	      local_rhs(i) += (Rayleigh_number *
			       gravity * phi_u[i] * old_temperature)*
			      fe_values.JxW(q);
          }
	}

				 // Next follows the assembly 
				 // of the face terms, result
				 // from Neumann boundary 
				 // conditions. Since these
				 // terms only enter the right
				 // hand side vector and not
				 // the matrix, there is no
				 // substantial benefit from
				 // extracting the data 
				 // before using it, so 
				 // we remain in the lines 
				 // of step-20 at this point.
      for (unsigned int face_no=0;
           face_no<GeometryInfo<dim>::faces_per_cell;
           ++face_no)
        if (cell->at_boundary(face_no))
          {
            fe_face_values.reinit (cell, face_no);

            pressure_boundary_values
              .value_list (fe_face_values.get_quadrature_points(),
                           boundary_values);

            for (unsigned int q=0; q<n_face_q_points; ++q)
              for (unsigned int i=0; i<dofs_per_cell; ++i)
                {
                  const Tensor<1,dim>
                    phi_i_u = fe_face_values[velocities].value (i, q);

                  local_rhs(i) += -(phi_i_u *
                                    fe_face_values.normal_vector(q) *
                                    boundary_values[q] *
                                    fe_face_values.JxW(q));
                }
          }

				 // The last step in the loop 
				 // over all cells is to
				 // enter the local contributions
				 // into the global matrix and 
				 // vector structures to the
				 // positions specified in 
				 // <code>local_dof_indices</code>.
				 // Again, we only add the 
				 // matrix data when it is 
				 // requested.
      cell->get_dof_indices (local_dof_indices);

      if (rebuild_matrices == true)
	{
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    for (unsigned int j=0; j<dofs_per_cell; ++j)
	      if (std::abs(local_matrix(i,j)) > 1e-20)
		system_matrix.add (local_dof_indices[i],
				   local_dof_indices[j],
				   local_matrix(i,j));
	}

      for (unsigned int i=0; i<dofs_per_cell; ++i)
        system_rhs(local_dof_indices[i]) += local_rhs(i);
    }

				 // Back at the outermost
				 // level of this function,
				 // we continue the work
				 // by condensing hanging
				 // node constraints to the
				 // right hand side and, 
				 // possibly, to the matrix.
  if (rebuild_matrices == true)
    hanging_node_constraints.condense (system_matrix);

  hanging_node_constraints.condense (system_rhs);

  if (rebuild_matrices == true)
    {
//       std::map<unsigned int,double> boundary_values;

//       typename DoFHandler<dim>::active_cell_iterator
// 	cell = dof_handler.begin_active(),
// 	emdc = dof_handler.end();
//       for (; cell!=endc; ++cell)
// 	for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
// 	  if (cell->vertex(v).distance(dim == 2
// 				       ?
// 				       Point<dim>(0,-1)
// 				       :
// 				       Point<dim>(0,0,-1)) < 1e-6)
// 	    {
// 	      std::cout << "Found cell and vertex: " << cell << ' '
// 			<< v << std::endl;

// 	      boundary_values[cell->vertex_dof_index(v,0)] = 0;
// 	      break;
// 	    }

//      std::vector<bool> component_mask (dim+2, true);
//       component_mask[dim] = component_mask[dim+1] = false;
//       VectorTools::interpolate_boundary_values (dof_handler,
// 						0,
// 						ZeroFunction<dim>(dim+2),
// 						boundary_values,
// 						component_mask);

//       MatrixTools::apply_boundary_values (boundary_values,
// 					  system_matrix,
// 					  solution,
// 					  system_rhs);
    }

				 // This last step of the assembly
				 // function sets up the preconditioners
				 // used for the solution of the
				 // system. We are going to use an
				 // ILU preconditioner for the
				 // velocity block (to be used
				 // by BlockSchurPreconditioner class)
				 // as well as an ILU preconditioner
				 // for the inversion of the 
				 // pressure mass matrix. Recall that
				 // the velocity-velocity block sits
				 // at position (0,0) in the 
				 // global system matrix, and
				 // the pressure mass matrix in
				 // (1,1). The 
				 // storage of these objects is
				 // as in step-22, that is, we
				 // include them using a 
				 // shared pointer structure from the
				 // boost library.
				 // 
				 // When all work is done, we 
				 // change the flags 
				 // <code>rebuild_preconditioner</code>
				 // and 
				 // <code>rebuild_matrices</code>
				 // to false.
  if (rebuild_preconditioner == true)
    {
      Assert (rebuild_matrices == true,
	      ExcMessage ("There is no point in rebuilding the preconditioner "
			  "without a rebuilt matrix!"));

      std::cout << "   Rebuilding preconditioner..." << std::flush;
      
      preconditioner_matrix.reinit (sparsity_pattern);
      assemble_preconditioner ();
      
      /*A_preconditioner
	= boost::shared_ptr<typename InnerPreconditioner<dim>::type>
		(new typename InnerPreconditioner<dim>::type());
      A_preconditioner->initialize (preconditioner_matrix.block(0,0),
		typename InnerPreconditioner<dim>::type::AdditionalData());*/
      
      Amg_preconditioner = 
	boost::shared_ptr<TrilinosAmgPreconditioner<dim> >
	  (new TrilinosAmgPreconditioner<dim>(dof_handler,
					      preconditioner_matrix.block(0,0),
					      true));

      Mp_preconditioner
	= boost::shared_ptr<SparseILU<double> >
		(new SparseILU<double>);
      Mp_preconditioner->initialize (system_matrix.block(1,1),
				     SparseILU<double>::AdditionalData());
      
				 // Throw away the preconditioner
				 // matrix since everything has been
				 // copied to the Epetra objects
				 // of the preconditioner.
      preconditioner_matrix.clear ();

      std::cout << std::endl;

      rebuild_preconditioner = false;
    }

  rebuild_matrices = false;
}





				 // @sect4{BoussinesqFlowProblem::assemble_rhs_T}
				 // 
				 // This function does the second
				 // part of the assembly work, the
				 // creation of the velocity-dependent
				 // right hand side of the
				 // temperature equation. The 
				 // declarations in this function
				 // are pretty much the same as the
				 // ones used in the other 
				 // assembly routine, except that we
				 // restrict ourselves to vectors
				 // this time. Though, we need to
				 // perform more face integrals 
				 // at this point, induced by the
				 // use of discontinuous elements for 
				 // the temperature (just
				 // as it was in the first DG 
				 // example in step-12) in combination
				 // with adaptive grid refinement
				 // and subfaces. The update 
				 // flags at face level are the 
				 // same as in step-12.
template <int dim>
void BoussinesqFlowProblem<dim>::assemble_rhs_T ()
{
  QGauss<dim>   quadrature_formula(degree+2);
  QGauss<dim-1> face_quadrature_formula(degree+2);
  FEValues<dim> fe_values (fe, quadrature_formula,
                           update_values    | update_gradients |
                           update_quadrature_points  | update_JxW_values);
  FEFaceValues<dim> fe_face_values (fe, face_quadrature_formula,
				    update_values    | update_normal_vectors |
				    update_quadrature_points  |
				    update_JxW_values);
  FESubfaceValues<dim> fe_subface_values (fe, face_quadrature_formula,
					  update_values | 
					  update_normal_vectors |
					  update_JxW_values);
  FEFaceValues<dim> fe_face_values_neighbor (fe, face_quadrature_formula,
                                             update_values);
  FESubfaceValues<dim> fe_subface_values_neighbor (fe, face_quadrature_formula,
						   update_values);

  const unsigned int   dofs_per_cell   = fe.dofs_per_cell;
  const unsigned int   n_q_points      = quadrature_formula.size();
  const unsigned int   n_face_q_points = face_quadrature_formula.size();

  Vector<double>       local_rhs (dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

				 // Here comes the declaration
				 // of vectors to hold the old
				 // and present solution values
				 // and gradients
				 // for both the cell as well as faces
				 // to the cell. Next comes the
				 // declaration of an object
				 // to hold the temperature 
				 // boundary values and a
				 // well-known extractor for
				 // accessing the temperature
				 // part of the FE system.
  std::vector<Vector<double> > old_solution_values(n_q_points,
						   Vector<double>(dim+2));

  std::vector<Vector<double> > old_solution_values_face(n_face_q_points, 
							Vector<double>(dim+2));
  std::vector<Vector<double> > old_solution_values_face_neighbor (
							n_face_q_points,
							Vector<double>(dim+2));
  std::vector<Vector<double> > present_solution_values (n_q_points, 
							Vector<double>(dim+2));
  std::vector<Vector<double> > present_solution_values_face(
							n_face_q_points, 
							Vector<double>(dim+2));

  std::vector<std::vector<Tensor<1,dim> > >  present_solution_grads(
				  n_q_points,
				  std::vector<Tensor<1,dim> >(dim+2));

  std::vector<double> neighbor_temperature (n_face_q_points);

  TemperatureBoundaryValues<dim> temperature_boundary_values;
  const FEValuesExtractors::Scalar temperature (dim+1);

				 // Now, let's start the loop
				 // over all cells in the
				 // triangulation. The first
				 // actions within the loop
				 // are, as usual, the evaluation
				 // of the FE basis functions 
				 // and the old and present
				 // solution at the quadrature 
				 // points.
  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

  for (; cell!=endc; ++cell)
    {
      local_rhs = 0;
      fe_values.reinit (cell);

      fe_values.get_function_values (old_solution, old_solution_values);
      fe_values.get_function_values (solution, present_solution_values);
      fe_values.get_function_gradients (solution, present_solution_grads);

      for (unsigned int q=0; q<n_q_points; ++q)
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            const double old_T = old_solution_values[q](dim+1);
            Tensor<1,dim> present_u;
            for (unsigned int d=0; d<dim; ++d)
              present_u[d] = present_solution_values[q](d);

	    double present_div_u = 0;
            for (unsigned int d=0; d<dim; ++d)
              present_div_u += present_solution_grads[q][d][d];

            const double        phi_i_T      = fe_values[temperature].value (i, q);
            const Tensor<1,dim> grad_phi_i_T = fe_values[temperature].gradient (i, q);

	    const Point<dim> p = fe_values.quadrature_point(q);

            local_rhs(i) += (time_step *
                             old_T *
                             (present_u *
			      grad_phi_i_T
			      +
			      present_div_u *
			      phi_i_T)
                             +
                             old_T * phi_i_T
			     +
			     time_step *
			     RightHandSide<dim>().value (p, dim+1)
			     * phi_i_T)
                            *
                            fe_values.JxW(q);
          }


//TODO: unify the code that actually does the assembly down below
      for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell;
           ++face_no)
        if (cell->at_boundary(face_no)
	    ||
	    ((cell->neighbor(face_no)->has_children() == false)
	     &&
	     (cell->neighbor(face_no)->level() == cell->level())))
	  {
					     // cell either at
					     // boundary or with a
					     // neighbor that has the
					     // same refinement level
					     // and is not further
					     // refined
	    fe_face_values.reinit (cell, face_no);

	    fe_face_values.get_function_values (old_solution,
						old_solution_values_face);
	    fe_face_values.get_function_values (solution,
						present_solution_values_face);

	    if (cell->at_boundary(face_no))
	      temperature_boundary_values
		.value_list (fe_face_values.get_quadrature_points(),
			     neighbor_temperature);
	    else
	      {
		const typename DoFHandler<dim>::active_cell_iterator
		  neighbor = cell->neighbor(face_no);

		fe_face_values_neighbor.reinit (neighbor,
						cell->neighbor_of_neighbor(face_no));

		fe_face_values_neighbor
		  .get_function_values (old_solution,
					old_solution_values_face_neighbor);

		for (unsigned int q=0; q<n_face_q_points; ++q)
		  neighbor_temperature[q] = old_solution_values_face_neighbor[q](dim+1);
	      }

	    for (unsigned int q=0; q<n_face_q_points; ++q)
	      {
		Tensor<1,dim> present_u_face;
		for (unsigned int d=0; d<dim; ++d)
		  present_u_face[d] = present_solution_values_face[q](d);

		const double normal_flux = present_u_face *
					   fe_face_values.normal_vector(q);

		const bool is_outflow_q_point = (normal_flux >= 0);

		for (unsigned int i=0; i<dofs_per_cell; ++i)
		  local_rhs(i) -= time_step *
				  normal_flux *
				  (is_outflow_q_point == true
				   ?
				   old_solution_values_face[q](dim+1)
				   :
				   neighbor_temperature[q]) *
				  fe_face_values[temperature].value (i,q) *
				  fe_face_values.JxW(q);
	      }
	  }
	else
	  if (cell->neighbor(face_no)->has_children())
	    {
					       // neighbor is further
					       // refined. loop over
					       // all sub faces
	      for (unsigned int subface_no=0;
		   subface_no<GeometryInfo<dim>::max_children_per_face;
		   ++subface_no)
		{
		  fe_subface_values.reinit (cell, face_no, subface_no);

		  fe_subface_values.get_function_values (old_solution,
							 old_solution_values_face);
		  fe_subface_values.get_function_values (solution,
							 present_solution_values_face);

		  const typename DoFHandler<dim>::active_cell_iterator
		    neighbor = cell->neighbor_child_on_subface (face_no, subface_no);

		  fe_face_values_neighbor.reinit (neighbor,
						  cell->neighbor_of_neighbor(face_no));

		  fe_face_values_neighbor
		    .get_function_values (old_solution,
					  old_solution_values_face_neighbor);

		  for (unsigned int q=0; q<n_face_q_points; ++q)
		    neighbor_temperature[q] = old_solution_values_face_neighbor[q](dim+1);

		  for (unsigned int q=0; q<n_face_q_points; ++q)
		    {
		      Tensor<1,dim> present_u_face;
		      for (unsigned int d=0; d<dim; ++d)
			present_u_face[d] = present_solution_values_face[q](d);

		      const double normal_flux = present_u_face *
						 fe_subface_values.normal_vector(q);

		      const bool is_outflow_q_point = (normal_flux >= 0);

		      for (unsigned int i=0; i<dofs_per_cell; ++i)
			local_rhs(i) -= time_step *
					normal_flux *
					(is_outflow_q_point == true
					 ?
					 old_solution_values_face[q](dim+1)
					 :
					 neighbor_temperature[q]) *
					fe_face_values[temperature].value (i,q) *
					fe_face_values.JxW(q);
		    }
		}
	    }
	  else
	    {
					       // neighbor is less
					       // refined. we need to
					       // use a subface values
					       // object for the
					       // neighbor's subface
	      fe_face_values.reinit (cell, face_no);

	      fe_face_values.get_function_values (old_solution, old_solution_values_face);
	      fe_face_values.get_function_values (solution, present_solution_values_face);

	      const typename DoFHandler<dim>::active_cell_iterator
		neighbor = cell->neighbor (face_no);

	      const std::pair<unsigned int, unsigned int> faceno_subfaceno=
		cell->neighbor_of_coarser_neighbor(face_no);
	      const unsigned int neighbor_face_no    = faceno_subfaceno.first,
				 neighbor_subface_no = faceno_subfaceno.second;

	      Assert (neighbor->neighbor_child_on_subface (neighbor_face_no,
							   neighbor_subface_no)
		      == cell,
		      ExcInternalError());

	      fe_subface_values_neighbor.reinit (neighbor,
						 neighbor_face_no,
						 neighbor_subface_no);

	      fe_subface_values_neighbor
		.get_function_values (old_solution,
				      old_solution_values_face_neighbor);

	      for (unsigned int q=0; q<n_face_q_points; ++q)
		neighbor_temperature[q] = old_solution_values_face_neighbor[q](dim+1);

	      for (unsigned int q=0; q<n_face_q_points; ++q)
		{
		  Tensor<1,dim> present_u_face;
		  for (unsigned int d=0; d<dim; ++d)
		    present_u_face[d] = present_solution_values_face[q](d);

		  const double normal_flux = present_u_face *
					     fe_face_values.normal_vector(q);

		  const bool is_outflow_q_point = (normal_flux >= 0);

		  for (unsigned int i=0; i<dofs_per_cell; ++i)
		    local_rhs(i) -= time_step *
				    normal_flux *
				    (is_outflow_q_point == true
				     ?
				     old_solution_values_face[q](dim+1)
				     :
				     neighbor_temperature[q]) *
				    fe_face_values[temperature].value (i,q) *
				    fe_face_values.JxW(q);
		}
	    }

      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        system_rhs(local_dof_indices[i]) += local_rhs(i);
    }
}




				 // @sect4{BoussinesqFlowProblem::solve}
template <int dim>
void BoussinesqFlowProblem<dim>::solve ()
{
  solution = old_solution;

				// Use the BlockMatrixArray structure
				// for extracting only the upper left
				// 2x2 blocks from the matrix that will
				// be used for the solution of the
				// blocked system.
  {
    GrowingVectorMemory<Vector<double> > simple_mem;
    BlockMatrixArray<double> stokes_submatrix(2, 2, simple_mem);

    stokes_submatrix.enter(system_matrix.block(0,0),0,0);
    stokes_submatrix.enter(system_matrix.block(0,1),0,1);
    stokes_submatrix.enter(system_matrix.block(1,0),1,0);

				// Define some temporary vectors
				// for the solution process.
				// TODO: Can we somehow avoid copying
				// these vectors back and forth? I.e.
				// accessing the block vectors in a
				// similar way as the matrix with the
				// BlockMatrixArray class?
    std::vector<unsigned int> block_sizes(2);
    block_sizes[0] = solution.block(0).size();
    block_sizes[1] = solution.block(1).size();

    BlockVector<double> up_rhs(block_sizes);
    BlockVector<double> up(block_sizes);

    up_rhs.block(0) = system_rhs.block(0);
    up_rhs.block(1) = system_rhs.block(1);

				// Set up inverse matrix for
				// pressure mass matrix
    InverseMatrix<SparseMatrix<double>,SparseILU<double> >
      mp_inverse (system_matrix.block(1,1), *Mp_preconditioner);

				// Set up block Schur preconditioner
    /*BlockSchurPreconditioner<typename InnerPreconditioner<dim>::type,
                           SparseILU<double> >
      preconditioner (system_matrix, mp_inverse, *A_preconditioner);*/
    BlockSchurPreconditioner<TrilinosAmgPreconditioner<dim>, SparseILU<double> >
      preconditioner (system_matrix, mp_inverse, *Amg_preconditioner);

				// Set up GMRES solver and
				// solve.
    SolverControl solver_control (system_matrix.m(),
				  1e-6*system_rhs.l2_norm());

    SolverGMRES<BlockVector<double> > gmres(solver_control,
			SolverGMRES<BlockVector<double> >::AdditionalData(100));

    gmres.solve(stokes_submatrix, up, up_rhs, preconditioner);

				// Produce a constistent solution field
    hanging_node_constraints.distribute (up);

    std::cout << "   "
              << solver_control.last_step()
              << " GMRES iterations for Stokes subsystem."
              << std::endl;
	      
    solution.block(0) = up.block(0);
    solution.block(1) = up.block(1);
  }
				   // for DGQ1 needs to be /15
  time_step = GridTools::minimal_cell_diameter(triangulation) /
              std::max (get_maximal_velocity(), .05) / 2;

  assemble_rhs_T ();
  {

    SolverControl solver_control (system_matrix.block(2,2).m(),
				  1e-8*system_rhs.block(2).l2_norm());
    SolverCG<>   cg (solver_control);
    PreconditionJacobi<> preconditioner;
    preconditioner.initialize (system_matrix.block(2,2));

    try
      {
	cg.solve (system_matrix.block(2,2), solution.block(2),
		  system_rhs.block(2), preconditioner);
      }
    catch (...)
      {
	abort ();
      }

				     // produce a consistent temperature field
    hanging_node_constraints.distribute (solution);

    std::cout << "   "
              << solver_control.last_step()
              << " CG iterations for temperature."
              << std::endl;
    std::cout << "   Max temperature: "
	      << *std::max_element (solution.block(2).begin(),
				    solution.block(2).end())
	      << std::endl;
  }
}



				 // @sect4{BoussinesqFlowProblem::output_results}
template <int dim>
void BoussinesqFlowProblem<dim>::output_results ()  const
{
  if (timestep_number % 10 != 0)
    return;

  std::vector<std::string> solution_names (dim, "velocity");
  solution_names.push_back ("p");
  solution_names.push_back ("T");

  DataOut<dim> data_out;

  data_out.attach_dof_handler (dof_handler);

  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation
    (dim+2, DataComponentInterpretation::component_is_scalar);
  for (unsigned int i=0; i<dim; ++i)
    data_component_interpretation[i]
      = DataComponentInterpretation::component_is_part_of_vector;

  data_out.add_data_vector (solution, solution_names,
			    DataOut<dim>::type_dof_data,
			    data_component_interpretation);

  data_out.build_patches ();

  std::ostringstream filename;
  filename << "solution-" << Utilities::int_to_string(timestep_number, 4) << ".vtk";

  std::ofstream output (filename.str().c_str());
  data_out.write_vtk (output);
}



				 // @sect4{BoussinesqFlowProblem::refine_mesh}
template <int dim>
void BoussinesqFlowProblem<dim>::refine_mesh ()
{
  Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

//TODO do this better
  DerivativeApproximation::approximate_gradient (dof_handler,
						 old_solution,
						 estimated_error_per_cell,
						 dim+1);

  typename Triangulation<dim>::active_cell_iterator
    cell = triangulation.begin_active(),
    endc = triangulation.end();
  for (unsigned int cell_index=0; cell!=endc; ++cell, ++cell_index)
    estimated_error_per_cell(cell_index) *= cell->diameter();

  GridRefinement::refine_and_coarsen_fixed_fraction (triangulation,
						     estimated_error_per_cell,
						     0.3, 0.03,
						     static_cast<unsigned int>
						     (triangulation.n_active_cells()*1.1));

  SolutionTransfer<dim, double> soltrans(dof_handler);

  triangulation.prepare_coarsening_and_refinement();

  Vector<double> x_old_solution (dof_handler.n_dofs());
  x_old_solution = old_solution;

  soltrans.prepare_for_coarsening_and_refinement(x_old_solution);

  triangulation.execute_coarsening_and_refinement ();
  setup_dofs (true);

  Vector<double> tmp (dof_handler.n_dofs());
  soltrans.interpolate(x_old_solution, tmp);

  rebuild_matrices       = true;
  rebuild_preconditioner = true;

  old_solution = tmp;
}



				 // @sect4{BoussinesqFlowProblem::get_maximal_velocity}
template <int dim>
double BoussinesqFlowProblem<dim>::get_maximal_velocity () const
{
  QGauss<dim>   quadrature_formula(degree+2);
  const unsigned int   n_q_points
    = quadrature_formula.size();

  FEValues<dim> fe_values (fe, quadrature_formula,
                           update_values);
  std::vector<Vector<double> > solution_values(n_q_points,
                                               Vector<double>(dim+2));
  double max_velocity = 0;

  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      fe_values.get_function_values (solution, solution_values);

      for (unsigned int q=0; q<n_q_points; ++q)
        {
          Tensor<1,dim> velocity;
          for (unsigned int i=0; i<dim; ++i)
            velocity[i] = solution_values[q](i);

          max_velocity = std::max (max_velocity,
                                   velocity.norm());
        }
    }

  return max_velocity;
}



				 // @sect4{BoussinesqFlowProblem::run}
template <int dim>
void BoussinesqFlowProblem<dim>::run ()
{
  switch (dim)
    {
      case 2:
      {
// 	GridGenerator::hyper_ball (triangulation);

// 	static const HyperBallBoundary<dim> boundary;
// 	triangulation.set_boundary (0, boundary);

	GridGenerator::hyper_cube (triangulation);

	triangulation.refine_global (6); 

	break;
      }

      case 3:
      {
	GridGenerator::hyper_shell (triangulation,
				    Point<dim>(), 0.5, 1.0);

	static HyperShellBoundary<dim> boundary;
	triangulation.set_boundary (0, boundary);

	triangulation.refine_global (2);

	break;
      }

      default:
	    Assert (false, ExcNotImplemented());
    }


  const bool do_adaptivity = false;

  if (do_adaptivity)
    {
      setup_dofs(false);

      VectorTools::project (dof_handler,
			    hanging_node_constraints,
			    QGauss<dim>(degree+2),
			    InitialValues<dim>(),
			    old_solution);

      for (unsigned int pre_refinement=0; pre_refinement<4-dim; ++pre_refinement)
	{
	  refine_mesh ();

	  VectorTools::project (dof_handler,
				hanging_node_constraints,
				QGauss<dim>(degree+2),
				InitialValues<dim>(),
				old_solution);
	}
    }
  else
    {
      setup_dofs(true);

      VectorTools::project (dof_handler,
			    hanging_node_constraints,
			    QGauss<dim>(degree+2),
			    InitialValues<dim>(),
			    old_solution);
    }

  timestep_number = 0;
  double time = 0;

  do
    {
      std::cout << "Timestep " << timestep_number
		<< ":  t=" << time
		<< ", dt=" << time_step
                << std::endl;

      std::cout << "   Assembling..." << std::endl;
      assemble_system ();

      std::cout << "   Solving..." << std::endl;
      solve ();

      output_results ();

      time += time_step;
      ++timestep_number;

      old_solution = solution;

      std::cout << std::endl;

      if (do_adaptivity)
	if (timestep_number % 10 == 0)
	  refine_mesh ();
    }
  while (time <= 50);
}



				 // @sect3{The <code>main</code> function}
int main ()
{
  try
    {
      deallog.depth_console (0);

      BoussinesqFlowProblem<2> flow_problem(1);
      flow_problem.run ();
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
    }

  return 0;
}
