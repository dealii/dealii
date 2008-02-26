/* $Id: step-31.cc 15679 2008-01-24 23:28:37Z bangerth $ */
/* Author: Wolfgang Bangerth, Texas A&M University, 2008 */

/*    $Id: step-31.cc 15679 2008-01-24 23:28:37Z bangerth $       */
/*    Version: $Name$                                          */
/*                                                                */
/*    Copyright (C) 2008 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */


				 // @sect3{Include files}
                        
				 // As usual, we start by including 
				 // some well-known files:
#include <base/quadrature_lib.h>
#include <base/logstream.h>
#include <base/function.h>
#include <base/utilities.h>

#include <lac/block_vector.h>
#include <lac/full_matrix.h>
#include <lac/block_sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>

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
#include <fe/fe_system.h>
#include <fe/fe_values.h>
#include <fe/mapping_q1.h>

#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <numerics/data_out.h>
#include <numerics/error_estimator.h>

				 // Then we need to include the header file
				 // for the sparse direct solver UMFPACK:
#include <lac/sparse_direct.h>

				 // This includes the libary for the 
				 // incomplete LU factorization that will 
				 // be used as a preconditioner in 3D:
#include <lac/sparse_ilu.h>

				 // This is C++:
#include <fstream>
#include <sstream>

				 // As in all programs, the namespace dealii
				 // in included:
using namespace dealii;

				 // @sect3{Defining the inner preconditioner type}
                        
				 // As explained in the introduction, we are
				 // going to use different preconditioners for
				 // two and three space dimensions,
				 // respectively. We differentiate between
				 // them by the use of the spatial dimension
				 // as a template parameter. See step-4 for
				 // details on templates.  We are not going to
				 // create any preconditioner object here, all
				 // we do is to create class that holds a
				 // local typedef determining the
				 // preconditioner class so we can write our
				 // program in a dimension-independent way.
template <int dim>
struct InnerPreconditioner;

				 // In 2D, we are going to use a sparse direct
				 // solver as preconditioner:
template <>
struct InnerPreconditioner<2> 
{
    typedef SparseDirectUMFPACK type;
};

				 // And the ILU preconditioning in 3D, called
				 // by SparseILU:
template <>
struct InnerPreconditioner<3> 
{
    typedef SparseILU<double> type;
};


				 // @sect3{The <code>StokesProblem</code> class template}
                    
				 // This is an adaptation of step-20, so the
				 // main class and the data types are the same
				 // as used there. In this example we also use
				 // adaptive grid refinement, which is handled
				 // in complete analogy to step-6:
template <int dim>
class StokesProblem 
{
  public:
    StokesProblem (const unsigned int degree);
    void run ();
    
  private:
    void setup_dofs ();
    void assemble_system ();
    void solve ();
    void output_results (const unsigned int refinement_cycle) const;
    void refine_mesh ();
    
    const unsigned int   degree;
    
    Triangulation<dim>   triangulation;
    FESystem<dim>        fe;
    DoFHandler<dim>      dof_handler;

    ConstraintMatrix     hanging_node_constraints;
    
    BlockSparsityPattern      sparsity_pattern;
    BlockSparseMatrix<double> system_matrix;

    BlockVector<double> solution;
    BlockVector<double> system_rhs;

				     // This one is new: We shall use a
				     // so-called shared pointer structure to
				     // access the preconditioner. Shared
				     // pointers are essentially just a
				     // convenient form of pointers. Several
				     // shared pointers can point to the same
				     // object (just like regular pointers),
				     // but when the last shared pointer
				     // object to point to a preconditioner
				     // object is deleted (for example if a
				     // shared pointer object goes out of
				     // scope, if the class of which it is a
				     // member is destroyed, or if the pointer
				     // is assigned a different preconditioner
				     // object) then the preconditioner object
				     // pointed to is also destroyed. This
				     // ensures that we don't have to manually
				     // track in how many places a
				     // preconditioner object is still
				     // referenced, it can never create a
				     // memory leak, and can never produce a
				     // dangling pointer to an already
				     // destroyed object:
    boost::shared_ptr<typename InnerPreconditioner<dim>::type> A_preconditioner;
};

				 // @sect3{%Boundary values and right hand side}

				 // As in step-20 and most other example
				 // programs, the next task is to define the
				 // data for the PDE: For the
				 // Stokes problem, we are going to use
				 // natural boundary values at some portion
				 // of the boundary (Neumann-type), and
				 // boundary conditions on the velocity
				 // (Dirichlet type) on the rest of the
				 // boundary.  The pressure boundary condition
				 // is scalar, and so is the respective
				 // function, whereas the Dirichlet (velocity)
				 // condition is vector-valued. Due to the
				 // structure of deal.II's libraries, we have
				 // to define the function on the (u,p)-space,
				 // but we are going to filter out the
				 // pressure component when condensating the
				 // Dirichlet data in
				 // <code>assemble_system</code>.
                    
				 // Given the problem described in the 
				 // introduction, we know which values to 
				 // set for the respective functions.
template <int dim>
class BoundaryValues : public Function<dim> 
{
  public:
    BoundaryValues () : Function<dim>(dim+1) {}
    
    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;

    virtual void vector_value (const Point<dim> &p, 
                               Vector<double>   &value) const;
};


template <int dim>
double
BoundaryValues<dim>::value (const Point<dim>  &p,
			    const unsigned int component) const 
{
  if (component == 0)
    return (p[0] < 0 ? -1 : (p[0] > 0 ? 1 : 0));
  return 0;
}


template <int dim>
void
BoundaryValues<dim>::vector_value (const Point<dim> &p,
				   Vector<double>   &values) const 
{
  for (unsigned int c=0; c<this->n_components; ++c)
    values(c) = BoundaryValues<dim>::value (p, c);
}



				 // We implement similar functions
				 // for the right hand side.
template <int dim>
class RightHandSide : public Function<dim> 
{
  public:
    RightHandSide () : Function<dim>(dim+1) {}
    
    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;

    virtual void vector_value (const Point<dim> &p, 
                               Vector<double>   &value) const;
    
};


template <int dim>
double
RightHandSide<dim>::value (const Point<dim>  &/*p*/,
                           const unsigned int /*component*/) const 
{
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
                        
				 // The linear solvers and preconditioners are
				 // discussed extensively in the
				 // introduction. Here, we create the
				 // respective objects that will be used.
                        
				 // @sect4{The <code>InverseMatrix</code> class template}
                        
				 // This is going to represent the data
				 // structure for an inverse matrix. This
				 // class is derived from the one in
				 // step-20. The only difference is that we
				 // now do include a preconditioner to the
				 // matrix.  This is going to happen via a
				 // template parameter <code>class
				 // Preconditioner</code>, so the
				 // preconditioner type will be set when an
				 // <code>InverseMatrix</code> object is
				 // created. The member function
				 // <code>vmult</code> is, as in step-20, a
				 // multiplication with a vector, obtained by
				 // solving a linear system.
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

    mutable GrowingVectorMemory<> vector_memory;    
};


template <class Matrix, class Preconditioner>
InverseMatrix<Matrix,Preconditioner>::InverseMatrix (const Matrix &m,
						     const Preconditioner &preconditioner)
		:
		matrix (&m),
		preconditioner (preconditioner)
{}


				 // This is the implementation of the 
				 // <code>vmult</code> function. We note 
				 // two things: 
                    
				 // Firstly, we use a rather large tolerance
				 // for the solver control. The reason for
				 // this is that the function is used very
				 // frequently, and hence, any additional
				 // effort to make the residual in the CG
				 // solve smaller makes the solution more
				 // expensive. Note that we do not only use
				 // this class as a preconditioner for the
				 // Schur complement, but also when forming
				 // the inverse of the Laplace matrix - which
				 // has to be accurate in order to obtain a
				 // solution to the right problem.
                    
				 // Secondly, we catch exceptions from the
				 // solver at this stage. While this is not of
				 // greater interest our general setting with
				 // the requirement of accurate inverses (and
				 // we indeed abort the program when any
				 // exception occurs), the situation would
				 // change if an object of the class
				 // <code>InverseMatrix</code> is only used
				 // for preconditioning. In such a setting,
				 // one could imagine to use a few CG sweeps
				 // as a preconditioner - which is done
				 // e.g. for mass matrices, see the results
				 // section below. Using <code>catch
				 // (SolverControl::NoConvergence) {}</code>
				 // in conjunction with only a few iterations,
				 // say 10, would result in that effect - the
				 // program would continue to run even though
				 // the solver has not converged.  Note,
				 // though, that applying the CG method is not
				 // a linear operation (see the actual CG
				 // algorithm for details on that), so
				 // unconverged preconditioners are to be used
				 // with care in order to not yield a wrong
				 // solution.
template <class Matrix, class Preconditioner>
void InverseMatrix<Matrix,Preconditioner>::vmult (Vector<double>       &dst,
						  const Vector<double> &src) const
{
  SolverControl solver_control (src.size(), 1e-6*src.l2_norm());
  SolverCG<> cg (solver_control, vector_memory);

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


				 // @sect4{The <code>SchurComplement</code> class template}

				 // This class implements the Schur complement
				 // discussed in the introduction.  It is in
				 // analogy to step-20.  Though, we now call
				 // it with a template parameter
				 // <code>Preconditioner</code> in order to
				 // access that when specifying the respective
				 // type of the inverse matrix class. As a
				 // consequence of the definition above, the
				 // declaration <code>InverseMatrix</code> now
				 // contains the second template parameter
				 // from preconditioning as above, which
				 // effects the <code>SmartPointer</code>
				 // object <code>m_inverse</code> as well.
template <class Preconditioner>
class SchurComplement : public Subscriptor
{
  public:
    SchurComplement (const BlockSparseMatrix<double> &A,
		     const InverseMatrix<SparseMatrix<double>, Preconditioner> &Minv);

    void vmult (Vector<double>       &dst,
		const Vector<double> &src) const;

  private:
    const SmartPointer<const BlockSparseMatrix<double> > system_matrix;
    const SmartPointer<const InverseMatrix<SparseMatrix<double>, Preconditioner> > m_inverse;
    
    mutable Vector<double> tmp1, tmp2;
};



template <class Preconditioner>
SchurComplement<Preconditioner>::
SchurComplement (const BlockSparseMatrix<double> &A,
		 const InverseMatrix<SparseMatrix<double>,Preconditioner> &Minv)
		:
		system_matrix (&A),
		m_inverse (&Minv),
		tmp1 (A.block(0,0).m()),
		tmp2 (A.block(0,0).m())
{}


template <class Preconditioner>
void SchurComplement<Preconditioner>::vmult (Vector<double>       &dst,
					     const Vector<double> &src) const
{
  system_matrix->block(0,1).vmult (tmp1, src);
  m_inverse->vmult (tmp2, tmp1);
  system_matrix->block(1,0).vmult (dst, tmp2);
}


				 // @sect3{StokesProblem class implementation}
                        
				 // @sect4{StokesProblem::StokesProblem}

				 // The constructor of this class looks very
				 // similar to the one of step-20. The
				 // constructor initializes the variables for
				 // the polynomial degree, triangulation,
				 // finite element system and the dof
				 // handler. The underlying polynomial
				 // functions are of order
				 // <code>degree+1</code> for the
				 // vector-valued velocity components and of
				 // order <code>degree</code> in pressure.
				 // This gives the LBB-stable element pair
				 // Q(degree+1)Q(degree).
				 //                    
				 // Note that we initialize the triangulation
				 // with a MeshSmoothing argument, which
				 // ensures that the refinement of cells is
				 // done in a way that the approximation of
				 // the PDE solution remains well-behaved
				 // (problems arise if grids are too
				 // unstructered), see the documentation of
				 // <code>Triangulation::MeshSmoothing</code>
				 // for details.
template <int dim>
StokesProblem<dim>::StokesProblem (const unsigned int degree)
                :
                degree (degree),
                triangulation (Triangulation<dim>::maximum_smoothing),
                fe (FE_Q<dim>(degree+1), dim,
                    FE_Q<dim>(degree), 1),
                dof_handler (triangulation)
{}


				 // @sect4{StokesProblem::setup_dofs}
                        
				 // Given a mesh, this function associates
				 // the degrees of freedom with it and
				 // creates the corresponding matrices and
				 // vectors.
template <int dim>
void StokesProblem<dim>::setup_dofs ()
{
				   // Release preconditioner from
				   // previous steps since it
				   // will definitely not be needed
				   // any more after this point.
  A_preconditioner.reset ();
  
  dof_handler.distribute_dofs (fe); 
  
				   // In order to make the ILU preconditioner
				   // (in 3D) to work efficiently, the dofs
				   // are renumbered using the Cuthill-McKee
				   // algorithm. Though, the block structure
				   // of velocity and pressure shall be as in
				   // step-20. This is done in two
				   // steps. First, all dofs are renumbered by
				   // <code>DoFRenumbering::Cuthill_McKee</code>,
				   // and then we renumber once again by
				   // components. Since
				   // <code>DoFRenumbering::component_wise</code>
				   // does not touch the renumbering within
				   // the individual blocks, the basic
				   // renumbering from Cuthill-McKee remains.
  DoFRenumbering::Cuthill_McKee (dof_handler);

				   // There is one more change: There is no
				   // reason in creating <code>dim</code>
				   // blocks for the velocity components, so
				   // they can all be grouped in only one
				   // block. The vector
				   // <code>block_component</code> does
				   // precisely this: velocity values
				   // correspond to block 0, and pressure
				   // values will sit in block 1.
  std::vector<unsigned int> block_component (dim+1,0);
  block_component[dim] = 1;
  DoFRenumbering::component_wise (dof_handler, block_component);

				   // Since we use adaptively refined grids
				   // the constraint matrix for hanging node
				   // constraints is generated from the dof
				   // handler.
  hanging_node_constraints.clear ();
  DoFTools::make_hanging_node_constraints (dof_handler,
                                           hanging_node_constraints);
  hanging_node_constraints.close ();

				   // In analogy to step-20, we count the dofs
				   // in the individual components.  We could
				   // do this in the same way as there, but we
				   // want to operate on the block structure
				   // we used already for the renumbering: The
				   // function
				   // <code>DoFTools::count_dofs_per_block</code>
				   // does the same as
				   // <code>DoFTools::count_dofs_per_component</code>,
				   // but now grouped as velocity and pressure
				   // block via <code>block_component</code>.
  std::vector<unsigned int> dofs_per_block (2);
  DoFTools::count_dofs_per_block (dof_handler, dofs_per_block, block_component);  
  const unsigned int n_u = dofs_per_block[0],
                     n_p = dofs_per_block[1];

  std::cout << "   Number of active cells: "
            << triangulation.n_active_cells()
            << std::endl
            << "   Number of degrees of freedom: "
            << dof_handler.n_dofs()
            << " (" << n_u << '+' << n_p << ')'
            << std::endl;

				   // Release the memory previously attached
				   // to the system matrix and untie it from
				   // the old sparsity pattern prior to
				   // generating the current data structure.
  system_matrix.clear ();
      
				   // The next task is to allocate a sparsity
				   // pattern for the system matrix we will
				   // create. We could do this in the same way
				   // as in step-20, though, there is a major
				   // reason not to do so. In 3D, the function
				   // <code>DoFTools::max_couplings_between_dofs</code>
				   // yields a very large number for the
				   // coupling between the individual dofs, so
				   // that the memory initially provided for
				   // the creation of the sparsity pattern of
				   // the matrix is far too much - so much
				   // actually that it won't even fit into the
				   // physical memory of most systems already
				   // for moderately-sized 3D problems. See
				   // also the discussion in step-18.
				   // Instead, we use a temporary object of
				   // the class
				   // BlockCompressedSparsityPattern, which is
				   // a block version of the compressed
				   // sparsity patterns from step-11 and
				   // step-18. All this is done inside a new
				   // scope, which means that the memory of
				   // <code>csp</code> will be released once
				   // the information has been copied to
				   // <code>sparsity_pattern</code>.
  {
    BlockCompressedSparsityPattern csp;

    csp.reinit (2,2);
    csp.block(0,0).reinit (n_u, n_u);
    csp.block(1,0).reinit (n_p, n_u);
    csp.block(0,1).reinit (n_u, n_p);
    csp.block(1,1).reinit (n_p, n_p);
  
    csp.collect_sizes();    
  
    DoFTools::make_sparsity_pattern (dof_handler, csp);
    hanging_node_constraints.condense (csp);
    sparsity_pattern.copy_from (csp);
  }
  
				   // Finally, the system matrix,
				   // solution and right hand side are 
				   // created from the block
				   // structure as in step-20.
  system_matrix.reinit (sparsity_pattern);
                                   
  solution.reinit (2);
  solution.block(0).reinit (n_u);
  solution.block(1).reinit (n_p);
  solution.collect_sizes ();
  
  system_rhs.reinit (2);
  system_rhs.block(0).reinit (n_u);
  system_rhs.block(1).reinit (n_p);
  system_rhs.collect_sizes ();
}


				 // @sect4{StokesProblem::assemble_system}
                        
				 // The assembly process follows the
				 // discussion in step-20 and in the
				 // introduction. We use the well-known
				 // abbreviations for the data structures
				 // that hold the local matrix, right
				 // hand side, and global
				 // numbers of the degrees of freedom
				 // for the present cell.
template <int dim>
void StokesProblem<dim>::assemble_system () 
{
  system_matrix=0;
  system_rhs=0;
  
  QGauss<dim>   quadrature_formula(degree+2); 
  QGauss<dim-1> face_quadrature_formula(degree+2);

  FEValues<dim> fe_values (fe, quadrature_formula,
                           update_values    |
                           update_quadrature_points  |
                           update_JxW_values |
                           update_gradients);
  FEFaceValues<dim> fe_face_values (fe, face_quadrature_formula, 
                                    update_values    | 
                                    update_normal_vectors |
                                    update_quadrature_points  |
                                    update_JxW_values);

  const unsigned int   dofs_per_cell   = fe.dofs_per_cell;
  
  const unsigned int   n_q_points      = quadrature_formula.size();

  FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       local_rhs (dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);
  
  const RightHandSide<dim>          right_hand_side;
  std::vector<Vector<double> >      rhs_values (n_q_points,
                                                Vector<double>(dim+1));

				   // Next, we need two objects that work as
				   // extractors for the FEValues
				   // object. Their use is explained in detail
				   // in the report on @ref vector_valued :
  const FEValuesExtractors::Vector velocities (0);
  const FEValuesExtractors::Scalar pressure (dim);

				   // This starts the loop over all
				   // cells. With the abbreviations
				   // <code>extract_u</code> etc. 
				   // introduced above, it is 
				   // evident what is going on.
  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    { 
      fe_values.reinit (cell);
      local_matrix = 0;
      local_rhs = 0;
      
      right_hand_side.vector_value_list(fe_values.get_quadrature_points(),
                                        rhs_values);
      
      for (unsigned int q=0; q<n_q_points; ++q)
	{
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    {
	      const SymmetricTensor<2,dim>
		phi_i_grads_u = fe_values[velocities].symmetric_gradient (i, q);
	      const Tensor<1,dim> phi_i_u       = fe_values[velocities].value (i, q);
	      const double        div_phi_i_u   = fe_values[velocities].divergence (i, q);
	      const double        phi_i_p       = fe_values[pressure].value (i, q);
          
	      for (unsigned int j=0; j<dofs_per_cell; ++j)
		{
		  const SymmetricTensor<2,dim>
		    phi_j_grads_u = fe_values[velocities].symmetric_gradient (j, q);
		  const double        div_phi_j_u   = fe_values[velocities].divergence (j, q);
		  const double        phi_j_p       = fe_values[pressure].value (j, q);
						   // Note the way we write
						   // the contributions <code>
						   // phi_i_p * phi_j_p
						   // </code>, yielding a
						   // pressure mass matrix,
						   // into the same data
						   // structure as the terms
						   // for the actual Stokes
						   // system - in accordance
						   // with the description in
						   // the introduction.  They
						   // won't be mixed up, since
						   // <code>phi_i_p *
						   // phi_j_p</code> is only
						   // non-zero when all the
						   // other terms vanish (and
						   // the other way around).
						   //
						   // Note also that operator*
						   // is overloaded for
						   // symmetric tensors,
						   // yielding the scalar
						   // product between the two
						   // tensors in the first
						   // line:
		  local_matrix(i,j) += (phi_i_grads_u * phi_j_grads_u
					- div_phi_i_u * phi_j_p
					- phi_i_p * div_phi_j_u
					+ phi_i_p * phi_j_p)
				       * fe_values.JxW(q);     

		}
	      const unsigned int component_i =
		fe.system_to_component_index(i).first;
	      local_rhs(i) += fe_values.shape_value(i,q) * 
			      rhs_values[q](component_i) *
			      fe_values.JxW(q);
	    }
	}      

				       // The final step is, as usual, the
				       // transfer of the local contributions
				       // to the global system matrix. This
				       // works also in the case of block
				       // vectors and matrices, and also the
				       // terms constituting the pressure mass
				       // matrix are written into the correct
				       // position without any further
				       // interaction.
      cell->get_dof_indices (local_dof_indices);

      for (unsigned int i=0; i<dofs_per_cell; ++i)
        for (unsigned int j=0; j<dofs_per_cell; ++j)
          system_matrix.add (local_dof_indices[i],
                             local_dof_indices[j],
                             local_matrix(i,j));
      
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        system_rhs(local_dof_indices[i]) += local_rhs(i);
    }

				   // After the addition of the local
				   // contributions, we have to condense the
				   // hanging node constraints and interpolate
				   // Dirichlet boundary conditions.  Note
				   // that Dirichlet boundary conditions are
				   // only condensed in boundary points that
				   // are labeled with "1", indicating that
				   // Dirichlet data is to be set.  There is
				   // one more thing, though.  The function
				   // describing the Dirichlet conditions was
				   // defined for all components, both
				   // velocity and pressure. However, the
				   // Dirichlet conditions are to be set for
				   // the velocity only.  To this end, we use
				   // a <code>component_mask</code> that
				   // filters away the pressure component, so
				   // that the condensation is performed on
				   // velocity dofs.
  hanging_node_constraints.condense (system_matrix);
  hanging_node_constraints.condense (system_rhs);  

  {
    std::map<unsigned int,double> boundary_values;
    std::vector<bool> component_mask (dim+1, true);
    component_mask[dim] = false;
    VectorTools::interpolate_boundary_values (dof_handler,
                                              1,
                                              BoundaryValues<dim>(),
                                              boundary_values,
                                              component_mask);

    MatrixTools::apply_boundary_values (boundary_values,
                                        system_matrix,
                                        solution,
                                        system_rhs);
  }
  
				   // Before we're going to solve this linear
				   // system, we generate a preconditioner for
				   // the velocity-velocity matrix, i.e.,
				   // <code>block(0,0)</code> in the system
				   // matrix. As mentioned above, this depends
				   // on the spatial dimension. Since this
				   // handled automatically by the template
				   // <code>dim</code> in
				   // <code>InnerPreconditioner</code>, we
				   // don't have to manually intervene at this
				   // point any further.
  std::cout << "   Computing preconditioner..." << std::endl << std::flush;
      
  A_preconditioner
    = boost::shared_ptr<typename InnerPreconditioner<dim>::type>(new typename InnerPreconditioner<dim>::type());
  A_preconditioner->initialize (system_matrix.block(0,0),
				typename InnerPreconditioner<dim>::type::AdditionalData());

}



				 // @sect4{StokesProblem::solve}
                        
				 // After the discussion in the introduction
				 // and the definition of the respective
				 // classes above, the implementation of the
				 // <code>solve</code> function is rather
				 // straigt-forward and done in a similar way
				 // as in step-20. To start with, we need an
				 // object of the <code>InverseMatrix</code>
				 // class that represents the inverse of the
				 // matrix A. As described in the
				 // introduction, the inverse is generated
				 // with the help of an inner preconditioner
				 // of type
				 // <code>InnerPreconditioner<dim></code>.
template <int dim>
void StokesProblem<dim>::solve () 
{
  const InverseMatrix<SparseMatrix<double>,
                      typename InnerPreconditioner<dim>::type>
    A_inverse (system_matrix.block(0,0), *A_preconditioner);
  Vector<double> tmp (solution.block(0).size());
  
				   // This is as in step-20. We generate the
				   // right hand side $B A^{-1} F - G$ for the
				   // Schur complement and an object that
				   // represents the respective linear
				   // operation $B A^{-1} B^T$, now with a
				   // template parameter indicating the
				   // preconditioner - in accordance with the
				   // definition of the class.
  {
    Vector<double> schur_rhs (solution.block(1).size());
    A_inverse.vmult (tmp, system_rhs.block(0));
    system_matrix.block(1,0).vmult (schur_rhs, tmp);
    schur_rhs -= system_rhs.block(1);
  
    SchurComplement<typename InnerPreconditioner<dim>::type>
      schur_complement (system_matrix, A_inverse);
    
				     // The usual control structures for
				     // the solver call are created...
    SolverControl solver_control (system_matrix.block(0,0).m(),
				  1e-6*schur_rhs.l2_norm());
    SolverCG<>    cg (solver_control);
    
				     // Now to the preconditioner to the Schur
				     // complement. As explained in the
				     // introduction, the preconditioning is
				     // done by a mass matrix in the pressure
				     // variable.  It is stored in the $(1,1)$
				     // block of the system matrix (that is
				     // not used elsewhere in this function).
				     //
				     // Actually, the solver needs to have the
				     // preconditioner in the form $P^{-1}$, so
				     // we need to create an inverse
				     // operation. Once again, we use an
				     // object of the class
				     // <code>InverseMatrix</code>, which
				     // implements the <code>vmult</code>
				     // operation that is needed by the
				     // solver.  In this case, we have to
				     // invert the pressure mass matrix. As it
				     // already turned out in earlier tutorial
				     // programs, the inversion of a mass
				     // matrix is a rather cheap and
				     // straight-forward operation (compared
				     // to, e.g., a Laplace matrix). The CG
				     // method with SSOR preconditioning
				     // converges in 10-20 steps,
				     // independently on the mesh size.  This
				     // is precisely what we do here: We
				     // choose an SSOR preconditioner with
				     // parameter 1.2 and take it along to the
				     // InverseMatrix object via the
				     // corresponding template parameter.  A
				     // CG solver is then called within the
				     // vmult operation of the inverse matrix.
    PreconditionSSOR<> preconditioner;
    preconditioner.initialize (system_matrix.block(1,1), 1.2);
  
    InverseMatrix<SparseMatrix<double>,PreconditionSSOR<> >
      m_inverse (system_matrix.block(1,1), preconditioner);
    
				     // With the Schur complement and an
				     // efficient preconditioner at hand,
				     // we can solve the respective
				     // equation in the usual way.
    cg.solve (schur_complement, solution.block(1), schur_rhs,
	      m_inverse);
  
				     // After this first solution step,
				     // the hanging node constraints have
				     // to be distributed to the solution -
				     // in order to achieve a consistent 
				     // pressure field.
    hanging_node_constraints.distribute (solution);
  
    std::cout << "  "
	      << solver_control.last_step()
	      << " outer CG Schur complement iterations for pressure"
	      << std::flush
	      << std::endl;    
  }
    
				   // As in step-20, we finally need to solve
				   // for the velocity equation where we plug
				   // in the the solution to the pressure
				   // equation. This involves only objects we
				   // already know - so we simply multiply p
				   // by $B^T$, subtract the right hand side and
				   // multiply by the inverse of A.
  {
    system_matrix.block(0,1).vmult (tmp, solution.block(1));
    tmp *= -1;
    tmp += system_rhs.block(0);
  
    A_inverse.vmult (solution.block(0), tmp);
  
				     // Again, we need to distribute the
				     // constraints from hanging nodes in
				     // order to obtain a constistent flow
				     // field.
    hanging_node_constraints.distribute (solution);
  }
}


				 // @sect4{StokesProblem::output_results}
                        
				 // The next function generates graphical
				 // output. In this example, we are going to
				 // use the VTK file format.  We attach names
				 // to the individual variables in the problem
				 // - <code>velocity</code> to the dim
				 // components of velocity and <code>p</code>
				 // to the pressure.  In order to tell the VTK
				 // file which components are vectors and
				 // which scalars, we need to add that
				 // information as well - achieved by the
				 // <code>DataComponentInterpretation</code>
				 // class.  The rest of the function is then
				 // the same as in step-20.
template <int dim>
void
StokesProblem<dim>::output_results (const unsigned int refinement_cycle)  const
{
  std::vector<std::string> solution_names (dim, "velocity");
  solution_names.push_back ("p");
  
  DataOut<dim> data_out;          

  data_out.attach_dof_handler (dof_handler);

  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation
    (dim+1, DataComponentInterpretation::component_is_scalar);
  for (unsigned int i=0; i<dim; ++i)
    data_component_interpretation[i]
      = DataComponentInterpretation::component_is_part_of_vector;
  
  data_out.add_data_vector (solution, solution_names,
			    DataOut<dim>::type_dof_data,
			    data_component_interpretation);
  
  data_out.build_patches ();
  
  std::ostringstream filename;
  filename << "solution-"
           << Utilities::int_to_string (refinement_cycle, 2)
           << ".vtk";

  std::ofstream output (filename.str().c_str());
  data_out.write_vtk (output);
}


				 // @sect4{StokesProblem::refine_mesh}
                        
				 // This is the last interesting function
				 // of the <code>StokesProblem</code> class.
				 // As indicated by its name, it takes the
				 // solution to the problem and
				 // refines the mesh where this is
				 // needed. The procedure is the same
				 // as in the respective step in
				 // step-6, with the exception that
				 // we base the refinement only on the
				 // change in pressure, i.e., we call
				 // the Kelly error estimator with a
				 // mask object. Additionally, we do
				 // not coarsen the grid again.
template <int dim>
void
StokesProblem<dim>::refine_mesh () 
{
  Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

  std::vector<bool> component_mask (dim+1, false);
  component_mask[dim] = true;
  KellyErrorEstimator<dim>::estimate (dof_handler,
                                      QGauss<dim-1>(degree+1),
                                      typename FunctionMap<dim>::type(),
                                      solution,
                                      estimated_error_per_cell,
                                      component_mask);

  GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                   estimated_error_per_cell,
                                                   0.3, 0.0);
  triangulation.execute_coarsening_and_refinement ();
}


				 // @sect4{StokesProblem::run}
                        
				 // The last step in the Stokes class
				 // is, as usual, the program that generates
				 // the initial grid and calls the other
				 // functions in the respective order.
template <int dim>
void StokesProblem<dim>::run () 
{
				   // We start off with a rectangle of size $4
				   // \times 1$ (in 2d) or $4 \times 1 times
				   // 1$ (in 3d), placed in $R^2/R^3$ as
				   // $(-2,2)times(-1,0)$ or
				   // $(-2,2)\times(0,1)\times(-1,1)$,
				   // respectively. It is natural to start
				   // with equal mesh size in each direction,
				   // so we subdivide the initial rectangle
				   // four times in the first coordinate
				   // direction.
  std::vector<unsigned int> subdivisions (dim, 1);
  subdivisions[0] = 4;
    
  GridGenerator::subdivided_hyper_rectangle (triangulation,
                                             subdivisions,
                                             (dim == 2 ?
                                              Point<dim>(-2,-1) :
                                              Point<dim>(-2,0,-1)),
                                             (dim == 2 ?
                                              Point<dim>(2,0) :
                                              Point<dim>(2,1,0)));
  
				   // A boundary indicator is set to all
				   // boundaries that are subject to Dirichlet
				   // boundary conditions, i.e.  to faces that
				   // are located at 0 in the last coordinate
				   // direction. See the example description
				   // above for details.
  for (typename Triangulation<dim>::active_cell_iterator
	 cell = triangulation.begin_active();
       cell != triangulation.end(); ++cell)
    for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
      if (cell->face(f)->center()[dim-1] == 0)
	{
	  cell->face(f)->set_boundary_indicator(1);
					   /*
					     for (unsigned int e=0; e<GeometryInfo<dim>::lines_per_face; ++e)
					     cell->face(f)->line(e)->set_boundary_indicator (1);
					   */
	}
  
  
				   // We employ an initial refinement before
				   // solving for the first time. In 3D, there
				   // are going to be more dofs, so we refine
				   // less there.
  triangulation.refine_global (4-dim);

				   // As first seen in step-6, we cycle over
				   // the different refinement levels and
				   // refine (if not the first step), setup
				   // the dofs and matrices, assemble, solve
				   // and create an output.
  for (unsigned int refinement_cycle = 0; refinement_cycle<7;
       ++refinement_cycle)
    {
      std::cout << "Refinement cycle " << refinement_cycle << std::endl;
      
      if (refinement_cycle > 0)
        refine_mesh ();
      
      setup_dofs ();

      std::cout << "   Assembling..." << std::endl << std::flush;
      assemble_system ();      

      std::cout << "   Solving..." << std::flush;
      solve ();
      
      output_results (refinement_cycle);

      std::cout << std::endl;
    }
}


				 // @sect3{The <code>main</code> function}

				 // The main function is the same as in
				 // step-20. We pass the element degree as a
				 // parameter and choose the space dimension
				 // at the well-known template slot.
int main () 
{
  try
    {
      deallog.depth_console (0);

      StokesProblem<3> flow_problem(1);
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
