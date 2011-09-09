/* $Id$ */
/* Author: Wolfgang Bangerth, Texas A&M University, 2008 */

/*    $Id$       */
/*                                                                */
/*    Copyright (C) 2008, 2009, 2010, 2011 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */


				 // @sect3{Include files}

				 // As usual, we start by including
				 // some well-known files:
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/numerics/vectors.h>
#include <deal.II/numerics/matrices.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>


#include <deal.II/lac/sparse_direct.h>

#include <deal.II/lac/sparse_ilu.h>

#include <deal.II/multigrid/multigrid.h>
#include <deal.II/multigrid/mg_dof_handler.h>
#include <deal.II/multigrid/mg_dof_accessor.h>
#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_matrix.h>

#include <fstream>
#include <sstream>


namespace Step42
{
  using namespace dealii;


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

  template <typename MATRIX>
  void copy(const MATRIX &matrix,
	    FullMatrix<double> &full_matrix)
  {
    const unsigned int m = matrix.m();
    const unsigned int n = matrix.n();
    full_matrix.reinit(n,m);

    Vector<double> unit (n);
    Vector<double> result (m);
    for(unsigned int i=0; i<n; ++i)
      {
	unit(i) = 1;
	for(unsigned int j=0; j<m; ++j)
	  {
	    matrix.vmult(result,unit);
	    full_matrix(i,j) = result(j);
	  }
	unit(i) = 0;
      }
  }

  template <int dim>
  class StokesProblem
  {
    public:
      StokesProblem (const unsigned int degree);
      void run ();

    private:
      void setup_dofs ();
      void assemble_system ();
      void assemble_multigrid ();
      void solve ();
      void solve_block ();

      void find_dofs_on_lower_level (std::vector<std::vector<bool> > &lower_dofs,
				     std::vector<std::vector<bool> > &boundary_dofs);

      void output_results (const unsigned int refinement_cycle) const;
      void refine_mesh ();

      const unsigned int   degree;

      Triangulation<dim>   triangulation;
      FESystem<dim>        fe;
      MGDoFHandler<dim>    dof_handler;

      ConstraintMatrix     constraints;

      BlockSparsityPattern      sparsity_pattern;
      BlockSparseMatrix<double> system_matrix;

      BlockVector<double> solution;
      BlockVector<double> system_rhs;

      MGLevelObject<ConstraintMatrix>           mg_constraints;
      MGLevelObject<BlockSparsityPattern>       mg_sparsity;
      MGLevelObject<BlockSparseMatrix<double> > mg_matrices;

      MGLevelObject<BlockSparseMatrix<double> > mg_interface_matrices;
      MGConstrainedDoFs                         mg_constrained_dofs;
      std::vector<std::vector<unsigned int> >   mg_dofs_per_component;

      std::vector<std_cxx1x::shared_ptr<typename InnerPreconditioner<dim>::type> > mg_A_preconditioner;
      std_cxx1x::shared_ptr<typename InnerPreconditioner<dim>::type> A_preconditioner;
  };


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
    Assert (component < this->n_components,
	    ExcIndexRange (component, 0, this->n_components));

    if (component == 0 && p[0] == 0)
      return (dim == 2 ? - p[1]*(p[1]-1.) : p[1]*(p[1]-1.) * p[2]*(p[2]-1.));
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
			     const unsigned int component) const
  {
    return (component == 1 ? 1 : 0);
  }


  template <int dim>
  void
  RightHandSide<dim>::vector_value (const Point<dim> &p,
				    Vector<double>   &values) const
  {
    for (unsigned int c=0; c<this->n_components; ++c)
      values(c) = RightHandSide<dim>::value (p, c);
  }




  template <class Matrix, class Preconditioner>
  class InverseMatrix : public Subscriptor
  {
    public:
      InverseMatrix (const Matrix         &m,
		     const Preconditioner &preconditioner);

      void vmult (Vector<double>       &dst,
		  const Vector<double> &src) const;

      mutable std::string name;
    private:
      const SmartPointer<const Matrix> matrix;
      const SmartPointer<const Preconditioner> preconditioner;
  };


  template <class Matrix, class Preconditioner>
  InverseMatrix<Matrix,Preconditioner>::InverseMatrix (const Matrix &m,
						       const Preconditioner &preconditioner)
		  :
		  matrix (&m),
		  preconditioner (&preconditioner)
  {}


  template <class Matrix, class Preconditioner>
  void InverseMatrix<Matrix,Preconditioner>::vmult (Vector<double>       &dst,
						    const Vector<double> &src) const
  {
    SolverControl solver_control (src.size(), 1.0e-12*src.l2_norm());
    SolverCG<>    cg (solver_control);

    dst = 0;

    try
      {
	cg.solve (*matrix, dst, src, *preconditioner);
      }
    catch (...)
      {
	std::cout << "Failure in " << __PRETTY_FUNCTION__ << std::endl;
	abort ();
      }

#ifdef STEP_42_TEST
    if (name == "in schur")
      std::cout << "      " << solver_control.last_step()
		<< " inner CG steps inside the Schur complement ";
    else if (name == "top left")
      std::cout << "    " << solver_control.last_step()
		<< " CG steps on the top left block ";
    else if (name == "rhs")
      std::cout << "    " << solver_control.last_step()
		<< " CG steps for computing the r.h.s. ";
    else
      abort ();

    std::cout << solver_control.initial_value() << "->" << solver_control.last_value()
	      << std::endl;
#endif
  }


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
  {}

				   // Now the interesting function, the multiplication of
				   // the preconditioner with a BlockVector.
  template <class PreconditionerA, class PreconditionerMp>
  void BlockSchurPreconditioner<PreconditionerA, PreconditionerMp>::vmult (
    BlockVector<double>       &dst,
    const BlockVector<double> &src) const
  {
				     // Form u_new = A^{-1} u
    a_preconditioner.vmult (dst.block(0), src.block(0));
				     // Form tmp = - B u_new + p
				     // (<code>SparseMatrix::residual</code>
				     // does precisely this)
    system_matrix->block(1,0).residual(tmp, dst.block(0), src.block(1));
				     // Change sign in tmp
    tmp *= -1;
				     // Multiply by approximate Schur complement
				     // (i.e. a pressure mass matrix)
    m_inverse->vmult (dst.block(1), tmp);
  }

  template <class Preconditioner>
  class SchurComplement : public Subscriptor
  {
    public:
      SchurComplement (const BlockSparseMatrix<double> &system_matrix,
		       const InverseMatrix<SparseMatrix<double>, Preconditioner> &A_inverse);

      void vmult (Vector<double>       &dst,
		  const Vector<double> &src) const;

      unsigned int m() const
	{
	  return system_matrix->block(1,1).m();
	}

      unsigned int n() const
	{
	  return system_matrix->block(1,1).n();
	}

    private:
      const SmartPointer<const BlockSparseMatrix<double> > system_matrix;
      const SmartPointer<const InverseMatrix<SparseMatrix<double>, Preconditioner> > A_inverse;

      mutable Vector<double> tmp1, tmp2;
  };



  template <class Preconditioner>
  SchurComplement<Preconditioner>::
  SchurComplement (const BlockSparseMatrix<double> &system_matrix,
		   const InverseMatrix<SparseMatrix<double>,Preconditioner> &A_inverse)
		  :
		  system_matrix (&system_matrix),
		  A_inverse (&A_inverse),
		  tmp1 (system_matrix.block(0,0).m()),
		  tmp2 (system_matrix.block(0,0).m())
  {}


  template <class Preconditioner>
  void SchurComplement<Preconditioner>::vmult (Vector<double>       &dst,
					       const Vector<double> &src) const
  {
    system_matrix->block(0,1).vmult (tmp1, src);
    A_inverse->name = "in schur";
    A_inverse->vmult (tmp2, tmp1);
    system_matrix->block(1,0).vmult (dst, tmp2);
    dst *= -1;
    system_matrix->block(1,1).vmult_add (dst, src);
    dst *= -1;
  }



  template <int dim>
  StokesProblem<dim>::StokesProblem (const unsigned int degree)
		  :
		  degree (degree),
		  triangulation (Triangulation<dim>::limit_level_difference_at_vertices),
		  fe (FE_Q<dim>(degree+1), dim,
		      FE_Q<dim>(degree), 1),
		  dof_handler (triangulation)
  {}





  template <int dim>
  void StokesProblem<dim>::setup_dofs ()
  {
    A_preconditioner.reset ();
    mg_A_preconditioner.resize (0);
    system_matrix.clear ();

    dof_handler.distribute_dofs (fe);
//  DoFRenumbering::Cuthill_McKee (dof_handler);

    std::vector<unsigned int> block_component (dim+1,0);
    block_component[dim] = 1;
    DoFRenumbering::component_wise (dof_handler, block_component);


    {
      constraints.clear ();
      typename FunctionMap<dim>::type      dirichlet_boundary;
      ZeroFunction<dim>                    homogeneous_dirichlet_bc (dim+1); //TODO: go back to BoundaryValues

      dirichlet_boundary[0] = &homogeneous_dirichlet_bc;
      MappingQ1<dim> mapping;

      std::vector<bool> component_mask (dim+1, true);
      component_mask[dim] = false;
      VectorTools::interpolate_boundary_values (mapping,
						dof_handler,
						dirichlet_boundary,
						constraints,
						component_mask);

      DoFTools::make_hanging_node_constraints (dof_handler,
					       constraints);

      mg_constrained_dofs.clear();
      mg_constrained_dofs.initialize(dof_handler, dirichlet_boundary);
    }

    constraints.close ();


    std::vector<unsigned int> dofs_per_block (2);
    DoFTools::count_dofs_per_block (dof_handler, dofs_per_block,
				    block_component);
    const unsigned int n_u = dofs_per_block[0],
		       n_p = dofs_per_block[1];

    std::cout << "   Number of active cells: "
	      << triangulation.n_active_cells()
	      << std::endl
	      << "   Number of degrees of freedom: "
	      << dof_handler.n_dofs()
	      << " (" << n_u << '+' << n_p << ')'
	      << std::endl;

    {
      BlockCompressedSimpleSparsityPattern csp (2,2);

      csp.block(0,0).reinit (n_u, n_u);
      csp.block(1,0).reinit (n_p, n_u);
      csp.block(0,1).reinit (n_u, n_p);
      csp.block(1,1).reinit (n_p, n_p);

      csp.collect_sizes();

      DoFTools::make_sparsity_pattern (
        static_cast<const DoFHandler<dim>&>(dof_handler),
        csp, constraints, false);
      sparsity_pattern.copy_from (csp);
    }


    system_matrix.reinit (sparsity_pattern);

    solution.reinit (2);
    solution.block(0).reinit (n_u);
    solution.block(1).reinit (n_p);
    solution.collect_sizes ();

    system_rhs.reinit (2);
    system_rhs.block(0).reinit (n_u);
    system_rhs.block(1).reinit (n_p);
    system_rhs.collect_sizes ();

				     //now setup stuff for mg
    const unsigned int nlevels = triangulation.n_levels();

    mg_matrices.resize(0, nlevels-1);
    mg_matrices.clear ();
    mg_interface_matrices.resize(0, nlevels-1);
    mg_interface_matrices.clear ();
    mg_sparsity.resize(0, nlevels-1);

    mg_dofs_per_component.resize (nlevels);
    for (unsigned int level=0; level<nlevels; ++level)
      mg_dofs_per_component[level].resize (2);

    MGTools::count_dofs_per_block (dof_handler, mg_dofs_per_component,
				   block_component);
    for (unsigned int level=0; level<nlevels; ++level)
      std::cout << "                        Level " << level << ": "
		<< dof_handler.n_dofs (level) << " ("
		<< mg_dofs_per_component[level][0] << '+'
		<< mg_dofs_per_component[level][1] << ')'
		<< std::endl;

    for (unsigned int level=0; level<nlevels; ++level)
      {
	DoFRenumbering::component_wise (dof_handler, level, block_component);

	BlockCompressedSparsityPattern bcsp (mg_dofs_per_component[level],
					     mg_dofs_per_component[level]);
	MGTools::make_sparsity_pattern(dof_handler, bcsp, level);
	mg_sparsity[level].copy_from (bcsp);
	mg_matrices[level].reinit (mg_sparsity[level]);
	mg_interface_matrices[level].reinit (mg_sparsity[level]);
      }
  }



  template <int dim>
  void StokesProblem<dim>::assemble_system ()
  {
    system_matrix=0;
    system_rhs=0;

    QGauss<dim>   quadrature_formula(degree+2);

    FEValues<dim> fe_values (fe, quadrature_formula,
			     update_values    |
			     update_quadrature_points  |
			     update_JxW_values |
			     update_gradients);

    const unsigned int   dofs_per_cell   = fe.dofs_per_cell;

    const unsigned int   n_q_points      = quadrature_formula.size();

    FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>       local_rhs (dofs_per_cell);

    std::vector<unsigned int> local_dof_indices (dofs_per_cell);

    const RightHandSide<dim>          right_hand_side;
    std::vector<Vector<double> >      rhs_values (n_q_points,
						  Vector<double>(dim+1));


    const FEValuesExtractors::Vector velocities (0);
    const FEValuesExtractors::Scalar pressure (dim);



    std::vector<Tensor<2,dim> >          phi_grads_u (dofs_per_cell);
    std::vector<double>                  div_phi_u   (dofs_per_cell);
    std::vector<double>                  phi_p       (dofs_per_cell);

    typename MGDoFHandler<dim>::active_cell_iterator
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
	    for (unsigned int k=0; k<dofs_per_cell; ++k)
	      {
		phi_grads_u[k] = fe_values[velocities].gradient (k, q);
		div_phi_u[k]   = fe_values[velocities].divergence (k, q);
		phi_p[k]       = fe_values[pressure].value (k, q);
	      }

	    for (unsigned int i=0; i<dofs_per_cell; ++i)
	      {
		for (unsigned int j=0; j<dofs_per_cell; ++j)
		  {
		    local_matrix(i,j) += (scalar_product(phi_grads_u[i], phi_grads_u[j])
					  - div_phi_u[i] * phi_p[j]
					  - phi_p[i] * div_phi_u[j]
					  - phi_p[i] * phi_p[j]
		    )
					 * fe_values.JxW(q);
		  }

		const unsigned int component_i =
		  fe.system_to_component_index(i).first;
		local_rhs(i) += fe_values.shape_value(i,q) *
				rhs_values[q](component_i) *
				fe_values.JxW(q);
	      }
	  }




	cell->get_dof_indices (local_dof_indices);
	constraints.distribute_local_to_global (local_matrix, local_rhs,
						local_dof_indices,
						system_matrix, system_rhs);
      }
  }


  template <int dim>
  void StokesProblem<dim>::assemble_multigrid ()
  {
    QGauss<dim>   quadrature_formula(degree+2);
    FEValues<dim> fe_values (fe, quadrature_formula,
			     update_values    |
			     update_quadrature_points  |
			     update_JxW_values |
			     update_gradients);

    const unsigned int   dofs_per_cell   = fe.dofs_per_cell;
    const unsigned int   n_q_points      = quadrature_formula.size();

    FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);

    std::vector<unsigned int> local_dof_indices (dofs_per_cell);

    const FEValuesExtractors::Vector velocities (0);
    const FEValuesExtractors::Scalar pressure (dim);


    std::vector<Tensor<2,dim> >          phi_grads_u (dofs_per_cell);
    std::vector<double>                  div_phi_u   (dofs_per_cell);
    std::vector<double>                  phi_p       (dofs_per_cell);

    std::vector<std::vector<bool> > interface_dofs
      = mg_constrained_dofs.get_refinement_edge_indices ();
    std::vector<std::vector<bool> > boundary_interface_dofs
      = mg_constrained_dofs.get_refinement_edge_boundary_indices ();

    std::vector<ConstraintMatrix> boundary_constraints (triangulation.n_levels());
    std::vector<ConstraintMatrix> boundary_interface_constraints (triangulation.n_levels());
    for (unsigned int level=0; level<triangulation.n_levels(); ++level)
      {
	boundary_constraints[level].add_lines (interface_dofs[level]);
	boundary_constraints[level].add_lines (mg_constrained_dofs.get_boundary_indices()[level]);
	boundary_constraints[level].close ();

	boundary_interface_constraints[level]
	  .add_lines (boundary_interface_dofs[level]);
	boundary_interface_constraints[level].close ();
      }

    typename MGDoFHandler<dim>::cell_iterator
      cell = dof_handler.begin(),
      endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
					 // Remember the level of the
					 // current cell.
	const unsigned int level = cell->level();
					 // Compute the values specified
					 // by update flags above.
	fe_values.reinit (cell);
	local_matrix = 0;

	for (unsigned int q=0; q<n_q_points; ++q)
	  {
	    for (unsigned int k=0; k<dofs_per_cell; ++k)
	      {
		phi_grads_u[k] = fe_values[velocities].gradient (k, q);
		div_phi_u[k]   = fe_values[velocities].divergence (k, q);
		phi_p[k]       = fe_values[pressure].value (k, q);
	      }

	    for (unsigned int i=0; i<dofs_per_cell; ++i)
	      for (unsigned int j=0; j<dofs_per_cell; ++j)
		local_matrix(i,j) += (
		  scalar_product(phi_grads_u[i], phi_grads_u[j])
		  - div_phi_u[i] * phi_p[j]
		  - phi_p[i] * div_phi_u[j]
//                                        - phi_p[i] * phi_p[j]
		)
				     * fe_values.JxW(q);
	  }

	cell->get_mg_dof_indices (local_dof_indices);
	boundary_constraints[level]
	  .distribute_local_to_global (local_matrix,
				       local_dof_indices,
				       mg_matrices[level]);

	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  for (unsigned int j=0; j<dofs_per_cell; ++j)
	    if( !(interface_dofs[level][local_dof_indices[i]]==true &&
		  interface_dofs[level][local_dof_indices[j]]==false))
	      local_matrix(i,j) = 0;

	boundary_interface_constraints[level]
	  .distribute_local_to_global (local_matrix,
				       local_dof_indices,
				       mg_interface_matrices[level]);
      }

    mg_A_preconditioner.resize (triangulation.n_levels());
    for (unsigned int level=0; level<triangulation.n_levels(); ++level)
      {
	mg_A_preconditioner[level]
	  = std_cxx1x::shared_ptr<typename InnerPreconditioner<dim>::type>(new typename InnerPreconditioner<dim>::type());
	mg_A_preconditioner[level]
	  ->initialize (mg_matrices[level].block(0,0),
			typename InnerPreconditioner<dim>::type::AdditionalData());
      }
  }


  template <typename InnerPreconditioner>
  class SchurComplementSmoother
  {
    public:
      struct AdditionalData
      {
	  const InnerPreconditioner *A_preconditioner;
      };

      void initialize (const BlockSparseMatrix<double> &system_matrix,
		       const AdditionalData            &data);

      void vmult (BlockVector<double> &dst,
		  const BlockVector<double> &src) const;

      void Tvmult (BlockVector<double> &dst,
		   const BlockVector<double> &src) const;

      void clear ();

    private:
      SmartPointer<const BlockSparseMatrix<double> > system_matrix;
      SmartPointer<const InnerPreconditioner>        A_preconditioner;
  };


  template <typename InnerPreconditioner>
  void
  SchurComplementSmoother<InnerPreconditioner>::
  initialize (const BlockSparseMatrix<double> &system_matrix,
	      const AdditionalData            &data)
  {
    this->system_matrix    = &system_matrix;
    this->A_preconditioner = data.A_preconditioner;
  }




  template <typename InnerPreconditioner>
  void
  SchurComplementSmoother<InnerPreconditioner>::
  vmult (BlockVector<double> &dst,
	 const BlockVector<double> &src) const
  {
#ifdef STEP_42_TEST
    std::cout << "Entering smoother with " << dst.size() << " unknowns" << std::endl;
#endif

    SparseDirectUMFPACK direct_solver;
    direct_solver.initialize(*system_matrix);
    Vector<double> solution, rhs;
    solution = dst;
    rhs = src;
    direct_solver.vmult(solution, rhs);
    dst = solution;
/*
  const InverseMatrix<SparseMatrix<double>,InnerPreconditioner>
  A_inverse (system_matrix->block(0,0), *A_preconditioner);
  Vector<double> tmp (dst.block(0).size());


  {
  Vector<double> schur_rhs (dst.block(1).size());
  A_inverse.name = "rhs";
  A_inverse.vmult (tmp, src.block(0));
//    std::cout << " TMP " << tmp.l2_norm() << std::endl;
system_matrix->block(1,0).vmult (schur_rhs, tmp);
schur_rhs -= src.block(1);
//    std::cout << " BLOCK 1 " << src.block(1).l2_norm() << std::endl;
//    std::cout << " SCHUR RHS " << schur_rhs.l2_norm() << std::endl;

SchurComplement<InnerPreconditioner>
schur_complement (*system_matrix, A_inverse);

				 // The usual control structures for
								  // the solver call are created...
								  SolverControl solver_control (dst.block(1).size(),
								  1e-1*schur_rhs.l2_norm());
								  SolverGMRES<>    cg (solver_control);

								  #ifdef STEP_42_TEST
								  std::cout << "    Starting Schur complement solver -- "
								  << schur_complement.m() << " unknowns"
								  << std::endl;
								  #endif
								  try
								  {
								  cg.solve (schur_complement, dst.block(1), schur_rhs,
								  PreconditionIdentity());
								  }
								  catch (...)
								  {
								  std::cout << "Failure in " << __PRETTY_FUNCTION__ << std::endl;
								  std::cout << schur_rhs.l2_norm () << std::endl;
								  abort ();
								  }

// no constraints to be taken care of here
#ifdef STEP_42_TEST
std::cout << "    "
<< solver_control.last_step()
<< " CG Schur complement iterations in smoother "
<< solver_control.initial_value() << "->" << solver_control.last_value()
<< std::endl;
#endif
}


{
system_matrix->block(0,1).vmult (tmp, dst.block(1));
tmp *= -1;
tmp += src.block(0);

A_inverse.name = "top left";
A_inverse.vmult (dst.block(0), tmp);
// no constraints here either
}
#ifdef STEP_42_TEST
std::cout << "Exiting smoother with " << dst.size() << " unknowns" << std::endl;
#endif
*/
  }



  template <typename InnerPreconditioner>
  void
  SchurComplementSmoother<InnerPreconditioner>::clear ()
  {}



  template <typename InnerPreconditioner>
  void
  SchurComplementSmoother<InnerPreconditioner>::
  Tvmult (BlockVector<double> &,
	  const BlockVector<double> &) const
  {
    Assert (false, ExcNotImplemented());
  }



  template <int dim>
  void StokesProblem<dim>::solve ()
  {
    system_matrix.block(1,1) = 0;
    assemble_multigrid ();
    typedef PreconditionMG<dim, BlockVector<double>, MGTransferPrebuilt<BlockVector<double> > >
      MGPREC;

    GrowingVectorMemory<BlockVector<double> >  mg_vector_memory;

    MGTransferPrebuilt<BlockVector<double> > mg_transfer(constraints, mg_constrained_dofs);
    std::vector<unsigned int> block_component (dim+1,0);
    block_component[dim] = 1;
    mg_transfer.set_component_to_block_map (block_component);
    mg_transfer.build_matrices(dof_handler);

    FullMatrix<float> mg_coarse_matrix;
    mg_coarse_matrix.copy_from (mg_matrices[0]);
    MGCoarseGridHouseholder<float, BlockVector<double> > mg_coarse;
    mg_coarse.initialize(mg_coarse_matrix);

    MGMatrix<BlockSparseMatrix<double>, BlockVector<double> >
      mg_matrix(&mg_matrices);
    MGMatrix<BlockSparseMatrix<double>, BlockVector<double> >
      mg_interface_up(&mg_interface_matrices);
    MGMatrix<BlockSparseMatrix<double>, BlockVector<double> >
      mg_interface_down(&mg_interface_matrices);

    typedef
      SchurComplementSmoother<typename InnerPreconditioner<dim>::type>
      Smoother;

    MGSmootherPrecondition<BlockSparseMatrix<double>,
      Smoother,
      BlockVector<double> >
      mg_smoother(mg_vector_memory);

    MGLevelObject<typename Smoother::AdditionalData>
      smoother_data (0, triangulation.n_levels()-1);

    for (unsigned int level=0; level<triangulation.n_levels(); ++level)
      smoother_data[level].A_preconditioner = mg_A_preconditioner[level].get();

    mg_smoother.initialize(mg_matrices, smoother_data);
    mg_smoother.set_steps(2);

    Multigrid<BlockVector<double> > mg(dof_handler,
				       mg_matrix,
				       mg_coarse,
				       mg_transfer,
				       mg_smoother,
				       mg_smoother);
    mg.set_debug(3);
    mg.set_edge_matrices(mg_interface_down, mg_interface_up);

    MGPREC  preconditioner(dof_handler, mg, mg_transfer);

    SolverControl solver_control (system_matrix.m(),
				  1e-6*system_rhs.l2_norm());
    GrowingVectorMemory<BlockVector<double> > vector_memory;
    SolverGMRES<BlockVector<double> >::AdditionalData gmres_data;
    gmres_data.max_n_tmp_vectors = 100;

    SolverGMRES<BlockVector<double> > gmres(solver_control, vector_memory,
					    gmres_data);

//  PreconditionIdentity precondition_identity;
#ifdef STEP_42_TEST
    std::cout << "Starting outer GMRES complement solver" << std::endl;
#endif
    try
      {
	gmres.solve(system_matrix, solution, system_rhs,
		    preconditioner);
      }
    catch (...)
      {
	std::cout << "Failure in " << __PRETTY_FUNCTION__ << std::endl;
	abort ();
      }

    constraints.distribute (solution);

    std::cout << solver_control.last_step()
	      << " outer GMRES iterations ";
  }


  template <int dim>
  void StokesProblem<dim>::solve_block ()
  {
    std::cout << "   Computing preconditioner..." << std::endl << std::flush;

    A_preconditioner
      = std_cxx1x::shared_ptr<typename InnerPreconditioner<dim>::type>(new typename InnerPreconditioner<dim>::type());
    A_preconditioner->initialize (system_matrix.block(0,0),
				  typename InnerPreconditioner<dim>::type::AdditionalData());

    SparseMatrix<double> pressure_mass_matrix;
    pressure_mass_matrix.reinit(sparsity_pattern.block(1,1));
    pressure_mass_matrix.copy_from(system_matrix.block(1,1));
    system_matrix.block(1,1) = 0;

    SparseILU<double> pmass_preconditioner;
    pmass_preconditioner.initialize (pressure_mass_matrix,
				     SparseILU<double>::AdditionalData());

    InverseMatrix<SparseMatrix<double>,SparseILU<double> >
      m_inverse (pressure_mass_matrix, pmass_preconditioner);

    BlockSchurPreconditioner<typename InnerPreconditioner<dim>::type,
      SparseILU<double> >
      preconditioner (system_matrix, m_inverse, *A_preconditioner);

    SolverControl solver_control (system_matrix.m(),
				  1e-6*system_rhs.l2_norm());
    GrowingVectorMemory<BlockVector<double> > vector_memory;
    SolverGMRES<BlockVector<double> >::AdditionalData gmres_data;
    gmres_data.max_n_tmp_vectors = 100;

    SolverGMRES<BlockVector<double> > gmres(solver_control, vector_memory,
					    gmres_data);

    gmres.solve(system_matrix, solution, system_rhs,
		preconditioner);

    constraints.distribute (solution);

    std::cout << " "
	      << solver_control.last_step()
	      << " block GMRES iterations ";
  }


  template <int dim>
  void
  StokesProblem<dim>::output_results (const unsigned int refinement_cycle)  const
  {
    std::vector<std::string> solution_names (dim, "velocity");
    solution_names.push_back ("pressure");

    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation
      (dim, DataComponentInterpretation::component_is_part_of_vector);
    data_component_interpretation
      .push_back (DataComponentInterpretation::component_is_scalar);

    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
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



  template <int dim>
  void
  StokesProblem<dim>::refine_mesh ()
  {
    Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

    std::vector<bool> component_mask (dim+1, false);
    component_mask[dim] = true;
    KellyErrorEstimator<dim>::estimate (static_cast<const DoFHandler<dim>&>(dof_handler),
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



  template <int dim>
  void StokesProblem<dim>::run ()
  {
    {
      std::vector<unsigned int> subdivisions (dim, 1);
      subdivisions[0] = 1;

      const Point<dim> bottom_left = (dim == 2 ?
				      Point<dim>(0,0) :
				      Point<dim>(0,0,0));
      const Point<dim> top_right   = (dim == 2 ?
				      Point<dim>(1,1) :
				      Point<dim>(1,1,1));

      GridGenerator::subdivided_hyper_rectangle (triangulation,
						 subdivisions,
						 bottom_left,
						 top_right);
    }


    for (typename Triangulation<dim>::active_cell_iterator
	   cell = triangulation.begin_active();
	 cell != triangulation.end(); ++cell)
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
	if (cell->face(f)->center()[0] == 1)
	  cell->face(f)->set_all_boundary_indicators(1);



    triangulation.refine_global (1);


    for (unsigned int refinement_cycle = 0; refinement_cycle<10;
	 ++refinement_cycle)
      {
	std::cout << "Refinement cycle " << refinement_cycle << std::endl;

	if (refinement_cycle > 0)
	  refine_mesh ();

	std::ostringstream out_filename;
	out_filename << "gitter"
		     << refinement_cycle
		     << ".eps";

	std::ofstream grid_output (out_filename.str().c_str());
	GridOut grid_out;
	grid_out.write_eps (triangulation, grid_output);

	setup_dofs ();

	std::cout << "   Assembling..." << std::endl << std::flush;
	assemble_system ();

	std::cout << "   Solving..." << std::flush;

	solve_block ();
	output_results (refinement_cycle);
	system ("mv solution-* block");

	solution = 0;

	solve ();
	output_results (refinement_cycle);
	system ("mv solution-* mg");

	std::cout << std::endl;
      }
  }
}



int main ()
{
  try
    {
      using namespace dealii;
      using namespace Step42;

      deallog.depth_console (0);

      StokesProblem<2> flow_problem(1);
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
