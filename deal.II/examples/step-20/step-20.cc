/* $Id$ */
/* Author: Wolfgang Bangerth, Texas A&M University, 2005, 2006 */

/*    $Id$       */
/*    Version: $Name$                                          */
/*                                                                */
/*    Copyright (C) 2005, 2006 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */

                                 // @sect3{Include files}

				 // Since this program is only an
				 // adaptation of step-4, there is not
				 // much new stuff in terms of header
				 // files. In deal.II, we usually list
				 // include files in the order
				 // base-lac-grid-dofs-fe-numerics,
				 // followed by C++ standard include
				 // files:
#include <base/quadrature_lib.h>
#include <base/logstream.h>
#include <base/function.h>
#include <lac/block_vector.h>
#include <lac/full_matrix.h>
#include <lac/block_sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/solver_gmres.h>
#include <lac/precondition.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_renumbering.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <fe/fe_dgq.h>
#include <fe/fe_system.h>
#include <fe/fe_values.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <numerics/data_out.h>

#include <fstream>
#include <iostream>

				 // This is the only new header,
				 // namely the one in which the
				 // Raviart-Thomas finite element is
				 // declared:
#include <fe/fe_raviart_thomas.h>


                                 // @sect3{The ``MixedLaplaceProblem'' class template}

				 // Again, since this is an adaptation
				 // of step-6, the main class is
				 // almost the same as the one in that
				 // tutorial program. In terms of
				 // member functions, the main
				 // differences are that the
				 // constructor takes the degree of
				 // the Raviart-Thomas element as an
				 // argument (and that there is a
				 // corresponding member variable to
				 // store this value) and the addition
				 // of the ``compute_error'' function
				 // in which, no surprise, we will
				 // compute the difference between the
				 // exact and the numerical solution
				 // to determine convergence of our
				 // computations:
template <int dim>
class MixedLaplaceProblem 
{
  public:
    MixedLaplaceProblem (const unsigned int degree);
    void run ();
    
  private:
    void make_grid_and_dofs ();
    void assemble_system ();
    void solve ();
    void compute_errors () const;
    void output_results () const;

    const unsigned int   degree;
    
    Triangulation<dim>   triangulation;
    FESystem<dim>        fe;
    DoFHandler<dim>      dof_handler;

				     // The second difference is that
				     // the sparsity pattern, the
				     // system matrix, and solution
				     // and right hand side vectors
				     // are now blocked. What this
				     // means and what one can do with
				     // such objects is explained in
				     // the introduction to this
				     // program as well as further
				     // down below when we explain the
				     // linear solvers and
				     // preconditioners for this
				     // problem:
    BlockSparsityPattern      sparsity_pattern;
    BlockSparseMatrix<double> system_matrix;

    BlockVector<double>       solution;
    BlockVector<double>       system_rhs;
};


				 // @sect3{Right hand side, coefficient, and exact solution}

				 // Our next task is to define the
				 // right hand side of our problem
				 // (i.e., the scalar right hand side
				 // for the pressure in the original
				 // Laplace equation), boundary values
				 // for the pressure, as well as a
				 // function that describes both the
				 // pressure and the velocity of the
				 // exact solution for later
				 // computations of the error. Note
				 // that these functions have one,
				 // one, and ``dim+1'' components,
				 // respectively, and that we pass the
				 // number of components down to the
				 // ``Function<dim>'' base class. For
				 // the exact solution, we only
				 // declare the function that actually
				 // returns the entire solution vector
				 // (i.e. all components of it) at
				 // once. Here are the respective
				 // declarations:
template <int dim>
class RightHandSide : public Function<dim> 
{
  public:
    RightHandSide () : Function<dim>(1) {};
    
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;
};



template <int dim>
class PressureBoundaryValues : public Function<dim> 
{
  public:
    PressureBoundaryValues () : Function<dim>(1) {};
    
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;
};


template <int dim>
class ExactSolution : public Function<dim> 
{
  public:
    ExactSolution () : Function<dim>(dim+1) {};
    
    virtual void vector_value (const Point<dim> &p, 
			       Vector<double>   &value) const;
};


				 // And then we also have to define
				 // these respective functions, of
				 // course. Given the ones that we
				 // discussed in the introduction, the
				 // following computations should be
				 // straightforward:
template <int dim>
double RightHandSide<dim>::value (const Point<dim>  &/*p*/,
				  const unsigned int /*component*/) const 
{
  return 0;
}



template <int dim>
double PressureBoundaryValues<dim>::value (const Point<dim>  &p,
					   const unsigned int /*component*/) const 
{
  const double alpha = 0.1;
  const double beta = 1;
  return -(alpha*p[0]*p[1]*p[1]/2 + beta*p[0] - alpha*p[0]*p[0]*p[0]/6);
}



template <int dim>
void
ExactSolution<dim>::vector_value (const Point<dim> &p,
				  Vector<double>   &values) const 
{
  Assert (values.size() == dim+1,
	  ExcDimensionMismatch (values.size(), dim+1));

  const double alpha = 0.1;
  const double beta = 1;

  values(0) = alpha*p[1]*p[1]/2 + beta - alpha*p[0]*p[0]/2;
  values(1) = alpha*p[0]*p[1];
  values(2) = -(alpha*p[0]*p[1]*p[1]/2 + beta*p[0] - alpha*p[0]*p[0]*p[0]/6);
}



template <int dim>
MixedLaplaceProblem<dim>::MixedLaplaceProblem (const unsigned int degree)
		:
		degree (degree),
                fe (FE_RaviartThomas<dim>(degree),1,FE_DGQ<dim>(degree),1),
		dof_handler (triangulation)
{}


template <int dim>
void MixedLaplaceProblem<dim>::make_grid_and_dofs ()
{
  GridGenerator::hyper_cube (triangulation, -1, 1);
  triangulation.refine_global (3);
  
  std::cout << "   Number of active cells: "
	    << triangulation.n_active_cells()
	    << std::endl
	    << "   Total number of cells: "
	    << triangulation.n_cells()
	    << std::endl;

  dof_handler.distribute_dofs (fe);
  DoFRenumbering::component_wise (dof_handler);
  
  std::vector<unsigned int> dofs_per_component (dim+1);
  DoFTools::count_dofs_per_component (dof_handler, dofs_per_component);
  const unsigned int n_u = dofs_per_component[0],
                     n_p = dofs_per_component[dim];

  std::cout << "   Number of degrees of freedom: "
	    << dof_handler.n_dofs()
            << " (" << n_u << '+' << n_p << ')'
	    << std::endl;
  
  sparsity_pattern.reinit (2,2);
  sparsity_pattern.block(0,0).reinit (n_u, n_u,
                                      dof_handler.max_couplings_between_dofs());
  sparsity_pattern.block(1,0).reinit (n_p, n_u,
                                      dof_handler.max_couplings_between_dofs());
  sparsity_pattern.block(0,1).reinit (n_u, n_p,
                                      dof_handler.max_couplings_between_dofs());
  sparsity_pattern.block(1,1).reinit (n_p, n_p,
				      dof_handler.max_couplings_between_dofs());
  sparsity_pattern.collect_sizes();
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
  sparsity_pattern.compress();

  system_matrix.reinit (sparsity_pattern);

  std::vector<unsigned int> block_components (2);
  block_components[0] = n_u;
  block_components[1] = n_p;
  solution.reinit (block_components);
  system_rhs.reinit (block_components);
}




Tensor<1,2> extract_u (const FEValuesBase<2> &fe_values,
                       const unsigned int j,
                       const unsigned int q)
{
  Tensor<1,2> tmp;
  tmp[0] = fe_values.shape_value_component (j,q,0);
  tmp[1] = fe_values.shape_value_component (j,q,1);
  return tmp;
}



Tensor<1,3> extract_u (const FEValuesBase<3> &fe_values,
                       const unsigned int j,
                       const unsigned int q)
{
  Tensor<1,3> tmp;
  tmp[0] = fe_values.shape_value_component (j,q,0);
  tmp[1] = fe_values.shape_value_component (j,q,1);
  tmp[2] = fe_values.shape_value_component (j,q,2);
  return tmp;
}





double extract_div_u (const FEValuesBase<2> &fe_values,
                      const unsigned int j,
                      const unsigned int q)
{
  return fe_values.shape_grad_component (j,q,0)[0] +
    fe_values.shape_grad_component (j,q,1)[1];
}


double extract_div_u (const FEValuesBase<3> &fe_values,
                      const unsigned int j,
                      const unsigned int q)
{
  return fe_values.shape_grad_component (j,q,0)[0] +
    fe_values.shape_grad_component (j,q,1)[1] +
    fe_values.shape_grad_component (j,q,2)[2];
}

  
template <int dim>
double extract_p (const FEValuesBase<dim> &fe_values,
                  const unsigned int j,
                  const unsigned int q)
{
  return fe_values.shape_value_component (j,q,dim);
}


template <int dim>
void MixedLaplaceProblem<dim>::assemble_system () 
{  
  QGauss<dim>   quadrature_formula(degree+2);
  QGauss<dim-1> face_quadrature_formula(degree+2);

  FEValues<dim> fe_values (fe, quadrature_formula, 
			   update_values    | update_gradients |
                           update_q_points  | update_JxW_values);
  FEFaceValues<dim> fe_face_values (fe, face_quadrature_formula, 
				    update_values    | update_normal_vectors |
				    update_q_points  | update_JxW_values);

  const unsigned int   dofs_per_cell   = fe.dofs_per_cell;
  const unsigned int   n_q_points      = quadrature_formula.n_quadrature_points;
  const unsigned int   n_face_q_points = face_quadrature_formula.n_quadrature_points;

  FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       local_rhs (dofs_per_cell);


  const RightHandSide<dim> right_hand_side;
  const PressureBoundaryValues<dim> pressure_boundary_values;

  std::vector<double> rhs_values (n_q_points);
  std::vector<double> boundary_values (n_face_q_points);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      local_matrix = 0;
      local_rhs = 0;

      right_hand_side.value_list (fe_values.get_quadrature_points(),
                                  rhs_values);
      
      for (unsigned int q=0; q<n_q_points; ++q) 
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            const Tensor<1,dim> phi_i_u = extract_u (fe_values, i, q);
            const double div_phi_i_u = extract_div_u (fe_values, i, q);
            const double phi_i_p = extract_p (fe_values, i, q);
            
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              {
                const Tensor<1,dim> phi_j_u = extract_u (fe_values, j, q);
                const double div_phi_j_u = extract_div_u (fe_values, j, q);
                const double phi_j_p = extract_p (fe_values, j, q);
                
                local_matrix(i,j) += (phi_i_u * phi_j_u
                                      - div_phi_i_u * phi_j_p
                                      - phi_i_p * div_phi_j_u)
                                     * fe_values.JxW(q);
              }

            local_rhs(i) += -phi_i_p *
                            rhs_values[q] *
                            fe_values.JxW(q);
          }

      for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
	if (cell->at_boundary(face_no))
	  {
	    fe_face_values.reinit (cell, face_no);
	    
	    pressure_boundary_values.value_list (fe_face_values.get_quadrature_points(),
						 boundary_values);

	    for (unsigned int q=0; q<n_face_q_points; ++q) 
	      for (unsigned int i=0; i<dofs_per_cell; ++i)
		{
		  const Tensor<1,dim> phi_i_u = extract_u (fe_face_values, i, q);

		  local_rhs(i) += -(phi_i_u *
				    fe_face_values.normal_vector(q) *
				    boundary_values[q] *
				    fe_face_values.JxW(q));
		}
	  }
      
      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        for (unsigned int j=0; j<dofs_per_cell; ++j)
          system_matrix.add (local_dof_indices[i],
                             local_dof_indices[j],
                             local_matrix(i,j));
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        system_rhs(local_dof_indices[i]) += local_rhs(i);
    }
}


class SchurComplement 
{
  public:
    SchurComplement (const BlockSparseMatrix<double> &A)
                    :
                    A (A),
                    tmp1 (A.block(0,0).m()),
                    tmp2 (A.block(0,0).m())
      {}

    void vmult (Vector<double>       &dst,
                const Vector<double> &src) const
      {
        A.block(0,1).vmult (tmp1, src);

        SolverControl           solver_control (tmp1.size(),
                                                1e-8*tmp1.l2_norm());
        SolverCG<>              cg (solver_control, vector_memory);

        PreconditionSSOR<> precondition;
        precondition.initialize(A.block(0,0));
        cg.solve (A.block(0,0), tmp2, tmp1, precondition);

        std::cout << "     " << solver_control.last_step()
                  << " inner iterations needed to obtain convergence."
                  << std::endl;
        
        A.block(1,0).vmult (dst, tmp2);
      }

  private:
    const BlockSparseMatrix<double> &A;

    mutable GrowingVectorMemory<> vector_memory;
    
    mutable Vector<double> tmp1, tmp2;
};



template <int dim>
void MixedLaplaceProblem<dim>::solve () 
{
  {
    Vector<double> schur_rhs (solution.block(1).size());
    {
      Vector<double> tmp (solution.block(0).size());
      
      SolverControl solver_control (system_matrix.block(0,0).m(),
				    1e-6*system_rhs.l2_norm());
      SolverCG<> cg (solver_control);

      cg.solve (system_matrix.block(0,0), tmp,
		system_rhs.block(0), PreconditionIdentity());
  
      std::cout << "   " << solver_control.last_step()
		<< " CG mass matrix iterations needed to obtain convergence."
		<< std::endl;

      system_matrix.block(1,0).vmult (schur_rhs, tmp);
      schur_rhs -= system_rhs.block(1);
    }

    SolverControl solver_control (system_matrix.block(0,0).m(),
				  1e-6*schur_rhs.l2_norm());
    SolverCG<>    cg (solver_control);

    cg.solve (SchurComplement(system_matrix), solution.block(1),
              schur_rhs,
              PreconditionIdentity());
  
    std::cout << "   " << solver_control.last_step()
              << " CG Schur complement iterations needed to obtain convergence."
              << std::endl;
  }
  {
    Vector<double> tmp (solution.block(0).size());
    system_matrix.block(0,1).vmult (tmp, solution.block(1));
    tmp *= -1;
    tmp += system_rhs.block(0);
    
    SolverControl solver_control (system_matrix.block(0,0).m(),
				  1e-6*tmp.l2_norm());
    SolverCG<> cg (solver_control);

    cg.solve (system_matrix.block(0,0), solution.block(0),
              tmp, PreconditionIdentity());
  
    std::cout << "   " << solver_control.last_step()
              << " CG mass matrix iterations needed to obtain convergence."
              << std::endl;
  }
}



template <int dim>
void MixedLaplaceProblem<dim>::compute_errors () const
{
  Vector<double> tmp (triangulation.n_active_cells());
  ExactSolution<dim> exact_solution;

				     // do NOT use QGauss here!
  QTrapez<1> q_trapez;
  QIterated<dim> quadrature (q_trapez, 5);
  {
    const ComponentSelectFunction<dim> mask (dim, 1., dim+1);
    VectorTools::integrate_difference (dof_handler, solution, exact_solution,
				       tmp, quadrature,
				       VectorTools::L2_norm,
				       &mask);
  }
  const double p_l2_error = tmp.l2_norm();
  
  double u_l2_error = 0;
  for (unsigned int d=0; d<dim; ++d)
    {
      const ComponentSelectFunction<dim> mask(d, 1., dim+1);
      VectorTools::integrate_difference (dof_handler, solution, exact_solution,
					 tmp, quadrature,
					 VectorTools::L2_norm,
					 &mask);
      u_l2_error = std::sqrt (u_l2_error*u_l2_error +
			      tmp.l2_norm() * tmp.l2_norm());
    }
  
//   double u_h1_error = 0;
//   for (unsigned int d=0; d<dim; ++d)
//     {
//       const ComponentSelectFunction<dim> mask(d, 1., dim+1);
//       VectorTools::integrate_difference (dof_handler, solution, exact_solution,
// 					 tmp, QGauss<dim>(degree+1),
// 					 VectorTools::H1_seminorm,
// 					 &mask);
//       u_h1_error = std::sqrt (u_h1_error*u_h1_error +
// 			      tmp.l2_norm() * tmp.l2_norm());
//     }


  std::cout << "Errors: ||e_p||_L2 = " << p_l2_error
	    << ",   ||e_u||_L2 = " << u_l2_error
//    	    << ",   |e_u|_H1 = " << u_h1_error
	    << std::endl;
}


template <int dim>
void MixedLaplaceProblem<dim>::output_results () const
{
  DataOut<dim> data_out;

  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution, "solution");
  data_out.add_data_vector (system_rhs, "rhs");

  data_out.build_patches (degree+1);

  std::ofstream output (dim == 2 ?
			"solution-2d.gmv" :
			"solution-3d.gmv");
  data_out.write_gmv (output);
}



template <int dim>
void MixedLaplaceProblem<dim>::run () 
{
  std::cout << "Solving problem in " << dim << " space dimensions." << std::endl;
  
  make_grid_and_dofs();
  assemble_system ();
  solve ();
  compute_errors ();
  output_results ();
}

    
int main () 
{
  deallog.depth_console (0);

  MixedLaplaceProblem<2> mixed_laplace_problem (1);
  mixed_laplace_problem.run ();

  return 0;
}
