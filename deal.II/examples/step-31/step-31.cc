/* $Id: step-22.cc 15679 2008-01-24 23:28:37Z bangerth $ */
/* Author: Wolfgang Bangerth, Texas A&M University, 2008 */

/*    $Id: step-22.cc 15679 2008-01-24 23:28:37Z bangerth $       */
/*    Version: $Name$                                          */
/*                                                                */
/*    Copyright (C) 2008 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */



#include <base/quadrature_lib.h>
#include <base/logstream.h>
#include <base/function.h>
#include <base/utilities.h>

#include <lac/block_vector.h>
#include <lac/full_matrix.h>
#include <lac/block_sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>
#include <lac/sparse_direct.h>
#include <lac/sparse_ilu.h>

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
#include <fe/mapping_c1.h>

#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <numerics/data_out.h>
#include <numerics/error_estimator.h>
#include <numerics/solution_transfer.h>

#include <fstream>
#include <sstream>

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

    boost::shared_ptr<typename InnerPreconditioner<dim>::type> A_preconditioner;
};





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
RightHandSide<dim>::value (const Point<dim>  &p,
                           const unsigned int component) const 
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






template <int dim>
Tensor<1,dim>
extract_u (const FEValuesBase<dim> &fe_values,
           const unsigned int i,
           const unsigned int q)
{
  Tensor<1,dim> tmp;

  const unsigned int component
    = fe_values.get_fe().system_to_component_index(i).first;

  if (component < dim)
    tmp[component] = fe_values.shape_value (i,q);

  return tmp;
}



template <int dim>
double
extract_div_u (const FEValuesBase<dim> &fe_values,
               const unsigned int i,
               const unsigned int q)
{
  const unsigned int component
    = fe_values.get_fe().system_to_component_index(i).first;

  if (component < dim)
    return fe_values.shape_grad (i,q)[component];
  else
    return 0;
}



template <int dim>
Tensor<2,dim>
extract_grad_s_u (const FEValuesBase<dim> &fe_values,
		  const unsigned int i,
		  const unsigned int q)
{
  Tensor<2,dim> tmp;

  const unsigned int component
    = fe_values.get_fe().system_to_component_index(i).first;
  
  if (component < dim)
    {
      const Tensor<1,dim> grad_phi_over_2 = fe_values.shape_grad (i,q) / 2;
      
      for (unsigned int e=0; e<dim; ++e)
	tmp[component][e] += grad_phi_over_2[e];
      for (unsigned int d=0; d<dim; ++d)
	tmp[d][component] += grad_phi_over_2[d];
    }
  
  return tmp;
}


  
template <int dim>
double extract_p (const FEValuesBase<dim> &fe_values,
                  const unsigned int i,
                  const unsigned int q)
{
  const unsigned int component
    = fe_values.get_fe().system_to_component_index(i).first;

  if (component == dim)
    return fe_values.shape_value (i,q);
  else
    return 0;
}



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



template <class Preconditioner>
class SchurComplement : public Subscriptor
{
  public:
    SchurComplement (const BlockSparseMatrix<double> &A,
		     const InverseMatrix<SparseMatrix<double>,Preconditioner> &Minv);

    void vmult (Vector<double>       &dst,
		const Vector<double> &src) const;

  private:
    const SmartPointer<const BlockSparseMatrix<double> > system_matrix;
    const SmartPointer<const InverseMatrix<SparseMatrix<double>,Preconditioner> > m_inverse;
    
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



template <int dim>
StokesProblem<dim>::StokesProblem (const unsigned int degree)
                :
                degree (degree),
		triangulation (Triangulation<dim>::maximum_smoothing),
                fe (FE_Q<dim>(degree+1), dim,
                    FE_Q<dim>(degree), 1),
                dof_handler (triangulation)
{}




template <int dim>
void StokesProblem<dim>::setup_dofs ()
{
				   // release preconditioner since it
				   // will definitely not be needed
				   // any more after this point
  A_preconditioner.reset ();
  
  dof_handler.distribute_dofs (fe); 
  DoFRenumbering::Cuthill_McKee (dof_handler);
  DoFRenumbering::component_wise (dof_handler);

  hanging_node_constraints.clear ();
  DoFTools::make_hanging_node_constraints (dof_handler,
					   hanging_node_constraints);
  hanging_node_constraints.close ();

  std::vector<unsigned int> dofs_per_component (dim+1);
  DoFTools::count_dofs_per_component (dof_handler, dofs_per_component);  
  const unsigned int n_u = dofs_per_component[0] * dim,
                     n_p = dofs_per_component[dim];

  std::cout << "   Number of active cells: "
            << triangulation.n_active_cells()
            << std::endl
            << "   Number of degrees of freedom: "
            << dof_handler.n_dofs()
            << " (" << n_u << '+' << n_p <<')'
            << std::endl;
  
  system_matrix.clear ();
      
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


template <int dim>
double
scalar_product (const Tensor<2,dim> &a,
		const Tensor<2,dim> &b)
{
  double tmp = 0;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      tmp += a[i][j] * b[i][j];
  return tmp;
}



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
                                    update_values    | update_normal_vectors |
                                    update_quadrature_points  | update_JxW_values);

  const unsigned int   dofs_per_cell   = fe.dofs_per_cell;
  
  const unsigned int   n_q_points      = quadrature_formula.size();
  const unsigned int   n_face_q_points = face_quadrature_formula.size();

  FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       local_rhs (dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);
  
  const PressureBoundaryValues<dim> pressure_boundary_values;
  
  std::vector<double>               boundary_values (n_face_q_points);
  
  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    { 
      fe_values.reinit (cell);
      local_matrix = 0;
      local_rhs = 0;

      for (unsigned int q=0; q<n_q_points; ++q)
	{
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    {

	      const Tensor<1,dim> phi_i_u      = extract_u (fe_values, i, q);

	      const Tensor<2,dim> phi_i_grads_u= extract_grad_s_u (fe_values, i, q);
	      const double        div_phi_i_u  = extract_div_u (fe_values, i, q);
	      const double        phi_i_p      = extract_p (fe_values, i, q);
            
	      for (unsigned int j=0; j<dofs_per_cell; ++j)
		{
		  const Tensor<2,dim> phi_j_grads_u     = extract_grad_s_u (fe_values, j, q);
		  const double        div_phi_j_u = extract_div_u (fe_values, j, q);
		  const double        phi_j_p     = extract_p (fe_values, j, q);
                
		  local_matrix(i,j) += (scalar_product(phi_i_grads_u, phi_j_grads_u)
					- div_phi_i_u * phi_j_p
					- phi_i_p * div_phi_j_u
					+ phi_i_p * phi_j_p)
				       * fe_values.JxW(q);     
		}
          }
	}
      

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
                    phi_i_u = extract_u (fe_face_values, i, q);

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
  
  std::cout << "   Computing preconditioner..." << std::flush;
      
  A_preconditioner
    = boost::shared_ptr<typename InnerPreconditioner<dim>::type>(new typename InnerPreconditioner<dim>::type());
  A_preconditioner->initialize (system_matrix.block(0,0),
				typename InnerPreconditioner<dim>::type::AdditionalData());

  std::cout << std::endl;
}



template <int dim>
void StokesProblem<dim>::solve () 
{
  const InverseMatrix<SparseMatrix<double>,typename InnerPreconditioner<dim>::type>
    A_inverse (system_matrix.block(0,0), *A_preconditioner);
  Vector<double> tmp (solution.block(0).size());
  Vector<double> schur_rhs (solution.block(1).size());
  
  {
    A_inverse.vmult (tmp, system_rhs.block(0));
    system_matrix.block(1,0).vmult (schur_rhs, tmp);
    schur_rhs -= system_rhs.block(1);

    SchurComplement<typename InnerPreconditioner<dim>::type>
      schur_complement (system_matrix, A_inverse);
    
    SolverControl solver_control (system_matrix.block(0,0).m(),
                                  1e-6*schur_rhs.l2_norm());
    SolverCG<>    cg (solver_control);
    
    PreconditionSSOR<> preconditioner;
    preconditioner.initialize (system_matrix.block(1,1), 1.2);

    InverseMatrix<SparseMatrix<double>,PreconditionSSOR<> >
      m_inverse (system_matrix.block(1,1), preconditioner);
    
    try
      {
	cg.solve (schur_complement, solution.block(1), schur_rhs,
		  m_inverse);
      }
    catch (...)
      {
	abort ();
      }

				     // produce a consistent flow field
    hanging_node_constraints.distribute (solution);
  
    std::cout << "   "
              << solver_control.last_step()
              << " CG Schur complement iterations for pressure."
              << std::endl;    
  }

  {
    system_matrix.block(0,1).vmult (tmp, solution.block(1));
    tmp *= -1;
    tmp += system_rhs.block(0);

    A_inverse.vmult (solution.block(0), tmp);

				     // produce a consistent pressure field
    hanging_node_constraints.distribute (solution);
  }
}
                                 


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


template <int dim>
void StokesProblem<dim>::run () 
{
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
  for (typename Triangulation<dim>::active_cell_iterator
	 cell = triangulation.begin_active();
       cell != triangulation.end(); ++cell)
    for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
      if (cell->face(f)->center()[dim-1] == 0)
	{
	  cell->face(f)->set_boundary_indicator(1);
	  
// 	  for (unsigned int e=0; e<GeometryInfo<dim>::lines_per_face; ++e)
// 	    cell->face(f)->line(e)->set_boundary_indicator (1);
	}
  
  
  triangulation.refine_global (4-dim);

  for (unsigned int refinement_cycle = 0; refinement_cycle<7;
       ++refinement_cycle)
    {
      std::cout << "Refinement cycle " << refinement_cycle << std::endl;
      
      if (refinement_cycle > 0)
	refine_mesh ();
      
      setup_dofs ();

      std::cout << "   Assembling..." << std::endl;
      assemble_system ();      

      std::cout << "   Solving..." << std::endl;
      solve ();
      
      output_results (refinement_cycle);

      std::cout << std::endl;
    }
}

    

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
